##############################################################################
#  SSTMap: A Python library for the calculation of water structure and
#         thermodynamics on solute surfaces from molecular dynamics
#         trajectories.
# MIT License
# Copyright 2016-2017 Lehman College City University of New York and the Authors
#
# Authors: Kamran Haider, Steven Ramsay, Anthony Cruz Balberdy
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
###############################################################################
"""
This module contains implementations of a parent class for water analysis in
MD trajectories.
"""
##############################################################################
# Imports
##############################################################################


import subprocess
import shutil, sys
import numpy as np
import mdtraj as md
from scipy import spatial
import mdtraj as md
from sstmap.utils import *
from sstmap.water_analysis import WaterAnalysis
import _sstmap_ext as calc
import _sstmap_entropy as ext1
import _sstmap_probableconfig as ext2



##############################################################################
# SiteWaterAnalysis class definition
##############################################################################


class SiteWaterAnalysis(WaterAnalysis):
    """
    Hydration site analysis calculations for water molecules on solute surfaces
    in a molecular dynamics simulation.
    """

    @function_timer
    def __init__(self, topology_file, trajectory, start_frame=0, num_frames=None,
                 supporting_file=None, ligand_file=None, hsa_region_radius=5.0, clustercenter_file=None,
                 rho_bulk=0.0334, prefix="hsa"):
        """Initialize a SiteWaterAnalysis object for site-based solvation structure
        and thermodynamic calculations.

        Parameters
        ----------
        topology_file : string
            Filename for the system topology file.
        trajectory : string
            Filename for the molecular dynamics trajectory.
        start_frame : int, optional
            The frame index from which the calculations will begin. Default: 0
        num_frames : int, optional
            The total number of frames or the length of simulation over which
            calculations will be performed. Default: 10000
        supporting_file : None, optional
            Filename of additional file containing non-bonded parameters for
            every particle in the system. Default: None
        ligand_file : None, optional
            Filename for a PDB file containing coordinates of a co-crystallized
            ligand. This is used by the clustering method to identify binding
            site water molecules in the trjectory.
        hsa_region_radius : float, optional
            Distance cutoff (in Angstrom) used to identify hsa region. All waters within this
            distance from any of the ligand atom are included in the analysis. Default: 5.0
        clustercenter_file : None, optional
            Filename for a PDB file containing a set of pre-generated cluster
            centers. If not provided, default clustering method would generate
            a set of clusters.
        prefix : string, optional
            A prefix for all results and intermediate files generated during
            calculations.
        rho_bulk : float
            Reference bulk water density to be used in calculations. Default: None
        """
        print("Initializing ...")
        self.start_frame = start_frame
        self.num_frames = num_frames
        super(SiteWaterAnalysis, self).__init__(topology_file, trajectory, supporting_file)

        self.prefix = prefix
        self.site_waters = None
        if clustercenter_file is None and ligand_file is None:
            sys.exit("Please provide either a ligand file for clustering or\
                        a cluster center file to generate hydration sites.")
        if self.num_frames == 0:
            sys.exit("Number of frames = %d, no calculations will be performed" % self.num_frames)

        self.rho_bulk = float(rho_bulk)
        self.ligand = ligand_file
        self.clustercenter_file = clustercenter_file
        self.hsa_region_radius = hsa_region_radius * 0.1
        if hsa_region_radius > 10.0:
            print("Warning: Currently, clustering region is restricted to a 10.0A sphere around the ligand molecule.")
            self.hsa_region_radius = 10.0
        self.hsa_data = None
        self.hsa_dict = None
        self.site_waters = None
        self.is_site_waters_populated = False
        self.hsa_region_O_ids = []
        self.hsa_region_flat_ids = []
        self.hsa_region_water_coords = None
        self.data_titles = ["index", "x", "y", "z",
                            "nwat", "occupancy",
                            "Esw", "EswLJ", "EswElec",
                            "Eww", "EwwLJ", "EwwElec", "Etot", "Ewwnbr",
                            "TSsw_trans", "TSsw_orient", "TStot",
                            "Nnbrs", "Nhbww", "Nhbsw", "Nhbtot",
                            "f_hb_ww", "f_enc",
                            "Acc_ww", "Don_ww", "Acc_sw", "Don_sw",
                            "solute_acceptors", "solute_donors"]
        self.energy_ww_lr_breakdown = None
        self.angular_st_distribution = None

    @function_timer
    def initialize_hydration_sites(self, clustering_density_cutoff=2.0):
        """
        Generates hydration sites and initialize data structures for storing hydration site
        information for subsequent calculations. Hydration sites are generated either from
        a pre-spcified set of cluster centers or from a clustering algorithm.

        Notes
        -----
        This function initializes several attributes of SiteWaterAnalysis object, which are
        previously assumed to be set to None.
        """

        cluster_coords, self.site_waters = self.generate_clusters(clustering_density_cutoff, self.ligand,
                                                                  self.clustercenter_file)
        self.hsa_data, self.hsa_dict = self.initialize_site_data(cluster_coords)
        self.is_site_waters_populated = True

    @function_timer
    def initialize_site_data(self, cluster_coords):
        """Initializes data elements to store the results of site-based calculations

        Parameters
        ----------
        cluster_coords : np.ndarray, float, shape(N_sites, 3)
            Array of Cartesian coordinates for cluster centers where N is the
            total number of cluster centers.

        Returns
        -------
        site_array : np.ndarray, float, shape(N,len(self.data_titles))
            Array storing averaged strucure and thermodynamic quantities for each of
            the N cluster centers. The number and types of quantities calculatied are
            based on data_titles attribute of this class.
        site_dict : dictionary
            A dictionary containing full results for each quantity for each of the N
            cluster centers. The full results for a quantity consist of a list of all
            measurements obtained during the calculation.
        """
        n_sites = cluster_coords.shape[0]
        site_array = np.zeros((n_sites, len(self.data_titles)))
        site_dict = {}
        for site_i in range(n_sites):
            site_array[site_i, 0] = site_i
            site_array[site_i, 1] = cluster_coords[site_i, 0]
            site_array[site_i, 2] = cluster_coords[site_i, 1]
            site_array[site_i, 3] = cluster_coords[site_i, 2]
            site_dict[site_i] = [[] for i in range(len(self.data_titles))]
            site_dict[site_i].append(np.zeros((self.num_frames * 3, 3)))
        return site_array, site_dict

    @function_timer
    def generate_clusters(self, density_factor, ligand_file, clustercenter_file):
        """Generate hydration sites from water molecules found in the binding site
        during the simulation. Clustering is done in two steps; i). An initial clustering over a 10%
        of frames, and ii). A refinement step where all frames are used.

        Parameters

        ----------
        ligand_file : string
            Name of the PDB file containing atomic coordinates of the ligand,
            assumed to be co-crystallized with the protein.

        Returns
        -------
        final_cluster_coords : numpy.ndarray
            Coordinates of hydration sites, represented by a 2-D array with shape N x 3,
            where N is the number of hydration sites identified during clustering.

        site_waters : list
            List of N sub-lists where N is the number of identified hydration sites, each sublist
            consist of a 3-element tuple for every water identified in that site. First element of
            the tuple is frame number, second element is correct index of the oxygen atom in the
            the original topology and third element is the offset index as read from a version of
            a trimmed version trajectory for clustering.

        Notes
        -----
        The following attributes of the object are updated when the clustering is successfully completed.
        self.hsa_region_O_ids:
            The indices of water oxygen atoms in HSA region for each frame are stored
            in the corresponding lists.
        self.hsa_region_flat_ids:
            Same as above except that indices are not atom indices from the topology
            but in a sequence from 0 to N, where N is the total number of water oxygen atoms found in the
            HSA region throughout the simulation.
        self.hsa_region_water_coords:
            An N x 3 numpy array is initialized, where N is the total number of water water oxygen atoms found in the
            HSA region throughout the simulation. The array gets populated during individual frame processing.
        """
        sphere_radius = md.utils.in_units_of(1.0, "angstroms", "nanometers")
        if not self.topology_file.endswith(".h5"):
            topology = md.load_topology(self.topology_file)
        if self.non_water_atom_ids.shape[0] == 0:
            raise Exception(ValueError,
                            "Clustering is supported only for solute-solvent systems, no solute atoms found.")

        ligand = md.load_pdb(ligand_file, no_boxchk=True)
        ligand_coords = ligand.xyz[0, :, :]
        binding_site_atom_indices = np.asarray(list(range(ligand_coords.shape[0])))
        init_cluster_coords = None
        # Step 1: Initial Clustering if user didn't provide cluster centers
        if clustercenter_file is None:
            clustering_stride = 10
            print("Reading trajectory for clustering.")
            with md.open(self.trajectory) as f:
                f.seek(self.start_frame)
                # read all frames if no frames specified by user
                if self.num_frames is None:
                    if not self.trajectory.endswith(".h5"):
                        trj_short = f.read_as_traj(topology,
                            atom_indices=np.concatenate(
                            (binding_site_atom_indices, self.wat_oxygen_atom_ids)), stride=clustering_stride)
                    else:
                        trj_short = f.read_as_traj(atom_indices=np.concatenate(
                            (binding_site_atom_indices, self.wat_oxygen_atom_ids)), stride=clustering_stride)
                else:
                    if not self.trajectory.endswith(".h5"):
                        trj_short = f.read_as_traj(topology,
                            atom_indices=np.concatenate(
                            (binding_site_atom_indices, self.wat_oxygen_atom_ids)), n_frames=self.num_frames, stride=clustering_stride)
                    else:
                        trj_short = f.read_as_traj(atom_indices=np.concatenate(
                            (binding_site_atom_indices, self.wat_oxygen_atom_ids)), n_frames=self.num_frames, stride=clustering_stride)
                if trj_short.n_frames < 10:
                    sys.exit("Clustering requires at least 100 frames, current trajectory contains {0:d} frames.".format(
                            trj_short.n_frames))
                print("Performing an initial clustering over {0:d} frames.".format(trj_short.n_frames))
                # Obtain water molecules solvating the binding site
                # FIXME: This is a workaround to use MDTraj compute_neighbor function xyz coordinates of the trajectory are
                # modified such that first n atoms coordinates are switched to n atoms of ligand coordinates.
                # Unexpected things will happen if the number of solute atoms less than the number of ligand atoms, which is
                # highly unlikely.
                coords = trj_short.xyz
                for i_frame in range(trj_short.n_frames):
                    for pseudo_index in range(binding_site_atom_indices.shape[0]):
                        coords[i_frame, pseudo_index, :] = ligand_coords[pseudo_index, :]

                haystack = np.setdiff1d(trj_short.topology.select("all"), binding_site_atom_indices)
                binding_site_waters = md.compute_neighbors(trj_short, self.hsa_region_radius,
                                                           binding_site_atom_indices, haystack_indices=haystack)
                # generate a list of tuples, each tuple is a water and corresponding frame number in trj_short
                water_id_frame_list = [(i, nbr) for i in range(len(binding_site_waters)) for nbr in
                                       binding_site_waters[i]]

                # Start initial clustering by building a KDTree and get initial neighbor count for all waters
                water_coordinates = np.ma.array([coords[wat[0], wat[1], :] for wat in water_id_frame_list], mask=False)
                tree = spatial.cKDTree(water_coordinates)
                nbr_list = tree.query_ball_point(water_coordinates, sphere_radius)
                nbr_count_list = np.ma.array([len(nbrs) for nbrs in nbr_list], mask=False)
                cutoff = trj_short.n_frames * density_factor * 0.1401
                if np.ceil(cutoff) - cutoff <= 0.5:
                    cutoff = np.ceil(cutoff)
                else:
                    cutoff = np.floor(cutoff)
                n_wat = 3 * cutoff

                # Set up clustering loop
                cluster_list = []
                cluster_iter = 0
                while n_wat > cutoff:
                    # Get water with max nbrs and retrieve its neighbors and marked for exclusion in next iteration
                    max_index = np.argmax(nbr_count_list)
                    to_exclude = np.array(nbr_list[max_index])
                    # Set current water count to current neighbors plus one for the water itself
                    n_wat = len(to_exclude) + 1

                    # Mask current water, its neighbors so that they are not considered in the next iteration
                    nbr_count_list.mask[to_exclude] = True
                    nbr_count_list.mask[max_index] = True
                    # Mask current waters' and its neighbors' coords so that they are not considered in the next iteration
                    water_coordinates.mask[to_exclude] = True
                    water_coordinates.mask[max_index] = True

                    # Accumulate neighbors for each water in current cluster, removing common neighbors
                    nbrs_of_to_exclude = np.unique(np.array([n_excluded for excluded_nbrs in
                                                             nbr_list[to_exclude] for n_excluded in excluded_nbrs]))

                    # Obtain the list of waters whose neighbors need to be updated due to exclusion of the waters above
                    to_update = np.setxor1d(to_exclude, nbrs_of_to_exclude)
                    to_update = np.setdiff1d(to_update, np.asarray(max_index))

                    # Update the neighbor count for each water from the list generated above
                    if to_update.shape[0] != 0:
                        tree = spatial.cKDTree(water_coordinates)
                        updated_nbr_list = tree.query_ball_point(water_coordinates[to_update], sphere_radius)
                        # for each updated member, get its original index and update the original neighbor search list
                        for index, nbrs in enumerate(updated_nbr_list):
                            if not nbr_count_list.mask[to_update[index]]:
                                nbr_count_list[to_update[index]] = len(nbrs)

                    # Check distances with previously identified clusters and do not consider if within 1.2 A
                    # of an existing cluster
                    current_wat = water_id_frame_list[max_index]
                    current_wat_coords = md.utils.in_units_of(coords[current_wat[0], current_wat[1], :],
                                                              "nanometers", "angstroms")
                    near_flag = 0
                    if len(cluster_list) != 0:
                        for clust in cluster_list:
                            clust_coords = coords[clust[0], clust[1], :]
                            dist = np.linalg.norm(current_wat_coords - clust_coords)
                            if dist < 1.20:
                                near_flag += 1
                    if near_flag == 0:
                        cluster_iter += 1
                        cluster_list.append(water_id_frame_list[max_index])
                init_cluster_coords = [coords[cluster[0], cluster[1], :] for cluster in cluster_list]
        else:
            clusters_pdb_file = md.load_pdb(clustercenter_file, no_boxchk=True)
            init_cluster_coords = clusters_pdb_file.xyz[0, :, :]

        # Read full trajectory
        print("Reading trajectory to obtain water molecules for each cluster.")
        with md.open(self.trajectory) as f:
            f.seek(self.start_frame)
            if self.num_frames is None:
                if not self.trajectory.endswith(".h5"):
                    trj = f.read_as_traj(topology, stride=1,
                                     atom_indices=np.concatenate((binding_site_atom_indices, self.wat_oxygen_atom_ids)))
                    self.num_frames = trj.n_frames
                else:
                    trj = f.read_as_traj(stride=1,
                            atom_indices=np.concatenate((binding_site_atom_indices, self.wat_oxygen_atom_ids)))
                    self.num_frames = trj.n_frames
            else:
                if not self.trajectory.endswith(".h5"):
                    trj = f.read_as_traj(topology, n_frames=self.num_frames, stride=1,
                                     atom_indices=np.concatenate((binding_site_atom_indices, self.wat_oxygen_atom_ids)))
                else:
                    trj = f.read_as_traj(n_frames=self.num_frames, stride=1,
                                     atom_indices=np.concatenate((binding_site_atom_indices, self.wat_oxygen_atom_ids)))
                if trj.n_frames < self.num_frames:
                    print(("Warning: {0:d} frames found in the trajectory, resetting self.num_frames.".format(
                        trj.n_frames)))
                    self.num_frames = trj.n_frames
            for i_frame in range(trj.n_frames):
                for pseudo_index in range(binding_site_atom_indices.shape[0]):
                    trj.xyz[i_frame, pseudo_index, :] = ligand_coords[pseudo_index, :]
            haystack = np.setdiff1d(trj.topology.select("all"), binding_site_atom_indices)
            start_point = haystack[0]
            binding_site_waters = md.compute_neighbors(trj, self.hsa_region_radius,
                                                       binding_site_atom_indices, haystack_indices=haystack)
            # From the full frame-wise set of waters in the binding site, build two more frame-wise lists
            # one where each frame has a correct index of waters and another with a new index which ranges from
            # 0 to M, where M is the total number of hsa region waters - 1
            start = 0
            for i in range(len(binding_site_waters)):
                self.hsa_region_O_ids.append([])
                self.hsa_region_flat_ids.append([])
                for wat in binding_site_waters[i]:
                    wat_0 = wat - start_point
                    wat_offset = (wat_0 * self.water_sites) + self.wat_oxygen_atom_ids[0]
                    self.hsa_region_O_ids[i].append(wat_offset)
                    self.hsa_region_flat_ids[i].append(start)
                    start += 3

            water_id_frame_list = [(i, nbr) for i in range(len(binding_site_waters)) for nbr in binding_site_waters[i]]
            water_coordinates = np.array([trj.xyz[wat[0], wat[1], :] for wat in water_id_frame_list])

        # Initialize array that stores coordinates all water molecules in HSA region, used for entropy calcs
        self.hsa_region_water_coords = np.zeros((len(water_id_frame_list) * 3, 3), dtype=float)
        tree = spatial.cKDTree(water_coordinates)
        nbr_list = tree.query_ball_point(init_cluster_coords, sphere_radius)
        final_cluster_coords = []
        cutoff = int(self.num_frames * density_factor * 0.1401)
        if np.ceil(cutoff) - cutoff <= 0.5:
            cutoff = np.ceil(cutoff)
        else:
            cutoff = np.floor(cutoff)

        # apply refinement if user defined clusters not provided
        if clustercenter_file is None:
            # Step 2: Refinement
            # Initialize variables and data structures
            # Read in the trajectory but only first N solute atoms where N equals the number of ligand atoms
            # plus all water oxygen atoms
            # WARNING: This shifts indices of waters and once they are assigned to clusters, the indices need to
            # be corrected.

            print(("Refining initial cluster positions by considering {0:d} frames.".format(self.num_frames)))
            # For each cluster, set cluster center equal to geometric center of all waters in the cluster
            site_waters = []
            cluster_index = 1
            for cluster in nbr_list:
                cluster_water_coords = water_coordinates[cluster]
                if len(cluster) > cutoff:
                    near_flag = 0
                    waters_offset = [(water_id_frame_list[wat][0] + self.start_frame,
                                      ((water_id_frame_list[wat][1] - start_point) * self.water_sites)
                                      + self.wat_oxygen_atom_ids[0]) for wat in cluster]

                    com = np.zeros(3)
                    masses = np.ones(cluster_water_coords.shape[0])
                    masses /= masses.sum()
                    com[:] = water_coordinates[cluster].T.dot(masses)
                    cluster_center = com[:]
                    # Raise flag if the current cluster center is within 1.2 A of existing cluster center
                    for other, coord in enumerate(final_cluster_coords[:-1]):
                        dist = np.linalg.norm(md.utils.in_units_of(cluster_center, "nanometers", "angstroms") - coord)
                        if dist < 1.20:
                            near_flag += 1
                    # Only add cluster center if it is at a safe distance from others
                    if near_flag == 0:
                        final_cluster_coords.append(md.utils.in_units_of(cluster_center, "nanometers", "angstroms"))
                        site_waters.append(waters_offset)
                        cluster_index += 1
        # otherwise store data for each user defined cluster
        else:
            # For each cluster, set cluster center equal to geometric center of all waters in the cluster
            final_cluster_coords = md.utils.in_units_of(init_cluster_coords, "nanometers", "angstroms")
            site_waters = []
            cluster_index = 1
            for cluster in nbr_list:
                waters_offset = [(water_id_frame_list[wat][0] + self.start_frame,
                                  ((water_id_frame_list[wat][1] - start_point) * self.water_sites)
                                  + self.wat_oxygen_atom_ids[0]) for wat in cluster]
                site_waters.append(waters_offset)
                cluster_index += 1

        # Write clustercenter file
        write_watpdb_from_coords("clustercenterfile", final_cluster_coords)
        self.clustercenter_file = "clustercenterfile.pdb"
        print(("Final number of clusters: {0:d}".format(len(final_cluster_coords))))
        return np.asarray(final_cluster_coords), site_waters

    def _process_frame(self, trj, frame_i, energy, hbonds, entropy,
                       energy_lr_breakdown, angular_structure,
                       shell_radii, r_theta_cutoff):
        """Calculates hydration site properties for a given frame.

        Parameters
        ----------
        trj : mdtraj.trajectory
            A trajectory object containing only one frame.
        frame_i : int
            Index of the frame to be processed
        energy : bool
            Flag for energy calculations
        hbonds : bool
            Flag for hydrogen bond calculations
        entropy :bool
            Flag for entropy calculations

        Returns
        -------
        None : NoneType

        """

        site_waters_copy = list(self.site_waters)
        nbr_cutoff_sq = 3.5 ** 2
        nbr_cutoff_index = -1
        nbr_cutoff_found = False
        shell_radii_copy = list()
        # We assume that shell_radii contains squared radii
        if type(shell_radii) == type(None):
            shell_radii_copy.append([nbr_cutoff_sq])
            nbr_cutoff_index = 0
        else:
            for shell in shell_radii:
                shell_radii_copy.append(shell)
            shell_radii_copy.sort()
            # The shell nbr_cutoff_sq will be found at shell_radii[nbr_cutoff_index].
            # Make sure nbr_cutoff_sq is in shell_radii and we know nbr_cutoff_index
            for radius_i, radius in enumerate(shell_radii_copy):
                if radius > nbr_cutoff_sq-0.0001 \
                and radius < nbr_cutoff_sq+0.0001:
                    nbr_cutoff_index = radius_i
                    nbr_cutoff_found = True
            if nbr_cutoff_index == -1:
                shell_radii_copy.append(nbr_cutoff_sq)
                shell_radii_copy.sort()
                for radius_i, radius in enumerate(shell_radii_copy):
                    if radius > nbr_cutoff_sq-0.0001 \
                    and radius < nbr_cutoff_sq+0.0001:
                        nbr_cutoff_index = radius_i
        # TODO: Use a robust unit conversion approach
        trj.xyz *= 10.0
        coords = trj.xyz
        trj.unitcell_lengths *= 10.0
        uc = trj.unitcell_vectors[0]*10.

        # Iterate over each site in the current frame if it has a water present
        distance_matrix = np.zeros((self.water_sites, self.all_atom_ids.shape[0]))
        for site_i in range(self.hsa_data.shape[0]):
            wat_O = None
            if self.is_site_waters_populated:
                if len(site_waters_copy[site_i]) != 0:
                    if site_waters_copy[site_i][0][0] == frame_i:
                        wat_O = site_waters_copy[site_i].pop(0)[1]
                        index = int(self.hsa_data[site_i, 4]) * 3
                        index_pairs = list(zip(list(range(wat_O, wat_O + 3)), list(range(index, index + 3))))
                        for index_pair in index_pairs:
                            self.hsa_dict[site_i][-1][index_pair[1]] += coords[0, index_pair[0], :]
                        self.hsa_data[site_i, 4] += 1

            if wat_O is not None and (energy or hbonds):
                ### The self.neighbor_ids array contains atom indices of all atoms that should
                ### be considered as potential neighbors. Valid_neighbors has same shape as 
                ### self.neighbor_ids and is False at the position where index wat occures in
                ### self.neighbor_ids, otherwise it is True. neighbor_ids stores the indices of
                ### the actual neighbor candidates that will be commited to get_pairwise_distances
                ### routine and has length of self.neighbor_ids-1. wat_nbrs_shell is of length neighbor_ids
                ### and holds the shell_index of each neighbor candidate atom (0:first shell, 1: second shell,
                ### ...,S: >highest shell). shell_radii_copy is an ordered (ascending) array
                ### with all S neighbor shell radii.
                valid_neighbors = np.ones(self.neighbor_ids.shape[0], dtype=bool)
                valid_neighbors[np.where(self.neighbor_ids==wat_O)] = False
                neighbor_ids   = self.neighbor_ids[valid_neighbors]
                wat_nbrs_shell = self.wat_nbrs_shell[valid_neighbors]
                calc.get_pairwise_distances(np.asarray([site_i, wat_O]), self.all_atom_ids,
                                            np.array(shell_radii_copy), neighbor_ids, wat_nbrs_shell,
                                            coords, uc, distance_matrix, 0)
                wat_nbrs = neighbor_ids[np.where(wat_nbrs_shell<(nbr_cutoff_index+1))]
                self.hsa_dict[site_i][17].append(wat_nbrs.shape[0])
                if energy:
                    e_lj_array, e_elec_array = np.copy(self.acoeff), np.copy(self.chg_product)
                    calc.calculate_energy(wat_O, distance_matrix, e_elec_array, e_lj_array, self.bcoeff)

                    e_lj_sw = np.sum(e_lj_array[:, self.non_water_atom_ids])
                    e_elec_sw = np.sum(e_elec_array[:, self.non_water_atom_ids])

                    e_lj_ww_left = e_lj_array[:, self.wat_oxygen_atom_ids[0]:wat_O]
                    e_lj_ww_right = e_lj_array[:, wat_O + self.water_sites:]
                    e_lj_ww = np.sum(e_lj_ww_left) + np.sum(e_lj_ww_right)
                    e_elec_ww_left = e_elec_array[:, self.wat_oxygen_atom_ids[0]:wat_O]
                    e_elec_ww_right = e_elec_array[:, wat_O + self.water_sites:]
                    e_elec_ww = np.sum(e_elec_ww_left) + np.sum(e_elec_ww_right)
                    e_nbr_list = [
                        np.sum(e_lj_array[:, nbr:nbr + self.water_sites] + e_elec_array[:, nbr:nbr + self.water_sites])
                        for nbr in wat_nbrs]

                    e_lj_ww = np.sum(
                        e_lj_array[:, self.wat_oxygen_atom_ids[0]:wat_O]) + np.sum(
                        e_lj_array[:, wat_O + self.water_sites:])
                    e_elec_ww = np.sum(
                        e_elec_array[:, self.wat_oxygen_atom_ids[0]:wat_O]) + np.sum(
                        e_elec_array[:, wat_O + self.water_sites:])

                    self.hsa_dict[site_i][7].append(e_lj_sw)
                    self.hsa_dict[site_i][8].append(e_elec_sw)
                    self.hsa_dict[site_i][10].append(e_lj_ww)
                    self.hsa_dict[site_i][11].append(e_elec_ww)
                    self.hsa_dict[site_i][6].append(e_lj_sw + e_elec_sw)
                    self.hsa_dict[site_i][9].append(e_lj_ww + e_elec_ww)
                    self.hsa_dict[site_i][12].append(e_lj_sw + e_elec_sw + e_lj_ww + e_elec_ww)
                    self.hsa_dict[site_i][13].extend(e_nbr_list)  # print(e_lj_sw/2.0)

                    if energy_lr_breakdown:
                        if nbr_cutoff_found:
                            self.energy_ww_lr_breakdown[site_i][nbr_cutoff_index] += sum(e_nbr_list)
                        shift = 0
                        for radius_i, radius in enumerate(shell_radii_copy):
                            if radius_i == nbr_cutoff_index:
                                continue
                            wat_nbrs_shell_i = neighbor_ids[np.where(wat_nbrs_shell==radius_i)]
                            if not nbr_cutoff_found \
                            and radius_i == nbr_cutoff_index+1:
                                wat_nbrs_shell_i = np.concatenate((wat_nbrs_shell_i, wat_nbrs))
                                shift = 1
                            if wat_nbrs_shell_i.shape[0] > 0:
                                e_nbr_list = [
                                    np.sum(e_lj_array[:, nbr:nbr + self.water_sites] + e_elec_array[:, nbr:nbr + self.water_sites])
                                    for nbr in wat_nbrs_shell_i]
                                self.energy_ww_lr_breakdown[site_i][radius_i-shift] += sum(e_nbr_list)

                    ### Might be usefull for API to have the neighbors and shell
                    ### indices available at all time.
                    self.wat_nbrs_shell[valid_neighbors] = wat_nbrs_shell
                    self.neighbor_ids[valid_neighbors]   = neighbor_ids

                if hbonds:
                    hbtot = 0
                    prot_nbrs_all = self.prot_atom_ids[
                        np.where(distance_matrix[0, :][self.prot_atom_ids] <= nbr_cutoff_sq)]
                    prot_nbrs_hb = prot_nbrs_all[np.where(self.prot_hb_types[prot_nbrs_all] != 0)]
                    if wat_nbrs.shape[0] > 0:
                        hb_ww = self.calculate_hydrogen_bonds(trj, wat_O, wat_nbrs)
                        acc_ww = hb_ww[:, 0][np.where(hb_ww[:, 0] == wat_O)].shape[0]
                        don_ww = hb_ww.shape[0] - acc_ww
                        self.hsa_dict[site_i][18].append(hb_ww.shape[0])
                        self.hsa_dict[site_i][23].append(acc_ww)
                        self.hsa_dict[site_i][24].append(don_ww)
                        hbtot += hb_ww.shape[0]
                        if wat_nbrs.shape[0] != 0 and hb_ww.shape[0] != 0:
                            self.hsa_dict[site_i][21].append(
                                hb_ww.shape[0] / wat_nbrs.shape[0])
                    if prot_nbrs_hb.shape[0] > 0:
                        hb_sw = self.calculate_hydrogen_bonds(trj, wat_O, prot_nbrs_hb, water_water=False)
                        acc_sw = hb_sw[:, 0][np.where(hb_sw[:, 0] == wat_O)].shape[0]
                        don_sw = hb_sw.shape[0] - acc_sw
                        don_sw_ids = hb_sw[:, 1][np.where(hb_sw[:, 0] == wat_O)]
                        acc_sw_ids = hb_sw[:, 0][np.where(hb_sw[:, 0] != wat_O)]
                        self.hsa_dict[site_i][19].append(hb_sw.shape[0])
                        self.hsa_dict[site_i][25].append(acc_sw)
                        self.hsa_dict[site_i][26].append(don_sw)
                        self.hsa_dict[site_i][27].extend(acc_sw_ids)
                        self.hsa_dict[site_i][28].extend(don_sw_ids)
                        hbtot += hb_sw.shape[0]
                    self.hsa_dict[site_i][20].append(hbtot)
                    if angular_structure:
                        wat_nbrs = self.wat_oxygen_atom_ids[np.where(
                            (distance_matrix[0, :][self.wat_oxygen_atom_ids] <= r_theta_cutoff**2) & (
                                    distance_matrix[0, :][self.wat_oxygen_atom_ids] > 0.0))]
                        angles = self.water_nbr_orientations(trj, wat_O, wat_nbrs)
                        dist = np.sqrt(distance_matrix[0, wat_nbrs])
                        self.angular_st_distribution[site_i].extend(zip(dist, angles))

        if entropy:
            # save coordinates of hsa region waters in current frame
            for index, wat_O in enumerate(self.hsa_region_O_ids[frame_i - self.start_frame]):
                flat_id = self.hsa_region_flat_ids[frame_i - self.start_frame][index]
                index_pairs = list(zip(list(range(wat_O, wat_O + 3)), list(range(flat_id, flat_id + 3))))
                for index_pair in index_pairs:
                    self.hsa_region_water_coords[index_pair[1], :] += coords[0, index_pair[0], :]

    @function_timer
    def calculate_site_quantities(self, energy=True, entropy=True, hbonds=True,
                                        energy_lr_breakdown=False, angular_structure=False,
                                        shell_radii=None, r_theta_cutoff=6.0):
        """
        Performs site-based solvation thermodynamics and structure calculations by iterating
        over frames in the trajectory. If water molecules in hydration sites are already determined
        (the case when clustering is already done), then the list of hydration site waters in
        each frame is used to iterate over each water and calculate its properties. If externally
        determined hydration sites are provided (when self.clustercenter_file is set to a pdb file of
        hydration sites) then for each site, corresponding water is found in each frame and is used
        for caclulations.

        Parameters
        ----------
        energy : bool, optional
            Description
        hbonds : bool, optional
            Description
        entropy : bool, optional
            Description

        Returns
        -------
        None : NoneType
            This function updates hydration site data structures to store the results of calculations.
        """
        print_progress_bar(0, self.num_frames)
        if not self.trajectory.endswith(".h5"):
            topology = md.load_topology(self.topology_file)
        read_num_frames = 0
        if energy_lr_breakdown:
            if shell_radii is None:
                shell_radii = [3.5, 5.5, 8.5]
            #else:
            #    assert len(shell_radii) == 3, "Water-water energy decomposition supported only upto 3 solvation shells." \
            #                                  "Please provide outer radii for three shells."
            ### We assume that this is sorted, right?
            ### Otherwise _process_frame should give strange restults.
            shell_radii = [i**2 for i in shell_radii]
            ### What about interactions > shell_radii[-1] ?
            self.energy_ww_lr_breakdown = [[0.0 for s in shell_radii] for i in range(self.hsa_data.shape[0])]

        if angular_structure:
            if r_theta_cutoff > 8.0:
                print("Warning: r_theta_cutoff > 8.0 can take a long time."
                      "Resetting angular structure distance cutoff to 8.0 Angstrom")
                r_theta_cutoff = 8.0
            self.angular_st_distribution = [[] for i in range(self.hsa_data.shape[0])]

        with md.open(self.trajectory) as f:
            for frame_i in range(self.start_frame, self.start_frame + self.num_frames):
                print_progress_bar(frame_i - self.start_frame, self.num_frames)
                f.seek(frame_i)
                if not self.trajectory.endswith(".h5"):
                    trj = f.read_as_traj(topology, n_frames=1, stride=1)
                else:
                    trj = f.read_as_traj(n_frames=1, stride=1)

                if trj.n_frames == 0:
                    print("No more frames to read.")
                    break
                else:
                    self._process_frame(trj, frame_i, energy, hbonds, entropy,
                                        energy_lr_breakdown, angular_structure,
                                        shell_radii, r_theta_cutoff)
                    read_num_frames += 1
            if read_num_frames < self.num_frames:
                print(("{0:d} frames found in the trajectory, resetting self.num_frames.".format(read_num_frames)))
                self.num_frames = read_num_frames

        if entropy:
            self.generate_data_for_entropycalcs(self.start_frame, self.num_frames)
            self.run_entropy_scripts()
        self.normalize_site_quantities(self.num_frames)



    @function_timer
    def generate_data_for_entropycalcs(self, start_frame, num_frames, user_defined_clusters=False):
        """
        """
        print("Writing PDB file containing all HSA region water molecules for entropy calculations.")
        write_watpdb_from_coords("within5Aofligand", self.hsa_region_water_coords, full_water_res=True)
        print("Done.")
        print("Writing PDB files for all water molecules in each hydration site.")
        for site_i in range(self.hsa_data.shape[0]):
            # print site_i, len(self.hsa_dict[site_i][-1])/3.0, self.hsa_data[site_i, 4]
            num_wat = int(self.hsa_data[site_i, 4]) * 3
            # print num_wat, self.hsa_dict[site_i][-1].shape
            cluster_name = '{0:06d}'.format(site_i + 1)
            write_watpdb_from_coords("cluster." + cluster_name, self.hsa_dict[site_i][-1][:num_wat, :],
                                     full_water_res=True)
        print("Done.")

    @function_timer
    def run_entropy_scripts(self, output_dir=None):
        """Rerurns list of trans and orient entropies for each cluster.

        Parameters
        ----------
        output_dir: string
            Name of the output directory
        """
        # run bruteclust in a newly created output directory
        curr_dir = os.getcwd()
        if output_dir is not None:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            else:
                shutil.rmtree(output_dir)
                os.makedirs(output_dir)

        else:
            output_dir = "entropy_output"
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            else:
                shutil.rmtree(output_dir)
                os.makedirs(output_dir)

        # prepare arguments for bruteclust, which are file paths
        input_c_arg = os.path.abspath(self.clustercenter_file)
        input_w_arg = os.path.abspath("within5Aofligand.pdb")
        # run bruteclust inside the output directory
        os.chdir(curr_dir + "/" + output_dir)
        print("Generating expanded cluster water files...")
        try:
            # subprocess.check_call("bruteclust  -c " + input_c_arg + " -w " + input_w_arg, shell=True)
            ext1.run_bruteclust(input_c_arg, input_w_arg)
        except Exception as e:
            print(e)
        os.chdir(curr_dir)
        # run entropy code, make sure previous output files do not exist
        trans_dat, orient_dat = os.path.abspath("trans.dat"), os.path.abspath("orient.dat")
        if os.path.isfile(trans_dat):
            os.remove(trans_dat)
        if os.path.isfile(orient_dat):
            os.remove(orient_dat)

        # run entropy code and generate most probable config
        input_o_arg = os.path.abspath(output_dir + "/probable.pdb")
        print("Running entropy calculation from extension module.")
        for site_i in range(self.hsa_data.shape[0]):
            cluster_filename = "cluster.{0:06d}.pdb".format(site_i + 1)
            input_i_arg = os.path.abspath(cluster_filename)
            input_e_arg = os.path.abspath(output_dir + "/" + cluster_filename)
            try:
                ext1.run_kdhsa102(input_i_arg, input_e_arg)
                ext2.run_probconfig(input_i_arg, input_o_arg)
            except Exception as e:
                print(e)

        a = np.loadtxt(input_o_arg, usecols=(6, 7, 8))
        write_watpdb_from_coords("probable_configs", a, full_water_res=True)
        # extract entropy data and put into summary data
        if os.path.isfile(trans_dat) and os.path.isfile(orient_dat):
            trans_ent, orient_ent = np.loadtxt(trans_dat), np.loadtxt(orient_dat)
            if trans_ent.shape[0] == self.hsa_data.shape[0] and orient_ent.shape[0] == self.hsa_data.shape[0]:
                self.hsa_data[:, 14] += trans_ent
                self.hsa_data[:, 15] += orient_ent
                self.hsa_data[:, 16] += trans_ent + orient_ent
        shutil.rmtree(output_dir)
        os.remove(input_w_arg)

    @function_timer
    def normalize_site_quantities(self, num_frames):
        """
        """
        sphere_volume = (4 / 3) * np.pi
        bulk_water_per_site = self.rho_bulk * sphere_volume * num_frames
        skip_normalization = ["index", "x", "y", "z", "nwat", "occupancy", "gO",
                              "TSsw_trans", "TSsw_orient", "TStot", "solute_acceptors", "solute_donors"]
        for site_i in range(self.hsa_data.shape[0]):
            n_wat = self.hsa_data[site_i, 4]
            if n_wat != 0:
                self.hsa_data[site_i, 5] = n_wat / (self.num_frames)
                for quantity_i in range(len(self.data_titles)):
                    if self.data_titles[quantity_i] not in skip_normalization:
                        if self.data_titles[quantity_i] in ["Esw", "EswLJ", "EswElec", "Eww", "EwwLJ", "EwwElec",
                                                            "Etot"]:
                            self.hsa_data[site_i, quantity_i] = (np.sum(
                                self.hsa_dict[site_i][quantity_i]) / n_wat) * 0.5
                        elif self.data_titles[quantity_i] in ["Ewwnbr"]:
                            if len(self.hsa_dict[site_i][17]) != 0:
                                self.hsa_data[site_i, quantity_i] = (np.sum(self.hsa_dict[site_i][quantity_i]) / len(
                                    self.hsa_dict[site_i][quantity_i])) * 0.5
                        else:
                            self.hsa_data[site_i, quantity_i] = np.sum(self.hsa_dict[site_i][quantity_i]) / n_wat
                    if self.data_titles[quantity_i] in ["solute_acceptors", "solute_donors"]:
                        self.hsa_dict[site_i][quantity_i] = np.unique(self.hsa_dict[site_i][quantity_i])
                    if self.data_titles[quantity_i] in ["f_enc"]:
                        self.hsa_data[site_i, quantity_i] = None
                if self.energy_ww_lr_breakdown is not None:
                    self.energy_ww_lr_breakdown[site_i] = [(shell_e / n_wat) * 0.5 for shell_e in self.energy_ww_lr_breakdown[site_i]]


    def print_system_summary(self):
        """Summary

        Returns
        -------
        TYPE
            Description
        """
        print("System information:")
        print(("\tParameter file: %s\n" % self.topology_file))
        print(("\tTrajectory: %s\n" % self.trajectory))
        print(("\tTotal Atoms: %d, Waters: %d, Solute Atoms: %d\n"
              % (self.all_atom_ids.shape[0], self.wat_oxygen_atom_ids.shape[0],
                 self.non_water_atom_ids.shape[0])))
        if self.hsa_data is not None:
            print(("\tNumber of clusters: %d\n" % len(self.hsa_data)))

    @function_timer
    def write_calculation_summary(self):
        """Write a summary of calculations in the form of a table of hydration sites and the average of the
        calculated quantities.


        Returns
        -------
        TYPE
            Description
        """
        with open(self.prefix + "_hsa_summary.txt", "w") as f:
            header = " ".join(self.data_titles) + "\n"
            f.write(header)

            # format first six columns
            formatted_output = "{0[0]:.0f} {0[1]:.2f} {0[2]:.2f} {0[3]:.2f} {0[4]:.0f} {0[5]:.2f} "

            # format site energetic, entropic and structural data
            for quantity_i in range(6, len(self.data_titles) - 2):
                formatted_output += "{0[%d]:.6f} " % quantity_i
            # format solute acceptors and donors
            formatted_output += "{1} {2}\n"
            for site_i in range(self.hsa_data.shape[0]):
                solute_acceptors = [str(self.topology.atom(acceptor))
                                    for acceptor in self.hsa_dict[site_i][27]]
                solute_donors = [str(self.topology.atom(donor))
                                 for donor in self.hsa_dict[site_i][28]]
                site_data_line = formatted_output.format(self.hsa_data[site_i, :],
                                                         ",".join(
                                                             solute_acceptors),
                                                         ",".join(solute_donors))
                f.write(site_data_line)

    @function_timer
    def write_data(self):
        """
        """

        skip_write_data = ["x", "y", "z", "nwat", "occupancy", "gO",
                           "TSsw_trans", "TSsw_orient", "TStot", "f_enc", "solute_acceptors", "solute_donors"]
        # create directory to store detailed data for individual columns in HSA
        directory = self.prefix + "_hsa_data"
        if not os.path.exists(directory):
            os.makedirs(directory)
        # for each cluster, go through time series data
        for site_i in range(self.hsa_data.shape[0]):
            site_index = "/%03d_" % site_i
            for quantity_i in range(len(self.data_titles)):
                if self.data_titles[quantity_i] not in skip_write_data and len(
                        self.hsa_dict[site_i][quantity_i]) != 0:
                    data_file_name = directory + site_index + self.prefix + \
                                     "_" + self.data_titles[quantity_i] + ".txt"
                    with open(data_file_name, "w") as data_file:
                        data_file.writelines("%s\n" % item for item in self.hsa_dict[
                            site_i][quantity_i])


    @function_timer
    def write_energy_ww_breakdown(self):
        """
        Note: Output is mean water-water interaction energy
              with all water molecules in within a shell
              normalized on number of water molecules in
              hydration site.
              Whereas the Ewwnbr quantity in the HSA output
              is the first shell (within 3.5 Ang) mean water-
              water interaction energy normliazed on number
              of water molecules in hydration site AND number
              of water molecules in first shell.
        """
        n_shells = len(self.energy_ww_lr_breakdown[0])
        if self.energy_ww_lr_breakdown is not None:
            with open(self.prefix + "_energy_ww_by_shell.txt", "w") as f:
                f.write("index ")
                for shell_i in range(n_shells):
                    f.write("shell_%d " %shell_i)
                f.write("\n")
                for site_i in range(self.hsa_data.shape[0]):
                    f.write("%d " %site_i)
                    for shell in self.energy_ww_lr_breakdown[site_i]:
                        f.write("%6.3f " %shell)
                    f.write("\n")
                f.write("\n")


    def write_angular_structure_distribution(self):
        directory = self.prefix + "_hsa_data"
        if not os.path.exists(directory):
            os.makedirs(directory)
        if self.angular_st_distribution is not None:
            for site_i in range(self.hsa_data.shape[0]):
                site_index = "/%03d_" % site_i
                data_file_name = directory + site_index + self.prefix + "_r_theta.txt"
                with open(data_file_name, "w") as data_file:
                    lines = ["{0[0]:.3f} {0[1]:.3f}\n".format(item) for item in self.angular_st_distribution[site_i]]
                    data_file.writelines(lines)
