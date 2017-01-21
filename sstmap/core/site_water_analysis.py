"""Summary
"""
##############################################################################
# SSTMap: A Python library for the calculation of water structure and
#         thermodynamics on solute surfaces from molecular dynamics
#         trajectories.
# Copyright 2016-2017 Lehman College City University of New York
# and the Authors
#
# Authors: Kamran Haider
# Contributors: Steven Ramsay, Anthony Cruz Balberdy
#
# SSTMap is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################
"""
This module contains implementations of a parent class for water analysis in
MD trajectories.
"""

##############################################################################
# Imports
##############################################################################
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from builtins import str
from builtins import range
from past.utils import old_div
import sys
import os

import numpy as np
from scipy import spatial
import mdtraj as md
from progressbar import Bar, Percentage, ProgressBar, ETA

from water_analysis import WaterAnalysis, NeighborSearch, function_timer
import _sstmap_ext as calc

##############################################################################
# SiteWaterAnalysis class definition
##############################################################################


class SiteWaterAnalysis(WaterAnalysis):

    """
    Hydration site analysis calculations for water molecules on solute surfaces
    in a simulation.
    
    Attributes
    ----------
    data_titles : TYPE
        Description
    prefix : TYPE
        Description
    site_waters : TYPE
        Description
    """
    @function_timer
    def __init__(self, topology_file, trajectory, start_frame=0, num_frames=0,
                 ligand_file=None, cluster_center_file=None,
                 desmond_helper_file=None, prefix="hsa"):
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
        ligand_file : None, optional
            Filename for a PDB file containing coordinates of a co-crystallized
            ligand. This is used by the clustering method to identify binding
            site water molecules in the trjectory.
        cluster_center_file : None, optional
            Filename for a PDB file containing a set of pre-generated cluster
            centers. If not provided, default clustering method would generate
            a set of clusters.
        desmond_helper_file : None, optional
            Filename for a pre-generated text file containing non-bonded
            parametersfor every particle in the system, applicable only
            when a Desmond MD trajectory and topology is to be processed.
            Default: None
        prefix : string, optional
            A prefix for all results and intermediate files generated during
            calculations.
        """
        print("Initializing...")
        super(SiteWaterAnalysis, self).__init__(
            topology_file, trajectory,
            start_frame, num_frames, desmond_helper_file)
        self.prefix = prefix
        self.site_waters = None
        if cluster_center_file is None and ligand_file is None:
            sys.exit("Please provide either a ligand file for clustering or\
                        a cluster center file to generate hydration sites.")
        if ligand_file is not None:
            cluster_coords, self.site_waters = self.generate_clusters(
                ligand_file)
            self.hsa_data, self.hsa_dict = self.initialize_site_data(
                cluster_coords)

        if cluster_center_file is not None:
            clusters_pdb_file = md.load_pdb(cluster_center_file)
            cluster_coords = md.utils.in_units_of(
                clusters_pdb_file.xyz[0, :, :], "nanometers", "angstroms")
            self.hsa_data, self.hsa_dict = self.initialize_site_data(
                cluster_coords)

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
        self.data_titles = ["index", "x", "y", "z",
                            "nwat", "occupancy",
                            "Esw", "EswLJ", "EswElec",
                            "Eww", "EwwLJ", "EwwElec", "Etot", "Ewwnbr",
                            "TSsw", "TSww", "TStot",
                            "Nnbrs", "Nhbww", "Nhbsw", "Nhbtot",
                            "f_hb_ww", "f_enc",
                            "Acc_ww", "Don_ww", "Acc_sw", "Don_sw",
                            "solute_acceptors", "solute_donors"]
        n_sites = cluster_coords.shape[0]
        site_array = np.zeros((n_sites, len(self.data_titles)))
        site_dict = {}
        for site_i in range(n_sites):
            site_array[site_i, 0] = site_i
            site_array[site_i, 1] = cluster_coords[site_i, 0]
            site_array[site_i, 2] = cluster_coords[site_i, 1]
            site_array[site_i, 3] = cluster_coords[site_i, 2]
            site_dict[site_i] = [[] for i in range(len(self.data_titles))]
        return site_array, site_dict

    @function_timer
    def generate_clusters(self, ligand_file):
        """Generate hydration sites from water molecules found in the binding site
        during the simulation.
        
        Parameters
        ----------
        ligand_file : string
            Name of the PDB file containing atomic coordinates of the ligand,
            assumed to be co-crystallized with the protein.
        
        Returns
        -------
        TYPE
            Description
        """
        # Obtain binding site solute atoms using ligand atom coordinates
        ligand = md.load_pdb(ligand_file)
        ligand_coords = ligand.xyz[0, :, :]
        binding_site_atom_indices = list(range(ligand_coords.shape[0]))
        # Obtain water molecules solvating the binding site
        stride = 10
        print("Reading in trajectory for clustering.")
        trj = md.load(self.trajectory, top=self.topology)
        for i_frame in range(trj.n_frames):
            for pseudo_index in range(ligand_coords.shape[0]):
                trj.xyz[i_frame, pseudo_index,
                        :] = ligand_coords[pseudo_index, :]

        trj_short = trj[
            self.start_frame:self.start_frame + trj.n_frames:stride]

        print(
            "Obtaining a superconfiguration of all water molecules\
            found in the binding site throught the trajectory.")
        binding_site_waters = md.compute_neighbors(
            trj_short, 0.50, binding_site_atom_indices,
            haystack_indices=self.wat_oxygen_atom_ids)
        # generate a list of all waters with their frame ids
        water_id_frame_list = [(i, nbr) for i in range(
            len(binding_site_waters)) for nbr in binding_site_waters[i]]
        # Set up clustering loop
        print("Performing clustering on the superconfiguration.")
        cutoff = trj_short.n_frames * 2 * 0.1401
        if np.ceil(cutoff) - cutoff <= 0.5:
            cutoff = np.ceil(cutoff)
        else:
            cutoff = np.floor(cutoff)
        n_wat = 3 * cutoff
        cluster_list = []
        cluster_iter = 0
        sphere_radius = 0.1
        # Build KDTree and get initial neighbor count for all waters
        water_coordinates = np.ma.array(
            [trj_short.xyz[wat[0], wat[1], :] for wat in water_id_frame_list],
            mask=False)
        tree = spatial.cKDTree(water_coordinates)
        nbr_list = tree.query_ball_point(water_coordinates, sphere_radius)
        nbr_count_list = np.ma.array([len(nbrs)
                                      for nbrs in nbr_list], mask=False)
        # Clustering loop
        while n_wat > cutoff:
            # get water with max nbrs and retrieve its nbrs, which are marked
            # for exclusion
            max_index = np.argmax(nbr_count_list)
            to_exclude = np.array(nbr_list[max_index])
            # set current water count to current neighbors plus one for the
            # water itself
            n_wat = len(to_exclude) + 1
            # Mask current water, its nbrs so that they are not considered in
            # the next iteration
            nbr_count_list.mask[to_exclude] = True
            nbr_count_list.mask[max_index] = True
            water_coordinates.mask[to_exclude] = True
            water_coordinates.mask[max_index] = True
            # For each members of this cluster, get its neighbors, create a
            # unique set out of this list
            nbrs_of_to_exclude = np.unique(
                np.array([n_excluded for excluded_nbrs in
                          nbr_list[to_exclude] for n_excluded in
                          excluded_nbrs]))
            # Remove original members of the cluster from this list to get the
            # final list of waters to update
            to_update = np.setxor1d(to_exclude, nbrs_of_to_exclude)
            to_update = np.setdiff1d(to_update, np.asarray(max_index))
            # Update the neighbor count for each water in the update list, if
            # update list is not empty
            if to_update.shape[0] != 0:
                tree = spatial.cKDTree(water_coordinates)
                updated_nbr_list = tree.query_ball_point(
                    water_coordinates[to_update], sphere_radius)
                # for each updated member, get its original index and update
                # the original neighbor search list
                for index, nbrs in enumerate(updated_nbr_list):
                    if not nbr_count_list.mask[to_update[index]]:
                        nbr_count_list[to_update[index]] = len(nbrs)
            # check distance with previously identified clusters and do not
            # consider if within 1.2A of an existing cluster
            current_wat = water_id_frame_list[max_index]
            current_wat_coords = md.utils.in_units_of(
                trj_short.xyz[current_wat[0], current_wat[1], :],
                "nanometers", "angstroms")
            near_flag = 0
            if len(cluster_list) != 0:
                for clust in cluster_list:
                    clust_coords = md.utils.in_units_of(
                        trj_short.xyz[clust[0], clust[1], :],
                        "nanometers", "angstroms")
                    dist = np.linalg.norm(current_wat_coords - clust_coords)
                    if round(dist, 2) <= 1.2:
                        near_flag += 1
            if near_flag == 0:
                cluster_iter += 1
                print("\tCluster iteration: ", cluster_iter)
                cluster_list.append(water_id_frame_list[max_index])

        self.write_wat_pdb(trj_short,
                           self.prefix + "_initial_clustercenterfile",
                           water_id_list=cluster_list)
        init_cluster_coords = [trj_short.xyz[cluster[0], cluster[1], :]
                               for cluster in cluster_list]
        print("Refining initial cluster positions by considering %d frames." %
              trj.n_frames)
        binding_site_waters = md.compute_neighbors(
            trj, 0.50, binding_site_atom_indices,
            haystack_indices=self.wat_oxygen_atom_ids)
        water_id_frame_list = [(i, nbr) for i in
                               range(len(binding_site_waters))
                               for nbr in binding_site_waters[i]]
        print("Writing all waters within 5 Angstrom of ligand to pdb file.")
        self.write_wat_pdb(trj, self.prefix + "_within5Aofligand",
                           water_id_list=water_id_frame_list)
        water_coordinates = np.array(
            [trj.xyz[wat[0], wat[1], :] for wat in water_id_frame_list])
        tree = spatial.cKDTree(water_coordinates)
        nbr_list = tree.query_ball_point(init_cluster_coords, sphere_radius)
        final_cluster_coords = []
        cutoff = int(trj.n_frames * 2 * 0.1401)
        if np.ceil(cutoff) - cutoff <= 0.5:
            cutoff = np.ceil(cutoff)
        else:
            cutoff = np.floor(cutoff)
        # for each cluster, set cluster center equal to geometric center of all
        # waters in the cluster
        site_waters = []
        for cluster in nbr_list:
            cluster_water_coords = water_coordinates[cluster]
            if len(cluster) > cutoff:
                waters = [water_id_frame_list[water] for water in cluster]
                site_waters.append(waters)
                com = np.zeros(3)
                masses = np.ones(cluster_water_coords.shape[0])
                masses /= masses.sum()
                com[:] = water_coordinates[cluster].T.dot(masses)
                cluster_center = com[:]
                final_cluster_coords.append(md.utils.in_units_of(
                    cluster_center, "nanometers", "angstroms"))

        self.write_wat_pdb(trj_short, self.prefix +
                           "_clustercenterfile",
                           wat_coords=final_cluster_coords)
        print("Final number of clusters: ", len(final_cluster_coords))
        return np.asarray(final_cluster_coords), site_waters

    @function_timer
    def calculate_site_quantities(self, energy=True, hbonds=True, entropy=True):
        '''
        TODO: replace TIP3P nbr count constant with variable that depends on water model
        
        Parameters
        ----------
        energy : bool, optional
            Description
        hbonds : bool, optional
            Description
        entropy : bool, optional
            Description
        '''
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA(
        )], maxval=self.start_frame + self.num_frames).start()
        for frame_i in range(self.start_frame, self.start_frame + self.num_frames):
            frame = md.load_frame(self.trajectory, frame_i, top=self.topology)
            pos = md.utils.in_units_of(frame.xyz, "nanometers", "angstroms")
            pbc = md.utils.in_units_of(
                frame.unitcell_lengths, "nanometers", "angstroms")
            oxygen_pos = pos[0, self.wat_oxygen_atom_ids, :]
            cluster_search_space = NeighborSearch(oxygen_pos, 1.0)
            water_search_space = NeighborSearch(oxygen_pos, 3.5)
            for site_i in range(self.hsa_data.shape[0]):
                wat_O = None
                if self.site_waters is not None:
                    if len(self.site_waters[site_i]) != 0:
                        if self.site_waters[site_i][0][0] == frame_i:
                            wat_O = self.site_waters[site_i].pop(0)[1]
                            self.hsa_data[site_i, 4] += 1
                else:
                    cluster_center_coords = (self.hsa_data[site_i, 1], self.hsa_data[
                                             site_i, 2], self.hsa_data[site_i, 3])
                    nbr_indices = cluster_search_space.query_nbrs_single_point(
                        cluster_center_coords)
                    cluster_wat_oxygens = [self.wat_oxygen_atom_ids[
                        nbr_index] for nbr_index in nbr_indices]
                    if len(cluster_wat_oxygens) != 0:
                        wat_O = cluster_wat_oxygens[0]
                        self.hsa_data[site_i, 4] += 1
                if wat_O is not None and (energy or hbonds):
                    distance_matrix = np.zeros(
                        (self.water_sites, self.all_atom_ids.shape[0]), np.float_)
                    calc.get_pairwise_distances(np.asarray(
                        [site_i, wat_O]), self.all_atom_ids, pos, pbc, distance_matrix)
                    wat_nbrs = self.wat_oxygen_atom_ids[np.where((distance_matrix[0, :][
                                                                 self.wat_oxygen_atom_ids] <= 3.5) & (distance_matrix[0, :][self.wat_oxygen_atom_ids] > 0.0))]
                    prot_nbrs = self.non_water_atom_ids[
                        np.where(distance_matrix[0, :][self.non_water_atom_ids] <= 3.5)]
                    prot_nbrs = np.asarray([prot_nbr for prot_nbr in prot_nbrs if self.topology.atom(
                        prot_nbr).name[0] not in ["C", "H"]])
                    self.hsa_dict[site_i][17].append(wat_nbrs.shape[0])
                    # TODO: 5.25 should be replaced with a variable
                    # 5.25 comes from TIP3P water model, define dictionary
                    # based on water residue names
                    self.hsa_dict[site_i][22].append(
                        1.0 - (old_div(wat_nbrs.shape[0], 5.25)))
                    if energy:
                        energy_lj, energy_elec = self.calculate_energy(
                            distance_matrix)
                        e_lj_sw = np.sum(
                            energy_lj[:self.wat_oxygen_atom_ids[0]:])
                        e_elec_sw = np.sum(
                            energy_elec[:, self.non_water_atom_ids])
                        e_lj_ww = np.nansum(
                            energy_lj[self.wat_oxygen_atom_ids[0]:])
                        e_elec_ww = np.sum(energy_elec[:, self.wat_atom_ids[
                                           0]:wat_O]) + np.sum(energy_elec[:, wat_O + self.water_sites:])
                        self.hsa_dict[site_i][7].append(e_lj_sw)
                        self.hsa_dict[site_i][8].append(e_elec_sw)
                        self.hsa_dict[site_i][10].append(e_lj_ww)
                        self.hsa_dict[site_i][11].append(e_elec_ww)
                        self.hsa_dict[site_i][6].append(e_lj_sw + e_elec_sw)
                        self.hsa_dict[site_i][9].append(e_lj_ww + e_elec_ww)
                        self.hsa_dict[site_i][12].append(
                            e_lj_sw + e_elec_sw + e_lj_ww + e_elec_ww)
                        # print "Solute-water LJ Energy of this water: ", e_lj_sw
                        # print "Solute-water Elec Energy of this water: ", e_elec_sw
                        # print "Water-water LJ Energy of this water: ", e_lj_ww
                        # print "Water-water Elec Energy of this water: ", e_elec_ww
                        # calc Enbr and other water structure metrics here
                        for nbr_i in wat_nbrs:
                            e_nbr_i = 0.0
                            e_nbr_i += energy_lj[self.wat_oxygen_atom_ids[0]                                                 :][old_div((nbr_i - self.wat_oxygen_atom_ids[0]), 3)]
                            for i in range(self.water_sites):
                                e_nbr_i += np.sum(energy_elec[:, nbr_i + i])
                        self.hsa_dict[site_i][13].append(e_nbr_i)
                    if hbonds:
                        hb_ww, hb_sw = self.calculate_hydrogen_bonds(
                            frame, wat_O, wat_nbrs, prot_nbrs)
                        acc_ww = hb_ww[:, 0][
                            np.where(hb_ww[:, 0] == wat_O)].shape[0]
                        don_ww = hb_ww.shape[0] - acc_ww
                        acc_sw = hb_sw[:, 0][
                            np.where(hb_sw[:, 0] == wat_O)].shape[0]
                        don_sw = hb_sw.shape[0] - acc_sw
                        # FIXME: Spurious atom names showing up in summary
                        don_sw_ids = hb_sw[:, 1][
                            np.where(hb_sw[:, 0] == wat_O)]
                        acc_sw_ids = hb_sw[:, 0][
                            np.where(hb_sw[:, 0] != wat_O)]
                        self.hsa_dict[site_i][18].append(hb_ww.shape[0])
                        self.hsa_dict[site_i][19].append(hb_sw.shape[0])
                        self.hsa_dict[site_i][20].append(
                            hb_ww.shape[0] + hb_sw.shape[0])
                        self.hsa_dict[site_i][23].append(acc_ww)
                        self.hsa_dict[site_i][24].append(don_ww)
                        self.hsa_dict[site_i][25].append(acc_sw)
                        self.hsa_dict[site_i][26].append(don_sw)
                        self.hsa_dict[site_i][27].extend(acc_sw_ids)
                        self.hsa_dict[site_i][28].extend(don_sw_ids)
                        if wat_nbrs.shape[0] != 0 and hb_ww.shape[0] != 0:
                            self.hsa_dict[site_i][28].append(
                                old_div(wat_nbrs.shape[0], hb_ww.shape[0]))
            pbar.update(frame_i)
        pbar.finish()
        # normalize site quantities
        print("Start Normalization")
        sphere_volume = (old_div(4, 3)) * np.pi
        bulk_water_per_site = self.rho_bulk * sphere_volume * self.num_frames
        skip_normalization = ["index", "x", "y", "z", "nwat", "occupancy", "gO",
                              "TSsw", "TSww", "TStot", "solute_acceptors", "solute_donors"]
        for site_i in range(self.hsa_data.shape[0]):
            n_wat = self.hsa_data[site_i, 4]
            if n_wat != 0:
                self.hsa_data[site_i, 5] = old_div(n_wat, \
                    (self.start_frame + self.num_frames))
                for quantity_i in range(len(self.data_titles)):
                    if self.data_titles[quantity_i] not in skip_normalization:
                        if self.data_titles[quantity_i] in ["Esw", "EswLJ", "EswElec", "Eww", "EwwLJ", "EwwElec", "Etot"]:
                            self.hsa_data[site_i, quantity_i] = (
                                old_div(np.sum(self.hsa_dict[site_i][quantity_i]), n_wat)) * 0.5
                        elif self.data_titles[quantity_i] in ["Ewwnbr"]:
                            self.hsa_data[site_i, quantity_i] = (old_div(np.sum(self.hsa_dict[site_i][
                                                                 quantity_i]), len(self.hsa_dict[site_i][quantity_i]))) * 0.5
                        else:
                            self.hsa_data[site_i, quantity_i] = old_div(np.sum(
                                self.hsa_dict[site_i][quantity_i]), n_wat)
                    if self.data_titles[quantity_i] in ["solute_acceptors", "solute_donors"]:
                        self.hsa_dict[site_i][quantity_i] = np.unique(
                            self.hsa_dict[site_i][quantity_i])
        print("End Normalization")

    @function_timer
    def calculate_angular_structure(self, site_indices=[], dist_cutoff=6.0):
        '''
        Returns energetic quantities for each hydration site
        
        Parameters
        ----------
        site_indices : list, optional
            Description
        dist_cutoff : float, optional
            Description
        
        '''

        if len(site_indices) == 0:
            site_indices = [int(i) for i in self.hsa_data[:, 0]]
        else:
            for index in site_indices:
                if index > self.hsa_data[:, 0][-1]:
                    sys.exit(
                        "Site %d does not exits, please provide valid site indices." % index)
        n_sites = len(site_indices)
        r_theta_data = [[] for i in range(n_sites)]
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA(
        )], maxval=self.start_frame + self.num_frames).start()
        for frame_i in range(self.start_frame, self.start_frame + self.num_frames):
            frame = md.load_frame(self.trajectory, frame_i, top=self.topology)
            pos = md.utils.in_units_of(frame.xyz, "nanometers", "angstroms")
            pbc = md.utils.in_units_of(
                frame.unitcell_lengths, "nanometers", "angstroms")
            oxygen_pos = pos[0, self.wat_oxygen_atom_ids, :]
            cluster_search_space = NeighborSearch(oxygen_pos, 1.0)
            for index, site_i in enumerate(site_indices):
                wat_O = None
                if self.site_waters is not None:
                    if len(self.site_waters[site_i]) != 0:
                        if self.site_waters[site_i][0][0] == frame_i:
                            wat_O = self.site_waters[site_i].pop(0)[1]
                            self.hsa_data[site_i, 4] += 1
                else:
                    cluster_center_coords = (self.hsa_data[site_i, 1], self.hsa_data[
                                             site_i, 2], self.hsa_data[site_i, 3])
                    nbr_indices = cluster_search_space.query_nbrs_single_point(
                        cluster_center_coords)
                    cluster_wat_oxygens = [self.wat_oxygen_atom_ids[
                        nbr_index] for nbr_index in nbr_indices]
                    if len(cluster_wat_oxygens) != 0:
                        wat_O = cluster_wat_oxygens[0]
                        self.hsa_data[site_i, 4] += 1
                if wat_O is not None:
                    distance_matrix = np.zeros(
                        (self.water_sites, self.all_atom_ids.shape[0]), np.float_)
                    calc.get_pairwise_distances(np.asarray(
                        [site_i, wat_O]), self.all_atom_ids, pos, pbc, distance_matrix)
                    wat_nbrs = self.wat_oxygen_atom_ids[np.where((distance_matrix[0, :][
                                                                 self.wat_oxygen_atom_ids] <= dist_cutoff) & (distance_matrix[0, :][self.wat_oxygen_atom_ids] > 0.0))]
                    angle_triplets = []
                    for wat_nbr in wat_nbrs:
                        angle_triplets.extend([[wat_O, wat_nbr, wat_nbr + 1], [wat_O, wat_nbr, wat_nbr + 2],
                                               [wat_nbr, wat_O, wat_O + 1], [wat_nbr, wat_O, wat_O + 2]])
                    angle_triplets = np.asarray(angle_triplets)
                    angles = md.utils.in_units_of(md.compute_angles(
                        frame, angle_triplets), "radians", "degrees")
                    for angle_index in range(0, angles.shape[1], 4):
                        hb_angle = np.min(
                            angles[0, angle_index:angle_index + 4])
                        nbr_index = old_div(angle_index, 4)
                        r_theta_data[index].append(
                            (hb_angle, distance_matrix[0, wat_nbrs[nbr_index]]))
            pbar.update(frame_i)
        pbar.finish()
        cwd = os.getcwd()
        directory = cwd + "/" + self.prefix + "_angular_structure_data"
        if not os.path.exists(directory):
            os.makedirs(directory)
        os.chdir(directory)
        for index, site_i in enumerate(site_indices):
            with open("%03d_r_theta.txt" % site_i, "w") as f:
                for d in r_theta_data[index]:
                    line = "{0[0]:.3f} {0[1]:.3f}\n".format(d)
                    f.write(line)
        os.chdir("../")

    @function_timer
    def calculate_lonranged_ww_energy(self, site_indices=[], shell_radii=[3.5, 5.5, 8.5]):
        """Summary
        
        Parameters
        ----------
        site_indices : list, optional
            Description
        shell_radii : list, optional
            Description
        
        Returns
        -------
        TYPE
            Description
        """
        if len(site_indices) == 0:
            site_indices = [int(i) for i in self.hsa_data[:, 0]]
        else:
            for index in site_indices:
                if index > self.hsa_data[:, 0][-1]:
                    sys.exit(
                        "Site %d does not exits, please provide valid site indices." % index)

        n_sites = len(site_indices)
        n_columns = 2 * len(shell_radii)
        shells = [(0.0, shell_radii[0])]
        for i in range(1, len(shell_radii)):
            shells.append((shell_radii[i - 1], shell_radii[i]))
        shells.append((shell_radii[-1], 100.0))
        longranged_data = np.zeros((n_sites, 2 * len(shells)))
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA(
        )], maxval=self.start_frame + self.num_frames).start()
        for frame_i in range(self.start_frame, self.start_frame + self.num_frames):
            frame = md.load_frame(self.trajectory, frame_i, top=self.topology)
            pos = md.utils.in_units_of(frame.xyz, "nanometers", "angstroms")
            pbc = md.utils.in_units_of(
                frame.unitcell_lengths, "nanometers", "angstroms")
            oxygen_pos = pos[0, self.wat_oxygen_atom_ids, :]
            cluster_search_space = NeighborSearch(oxygen_pos, 1.0)
            water_search_space = NeighborSearch(oxygen_pos, 3.5)
            for index, site_i in enumerate(site_indices):
                wat_O = None
                if self.site_waters is not None:
                    if len(self.site_waters[site_i]) != 0:
                        if self.site_waters[site_i][0][0] == frame_i:
                            wat_O = self.site_waters[site_i].pop(0)[1]
                            self.hsa_data[site_i, 4] += 1
                else:
                    cluster_center_coords = (self.hsa_data[site_i, 1], self.hsa_data[
                                             site_i, 2], self.hsa_data[site_i, 3])
                    nbr_indices = cluster_search_space.query_nbrs_single_point(
                        cluster_center_coords)
                    cluster_wat_oxygens = [self.wat_oxygen_atom_ids[
                        nbr_index] for nbr_index in nbr_indices]
                    if len(cluster_wat_oxygens) != 0:
                        wat_O = cluster_wat_oxygens[0]
                        self.hsa_data[site_i, 4] += 1
                if wat_O is not None:
                    distance_matrix = np.zeros(
                        (self.water_sites, self.all_atom_ids.shape[0]), np.float_)
                    calc.get_pairwise_distances(np.asarray(
                        [site_i, wat_O]), self.all_atom_ids, pos, pbc, distance_matrix)
                    energy_lj, energy_elec = self.calculate_energy(
                        distance_matrix)
                    for shell_index, shell in enumerate(shells):
                        wat_nbrs = self.wat_oxygen_atom_ids[np.where((distance_matrix[0, :][self.wat_oxygen_atom_ids] <= shell[
                                                                     1]) & (distance_matrix[0, :][self.wat_oxygen_atom_ids] > shell[0]))]
                        longranged_data[index, 2 *
                                        shell_index] += wat_nbrs.shape[0]
                        for i in range(self.water_sites):
                            longranged_data[
                                index, (2 * shell_index) + 1] += np.sum(energy_elec[:, wat_nbrs + i])
            pbar.update(frame_i)
        pbar.finish()
        # normalize site quantities
        for index, site_i in enumerate(site_indices):
            n_wat = self.hsa_data[site_i, 4]
            if n_wat != 0:
                longranged_data[index, :] /= n_wat * 2.0

        # write data
        with open(self.prefix + "_longrange_Eww_summary.txt", "w") as f:
            header = "index "
            formatted_output = "{0:.0f} "
            for shell_index, shell in enumerate(shells):
                if shell_index + 1 == len(shells):
                    header += "Nnbr_>" + \
                        str(shell_index) + " E_shell_>" + \
                        str(shell_index) + "\n"
                    formatted_output += "{1[%d]:.6f} {1[%d]:.6f}\n" % (
                        2 * shell_index, (2 * shell_index) + 1)
                else:
                    header += "Nnbr_" + \
                        str(shell_index + 1) + " E_shell_" + \
                        str(shell_index + 1) + " "
                    formatted_output += "{1[%d]:.6f} {1[%d]:.6f} " % (
                        2 * shell_index, (2 * shell_index) + 1)
            f.write(header)
            for index, site_i in enumerate(site_indices):
                n_wat = self.hsa_data[site_i, 4]
                site_data_line = formatted_output.format(
                    self.hsa_data[site_i, 0], longranged_data[index, :])
                f.write(site_data_line)

    def print_system_summary(self):
        """Summary
        
        Returns
        -------
        TYPE
            Description
        """
        print("System information:")
        print("\tParameter file: %s\n" % self.topology_file)
        print("\tTrajectory: %s\n" % self.trajectory)
        print("\tNumber of clusters: %d\n" % len(self.hsa_data))
        print("\tPeriodic Box: %s\n" % self.box_type)
        print("\tFrames: %d, Total Atoms: %d, Waters: %d, Solute Atoms: %d\n"
              % (self.num_frames, self.all_atom_ids.shape[0], self.wat_oxygen_atom_ids.shape[0], self.non_water_atom_ids.shape[0]))
        print("\tH-bond Donors: %d, H-bond Acceptors: %d, H-bond Donors & Acceptors: %d\n"
              % (self.solute_don_ids.shape[0], self.solute_acc_ids.shape[0], self.solute_acc_don_ids.shape[0]))

    @function_timer
    def write_summary(self):
        """Summary
        
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
        TODO: output energy quantities in half
        """
        skip_write_data = ["x", "y", "z", "nwat", "occupancy", "gO",
                           "TSsw", "TSww", "TStot", "solute_acceptors", "solute_donors"]
        cwd = os.getcwd()
        # create directory to store detailed data for individual columns in HSA
        directory = cwd + "/" + self.prefix + "_hsa_data"
        if not os.path.exists(directory):
            os.makedirs(directory)
        os.chdir(directory)
        # for each cluster, go through time series data
        for site_i in range(self.hsa_data.shape[0]):
            site_index = "%03d_" % site_i
            for quantity_i in range(len(self.data_titles)):
                if self.data_titles[quantity_i] not in skip_write_data and len(
                        self.hsa_dict[site_i][quantity_i]) != 0:
                    data_file_name = site_index + self.prefix + \
                        "_" + self.data_titles[quantity_i] + ".txt"
                    with open(data_file_name, "w") as data_file:
                        data_file.writelines("%s\n" % item for item in self.hsa_dict[
                                             site_i][quantity_i])
        os.chdir("../")

    @function_timer
    def write_wat_pdb(self, traj, filename, water_id_list=None, wat_coords=None, full_water_res=False):
        """Summary
        
        Parameters
        ----------
        traj : TYPE
            Description
        filename : TYPE
            Description
        water_id_list : None, optional
            Description
        wat_coords : None, optional
            Description
        full_water_res : bool, optional
            Description
        
        Returns
        -------
        TYPE
            Description
        """
        pdb_line_format = "{0:6}{1:>5} {2:<4}{3:<1}{4:>3} {5:1}{6:>4}{7:1}   {8[0]:>8.3f}{8[1]:>8.3f}{8[2]:>8.3f}{9:>6.2f}{10:>6.2f}\n"
        ter_line_format = "{0:3}   {1:>5}      {2:>3} {3:1}{4:4} \n"
        pdb_lines = []
        # write form the list of (water, frame) tuples
        if wat_coords is None:
            coords = traj.xyz
            if not full_water_res:
                # at_index, wat in enumerate(water_id_list):
                for at in range(len(water_id_list)):
                    wat = water_id_list[at]
                    at_index = at % 10000
                    res_index = at % 10000
                    wat_coords = md.utils.in_units_of(
                        coords[wat[0], wat[1], :], "nanometers", "angstroms")
                    #chain_id = possible_chains[chain_id_index]
                    chain_id = "A"
                    pdb_line = pdb_line_format.format(
                        "ATOM", at_index, "O", " ", "WAT", chain_id, res_index, " ", wat_coords, 0.00, 0.00)
                    pdb_lines.append(pdb_line)
                    # if full_water_res:
                    #    pdb_line = pdb_line_format.format("ATOM", at_index, "O", " ", "WAT", chain_id, res_index, " ", wat_coords, 0.00, 0.00)
                    #    pdb_line = pdb_line_format.format("ATOM", at_index, "O", " ", "WAT", chain_id, res_index, " ", wat_coords, 0.00, 0.00)
                    if res_index == 9999:
                        ter_line = ter_line_format.format(
                            "TER", at_index, "WAT", chain_id, res_index)
                        pdb_lines.append(ter_line)
                pdb_lines.append("END")

        # write directly from the array of coords
        else:
            for at in range(len(wat_coords)):
                wat_coord = wat_coords[at]
                at_index = at % 10000
                res_index = at % 10000
                chain_id = "A"
                pdb_line = pdb_line_format.format(
                    "ATOM", at_index, "O", " ", "WAT", chain_id, res_index, " ", wat_coord, 0.00, 0.00)
                pdb_lines.append(pdb_line)
                # if full_water_res:
                #    pdb_line = pdb_line_format.format("ATOM", at_index, "O", " ", "WAT", chain_id, res_index, " ", wat_coords, 0.00, 0.00)
                #    pdb_line = pdb_line_format.format("ATOM", at_index, "O", " ", "WAT", chain_id, res_index, " ", wat_coords, 0.00, 0.00)
                if res_index == 9999:
                    ter_line = ter_line_format.format(
                        "TER", at_index, "WAT", chain_id, res_index)
                    pdb_lines.append(ter_line)
        with open(filename + ".pdb", "w") as f:
            f.write("".join(pdb_lines))
