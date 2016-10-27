
import sys
import time
import datetime
import os

import numpy as np
from scipy import stats, spatial
import mdtraj as md
import parmed as pmd
from progressbar import Bar, Percentage, ProgressBar, ETA

from water_analysis import WaterAnalysis, NeighborSearch, function_timer
#import _sstmap_ext_old as calc_old
import _sstmap_ext as calc



class SiteWaterAnalysis(WaterAnalysis):
    """
    """
    def __init__(self, topology_file, trajectory, start_frame=0, num_frames=0,
                 ligand_file=None, cluster_center_file=None, 
                 desmond_helper_file=None, prefix="test"):
        
        super(SiteWaterAnalysis, self).__init__(topology_file, trajectory,
            start_frame, num_frames, desmond_helper_file)        
        self.prefix = prefix
        self.site_waters = None
        if cluster_center_file == None and ligand_file == None:
            sys.exit("Please provide either a ligand file for clustering or\
                        a cluster center file to generate hydration sites.")
        if ligand_file != None:
            cluster_coords, self.site_waters = self.generate_clusters(ligand_file)
            self.hsa_data, self.hsa_dict = self.initialize_site_data(cluster_coords)

        if cluster_center_file != None:
            clusters_pdb_file = md.load_pdb(cluster_center_file)
            cluster_coords = md.utils.in_units_of(clusters_pdb_file.xyz[0, :, :], "nanometers", "angstroms")
            self.hsa_data, self.hsa_dict = self.initialize_site_data(cluster_coords)

    def initialize_site_data(self, cluster_coords):
        self.data_titles = ["index", "x", "y", "z",
                            "nwat", "occupancy",
                            "Esw", "EswLJ", "EswElec", "Eww", "EwwLJ", "EwwElec", "Etot", "Ewwnbr",
                            "TSsw", "TSww", "TStot",
                            "Nnbrs", "Nhbww", "Nhbsw", "Nhbtot",
                            "f_hb_ww", "f_enc",
                            "Acc_ww", "Don_ww", "Acc_sw", "Don_sw",
                            "solute_acceptors", "solute_donors"]
        n_sites = cluster_coords.shape[0]
        site_array = np.zeros((n_sites, len(self.data_titles)))
        site_dict = {}
        for site_i in xrange(n_sites):
            site_array[site_i, 0] = site_i            
            site_array[site_i, 1] = cluster_coords[site_i, 0]
            site_array[site_i, 2] = cluster_coords[site_i, 1]
            site_array[site_i, 3] = cluster_coords[site_i, 2]
            site_dict[site_i] = [[] for i in range(len(self.data_titles))]
        return site_array, site_dict
    @function_timer
    def generate_clusters(self, ligand_file):
        # Obtain binding site solute atoms using ligand atom coordinates
        ligand = md.load_pdb(ligand_file)
        ligand_coords = md.utils.in_units_of(ligand.xyz[0, :, :], "nanometers", "angstroms", inplace=True)
        first_frame = md.load_frame(self.trajectory, 0, top=self.topology_file)
        solute_pos = md.utils.in_units_of(first_frame.xyz[0, self.non_water_atom_ids, :], "nanometers", "angstroms")
        search_space = NeighborSearch(solute_pos, 5.0)
        near_indices = search_space.query_nbrs_multiple_points(ligand_coords)
        binding_site_atom_indices = [self.non_water_atom_ids[nbr_index] for nbr_index in near_indices]
        # Obtain water molecules solvating the binding site
        stride = 10
        print "Reading in trajectory for clustering."
        trj = md.load(self.trajectory, top=self.topology)
        trj_short = trj[self.start_frame:self.start_frame + trj.n_frames:stride]
        print "Obtaining a superconfiguration of all water molecules found in the binding site throught the trajectory."
        binding_site_waters = md.compute_neighbors(trj_short, 0.50, binding_site_atom_indices, haystack_indices=self.wat_oxygen_atom_ids)
        # generate a list of all waters with their frame ids
        water_id_frame_list = [(i, nbr) for i in range(len(binding_site_waters)) for nbr in binding_site_waters[i]]
        # Set up clustering loop
        print "Performing clustering on the superconfiguration."
        cutoff = trj_short.n_frames * 2 * 0.1401
        if np.ceil(cutoff) - cutoff <= 0.5:
            cutoff = np.ceil(cutoff)
        else:
            cutoff = np.floor(cutoff)
        n_wat = 3 * cutoff
        cluster_list = []
        cluster_iter = 0
        sphere_radius = 1.0
        # Build KDTree and get initial neighbor count for all waters
        water_coordinates = np.ma.array([trj_short.xyz[wat[0], wat[1], :] for wat in water_id_frame_list], mask=False)*10.0        
        tree = spatial.cKDTree(water_coordinates)
        nbr_list = tree.query_ball_point(water_coordinates, sphere_radius)
        nbr_count_list = np.ma.array([len(nbrs) for nbrs in nbr_list], mask=False)
        # Clustering loop
        initial_time = time.time()
        while n_wat > cutoff:
            cluster_iter += 1
            # get water with max nbrs and retrieve its nbrs, which are later marked for exclusion
            print "\tCluster iteration: ", cluster_iter
            max_index = np.argmax(nbr_count_list)
            to_exclude = np.array(nbr_list[max_index])
            # set current water count to current neighbors plus one for the water itself
            n_wat = len(to_exclude) + 1
            cluster_list.append(water_id_frame_list[max_index])
            # Mask current water and its nbrs so that they are not considered in the next iteration
            nbr_count_list.mask[to_exclude] = True
            nbr_count_list.mask[max_index] = True
            # Create a version of water coordinates, after deleting the coordinates of waters set for exclsuion 
            updated_wat_coords = np.delete(water_coordinates, to_exclude, 0)
            updated_wat_coords = np.delete(updated_wat_coords, max_index, 0)
            
            # For each members of this cluster, get its neighbors, create a unique set out of this list
            nbrs_of_to_exclude = np.unique(np.array([n_excluded for excluded_nbrs in nbr_list[to_exclude] for n_excluded in excluded_nbrs]))
            # Remove original members of the cluster from this list to get the final list of waters to update
            to_update = np.setdiff1d(nbrs_of_to_exclude, to_exclude)
            to_update = np.setdiff1d(to_update, np.asarray(max_index))
            # Update the neighbor count for each water in the update list
            tree = spatial.cKDTree(updated_wat_coords)
            updated_nbr_list = tree.query_ball_point(water_coordinates[to_update], sphere_radius)
            # for each updated member, get its original index and update the original neighbor search list
            for index, nbrs in enumerate(updated_nbr_list):
                if not nbr_count_list.mask[to_update[index]]:
                    nbr_count_list[to_update[index]] = len(nbrs)
        init_cluster_coords = [trj_short.xyz[cluster[0], cluster[1], :]*10.0 for cluster in cluster_list]
        print "Initial number of clusters: ", len(init_cluster_coords)
        with open("initial_clustercenterfile.pdb", "w") as f:
            header = "REMARK Initial number of clusters: {0}\n".format(len(init_cluster_coords))
            f.write(header)
            for cluster in init_cluster_coords:
                pdb_line = "ATOM      1  O   WAT A   1      {0[0]:2.3f}  {0[1]:2.3f}  {0[2]:2.3f}  0.00  0.00\n".format(cluster)
                f.write(pdb_line)
        print "Refining initial cluster positions by considering %d frames." % trj.n_frames
        binding_site_waters = md.compute_neighbors(trj, 0.50, binding_site_atom_indices, haystack_indices=self.wat_oxygen_atom_ids)
        water_id_frame_list = [(i, nbr) for i in range(len(binding_site_waters)) for nbr in binding_site_waters[i]]
        water_coordinates = np.array([trj.xyz[wat[0], wat[1], :] for wat in water_id_frame_list])*10.0
        tree = spatial.cKDTree(water_coordinates)        
        nbr_list = tree.query_ball_point(init_cluster_coords, sphere_radius)
        final_cluster_coords = []
        cutoff = int(trj.n_frames * 2 * 0.1401)
        if np.ceil(cutoff) - cutoff <= 0.5:
            cutoff = np.ceil(cutoff)
        else:
            cutoff = np.floor(cutoff)
        # for each cluster, set cluster center equal to geometric center of all waters in the cluster
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
                final_cluster_coords.append(cluster_center)

        with open("clustercenterfile.pdb", "w") as f:
            header = "REMARK Initial number of clusters: {0}\n".format(
                len(final_cluster_coords))
            f.write(header)
            for cluster in final_cluster_coords:
                pdb_line = "ATOM      1  O   WAT A   1      {0[0]:2.3f}  {0[1]:2.3f}  {0[2]:2.3f}  0.00  0.00\n".format(cluster)
                f.write(pdb_line)
        print "Final number of clusters: ", len(final_cluster_coords)
        return np.asarray(final_cluster_coords), site_waters

    @function_timer
    def calculate_site_quantities(self, energy=True, hbonds=True, entropy=True):
        '''
        TODO: replace TIP3P nbr count constant with variable that depends on water model
        '''
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=self.start_frame + self.num_frames).start()
        for frame_i in xrange(self.start_frame, self.start_frame + self.num_frames):
            frame = md.load_frame(self.trajectory, frame_i, top=self.topology)
            pos = md.utils.in_units_of(frame.xyz, "nanometers", "angstroms")
            pbc = md.utils.in_units_of(frame.unitcell_lengths, "nanometers", "angstroms")
            oxygen_pos = pos[0, self.wat_oxygen_atom_ids, :]
            cluster_search_space = NeighborSearch(oxygen_pos, 1.0)
            water_search_space = NeighborSearch(oxygen_pos, 3.5)
            for site_i in xrange(self.hsa_data.shape[0]):
                wat_O = None
                if self.site_waters is not None:
                    if len(self.site_waters[site_i]) != 0:
                        if self.site_waters[site_i][0][0] == frame_i:
                            wat_O = self.site_waters[site_i].pop(0)[1]
                            self.hsa_data[site_i, 4] += 1
                else:
                    cluster_center_coords = (self.hsa_data[site_i, 1], self.hsa_data[site_i, 2], self.hsa_data[site_i, 3])
                    nbr_indices = cluster_search_space.query_nbrs_single_point(cluster_center_coords)
                    cluster_wat_oxygens = [self.wat_oxygen_atom_ids[nbr_index] for nbr_index in nbr_indices]
                    if len(cluster_wat_oxygens) != 0:
                        wat_O = cluster_wat_oxygens[0]
                        self.hsa_data[site_i, 4] += 1
                if wat_O is not None and (energy or hbonds):
                    distance_matrix = np.zeros((self.water_sites, self.all_atom_ids.shape[0]), np.float_)
                    calc.get_pairwise_distances(np.asarray([site_i, wat_O]), self.all_atom_ids, pos, pbc, distance_matrix)
                    wat_nbrs = self.wat_oxygen_atom_ids[np.where((distance_matrix[0, :][self.wat_oxygen_atom_ids] <= 3.5) & (distance_matrix[0, :][self.wat_oxygen_atom_ids] > 0.0))]
                    prot_nbrs = self.non_water_atom_ids[np.where(distance_matrix[0, :][self.non_water_atom_ids] <= 3.5)]
                    prot_nbrs = np.asarray([prot_nbr for prot_nbr in prot_nbrs if self.topology.atom(prot_nbr).name[0] not in ["C", "H"]])
                    self.hsa_dict[site_i][17].append(wat_nbrs.shape[0])
                    #TODO: 5.25 should be replaced with a variable
                    self.hsa_dict[site_i][22].append(1.0 - (wat_nbrs.shape[0]/5.25))
                    if energy:
                        energy_lj, energy_elec = self.calculate_energy(distance_matrix)
                        e_lj_sw = np.sum(energy_lj[:self.wat_oxygen_atom_ids[0]:])
                        e_elec_sw = np.sum(energy_elec[:, self.non_water_atom_ids])
                        e_lj_ww = np.nansum(energy_lj[self.wat_oxygen_atom_ids[0]:])
                        e_elec_ww = np.sum(energy_elec[:, self.wat_atom_ids[0]:wat_O]) + np.sum(energy_elec[:, wat_O + self.water_sites:])
                        self.hsa_dict[site_i][7].append(e_lj_sw)
                        self.hsa_dict[site_i][8].append(e_elec_sw)
                        self.hsa_dict[site_i][10].append(e_lj_ww)
                        self.hsa_dict[site_i][11].append(e_elec_ww)
                        self.hsa_dict[site_i][6].append(e_lj_sw + e_elec_sw)
                        self.hsa_dict[site_i][9].append(e_lj_ww + e_elec_ww)
                        self.hsa_dict[site_i][12].append(e_lj_sw + e_elec_sw + e_lj_ww + e_elec_ww)
                        #print "Solute-water LJ Energy of this water: ", e_lj_sw
                        #print "Solute-water Elec Energy of this water: ", e_elec_sw
                        #print "Water-water LJ Energy of this water: ", e_lj_ww
                        #print "Water-water Elec Energy of this water: ", e_elec_ww
                        # calc Enbr and other water structure metrics here
                        for nbr_i in wat_nbrs:
                            e_nbr_i = 0.0
                            e_nbr_i += energy_lj[self.wat_oxygen_atom_ids[0]:][(nbr_i - self.wat_oxygen_atom_ids[0])/3]
                            for i in xrange(self.water_sites):
                                e_nbr_i += np.sum(energy_elec[:, nbr_i + i])
                        self.hsa_dict[site_i][13].append(e_nbr_i)
                    if hbonds:                        
                        hb_ww, hb_sw = self.calculate_hydrogen_bonds(frame, wat_O, wat_nbrs, prot_nbrs)
                        acc_ww = hb_ww[:, 0][np.where(hb_ww[:, 0] == wat_O)].shape[0]
                        don_ww = hb_ww.shape[0] - acc_ww
                        acc_sw = hb_sw[:, 0][np.where(hb_sw[:, 0] == wat_O)].shape[0]
                        don_sw = hb_sw.shape[0] - acc_sw
                        # FIXME: Spurious atom names showing up in summary
                        don_sw_ids = hb_sw[:, 1][np.where(hb_sw[:, 0] == wat_O)]
                        acc_sw_ids = hb_sw[:, 0][np.where(hb_sw[:, 0] != wat_O)]
                        self.hsa_dict[site_i][18].append(hb_ww.shape[0])
                        self.hsa_dict[site_i][19].append(hb_sw.shape[0])
                        self.hsa_dict[site_i][20].append(hb_ww.shape[0] + hb_sw.shape[0])
                        self.hsa_dict[site_i][23].append(acc_ww)
                        self.hsa_dict[site_i][24].append(don_ww)
                        self.hsa_dict[site_i][25].append(acc_sw)
                        self.hsa_dict[site_i][26].append(don_sw)
                        self.hsa_dict[site_i][27].extend(acc_sw_ids)
                        self.hsa_dict[site_i][28].extend(don_sw_ids)
                        if wat_nbrs.shape[0]!= 0 and hb_ww.shape[0] != 0:
                            self.hsa_dict[site_i][28].append(wat_nbrs.shape[0]/hb_ww.shape[0])
            pbar.update(frame_i)
        pbar.finish()
        # normalize site quantities
        print "Start Normalization"        
        sphere_volume = (4 / 3) * np.pi
        bulk_water_per_site = self.rho_bulk * sphere_volume * self.num_frames
        skip_normalization = ["index", "x", "y", "z", "nwat", "occupancy", "gO",
            "TSsw", "TSww", "TStot", "solute_acceptors", "solute_donors"]
        for site_i in xrange(self.hsa_data.shape[0]):
            n_wat = self.hsa_data[site_i, 4]
            if n_wat != 0:
                self.hsa_data[site_i, 5] = n_wat/(self.start_frame + self.num_frames)
                for quantity_i in xrange(len(self.data_titles)):
                    if self.data_titles[quantity_i] not in skip_normalization:
                        if self.data_titles[quantity_i] in ["Esw", "EswLJ", "EswElec", "Eww", "EwwLJ", "EwwElec", "Etot"]:                        
                            self.hsa_data[site_i, quantity_i] = (np.sum(self.hsa_dict[site_i][quantity_i])/n_wat)*0.5
                        elif self.data_titles[quantity_i] in ["Ewwnbr"]:
                            self.hsa_data[site_i, quantity_i] = (np.sum(self.hsa_dict[site_i][quantity_i])/len(self.hsa_dict[site_i][quantity_i]))*0.5
                        else:
                            self.hsa_data[site_i, quantity_i] = np.sum(self.hsa_dict[site_i][quantity_i])/n_wat
                    if self.data_titles[quantity_i] in ["solute_acceptors", "solute_donors"]:
                        self.hsa_dict[site_i][quantity_i] = np.unique(self.hsa_dict[site_i][quantity_i])
        print "End Normalization"        

    def calculate_angular_structure(self, dist_cutoff=6.0, site_indices=[]):
        '''
        Returns energetic quantities for each hydration site
        '''
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=self.start_frame + self.num_frames).start()
        angular_strucure_data = [[] for site in self.hsa_data]
        for i in xrange(self.start_frame, self.start_frame + self.num_frames):
            # print "Processing frame: ", i+1, "..."
            frame = md.load_frame(self.trajectory, i, top=self.topology)
            pos = frame.xyz[0, :, :] * 10.0
            pbc = frame.unitcell_lengths[0] * 10.0
            #pos = self.trj[i].xyz[0,:,:]*10.0
            # obtain coords of O-atoms
            oxygen_pos = pos[self.wat_oxygen_atom_ids]

            cluster_search_space = NeighborSearch(oxygen_pos, 1.0)
            water_search_space = NeighborSearch(oxygen_pos, dist_cutoff)

            for site_i, site_data in enumerate(self.hsa_data):
                #cluster_center_coords = (site_data["x"], site_data["y"], site_data["z"])
                cluster_center_coords = (site_data[1][0], site_data[1][1], site_data[1][2])
                nbr_indices = cluster_search_space.query_nbrs_single_point(cluster_center_coords)
                cluster_wat_oxygens = [self.wat_oxygen_atom_ids[nbr_index] for nbr_index in nbr_indices]
                # begin iterating over water oxygens found in this cluster in
                # current frame
                for wat_O in cluster_wat_oxygens:
                    cluster_water_all_atoms = self.topology.select("resid " + str(self.topology.atom(wat_O).residue.index))
                    nbr_indices = water_search_space.query_point_and_distance(pos[wat_O])
                    firstshell_wat_oxygens = [self.wat_oxygen_atom_ids[nbr_index[0]] for nbr_index in nbr_indices]
                    nbr_dist = [nbr_index[1] for nbr_index in nbr_indices]
                    #nbr_angle = []
                    if len(firstshell_wat_oxygens) != 0:
                        nbr_water_residues = [self.topology.select("resid " + str(self.topology.atom(
                            nbr_wat_O).residue.index)) for nbr_wat_O in firstshell_wat_oxygens]
                        for idx, water_res in enumerate(nbr_water_residues):
                            hbangle_combinations = np.asarray([[wat_O, water_res[0], water_res[1]], [wat_O, water_res[0], 
                                water_res[2]], [water_res[0], wat_O, wat_O + 1], [water_res[0], wat_O, wat_O + 2]])
                            theta_list = md.compute_angles(frame, hbangle_combinations) * DEG_PER_RAD
                            #nbr_angle.append(np.min(theta_list))  # min angle is a potential Hbond
                            angular_strucure_data[site_i].append([nbr_dist[idx], np.min(theta_list)])
            pbar.update(i + 1)
        pbar.finish()

        cwd = os.getcwd()
        # create directory to store detailed data for individual columns in HSA
        directory = cwd + "/" + self.prefix + "_angular_data"
        if not os.path.exists(directory):
            os.makedirs(directory)
        os.chdir(directory)

        for site_i, data in enumerate(angular_strucure_data):
            with open("%03d_r_theta.txt" % site_i, "w") as f:
                for d in data:
                    line = "{0[0]:.3f} {0[1]:.3f}\n".format(d)
                    f.write(line)


    def calculate_lonranged_ww_energy(self, site_indices=[], shell_radii=[3.5, 5.5, 8.5]):

        if len(site_indices) == 0:
            site_indices = [int(i) for i in self.hsa_data[:, 0]]
        else:
            for index in site_indices:
                if index > self.hsa_data[:, 0][-1]:
                    sys.exit("Site %d does not exits, please provide valid site indices." % index)

        n_sites = len(site_indices)
        n_columns = 2*len(shell_radii)
        shells = [(0.0, shell_radii[0])]
        for i in xrange(1, len(shell_radii)):
            shells.append((shell_radii[i - 1], shell_radii[i]))
        shells.append((shell_radii[-1], 100.0))
        longranged_data = np.zeros((n_sites, 2*len(shells)))
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=self.start_frame + self.num_frames).start()
        for frame_i in xrange(self.start_frame, self.start_frame + self.num_frames):
            frame = md.load_frame(self.trajectory, frame_i, top=self.topology)
            pos = md.utils.in_units_of(frame.xyz, "nanometers", "angstroms")
            pbc = md.utils.in_units_of(frame.unitcell_lengths, "nanometers", "angstroms")
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
                    cluster_center_coords = (self.hsa_data[site_i, 1], self.hsa_data[site_i, 2], self.hsa_data[site_i, 3])
                    nbr_indices = cluster_search_space.query_nbrs_single_point(cluster_center_coords)
                    cluster_wat_oxygens = [self.wat_oxygen_atom_ids[nbr_index] for nbr_index in nbr_indices]
                    if len(cluster_wat_oxygens) != 0:
                        wat_O = cluster_wat_oxygens[0]
                        self.hsa_data[site_i, 4] += 1
                if wat_O is not None:
                    distance_matrix = np.zeros((self.water_sites, self.all_atom_ids.shape[0]), np.float_)
                    calc.get_pairwise_distances(np.asarray([site_i, wat_O]), self.all_atom_ids, pos, pbc, distance_matrix)
                    energy_lj, energy_elec = self.calculate_energy(distance_matrix)
                    for shell_index, shell in enumerate(shells):
                        wat_nbrs = self.wat_oxygen_atom_ids[np.where((distance_matrix[0, :][self.wat_oxygen_atom_ids] <= shell[1]) & (distance_matrix[0, :][self.wat_oxygen_atom_ids] > shell[0]))]
                        longranged_data[index, 2*shell_index] += wat_nbrs.shape[0]
                        for i in xrange(self.water_sites):
                            longranged_data[index, (2*shell_index) + 1] += np.sum(energy_elec[:, wat_nbrs + i])
            pbar.update(frame_i)
        pbar.finish()
        # normalize site quantities
        for index, site_i in enumerate(site_indices):
            n_wat = self.hsa_data[site_i, 4]
            if n_wat != 0:
                longranged_data[index, :] /= n_wat*2.0
        
        # write data
        with open(self.prefix + "_longrange_Eww_summary.txt", "w") as f:
            header = "index "
            formatted_output = "{0:.0f} "
            for shell_index, shell in enumerate(shells):
                if shell_index + 1 == len(shells):
                    header += "Nnbr_>" + str(shell_index) + " E_shell_>" + str(shell_index) + "\n"
                    formatted_output += "{1[%d]:.6f} {1[%d]:.6f}\n" % (2*shell_index, (2*shell_index) + 1)
                else:
                    header +=  "Nnbr_" + str(shell_index + 1) + " E_shell_" + str(shell_index + 1) + " "
                    formatted_output += "{1[%d]:.6f} {1[%d]:.6f} " % (2*shell_index, (2*shell_index) + 1)
            f.write(header)
            for index, site_i in enumerate(site_indices):
                n_wat = self.hsa_data[site_i, 4]
                site_data_line = formatted_output.format(self.hsa_data[site_i, 0], longranged_data[index, :])
                f.write(site_data_line)

    def print_system_summary(self):
        print "System information:"
        print "\tParameter file: %s\n" % self.topology_file
        print "\tTrajectory: %s\n" % self.trajectory
        print "\tNumber of clusters: %d\n" % len(self.hsa_data)
        print "\tPeriodic Box: %s\n" % self.box_type
        print "\tFrames: %d, Total Atoms: %d, Waters: %d, Solute Atoms: %d\n" \
                % (self.num_frames, self.all_atom_ids.shape[0], self.wat_oxygen_atom_ids.shape[0], self.non_water_atom_ids.shape[0])
        print "\tH-bond Donors: %d, H-bond Acceptors: %d, H-bond Donors & Acceptors: %d\n" \
                % (self.solute_don_ids.shape[0], self.solute_acc_ids.shape[0], self.solute_acc_don_ids.shape[0])
    @function_timer
    def write_summary(self):
        with open(self.prefix + "_hsa_summary.txt", "w") as f:
            header = " ".join(self.data_titles) + "\n"
            f.write(header)
            
            # format first six columns
            formatted_output = "{0[0]:.0f} {0[1]:.2f} {0[2]:.2f} {0[3]:.2f} {0[4]:.0f} {0[5]:.2f} "
            # format site energetic, entropic and structural data
            for quantity_i in xrange(6, len(self.data_titles) - 2):
                formatted_output += "{0[%d]:.6f} " % quantity_i
            # format solute acceptors and donors
            formatted_output += "{1} {2}\n"
            for site_i in xrange(self.hsa_data.shape[0]):
                solute_acceptors = [str(self.topology.atom(acceptor)) for acceptor in self.hsa_dict[site_i][27]]
                solute_donors = [str(self.topology.atom(donor)) for donor in self.hsa_dict[site_i][28]]
                site_data_line = formatted_output.format(self.hsa_data[site_i, :],
                                ",".join(solute_acceptors), 
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
        for site_i in xrange(self.hsa_data.shape[0]):
            site_index = "%03d_" % site_i
            for quantity_i in xrange(len(self.data_titles)):
                if self.data_titles[quantity_i] not in skip_write_data and len(
                        self.hsa_dict[site_i][quantity_i]) != 0:
                    data_file_name = site_index + self.prefix + "_" + self.data_titles[quantity_i] + ".txt"
                    with open(data_file_name, "w") as data_file:
                        data_file.writelines("%s\n" % item for item in self.hsa_dict[site_i][quantity_i])
        os.chdir("../")
        
def plot_enbr_distribution(site_data, data_dir, sites=[], nbr_norm = False):
    """
    generate an Enbr plot for an arbitrary list of sites. First site should be the reference system.
    sites: a list of keys which represent site labels
    data: a dictionary of sites
    x_values: data points on x-axis
    nbr_norm: Normalize by number of neighbors
    outname: name of output file
    """
    colors = ["green", "blue", "red", "orange"]
    for k in data.keys():
        "generating Enbr plot for: ", k
        y_values = data[k][0]*data[k][1]
        #for val in range(len(x_values)):
        #    if x_values[val] <= -3.26:
        #         y_values[val] = 0.0
        fig, ax = plt.subplots(1)
        fig.set_size_inches(3, 3)
        # set some basic settings for fonts and rendering
        # set x and y axes limits, hard coded for now
        plt.xlim(-4.0, 2.5)
        plt.ylim(0.0, 1.0)
        # set x and y axes label titles and font sizes
        x_label = r'$\mathit{E_{n} (kcal/mol)}$'
        y_label = r'$\mathit{\rho(E_{n})}$'
        # set some features for the tick marks and labels
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(start, end, 2.0))
        #start, end = ax.get_ylim()
        #ax.yaxis.set_ticks(np.arange(start, end, 0.1))
        if nbr_norm:
            y_label = r'$\mathit{\rho(E_{n})/N^{nbr}}$'
            plt.ylim(0.0, 0.5)
            start, end = ax.get_ylim()
            ax.yaxis.set_ticks(np.arange(start, end, 0.1))
        ax.set_xlabel(x_label, size=14)
        ax.set_ylabel(y_label, size=14)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        plt.minorticks_on()
        plt.tick_params(which='major', width=1, length=4, direction='in')
        plt.tick_params(which='minor', width=1, length=2, direction='in')
        plt.tick_params(axis='x', labelsize=12)
        plt.tick_params(axis='y', labelsize=12)
        ax.yaxis.tick_left()
        ax.xaxis.tick_bottom()
        # set legend locations
        # save figure with optimal settings
        plt.tight_layout()
        #plt.plot(x_1, y_1, antialiased=True, linewidth=1.0, color="red", label="Methane")
        plt.plot(x_values, y_values, antialiased=True, linewidth=1.0, color="red", label=k)

        plt.legend(loc='upper right', prop={'size':10}, frameon=False)
        plt.savefig(save_dir + "/" + k + "_Enbr_plot.png", dpi=300)
        plt.close()

    """
    set sites to all sites if sites==[] (take all sites from site_data)
    for each site in sites
      get enbr_data into an array
      construct kernel density estimate
      get a range of x_values
      evaluate density at x_values
      set plot formatting
      plot data
      save plot
    """

#site_data, data_dir, sites=[], 
def plot_rtheta_distribution(rtheta_file, site_data, legend_label, nwat):
    
    integ_counts = 16.3624445886
    #integ_counts = 22560
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    theta_data = np.loadtxt(theta_file)
    r_data = np.loadtxt(r_file)
    Nnbr = len(r_data)/nwat
    print nwat, Nnbr
    # generate index matrices
    X, Y = np.mgrid[0:130:131j, 2.0:5.0:31j]
    # generate kernel density estimates
    values = np.vstack([theta_data, r_data])
    kernel =  stats.gaussian_kde(values)
    positions = np.vstack([X.ravel(), Y.ravel()])
    Z = np.reshape(kernel(positions).T, X.shape)
    Z *= integ_counts*0.1
    #Z /= integ_counts
    sum_counts_kernel = 0
    #print kernel.n
    # correct Z
    for i in xrange(0, Y.shape[1]):
        d = Y[0,i]
        # get shell_vol
        d_low = d - 0.1
        vol = (4.0/3.0)*np.pi*(d**3)
        vol_low = (4.0/3.0)*np.pi*(d_low**3)
        shell_vol =  vol - vol_low

        counts_bulk = 0.0329*shell_vol
        sum_counts_kernel += np.sum(Z[:,i])
        #Z[:,i] /= counts_bulk
        Z[:,i] = Z[:,i]/counts_bulk

    print sum_counts_kernel
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0.5, antialiased=True, alpha=1.0, cmap=cm.coolwarm, label=legend_label)
    x_label = r"$\theta^\circ$"
    y_label = r"$r (\AA)$"
    ax.set_xlabel(x_label)
    ax.set_xlim(0, 130)
    ax.set_ylabel(y_label)
    ax.set_ylim(2.0, 5.0)
    z_label = r'$\mathrm{P(\theta, \AA)}$'
    ax.set_zlabel(z_label)
    ax.legend(legend_label, loc='upper left', prop={'size':6})
    ax.set_zlim(0.0, 0.15)
    plt.savefig(legend_label + ".png", dpi=300)
    plt.close()

def read_hsa_summary(hsa_data_file):
    '''
    Returns a dictionary with hydration site index as keys and a list of various attributes as values.
    Parameters
    ----------
    hsa_data_file : string
        Text file containing 
    Returns
    -------

    '''
    f = open(hsa_data_file, 'r')
    data = f.readlines()
    hsa_header = data[0]
    data_keys = hsa_header.strip("\n").split()
    hsa_data = {}
    for l in data[1:]:
        float_converted_data = [float(x) for x in l.strip("\n").split()[1:27]]
        hsa_data[int(l.strip("\n").split()[0])] = float_converted_data
    f.close()
    return hsa_data
def main():

    test_1 = False

    if test_1:
        hsa = SiteWaterAnalysis("/Users/kamranhaider/arg_gist_calcs/mdtraj/arg.prmtop", "/Users/kamranhaider/arg_gist_calcs/mdtraj/md10ns.ncdf", start_frame=0, num_frames=1000, ligand_file="/Users/kamranhaider/arg_gist_calcs/mdtraj/ligand_arg.pdb", prefix="test_1")
        hsa.print_system_summary()
        hsa.calculate_site_quantities()
        hsa.write_summary()
        hsa.write_data()
    else:
        hsa = SiteWaterAnalysis("/Users/kamranhaider/arg_gist_calcs/mdtraj/arg.prmtop", "/Users/kamranhaider/arg_gist_calcs/mdtraj/md10ns.ncdf", start_frame=0, num_frames=1000, cluster_center_file="/Users/kamranhaider/arg_gist_calcs/03_organizewaters/clustercenterfile.pdb", prefix="test_2")
        hsa.print_system_summary()
        #hsa.calculate_lonranged_ww_energy([1, 3])
        hsa.calculate_site_quantities()
        hsa.write_summary()
        hsa.write_data()


def entry_point():
    main()

if __name__ == '__main__':
    entry_point()
