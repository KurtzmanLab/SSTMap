
import sys
import time
import datetime
import os

import numpy as np
from scipy import stats
import mdtraj as md
import parmed as pmd
from progressbar import Bar, Percentage, ProgressBar, ETA

from water_analysis import WaterAnalysis, NeighborSearch
import _sstmap_ext_old as calc_old

DEG_PER_RAD = 180.0 / np.pi
RHO_BULK = 0.0329  # molecules/A^3 # 0.0329
SPHERE_VOL = (4 / 3) * np.pi


class SiteWaterAnalysis(WaterAnalysis):
    """
    """
    def __init__(self, topology_file, trajectory, start_frame=0, num_frames=0,
                 ligand_file=None, cluster_center_file=None, 
                 desmond_helper_file=None, prefix="test"):
        
        super(SiteWaterAnalysis, self).__init__(topology_file, trajectory,
            start_frame, num_frames, desmond_helper_file)        
        self.prefix = prefix

        if cluster_center_file == None and ligand_file == None:
            sys.exit("Please provide either a ligand file for clustering or\
                        a cluster center file to generate hydration sites.")

        if ligand_file != None:
            cluster_coords = self.generate_clusters()
            self.hsa_data = self.initialize_site_data(cluster_coords)
        if cluster_center_file != None:
            clusters_pdb_file = md.load_pdb(cluster_center_file)
            cluster_coords = clusters_pdb_file.xyz * 10.0
            self.hsa_data = self.initialize_site_data(cluster_coords[0, :, :])

    def initialize_site_data(cluster_coords):
        self.data_titles = ["x", "y", "z",
                            "nwat", "occupancy", "gO",
                            "Esw", "EswLJ", "EswElec", "Eww", "EwwLJ", "EwwElec", "Etot", "Ewwnbr",
                            "TSsw", "TSww", "TStot",
                            "Nnbrs", "Nhbww", "Nhbsw", "Nhbtot",
                            "PercentHBww", "f_enc",
                            "Acc_ww", "Don_ww", "Acc_sw", "Don_sw",
                            "solute_acceptors", "solute_donors"]
        n_sites = cluster_coords.shape[0]
        hydration_sites = [{} for i in range(n_sites)]
        for site_i in xrange(n_sites):
            hydration_sites[site_i] = [{quantity: [] for quantity in self.data_titles[6:]}]
            hydration_sites[site_i].append(np.zeros(len(self.data_titles), dtype="float64"))
            hydration_sites[site_i][1][0] = cluster_coords[0, site_i, 0]
            hydration_sites[site_i][1][1] = cluster_coords[0, site_i, 1]
            hydration_sites[site_i][1][2] = cluster_coords[0, site_i, 2]
    return hydration_sites
            
    def get_clusters_from_file(self):
        '''
        Returns a dictionary with hydration site indices as keys and their properties as values.
        '''
        self.data_titles = ["x", "y", "z",
                            "nwat", "occupancy", "gO",
                            "Esw", "EswLJ", "EswElec", "Eww", "EwwLJ", "EwwElec", "Etot", "Ewwnbr",
                            "TSsw", "TSww", "TStot",
                            "Nnbrs", "Nhbww", "Nhbsw", "Nhbtot",
                            "PercentHBww", "f_enc",
                            "Acc_ww", "Don_ww", "Acc_sw", "Don_sw",
                            "solute_acceptors", "solute_donors"]
            n_sites = clusters_pdb_file.n_atoms
            # data for each sites consist of a dictionary
            # within which we have a list dictionaries for each quantity
            # then we have an array to store averages for each quantity
            hydration_sites = [{} for i in range(n_sites)]
            for site_i in xrange(n_sites):
                hydration_sites[site_i] = [{quantity: [] for quantity in self.data_titles[6:]}]
                hydration_sites[site_i].append(np.zeros(len(self.data_titles), dtype="float64"))
                hydration_sites[site_i][1][0] = cluster_coords[0, site_i, 0]
                hydration_sites[site_i][1][1] = cluster_coords[0, site_i, 1]
                hydration_sites[site_i][1][2] = cluster_coords[0, site_i, 2]
        return hydration_sites

    def generate_clusters(self):
        # Obtain binding site solute atoms
        ligand = md.load_pdb(self.ligand_file)
        ligand_coords = ligand.xyz[0, :, :] * 10.0
        first_frame = md.load_frame(self.trajectory, 0, top=self.topology_file)
        solute_pos = first_frame.xyz[0, self.non_water_atom_ids, :] * 10.0
        search_space = NeighborSearch(solute_pos, 5.0)
        near_indices = search_space.query_nbrs_multiple_points(ligand_coords)
        solute_nbrs = [self.non_water_atom_ids[nbr_index]
                       for nbr_index in near_indices]
        ###DEBUGG###
        ligand_arg_atom_ids = self.topology.select("resname ARG")
        ###DEBUGG###

        # Obtain waters solvating the binding site
        print "Reading in trajectory for clustering..."
        trj = md.load(self.trajectory, top=self.topology)[0:1000]
        print "Done."
        initial_time = time.time()
        print "Obtaining a superconfiguration of all waters found in the binding site throught the trajectory, can take long time for big systems!"
        #nbrs = md.compute_neighbors(trj, 0.50, solute_nbrs, haystack_indices=self.wat_oxygen_atom_ids)
        ###DEBUGG###
        nbrs = md.compute_neighbors(trj, 0.50, ligand_arg_atom_ids, haystack_indices=self.wat_oxygen_atom_ids)
        ###DEBUGG###
        all_nbrs = [(i, nbr) for i in range(len(nbrs)) for nbr in nbrs[i]]
        # all_nbrs is a list of all waters with their frame ids
        
        print "Total number of waters to porcess: ", len(all_nbrs)
        elapsed_time = time.time() - initial_time
        print "Time took for creating superconfiguration: ", str(datetime.timedelta(seconds=elapsed_time))
        # Perform Clustering
        print "Performing clustering on superconfiguration..."
        cutoff = self.num_frames * 2 * 0.1401
        n_wat = 3 * cutoff
        cluster_list = []
        cluster_iter = 0
        # Build KDTree and get initial nbr count for all waters
        all_wat_coords = np.ma.array([trj.xyz[wat[0], wat[1], :] for wat in all_nbrs], mask=False)*10.0
        tree = spatial.cKDTree(all_wat_coords)
        nbr_list = tree.query_ball_point(all_wat_coords, 1.0)
        #dist_matrix = np.zeros((all_wat_coords.shape[0], all_wat_coords.shape[0]), np.float_)
        #calc.get_dist_matrix(all_wat_coords.shape[0], dist_matrix, all_wat_coords)
        #print sys.getsizeof(dist_matrix), dist_matrix.nbytes
        #nbr_list = []
        #for point in all_wat_coords:
        #    nbr_list.append(tree.query_ball_point(point, 0.10))


        #print sys.getsizeof(nbr_list)
        nbr_count_list = np.ma.array([len(nbrs)
                                      for nbrs in nbr_list], mask=False)

        # Clustering loop
        initial_time = time.time()
        while n_wat >= cutoff:
            cluster_iter += 1
            # get water with max nbrs, reterive its nbrs
            print "\tCluster iteration: ", cluster_iter
            print "\tretrieving maximum..."
            max_index = np.argmax(nbr_count_list)
            to_exclude = np.array(nbr_list[max_index])
            print "\tindex of maximum: ", max_index, "population: ", len(to_exclude)
            n_wat = len(to_exclude) + 1
            cluster_list.append(all_nbrs[max_index])
            # Mask current water and its nbrs so that they are not considered
            # in the next iteration
            nbr_count_list.mask[to_exclude] = True
            nbr_count_list.mask[max_index] = True
            updated_wat_coords = np.delete(all_wat_coords, to_exclude, 0)
            #updated_wat_coords = np.delete(updated_wat_coords, max_index, 0)

            # Update nbr_count_list
            # Recalculate neighbors only for the neighbors of exlucded waters
            # First obtain a list of nbrs, whose nbr counts are supposed to be
            # updated
            nbrs_of_to_exclude = np.unique(np.array([n_excluded for excluded_nbrs in nbr_list[
                                           to_exclude] for n_excluded in excluded_nbrs]))
            to_update = np.setdiff1d(nbrs_of_to_exclude, to_exclude)
            to_update = np.setdiff1d(to_update, np.asarray(max_index))
            # print "Current nbr counts are: "
            # print nbr_count_list[to_update]
            tree = spatial.cKDTree(updated_wat_coords)
            updated_nbr_list = tree.query_ball_point(
                all_wat_coords[to_update], 0.10)
            updated_nbr_count_list = np.array(
                [len(nbrs) for nbrs in updated_nbr_list])
            for index, nbrs in enumerate(updated_nbr_list):
                if not nbr_count_list.mask[to_update[index]]:
                   # print "Updating nbrs of index: ", to_update[index], "
                   # from: ", nbr_count_list[to_update[index]], " to: ",
                   # len(nbrs)
                    nbr_count_list[to_update[index]] = len(nbrs)
            # print "Updated nbr counts are: "
            # print updated_nbr_count_list

        elapsed_time = time.time() - initial_time
        print "Time took for Clustering: ", str(datetime.timedelta(seconds=elapsed_time))

        # print cluster_list

        f = open("test_clustercenterfile.pdb", "w")
        header = "REMARK Initial number of clusters: {0}\n".format(
            len(cluster_list))
        f.write(header)
        for cluster in cluster_list:
            cluster_coords = trj.xyz[cluster[0], cluster[1], :]
            # print "Cluster coords: ", cluster_coords
            pdb_line = "ATOM      1  O   WAT A   1      {0[0]:2.3f}  {0[1]:2.3f}  {0[2]:2.3f}  0.00  0.00\n".format(
                cluster_coords * 10.0)
            f.write(pdb_line)
        f.close()
        print "Number of clusters: ", len(cluster_list)

    def calculate_site_quantities(self, energy=True, hbonds=True):
        '''
        Returns energetic quantities for each hydration site
        '''
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=self.start_frame + self.num_frames).start()
        for i in xrange(self.start_frame, self.start_frame + self.num_frames):
            # print "Processing frame: ", i+1, "..."
            frame = md.load_frame(self.trajectory, i, top=self.topology)
            pos = frame.xyz[0, :, :] * 10.0
            pbc = frame.unitcell_lengths[0] * 10.0
            #pos = self.trj[i].xyz[0,:,:]*10.0
            # obtain coords of O-atoms
            oxygen_pos = pos[self.wat_oxygen_atom_ids]

            cluster_search_space = NeighborSearch(oxygen_pos, 1.0)
            water_search_space = NeighborSearch(oxygen_pos, 3.5)

            solute_pos = pos[self.non_water_atom_ids]
            solute_search_space = NeighborSearch(solute_pos, 3.5)

            for site_i, site_data in enumerate(self.hsa_data):
                #cluster_center_coords = (site_data["x"], site_data["y"], site_data["z"])
                cluster_center_coords = (
                    site_data[1][0], site_data[1][1], site_data[1][2])
                nbr_indices = cluster_search_space.query_nbrs_single_point(
                    cluster_center_coords)
                cluster_wat_oxygens = [self.wat_oxygen_atom_ids[
                    nbr_index] for nbr_index in nbr_indices]
                # begin iterating over water oxygens found in this cluster in
                # current frame
                for wat_O in cluster_wat_oxygens:
                    cluster_water_all_atoms = self.topology.select(
                        "resid " + str(self.topology.atom(wat_O).residue.index))
                    # at indices for rest of solvent water atoms
                    rest_wat_at_ids = np.setxor1d(
                        cluster_water_all_atoms, self.wat_atom_ids)
                    # at indices for rest of solvent water O-atoms
                    rest_wat_oxygen_at_ids = np.setxor1d(
                        wat_O, self.wat_oxygen_atom_ids)
                    # Solute-Water energy and h-bond calculations
                    Elec_sw = 0.0
                    LJ_sw = 0.0
                    hbsw = 0
                    acc_sw = 0
                    don_sw = 0
                    hbww = 0
                    acc_ww = 0
                    don_ww = 0
                    if len(self.non_water_atom_ids) != 0:
                        Elec_sw = calc_old.elecE( cluster_water_all_atoms, self.non_water_atom_ids,
                            pos, self.chg, pbc) * 0.5
                        LJ_sw = calc_old.vdwE(np.asarray([wat_O]), self.non_water_atom_ids,
                            pos, self.vdw, pbc) * 0.5

                        combined_don_accdon_ids = np.concatenate(
                            (self.solute_don_ids, self.solute_acc_don_ids))
                        combined_acc_accdon_ids = np.concatenate(
                            (self.solute_acc_ids, self.solute_acc_don_ids))

                        solute_nbr_indices = solute_search_space.query_nbrs_single_point(pos[
                                                                                         wat_O])
                        solute_nbr_atom_ids = [self.non_water_atom_ids[
                            nbr_index] for nbr_index in solute_nbr_indices]
                        nbr_solute_donors = np.intersect1d(
                            solute_nbr_atom_ids, combined_don_accdon_ids)
                        for donor in nbr_solute_donors:
                            for don_H_pair in self.don_H_pair_dict[donor]:
                                hbangle = md.compute_angles(frame, np.asarray(
                                    [[wat_O, donor, don_H_pair[1]]])) * DEG_PER_RAD
                                if hbangle <= 30.0:
                                    hbsw += 1
                                    acc_sw += 1
                                    if donor not in site_data[0]["solute_donors"]:
                                        site_data[0]["solute_donors"].append(donor)

                        nbr_solute_acceptors = np.intersect1d(
                            solute_nbr_atom_ids, combined_acc_accdon_ids)
                        for acceptor in nbr_solute_acceptors:
                            hbangle_combinations = np.asarray(
                                [[acceptor, wat_O, wat_O + 1], [acceptor, wat_O, wat_O + 2]])
                            theta_list = md.compute_angles(
                                frame, hbangle_combinations) * DEG_PER_RAD
                            hbangle, hbangle_index = np.min(theta_list), np.argmin(
                                theta_list)  # min angle is a potential Hbond
                            if hbangle <= 30:
                                hbsw += 1
                                don_sw += 1
                                if acceptor not in site_data[0]["solute_acceptors"]:
                                    site_data[0]["solute_acceptors"].append(acceptor)
                        site_data[0]["Esw"].append(Elec_sw + LJ_sw)
                        site_data[1][self.data_titles.index("Esw")] += Elec_sw + LJ_sw
                        site_data[0]["EswLJ"].append(LJ_sw)
                        site_data[1][self.data_titles.index("EswLJ")] += LJ_sw
                        site_data[0]["EswElec"].append(Elec_sw)
                        site_data[1][self.data_titles.index("EswElec")] += Elec_sw
                        site_data[0]["Nhbsw"].append(hbsw)
                        site_data[1][self.data_titles.index("Nhbsw")] += hbsw
                        site_data[0]["Acc_sw"].append(acc_sw)
                        site_data[1][self.data_titles.index("Acc_sw")] += acc_sw
                        site_data[0]["Don_sw"].append(don_sw)
                        site_data[1][self.data_titles.index("Don_sw")] += don_sw
                    # Water-Water energy and h-bond calculations
                    Elec_ww = calc_old.elecE(cluster_water_all_atoms, rest_wat_at_ids,
                        pos, self.chg, pbc) * 0.5
                    LJ_ww = calc_old.vdwE(np.asarray([wat_O]), rest_wat_oxygen_at_ids,
                        pos, self.vdw, pbc) * 0.5
                    nbr_indices = water_search_space.query_nbrs_single_point(pos[
                                                                             wat_O])
                    firstshell_wat_oxygens = [
                        self.wat_oxygen_atom_ids[nbr_index] for nbr_index in nbr_indices]
                    Enbr = 0.0
                    if len(firstshell_wat_oxygens) != 0:
                        nbr_water_residues = [self.topology.select("resid " + str(self.topology.atom(
                            nbr_wat_O).residue.index)) for nbr_wat_O in firstshell_wat_oxygens]
                        for idx, water_res in enumerate(nbr_water_residues):
                            nbr_elec_ww = calc_old.elecE(
                                cluster_water_all_atoms, water_res, pos, self.chg, pbc) * 0.5
                            nbr_LJ_ww = calc_old.vdwE(np.asarray([wat_O]), np.asarray(
                                [firstshell_wat_oxygens[idx]]), pos, self.vdw, pbc) * 0.5
                            site_data[0]["Ewwnbr"].append(nbr_elec_ww + nbr_LJ_ww)
                            site_data[1][
                                self.data_titles.index("Ewwnbr")] += nbr_elec_ww + nbr_LJ_ww
                            hbangle_combinations = np.asarray([[wat_O, water_res[0], water_res[1]], [wat_O, water_res[
                                                              0], water_res[2]], [water_res[0], wat_O, wat_O + 1], [water_res[0], wat_O, wat_O + 2]])
                            theta_list = md.compute_angles(frame, hbangle_combinations) * DEG_PER_RAD
                            hbangle, hbangle_index = np.min(theta_list), np.argmin(
                                theta_list)  # min angle is a potential Hbond
                            
                            if hbangle <= 30:
                                hbww += 1
                            if hbangle_index in [0, 1]:
                                acc_ww += 1
                            else:
                                don_ww += 1
                        site_data[0]["Nnbrs"].append(
                            len(firstshell_wat_oxygens))
                        site_data[1][self.data_titles.index(
                            "Nnbrs")] += len(firstshell_wat_oxygens)
                        site_data[0]["Nhbww"].append(hbww)
                        site_data[1][self.data_titles.index("Nhbww")] += hbww
                        site_data[0]["Acc_ww"].append(acc_ww)
                        site_data[1][
                            self.data_titles.index("Acc_ww")] += acc_ww
                        site_data[0]["Don_ww"].append(don_ww)
                        site_data[1][
                            self.data_titles.index("Don_ww")] += don_ww
                        site_data[0]["f_enc"].append(
                            1.0 - (len(firstshell_wat_oxygens) / 5.25))
                        site_data[1][self.data_titles.index(
                            "f_enc")] += 1.0 - (len(firstshell_wat_oxygens) / 5.25)
                        site_data[0]["PercentHBww"].append(
                            (hbww / len(firstshell_wat_oxygens)) * 100.0)
                        site_data[1][self.data_titles.index(
                            "PercentHBww")] += (hbww / len(firstshell_wat_oxygens)) * 100.0
                    site_data[1][self.data_titles.index("nwat")] += 1
                    site_data[0]["Eww"].append(Elec_ww + LJ_ww)
                    site_data[1][
                        self.data_titles.index("Eww")] += Elec_ww + LJ_ww
                    site_data[0]["EwwLJ"].append(LJ_ww)
                    site_data[1][self.data_titles.index("EwwLJ")] += LJ_ww
                    site_data[0]["EwwElec"].append(Elec_ww)
                    site_data[1][self.data_titles.index("EwwElec")] += Elec_ww
                    site_data[0]["Etot"].append(
                        Elec_sw + LJ_sw + Elec_ww + LJ_ww)
                    site_data[1][self.data_titles.index(
                        "Etot")] += Elec_sw + LJ_sw + Elec_ww + LJ_ww
                    site_data[0]["Nhbtot"].append(hbww + hbsw)
                    site_data[1][
                        self.data_titles.index("Nhbtot")] += hbww + hbsw
            pbar.update(i + 1)
        # normalize site quantities
        pbar.finish()
        BULKWAT_PER_SITE = RHO_BULK * SPHERE_VOL * self.num_frames
        skip_normalization = ["x", "y", "z", "nwat", "occupancy", "gO",
            "TSsw", "TSww", "TStot", "solute_acceptors", "solute_donors"]
        for site_i, site_data in enumerate(self.hsa_data):
            if site_data[1][self.data_titles.index("nwat")] != 0:
                site_data[1][self.data_titles.index(
                    "occupancy")] += site_data[1][self.data_titles.index("nwat")] / self.num_frames
                site_data[1][self.data_titles.index(
                    "gO")] += site_data[1][self.data_titles.index("nwat")] / BULKWAT_PER_SITE
                for key in site_data[0].keys():
                    if key not in skip_normalization:
                        site_data[1][self.data_titles.index(
                            key)] /= len(site_data[0][key])
                site_data[0]["solute_acceptors"] = [str(self.topology.atom(
                    acceptor)) for acceptor in site_data[0]["solute_acceptors"]]
                site_data[0]["solute_donors"] = [str(self.topology.atom(
                    donor)) for donor in site_data[0]["solute_donors"]]
        

    def calculate_angular_structure(self, dist_cutoff=6.0):
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


    def calculate_lonrange_energy(self, debug=False):
        '''
        Returns energetic quantities for each hydration site
        '''
        pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=self.start_frame + self.num_frames).start()# set up data storage array, each site has a list of two lists, one for energies and one for neighbors
        # divided into four shells, first three solvation shells and beyond
        long_range_data = [[[0.0, 0.0, 0.0, 0.0], [0, 0, 0, 0]] for site in self.hsa_data]
        nwat = [0.0 for site in self.hsa_data]
        for i in xrange(self.start_frame, self.start_frame + self.num_frames):
            frame = md.load_frame(self.trajectory, i, top=self.topology)
            pos = frame.xyz[0, :, :] * 10.0
            pbc = frame.unitcell_lengths[0] * 10.0
            # obtain coords of O-atoms
            oxygen_pos = pos[self.wat_oxygen_atom_ids]
            cluster_search_space = NeighborSearch(oxygen_pos, 1.0)
            # list of outer radii for first three shells
            shell_outer_radii = [3.5, 5.5, 8.5]
            for site_i, site_data in enumerate(self.hsa_data):
                cluster_center_coords = (
                    site_data[1][0], site_data[1][1], site_data[1][2])
                nbr_indices = cluster_search_space.query_nbrs_single_point(
                    cluster_center_coords)
                cluster_wat_oxygens = [self.wat_oxygen_atom_ids[
                    nbr_index] for nbr_index in nbr_indices]
                # begin iterating over water oxygens found in this cluster in
                # current frame
                for wat_O in cluster_wat_oxygens:
                    nwat[site_i] += 1
                    cluster_water_all_atoms = self.topology.select(
                        "resid " + str(self.topology.atom(wat_O).residue.index))
                    # at indices for rest of solvent water O-atoms
                    rest_wat_oxygen_at_ids = np.setxor1d(
                        wat_O, self.wat_oxygen_atom_ids)
                    # at indices for rest of solvent water atoms
                    rest_wat_at_ids = np.setxor1d(
                        cluster_water_all_atoms, self.wat_atom_ids)
                    # lists to store cumulative and actual nbrs in four shells
                    shell_nbrs_cumulative = [[], [], [], []]
                    shell_nbrs = [[], [], [], []]
                    # get cumulative nbrs for first three solvation shells
                    for shell_index, radius in enumerate(shell_outer_radii):
                        nbr_search_space = NeighborSearch(oxygen_pos, radius)
                        nbr_indices = nbr_search_space.query_nbrs_single_point(pos[
                                                                               wat_O])
                        nbr_oxygens = [self.wat_oxygen_atom_ids[
                            nbr_index] for nbr_index in nbr_indices]
                        if len(nbr_oxygens) != 0:
                            shell_nbrs_cumulative[shell_index] = nbr_oxygens
                    # get cumulative nbrs for beyond third shell
                    shell_nbrs_cumulative[-1] = rest_wat_oxygen_at_ids
                    # set actual nbrs for first solvation shell
                    shell_nbrs[0] = shell_nbrs_cumulative[0]
                    # iterate over solvation shells, second, third and beyond
                    for shell_index in range(1, len(shell_nbrs)):
                        shell_nbrs[shell_index] = np.setxor1d(
                            np.asarray(
                                shell_nbrs_cumulative[shell_index]), np.asarray(
                                shell_nbrs_cumulative[
                                    shell_index - 1])).astype(int)
                    # make sure the total number of nbrs sum up to total number
                    # of waters in the system (except the current water being
                    # analyzed)
                    Nnbr_list = [len(nbrs) for nbrs in shell_nbrs]
                    assert(sum(Nnbr_list) == len(rest_wat_oxygen_at_ids))
                    # storage array for water-water energy for first three
                    # shell
                    shell_energies = [0.0 for shell in shell_outer_radii]
                    # store water atoms up to 8.5A, later used to get water atoms beyond third shell
                    # since we do not need to calculated 'nbrs' beyond third
                    # shell
                    cumulative_8_5_wat_res = []
                    for shell_index, shell in enumerate(shell_outer_radii):
                        # calculate water-water energy for first three shells
                        if len(shell_nbrs[shell_index]) != 0:
                            nbr_water_residues = [self.topology.select("resid " + str(self.topology.atom(
                                nbr_wat_O).residue.index)) for nbr_wat_O in shell_nbrs[shell_index]]
                            for index, water_res in enumerate(
                                    nbr_water_residues):
                                cumulative_8_5_wat_res.extend(water_res)
                                nbr_elec_ww = calc_old.elecE(
                                    cluster_water_all_atoms, water_res, pos, self.chg, pbc) * 0.5
                                nbr_LJ_ww = calc_old.vdwE(np.asarray([wat_O]), np.asarray(
                                    [shell_nbrs[shell_index][index]]), pos, self.vdw, pbc) * 0.5
                                shell_energies[
                                    shell_index] += nbr_elec_ww + nbr_LJ_ww

                    # calculate water-water energy beyond the third shell
                    beyond_wat_oxygens = shell_nbrs[-1]
                    # print len(beyond_wat_oxygens)
                    beyond_wat_all_atoms = np.setxor1d(
                        np.asarray(cumulative_8_5_wat_res), rest_wat_at_ids)
                    beyond_Elec_ww = calc_old.elecE(
                        cluster_water_all_atoms, beyond_wat_all_atoms, pos, self.chg, pbc) * 0.5
                    beyond_LJ_ww = calc_old.vdwE(
                        np.asarray(
                            [wat_O]),
                        beyond_wat_oxygens,
                        pos,
                        self.vdw,
                        pbc) * 0.5
                    beyond_Etot = beyond_Elec_ww + beyond_LJ_ww
                    shell_energies.append(beyond_Etot)
                    # in debug mode, calculate total water-water and check if
                    # its equal to the sum of four shells
                    if debug:
                        Elec_ww = calc_old.elecE(
                            cluster_water_all_atoms, rest_wat_at_ids, pos, self.chg, pbc) * 0.5
                        LJ_ww = calc_old.vdwE(
                            np.asarray(
                                [wat_O]),
                            rest_wat_oxygen_at_ids,
                            pos,
                            self.vdw,
                            pbc) * 0.5
                        Etot = Elec_ww + LJ_ww
                        beyond_Etot += sum(shell_energies)
                        np.testing.assert_almost_equal(
                            beyond_Etot, Etot, decimal=7, err_msg='', verbose=True)
                    # update shell energy data for this site
                    for index, e_shell in enumerate(shell_energies):
                        long_range_data[site_i][0][index] += e_shell
                        long_range_data[site_i][1][index] += Nnbr_list[index]
            pbar.update(i + 1)
        pbar.finish()
        f = open(self.prefix + "_longrange_Eww_summary.txt", "w")
        header = "index Nnbr_1 E_shell_1 Nnbr_2 E_shell_2 Nnbr_3 E_shell_3 Nnbr_>3 E_shell_>3\n"
        f.write(header)

        for site_i, data in enumerate(long_range_data):
            if nwat[site_i] != 0:
                site_values = str(site_i) + " "
                for i in range(4):
                    site_values += "%0.3f %0.3f " % (long_range_data[site_i][1][i] / nwat[
                                                     site_i], long_range_data[site_i][0][i] / nwat[site_i])
                site_values += "\n"
                f.write(site_values)
        f.close()

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


    def write_summary(self):
        # Format and write data
        f = open(self.prefix + "_hsa_summary.txt", "w")
        structural_quantities = self.data_titles[17:-2]
        header = "index " + " ".join(self.data_titles) + "\n"
        f.write(header)
        # format first six columns
        formatted_output = "{0:d} {1[0]:.2f} {1[1]:.2f} {1[2]:.2f} {1[3]:.0f} {1[4]:.2f} {1[5]:.2f} "
        # format site energetic, entropic and structural data
        for column in self.data_titles[6:-2]:
            if column in structural_quantities:
                formatted_output += "{1[%d]:.2f} " % self.data_titles.index(column)
            else:
                formatted_output += "{1[%d]:.6f} " % self.data_titles.index(column)
        # format solute acceptors and donors
        formatted_output += "{2:s} {3:s}\n"
        for site_i, site_data in enumerate(self.hsa_data):
            site_data_line = formatted_output.format(
                site_i, site_data[1], ",".join(
                    site_data[0]["solute_acceptors"]), ",".join(
                    site_data[0]["solute_donors"]))
            f.write(site_data_line)
        f.close()

    def write_data(self):
        skip_write_data = ["x", "y", "z", "nwat", "occupancy", "gO",
                            "TSsw", "TSww", "TStot", "solute_acceptors", "solute_donors"]
        cwd = os.getcwd()
        # create directory to store detailed data for individual columns in HSA
        directory = cwd + "/" + self.prefix + "_hsa_data"
        if not os.path.exists(directory):
            os.makedirs(directory)
        os.chdir(directory)
        # for each cluster, go through time series data
        for site_i, site_data in enumerate(self.hsa_data):
            site_index = "%03d_" % site_i
            for title in self.data_titles:
                if title not in skip_write_data and len(
                        site_data[0][title]) != 0:
                    data_file_name = site_index + self.prefix + "_" + title
                    data_file = open(data_file_name, "w")
                    data_file.writelines("%s\n" %
                                         item for item in site_data[0][title])
                    data_file.close()
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
    data = read_hsa_summary()

def entry_point():
    main()

if __name__ == '__main__':
    entry_point()
