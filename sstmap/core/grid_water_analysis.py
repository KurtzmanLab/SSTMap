"""

"""

# import statements
import sys
import time

import numpy as np
import mdtraj as md
import parmed as pmd
from progressbar import Bar, Percentage, ProgressBar, ETA

from water_analysis import WaterAnalysis, function_timer
import _sstmap_ext as calc

DEG_PER_RAD = 180.0 / np.pi
GASKCAL = 0.0019872041


class GridWaterAnalysis(WaterAnalysis):
    """
    
    """
    @function_timer
    def __init__(self, topology_file, trajectory, start_frame=0, num_frames=0,
                 ligand_file=None, desmond_helper_file=None, prefix="test",
                 grid_center=[0.0, 0.0, 0.0], grid_dimensions=[5.0, 5.0, 5.0],
                 grid_resolution=[0.5, 0.5, 0.5]):

        super(GridWaterAnalysis, self).__init__(topology_file, trajectory, start_frame, num_frames, desmond_helper_file)
        self.paramname = topology_file
        self.param_class = pmd.load_file(self.paramname)
        self.resolution = grid_resolution[0]
        self.prefix = prefix
        if ligand_file is not None:
            # TODO: change this as a utility function
            # get ligad center
            lig = md.load_pdb(ligand_file)
            com = np.zeros((lig.n_frames, 3))
            masses = np.ones(lig.n_atoms)
            masses /= masses.sum()
            com[0, :] = lig.xyz[0, :].astype('float64').T.dot(masses)
            grid_center = com[0, :] * 10.0
        self.voxel_vol = self.resolution**3.0
        self.rho_bulk = 0.0334
        # set 3D grid around the region of interest
        self.initialize_grid(grid_center, grid_resolution, grid_dimensions)
        # initialize data structures to store voxel data
        self.voxeldata, self.voxeldict = self.initialize_voxel_data()
        #print "Reading in trajectory ..."        
        #self.trj = md.load(self.trajectory, top=self.paramname)[self.start_frame: self.start_frame + self.num_frames]
        #print "Done!"


    def initialize_grid(self, center, resolution, dimensions):
        # set grid center, res and dimension
        #self.center = np.array(center,dtype=np.float_)
        #self.dims = np.array(dimensions)
        #self.spacing = np.array(resolution,dtype=np.float_)
        self.center = np.array(center, dtype=np.float_)
        self.dims = np.array(dimensions, dtype=np.int_)
        self.spacing = np.array(resolution, dtype=np.float_)
        self.gridmax = self.dims * self.spacing + 1.5
        # set origin
        o = self.center - (0.5 * self.dims * self.spacing)
        self.origin = np.around(o, decimals=3)
        # set grid size (in terms of total points alog each axis)
        length = np.array(self.dims / self.spacing, dtype=np.float_)
        self.grid_size = np.ceil(length / self.spacing + 1.0)
        self.grid_size = np.cast['uint32'](self.grid_size)
        # Finally allocate the space for the grid
        self.grid = np.zeros(self.dims, dtype=np.int_)

    def initialize_voxel_data(self):
        voxel_dict = []
        v_count = 0
        voxel_array = np.zeros((self.grid.size, 35), dtype="float64")
        # print voxel_dict_new.shape
        for index, value in np.ndenumerate(self.grid):
            # point = grid.pointForIndex(index) # get cartesian coords for the
            # grid point
            _index = np.array(index, dtype=np.int32)
            # point = self.spacing * _index + self._origin
            point = _index * self.spacing + self.origin + 0.5 * self.spacing
            voxel_array[v_count, 1] = point[0]
            voxel_array[v_count, 2] = point[1]
            voxel_array[v_count, 3] = point[2]
            voxel_array[v_count, 0] = v_count
            # print voxel_dict_new[v_count, 0], voxel_dict_new[v_count, 1],
            # voxel_dict_new[v_count, 2]
            # create a dictionary key-value pair with voxel index as key and
            # it's coords as
            voxel_dict.append([[],[]])
            # voxel_dict[v_count].append(np.zeros(14, dtype="float64"))
            v_count += 1
        return voxel_array, voxel_dict

    def calculate_euler_angles(self, water, coords):
        pi = np.pi
        twopi = 2 * np.pi
        # define the lab frame of reference
        xlab = np.asarray([1.0, 0.0, 0.0], dtype="float64")
        #ylab = np.asarray([0.0, 1.0, 0.0], dtype="float64")
        zlab = np.asarray([0.0, 0.0, 1.0], dtype="float64")
        # create array for water oxygen atom coords, and append to this voxel's
        # coord list
        voxel_id = water[0]
        owat = coords[water[1], :]
        # create array for water's hydrogen 1 and 2
        h1wat = coords[water[1] + 1, :] - owat
        h2wat = coords[water[1] + 2, :] - owat
        # print frame_index, wat_O, owat, h1wat, h2wat
        # define water molecule's frame
        # H1 is water's x-axis, should be normalized
        xwat = np.copy(h1wat)
        xwat *= 1 / (np.linalg.norm(h1wat))
        # z-axis is the cross-product of H1 and H2
        zwat = np.cross(xwat, h2wat)
        zwat *= 1 / (np.linalg.norm(zwat))
        # y-axis is just perpendicular to z- and x-axis
        ywat = np.cross(zwat, xwat)
        ywat *= 1 / (np.linalg.norm(ywat))
        # Now we proceed to Euler angle calculations between water and lab frame
        # we start with theta and we will use cosine formula for the dot
        # product`
        dp = np.dot(zlab, zwat)
        # first we get theta which is angle between z-axes of two frames
        theta = np.arccos(dp)
        phi = 0
        psi = 0
        # if theta is between 0 and pi
        if theta > 1E-5 and theta < pi - 1E-5:
            # define a new vector which is perpendicular to both z-axes
            node = np.cross(zlab, zwat)
            node *= 1 / (np.linalg.norm(node))
            # get angle phi which is the angle between node and xlab
            dp = np.dot(node, xlab)
            if dp <= -1.0:
                phi = pi
            elif dp >= 1.0:
                phi = pi
            else:
                phi = np.arccos(dp)
                # check if angle phi is between 0 and 2pi
            if phi > 0.0 and phi < twopi:
                # define new vector v which is perpendicular to xlab and node
                v = np.cross(xlab, node)
                dp = np.dot(v, zlab)
                if dp < 0:
                    phi = twopi - phi
                    # get angle psi
                dp = np.dot(xwat, node)
            if dp <= - 1.0:
                psi = pi
            elif dp > 1.0:
                psi = 0.0
            else:
                psi = np.arccos(dp)
            if psi > 0.0 and psi < twopi:
                v = np.cross(node, xwat)
                dp = np.dot(v, zwat)
                if dp < 0:
                    psi = twopi - psi
            if not theta <= pi and theta >= 0 and phi <= twopi and phi >= 0 and psi <= twopi and psi >= 0:
                print "Error: Euler angles don't fall into range!"
                # add angle information for this water
        #self._voxel_wat_euler_map[frame_index, wat_index, 0] = theta
        #self._voxel_wat_euler_map[frame_index, wat_index, 1] = phi
        #self._voxel_wat_euler_map[frame_index, wat_index, 2] = psi
        #print theta, phi, psi
        self.voxeldict[voxel_id][0].append([theta, phi, psi])
        self.voxeldict[voxel_id][1].append(owat)

    def calculate_energy(self, distance_matrix):

        wat_wat_dist = distance_matrix[0, :][self.wat_oxygen_atom_ids]
        wat_solute_dist = distance_matrix[0, :][self.non_water_atom_ids]
        with np.errstate(invalid='ignore', divide='ignore'):
            energy_sw_lj = (self.solute_water_acoeff*np.power(wat_solute_dist, -12)) + (self.solute_water_bcoeff*np.power(wat_solute_dist, -6))
            energy_ww_lj = (self.water_water_acoeff*np.power(wat_wat_dist, -12)) + (self.water_water_bcoeff*np.power(wat_wat_dist, -6))
            water_chg = self.chg[self.wat_atom_ids[0:self.water_sites]].reshape(self.water_sites, 1)
            energy_elec = water_chg*np.tile(self.chg[self.all_atom_ids], (3, 1))
            energy_elec *= np.power(distance_matrix, -1)

        #print "Solute-water LJ Energy of this water: ", np.nansum(energy_sw_lj)
        #print "Solute-water Elec Energy of this water: ", np.sum(energy_elec[:, self.non_water_atom_ids])
        #print "Water-water LJ Energy of this water: ", np.nansum(energy_ww_lj)/2.0
        #print "Water-water Elec Energy of this water: ", (np.sum(energy_elec[:, self.wat_atom_ids[0]:water_id])/2.0) + (np.sum(energy_elec[:, water_id + self.water_sites:])/2.0)
        energy_lj = np.concatenate((energy_sw_lj, energy_ww_lj), axis=0)
        return energy_lj, energy_elec


    def calculate_hydrogen_bonds(self, traj, water, water_nbrs, solute_nbrs):
        hbond_data = []
        angle_triplets = []
        for wat_nbr in water_nbrs:
            angle_triplets.extend([[water, wat_nbr, wat_nbr+1], [water, wat_nbr, wat_nbr+2], [wat_nbr, water, water+1], [wat_nbr, water, water+2]])
        for solute_nbr in solute_nbrs:
            if self.prot_hb_types[solute_nbr] == 1 or self.prot_hb_types[solute_nbr] == 3:
                angle_triplets.extend([[solute_nbr, water, water+1], [solute_nbr, water, water+2]])
            if self.prot_hb_types[solute_nbr] == 2 or self.prot_hb_types[solute_nbr] == 3:
                for don_H_pair in self.don_H_pair_dict[solute_nbr]:
                    angle_triplets.extend([[water, solute_nbr, don_H_pair[1]]])
        angle_triplets = np.asarray(angle_triplets)

        #FIXME: Why does using self.trj[frame] take too long?
        angles = md.compute_angles(traj, angle_triplets) * DEG_PER_RAD
        angles_ww = angles[0, 0:water_nbrs.shape[0]*4]
        angles_sw = angles[0, water_nbrs.shape[0]*4:]
        hbonds_ww = angle_triplets[np.where(angles_ww <= 30.0)]
        hbonds_sw = angle_triplets[np.where(angles_sw <= 30.0)]
        acc_ww = hbonds_ww[:, 0][np.where(hbonds_ww[:, 0] == water)].shape[0]
        don_ww = hbonds_ww.shape[0] - acc_ww
        acc_sw = hbonds_sw[:, 0][np.where(hbonds_sw[:, 0] == water)].shape[0]
        don_sw = hbonds_sw.shape[0] - acc_sw

        return (hbonds_ww.shape[0], hbonds_sw.shape[0], acc_ww, don_ww, acc_sw, don_sw)

    @function_timer
    def calculate_entropy(self):
        bulkwaterpervoxel = self.voxel_vol * self.rho_bulk * self.num_frames        
        for voxel, waters in enumerate(self.voxeldict):
            if len(waters[0]) > 1.0:
                self.voxeldata[voxel, 5]  = len(waters) / bulkwaterpervoxel
                #density-weighted trans entropy
                dTStr_dens = -GASKCAL * 300 * self.rho_bulk * self.voxeldata[voxel, 5] * np.log(self.voxeldata[voxel, 5])
                self.voxeldata[voxel, 7] = dTStr_dens
                self.voxeldata[voxel, 8] = self.voxeldata[voxel, 7] * self.num_frames * self.voxel_vol / (1.0 * self.voxeldata[voxel, 4])
                #print voxel, self.voxeldata[voxel, 7], self.voxeldata[voxel, 8]
                angle_array = np.asarray(waters[0])
                #TODO: unintuitive name for this variable 
                #density-weighted orinet entropy
                dTS_nn_or = calc.getNNOrEntropy(len(waters[0]), angle_array)
                # normalized orientational entropy
                self.voxeldata[voxel, 10] = GASKCAL * 300 * ((dTS_nn_or / self.voxeldata[voxel, 4]) + 0.5772)
                # density-weighted orientational entropy
                self.voxeldata[voxel, 9] = self.voxeldata[voxel, 10] * self.voxeldata[voxel, 4] / (self.num_frames * self.voxel_vol)
                #coord_array = np.asarray(waters[1])
                #dTS_nn_tr = calc.getNNTrEntropy(len(waters[0]), self.num_frames, coord_array)
                # normalized translational entropy
                #self.voxeldata[voxel, 8] = GASKCAL * 300 * ((dTS_nn_tr/self.voxeldata[voxel, 4]) + 0.5772156649)
                # density-weighted trnaslationa entropy
                #self.voxeldata[voxel, 7] = self.voxeldata[voxel, 8] * self.voxeldata[voxel, 4]/(self.num_frames * self.voxel_vol)

    @function_timer
    def process_grid(self, energy=True, hbonds=True, entropy=True):
        pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=self.num_frames).start()
        for frame in xrange(self.start_frame, self.start_frame + self.num_frames):
            trj = md.load_frame(self.trajectory, frame, top=self.topology_file)
            pos = trj.xyz * 10.0
            periodic_box = trj.unitcell_lengths * 10.0
            frame_data = [[] for i in range(trj.n_frames)]
            calc.assign_voxels(pos, self.dims, self.gridmax, self.origin, frame_data, self.wat_oxygen_atom_ids)
            waters = frame_data[0]
            for wat in waters:
                self.voxeldata[wat[0], 4] += 1
                
                if energy or hbonds:
                    distance_matrix = np.zeros((self.water_sites, self.all_atom_ids.shape[0]), np.float_)
                    calc.get_pairwise_distances(wat, self.all_atom_ids, pos, periodic_box, distance_matrix)
                    wat_nbrs = self.wat_oxygen_atom_ids[np.where((distance_matrix[0, :][self.wat_oxygen_atom_ids] <= 3.5) &(distance_matrix[0, :][self.wat_oxygen_atom_ids] > 0.0))]
                    prot_nbrs = self.non_water_atom_ids[np.where(distance_matrix[0, :][self.non_water_atom_ids] <= 3.5)]
                    self.voxeldata[wat[0], 18] += wat_nbrs.shape[0]
                    if energy:
                        energy_lj, energy_elec = self.calculate_energy(distance_matrix)
                        self.voxeldata[wat[0], 11] += np.nansum(energy_lj[self.wat_oxygen_atom_ids[0]:])
                        self.voxeldata[wat[0], 11] += np.sum(energy_elec[:, self.wat_atom_ids[0]:wat[1]]) + np.sum(energy_elec[:, wat[1] + self.water_sites:])
                        self.voxeldata[wat[0], 13] += np.sum(energy_lj[:self.wat_oxygen_atom_ids[0]:])
                        self.voxeldata[wat[0], 13] += np.sum(energy_elec[:, self.non_water_atom_ids])
                        #print "Water-water LJ Energy of this water: ", np.nansum(energy_lj[self.wat_oxygen_atom_ids[0]:])
                        #print "Water-water Elec Energy of this water: ", np.sum(energy_elec[:, self.wat_atom_ids[0]:wat[1]]) + np.sum(energy_elec[:, wat[1] + self.water_sites:])
                        #print "Solute-water LJ Energy of this water: ", np.sum(energy_lj[:self.wat_oxygen_atom_ids[0]:])
                        #print "Solute-water Elec Energy of this water: ", np.sum(energy_elec[:, self.non_water_atom_ids])
                        # calc Enbr and other water structure metrics here
                        self.voxeldata[wat[0], 16] += np.sum(energy_lj[self.wat_oxygen_atom_ids[0]:][(wat_nbrs - self.wat_oxygen_atom_ids[0])/3])
                        for i in xrange(self.water_sites):
                            #print energy_elec[:, wat_nbrs + i].shape
                            self.voxeldata[wat[0], 16] += np.sum(energy_elec[:, wat_nbrs + i])
                    # H-bond calcuations
                    #TODO: document the algorithm
                    if hbonds:
                        hb_data = self.calculate_hydrogen_bonds(trj, wat[1], wat_nbrs, prot_nbrs)
                        self.voxeldata[wat[0], 24] += hb_data[0]
                        self.voxeldata[wat[0], 26] += hb_data[1]
                        self.voxeldata[wat[0], 28] += hb_data[2]
                        self.voxeldata[wat[0], 30] += hb_data[3]
                        self.voxeldata[wat[0], 32] += hb_data[4]
                        self.voxeldata[wat[0], 34] += hb_data[5]
                if entropy:
                    # calculate Euler angles
                    self.calculate_euler_angles(wat, pos[0, :, :])
                
            pbar.update(frame)
        pbar.finish()
        # get normalized and density-weighted values
        for voxel, waters in enumerate(self.voxeldict):
            #print voxel, waters
            if len(waters) >= 1.0:
                # deal with water-water energies
                self.voxeldata[voxel, 12] = self.voxeldata[voxel, 11] / len(waters)
                self.voxeldata[voxel, 11] /= (self.num_frames * self.voxel_vol * 2.0) 
                # deal with solute-water energies
                self.voxeldata[voxel, 14] = self.voxeldata[voxel, 13] / len(waters)
                self.voxeldata[voxel, 13] /= (self.num_frames * self.voxel_vol)
                # deal with water-water nbr energies
                if self.voxeldata[voxel, 18] > 0.0:
                    self.voxeldata[voxel, 16] = self.voxeldata[voxel, 15] / self.voxeldata[voxel, 18]
                    self.voxeldata[voxel, 15] /= (self.num_frames * self.voxel_vol * 2.0)
                # deal with rest
                for i in xrange(17, 35, 2):
                    self.voxeldata[voxel, i + 1] = self.voxeldata[voxel, i] / len(waters)
                    self.voxeldata[voxel, i] /= (self.num_frames * self.voxel_vol)
        if entropy:
            self.calculate_entropy()

    @function_timer
    def write_data(self, prefix=None):
        if prefix == None:
            prefix = self.prefix
        print "Writing voxel data ..."
        with open(prefix + "_gist_data.txt", "w") as f:
            gist_header = "voxel x y z nwat gO gH dTStr-dens dTStr-norm dTSor-dens dTSor-norm Esw-dens Esw-norm Eww-dens Eww-norm Eww-nbr-dens Eww-nbr-norm Nnbr-dens Nnbr-norm fHB-dens fHB-norm fenc-dens fenc-norm Nhbww-dens Nhbww-norm Nhbsw-dens Nhbsw-norm Ndonsw-dens Ndonsw-norm Naccsw-dens Naccsw-norm Ndonww-dens Ndonww-norm Naccww-dens Naccww-norm\n"
            f.write(gist_header)
            formatted_output_occupied_voxels = "{0[0]:.0f} {0[1]:.3f} {0[2]:.3f} {0[3]:.3f} {0[4]:.0f} {0[5]:.6f} {0[6]:.0f} " 
            formatted_output_empty_voxels = formatted_output_occupied_voxels
            for q in xrange(7, 35):
                formatted_output_occupied_voxels += "{0[%d]:.6f} " % q
                formatted_output_empty_voxels += "{0[%d]:.0f} " % q
            formatted_output_occupied_voxels += "\n"
            formatted_output_empty_voxels += "\n"
            for k in xrange(self.voxeldata.shape[0]):
                if self.voxeldata[k, 4] == 0.0:
                    f.write(formatted_output_empty_voxels.format(self.voxeldata[k, :]))
                else:
                    f.write(formatted_output_occupied_voxels.format(self.voxeldata[k, :]))

    @function_timer
    def generate_dx_files(self, prefix=None):
        if prefix == None:
            prefix = self.prefix
        print "Generating dx files ..."
        gist_header = "voxel x y z nwat gO gH dTStr-dens dTStr-norm dTSor-dens dTSor-norm Esw-dens Esw-norm Eww-dens Eww-norm Eww-nbr-dens Eww-nbr-norm Nnbr-dens Nnbr-norm fHB-dens fHB-norm fenc-dens fenc-norm Nhbww-dens Nhbww-norm Nhbsw-dens Nhbsw-norm Ndonsw-dens Ndonsw-norm Naccsw-dens Naccsw-norm Ndonww-dens Ndonww-norm Naccww-dens Naccww-norm\n"        
        dx_header = ""
        dx_header += 'object 1 class gridpositions counts %d %d %d\n' % (
            self.grid.shape[0], self.grid.shape[1], self.grid.shape[2])
        dx_header += 'origin %.3f %.3f %.3f\n' % (
            self.origin[0], self.origin[1], self.origin[2])
        dx_header += 'delta %.1f 0 0\n' % (self.spacing[0])
        dx_header += 'delta 0 %.1f 0\n' % (self.spacing[1])
        dx_header += 'delta 0 0 %.1f\n' % (self.spacing[2])
        dx_header += 'object 2 class gridconnections counts %d %d %d\n' % (
            self.grid.shape[0], self.grid.shape[1], self.grid.shape[2])
        dx_header += 'object 3 class array type float rank 0 items %d data follows\n' % (
            self.grid.shape[0] * self.grid.shape[1] * self.grid.shape[2])
        dx_file_objects = []
        
        data_keys = gist_header.strip("\n").split()

        for data_field, title in enumerate(data_keys):
            # if data_field > 4:# and data_field < 6:
            # print "Writing dx file for: ", title
            if data_field > 4:
                f = open(prefix + "_" + title + ".dx", 'w')
                f.write(dx_header)
                dx_file_objects.append(f)
            else:
                dx_file_objects.append(None)

        for k in xrange(1, len(self.voxeldata) + 1):
            # print "writing data for voxel: ", k
            if self.voxeldata[k - 1][4] > 1.0:
                for column_i in xrange(5, len(data_keys)):
                    dx_file_objects[column_i].write(
                        "%0.6f " % (self.voxeldata[k - 1][column_i]))
                    if k % 3 == 0:
                        dx_file_objects[column_i].write("\n")
            else:
                for column_i in xrange(5, len(data_keys)):
                    dx_file_objects[column_i].write(
                        "%i " % (self.voxeldata[k - 1][column_i]))
                    if k % 3 == 0:
                        dx_file_objects[column_i].write("\n")
        for f in dx_file_objects:
            if f is not None:
                f.close()

    def print_system_summary(self):
        print "System information:"
        print "\tParameter file: %s\n" % self.topology_file
        print "\tTrajectory: %s\n" % self.trajectory
        print "\tPeriodic Box: %s\n" % self.box_type
        print "\tFrames: %d, Total Atoms: %d, Waters: %d, Solute Atoms: %d\n" \
                % (self.num_frames, self.all_atom_ids.shape[0], self.wat_oxygen_atom_ids.shape[0], self.non_water_atom_ids.shape[0])
        print "\tH-bond Donors: %d, H-bond Acceptors: %d, H-bond Donors & Acceptors: %d\n" \
                % (self.solute_don_ids.shape[0], self.solute_acc_ids.shape[0], self.solute_acc_don_ids.shape[0])
        #print "\tWater Model: %s\n"
        print "Grid information:"
        print "\tGIST grid center: %5.3f %5.3f %5.3f\n" % (self.center[0], self.center[1], self.center[2])
        print "\tGIST grid dimensions: %i %i %i\n" % (self.dims[0], self.dims[1], self.dims[2])
        print "\tGIST grid spacing: %5.3f A^3\n" % (self.spacing[0])
    def print_calcs_summary(self):        
        print "Summary of main calculations:"
        nwat_grid = 0.0
        Eswtot = 0.0
        Ewwtot = 0.0
        dTStr_tot = 0.0
        dTSor_tot = 0.0
        for k in self.voxeldata:
            if k[4] > 1.0:
                nwat_grid += k[4] / (self.num_frames * self.voxel_vol)
                #print k[11]
                Ewwtot += k[11]
                Eswtot += k[13]
                dTStr_tot += k[7]
                dTSor_tot += k[9]

        nwat_grid *= self.voxel_vol
        Eswtot *= self.voxel_vol
        Ewwtot *= self.voxel_vol
        dTStr_tot *= self.voxel_vol
        dTSor_tot *= self.voxel_vol
        print "\tAverage number of water molecules over the grid: %d" % nwat_grid
        print "\tTotal Solute-Water Energy over the grid: %.6f" % Eswtot
        print "\tTotal Water-Water Energy over the grid: %.6f" % Ewwtot
        print "\tTotal Solute-Water Orientational Entropy over the grid: %.6f" % dTSor_tot
        print "\tTotal Solute-Water Translational Entropy over the grid: %.6f" % dTStr_tot

def entry_point():
    main()

if __name__ == '__main__':
    entry_point()
