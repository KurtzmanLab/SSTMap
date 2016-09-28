"""
docstring goes here
"""

# import statements
import sys
import time

import numpy as np
import mdtraj as md
import parmed as pmd


from water_analysis import WaterAnalysisBase
import _WATAnalysis_calcs as calc

DEG_PER_RAD = 180.0 / np.pi


class GridWaterAnalysis(WaterAnalysisBase):
    """docstring for SiteWaterAnalysis"""

    def __init__(self, topology_file, trajectory, start_frame=0, num_frames=0,
                 ligand_file=None, desmond_helper_file=None, prefix="test",
                 grid_center=[0.0, 0.0, 0.0], grid_dimensions=[5.0, 5.0, 5.0],
                 grid_resolution=[0.5, 0.5, 0.5]):

        super(
            GridWaterAnalysis,
            self).__init__(
            topology_file,
            trajectory,
            start_frame,
            num_frames,
            desmond_helper_file)
        # a strange way of checking if dimensions are non-zero
        # if self.ligand_file != None:
        print "Parsing parameters ..."
        # load in parameters
        self.paramname = topology_file
        self.param_class = pmd.load_file(self.paramname)
        self.resolution = grid_resolution[0]
        if ligand_file is not None:
            # TODO: change this as a utility function
            # get ligad center
            lig = md.load_pdb(ligand_file)
            com = np.zeros((lig.n_frames, 3))
            masses = np.ones(lig.n_atoms)
            masses /= masses.sum()
            com[0, :] = lig.xyz[0, :].astype('float64').T.dot(masses)
            grid_center = com[0, :] * 10.0
        print grid_center
        self.voxel_vol = self.resolution**3.0
        self.rho_bulk = 0.0329
        print "Done!"
        print "Setting up GIST grid"
        # set 3D grid around the region of interest
        self._initializeGrid(grid_center, grid_resolution, grid_dimensions)
        # initialize data structures to store voxel data
        self.voxeldata, self.voxeldict = self._initializeVoxelDict()
        # create data structures for storing water info
        voxel_wat_dim = np.array(
            [self.num_frames, len(self.wat_oxygen_atom_ids), 20])
        self._voxel_wat_map = np.zeros(voxel_wat_dim, dtype=np.int_)
        # create data structures for storing water info
        voxel_wat_euler_dim = np.array(
            [self.num_frames, len(self.wat_oxygen_atom_ids), 3])
        self._voxel_wat_euler_map = np.zeros(
            voxel_wat_euler_dim, dtype="float64")
        print "Done!"

    def _initializeGrid(self, center, resolution, dimensions):
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
        print "\tGIST grid center: %5.3f %5.3f %5.3f\n" % (self.center[0], self.center[1], self.center[2])
        print "\tGIST grid dimensions: %i %i %i\n" % (self.dims[0], self.dims[1], self.dims[2])
        print "\tGIST grid spacing: %5.3f A^3\n" % (self.spacing[0])

    def _initializeVoxelDict(self):
        voxel_dict = {}
        v_count = 0
        voxel_array = np.zeros((self.grid.size, 33), dtype="float64")
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
            voxel_dict[v_count] = [[]]
            # voxel_dict[v_count].append(np.zeros(14, dtype="float64"))
            v_count += 1
        return voxel_array, voxel_dict

    def _calcEulerAngles(self, frame_index, wat_index, coords):
        pi = np.pi
        twopi = 2 * np.pi
        # define the lab frame of reference
        xlab = np.asarray([1.0, 0.0, 0.0], dtype="float64")
        # ylab = np.asarray([0.0, 1.0, 0.0], dtype="float64")
        zlab = np.asarray([0.0, 0.0, 1.0], dtype="float64")
        # create array for water oxygen atom coords, and append to this voxel's
        # coord list
        wat_O = self._voxel_wat_map[frame_index, wat_index, 0]
        voxel_id = self._voxel_wat_map[frame_index, wat_index, 1]
        owat = coords[wat_O, :]
        # create array for water's hydrogen 1 and 2
        h1wat = coords[wat_O + 1, :] - owat
        h2wat = coords[wat_O + 2, :] - owat
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
        self._voxel_wat_euler_map[frame_index, wat_index, 0] = theta
        self._voxel_wat_euler_map[frame_index, wat_index, 1] = phi
        self._voxel_wat_euler_map[frame_index, wat_index, 2] = psi
        self.voxeldict[voxel_id][0].append([wat_index, frame_index])

    def _processHbonds(self, frame_index, wat_index, frame):
        # print self._voxel_wat_map[frame_index, wat_index, :]
        wat_O = self._voxel_wat_map[frame_index, wat_index, 0]
        voxel_id = self._voxel_wat_map[frame_index, wat_index, 1]
        N_nbr_wat = self._voxel_wat_map[frame_index, wat_index, 6]
        N_nbr_prot = self._voxel_wat_map[frame_index, wat_index, 7]
        hbww = 0
        hbsw = 0
        for i_wat_nbr in xrange(N_nbr_wat):
            # print wat_O, self._voxel_wat_map[frame_index, wat_index,
            # i_wat_nbr+8]
            nbr_O = self._voxel_wat_map[frame_index, wat_index, i_wat_nbr + 8]
            hbangle_combinations = np.asarray([[wat_O, nbr_O, nbr_O + 1], [wat_O, nbr_O, nbr_O + 2], [
                                              nbr_O, wat_O, wat_O + 1], [nbr_O, wat_O, wat_O + 2]])
            theta_list = md.compute_angles(
                frame, hbangle_combinations, opt=True) * DEG_PER_RAD
            hbangle, hbangle_index = np.min(theta_list), np.argmin(
                theta_list)  # min angle is a potential Hbond
            # print hbangle, hbangle_index
            if hbangle <= 30:  # if Hbond is made
                # print hbangle
                hbww += 1
                self.voxeldata[voxel_id, 22] += 1.0
                # further check if voxel water is acting as an acceptor
                if hbangle_index <= 1:
                    # add to cumulative of summ of HBwwacc
                    self.voxeldata[voxel_id, 28] += 1.0
                # or as a donor
                else:
                    # add to cumulative of summ of HBwwdon
                    self.voxeldata[voxel_id, 26] += 1.0
        for i_prot_nbr in xrange(N_nbr_prot):
            hbsw = 0
            nbr_prot = self._voxel_wat_map[
                frame_index, wat_index, i_prot_nbr + 8 + N_nbr_wat]
            nbr_prot_hbtype = self.prot_hb_types[nbr_prot]
            # Check protein acceptors
            if nbr_prot_hbtype == 1 or nbr_prot_hbtype == 3:
                hbangle_combinations = np.asarray(
                    [[nbr_prot, wat_O, wat_O + 1], [nbr_prot, wat_O, wat_O + 2]])
                theta_list = md.compute_angles(
                    frame, hbangle_combinations) * DEG_PER_RAD
                hbangle, hbangle_index = np.min(
                    theta_list), np.argmin(theta_list)
                if hbangle <= 30:
                    hbsw += 1
                    self.voxeldata[voxel_id, 24] += 1.0
                    # add to cumulative of summ of HBswdon
                    self.voxeldata[voxel_id, 30] += 1.0
            # Check protein donors, process each D-H pair separately
            if nbr_prot_hbtype == 2 or nbr_prot_hbtype == 3:
                for don_H_pair in self.don_H_pair_dict[nbr_prot]:
                    hbangle = md.compute_angles(frame, np.asarray(
                        [[wat_O, nbr_prot, don_H_pair[1]]]), opt=True) * DEG_PER_RAD
                    if hbangle <= 30:
                        hbsw += 1
                        self.voxeldata[voxel_id, 24] += 1.0
                        # add to cumulative of summ of HBswacc
                        self.voxeldata[voxel_id, 32] += 1.0
        # print hbwwat_index, hbsw

    # Energy Function
    def calculate_grid_quantities(
            self, energy=True, entropy=True, hbonds=True):
        getEnergy = int(energy)
        getEntropy = int(entropy)
        getHbonds = int(hbonds)
        t = time.time()
        print "******************************"
        print "Starting Energy and H-bonds calculations: "

        for i in xrange(self.start_frame, self.start_frame + self.num_frames):
            print "Processing frame: ", i
            trj = md.load_frame(self.trajectory, i, top=self.paramname)
            pos = trj.xyz[0, :, :] * 10.0
            periodicBox = trj.unitcell_lengths * 10.0
            calc.processGrid(getEnergy, getEntropy, getHbonds, i - self.start_frame, self.wat_O_sites,
                             pos, periodicBox,
                             self.vdw, self.chg,
                             self.dims, self.gridmax, self.origin, self._voxel_wat_map, self.voxeldata,
                             self.wat_oxygen_atom_ids, self.non_water_atom_ids, self.all_atom_ids,
                             self.prot_hb_types)

            for wat_index in xrange(self._voxel_wat_map.shape[1]):
                if self._voxel_wat_map[
                        i - self.start_frame, wat_index, 2] == 1:
                    self._processHbonds(i - self.start_frame, wat_index, trj)
                    self._calcEulerAngles(i - self.start_frame, wat_index, pos)

        print "Energy and H-bond calcs took seconds.", time.time() - t
        print "*******************************"
        print "Starting Entropy calculations: "
        t = time.time()
        self._normalizeVoxelQuantities()
        self._getVoxelEntropies()
        print "Entropy calcs took seconds.", time.time() - t
        print "*******************************"

    def _getVoxelEntropies(self):
        # set constants
        gas_kcal = 0.0019872041
        # initialize variables to store entropy results
        dTStr_tot = 0.0
        dTSor_tot = 0.0
        for k in self.voxeldict.keys():
            if self.voxeldata[k, 4] > 1.0:
                # print self.voxeldata[k, 4], self.voxeldict[k][0]
                dTStr_dens = -gas_kcal * 300 * self.rho_bulk * \
                    self.voxeldata[k, 5] * np.log(self.voxeldata[k, 5])
                # voxel translational entropy density
                self.voxeldata[k, 7] = dTStr_dens
                # voxel translational entropy normalized
                self.voxeldata[k, 8] = dTStr_dens * self.num_frames * \
                    self.voxel_vol / (1.0 * self.voxeldata[k, 4])

                dTSor_dens = calc.getNNEntropy(k, int(self.voxeldata[k, 4]), np.asarray(
                    self.voxeldict[k][0]), self._voxel_wat_euler_map, self.voxeldata)
                # voxel orientational entropy normalized
                self.voxeldata[k, 10] = gas_kcal * 300 * \
                    ((dTSor_dens / self.voxeldata[k, 4]) + 0.5772)
                # voxel orientational entropy density
                self.voxeldata[k, 9] = self.voxeldata[
                    k, 10] * (self.voxeldata[k, 4] / (self.num_frames * self.voxel_vol))

                dTStr_tot += self.voxeldata[k, 7]
                dTSor_tot += self.voxeldata[k, 9]

        dTStr_tot *= self.voxel_vol
        dTSor_tot *= self.voxel_vol
        print "Total Solute-Water Orientational Entropy over the grid: %.8f" % dTSor_tot
        print "Total Solute-Water Translational Entropy over the grid: %.8f" % dTStr_tot

    def _normalizeVoxelQuantities(self):
        bulkwaterpervoxel = self.voxel_vol * self.rho_bulk * self.num_frames

        Eswtot = 0
        Ewwtot = 0
        Ewwnbrtot = 0
        nwat = 0
        maxwat = 0
        for k in self.voxeldata:
            if k[4] > 1.0:
                # number density of oxygen centers
                k[5] = k[4] / bulkwaterpervoxel
                # number density of hydrogens
                #k[6] = k[6]/(bulkwaterpervoxel*2)
                if k[4] > maxwat:
                    maxwat = k[4]
                # energy density-weighted and normalized
                k[11] = k[12] / (self.num_frames * self.voxel_vol)  # E_sw_dens
                k[12] = k[12] / k[4]  # E_sw_norm
                Eswtot += k[11]

                k[13] = k[14] / \
                    (self.num_frames * self.voxel_vol * 2.0)  # E_ww_dens
                k[14] = k[14] / (k[4] * 2.0)  # E_ww_norm
                Ewwtot += k[13]

                if k[18] > 1.0:
                    # E_nbr
                    k[15] = k[16] / (self.num_frames * self.voxel_vol * 2.0)
                    k[16] = k[16] / k[18]
                    Ewwnbrtot += k[15]

                    # neighbours
                    k[17] = k[18] / (self.num_frames *
                                     self.voxel_vol)  # nbr_dens
                    k[18] = k[18] / k[4]  # nbr_norm

                #enclosure = 1 - (k[18]/5.25)
                # if enclosure < 0.0:
                #    enclosure = 0.0
                #k[19] += enclosure
                # k[20] = k[19]/(self.num_frames*self.voxel_vol) #
                # Enclosure-dens

                # H-bonds
                k[21] = k[22] / (self.num_frames *
                                 self.voxel_vol)  # HB_ww_dens
                k[22] = k[22] / k[4]  # HB_ww_norm
                k[23] = k[24] / (self.num_frames *
                                 self.voxel_vol)  # HB_sw_dens
                k[24] = k[24] / k[4]  # HB_sw_norm
                # Don_ww
                k[25] = k[26] / (self.num_frames * self.voxel_vol)
                k[26] = k[26] / k[4]
                # Acc_ww
                k[27] = k[28] / (self.num_frames * self.voxel_vol)
                k[28] = k[28] / k[4]
                # Don_sw
                k[29] = k[30] / (self.num_frames * self.voxel_vol)
                k[30] = k[30] / k[4]
                # Acc_sw
                k[31] = k[32] / (self.num_frames * self.voxel_vol)
                k[32] = k[32] / k[4]
                nwat += k[4] / (self.num_frames * self.voxel_vol)

        Eswtot *= self.voxel_vol
        Ewwtot *= self.voxel_vol
        Ewwnbrtot *= self.voxel_vol
        nwat *= self.voxel_vol
        print "Average number of water molecules over the grid: ", nwat
        print "Maximum number of waters found in one voxel for %i frames = %i" % (self.num_frames, maxwat)
        print "Total Solute-Water Energy over the grid: ", Eswtot
        print "Total Water-Water Energy over the grid: ", Ewwtot
        # print "Total Water-Water Nbr Energy over the grid: ", Ewwnbrtot

    def writeGistData(self, outfile):
        f = open(outfile + "_gist_data.txt", "w")
        gist_header = "voxel x y z wat g_O g_H dTStr-dens dTStr-norm dTSor-dens dTSor-norm Esw-dens Esw-norm Eww-dens Eww-norm Eww-nbr-dens Eww-nbr-norm nbr-dens nbr-norm enc-dens enc-norm hbww-dens hbww-norm hbsw-dens hbsw-norm donww-dens donww-norm accww-dens accww-norm donsw-dens donsw-norm accsw-dens accsw-norm\n"
        gist_header_list = ["voxel", "x", "y", "z", "wat", "g_O", "g_H",
                            "dTStr-dens", "dTStr-norm", "dTSor-dens", "dTSor-norm",
                            "Esw-dens", "Esw-norm", "Eww-dens", "Eww-norm", "Eww-nbr-dens", "Eww-nbr-norm",
                            "nbr-dens", "nbr-norm", "enc-dens", "enc-norm",
                            "hbww-dens", "hbww-norm", "hbsw-dens", "hbsw-norm",
                            "donww-dens", "donww-norm", "accww-dens", "accww-norm",
                            "donsw-dens", "donsw-norm", "accsw-dens", "accsw-norm"]
        f.write(gist_header)
        for k in self.voxeldata:
            if k[4] >= 1.0:
                l = "%i %.3f %.3f %.3f %i %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n" % \
                    (k[0], k[1], k[2], k[3], k[4], k[5], k[6],
                     k[7], k[8], k[9], k[10],
                     k[11], k[12], k[13], k[14], k[15], k[
                         16], k[17], k[18], k[19], k[20],
                     k[21], k[22], k[23], k[24],
                     k[25], k[26], k[27], k[28],
                     k[29], k[30], k[31], k[32])
                # print l
                f.write(l)
            else:
                l = "%i %.3f %.3f %.3f %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n" % \
                    (k[0], k[1], k[2], k[3], k[4], k[5], k[6],
                     k[7], k[8], k[9], k[10],
                     k[11], k[12], k[13], k[14], k[15], k[
                         16], k[17], k[18], k[19], k[20],
                     k[21], k[22], k[23], k[24],
                     k[25], k[26], k[27], k[28],
                     k[29], k[30], k[31], k[32])
                # print l
                f.write(l)
        f.close()
        # write dx files

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
                f = open(outfile + "_" + title + ".dx", 'w')
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

import time
import resource
import sys
from optparse import OptionParser


def main():
    parser = OptionParser()
    parser.add_option(
        "-i",
        "--input_prmtop",
        dest="prmtop",
        type="string",
        help="Input amber prmtop file")
    parser.add_option(
        "-t",
        "--input_trajectory",
        dest="trjname",
        type="string",
        help="Input trajectory file")
    parser.add_option(
        "-l",
        "--ligand_file",
        dest="ligand",
        type="string",
        help="Ligand file")
    parser.add_option(
        "-f",
        "--frames",
        dest="frames",
        type="int",
        help="Number of frames")
    parser.add_option(
        "-s",
        "--starting frame",
        dest="start_frame",
        type="int",
        help="Starting frame")
    parser.add_option(
        "-o",
        "--output_name",
        dest="prefix",
        type="string",
        help="Output log file")
    (options, args) = parser.parse_args()
    print "Setting things up..."
    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        print "no argument given!"
        parser.print_help()
    else:
        print "Setting up calculations..."
        g = GridWaterAnalysis(
            options.prmtop,
            options.trjname,
            start_frame=options.start_frame,
            num_frames=options.frames,
            ligand_file=options.ligand,
            prefix=options.prefix,
            grid_dimensions=[
                40,
                40,
                40])
        print "\tMemory consumption at Initialization: %5.3f MB" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000)
        gist_logfile = open("desmond-gist-energy.log", "w")
        gist_logfile.write("#Grid setup for the system in DX header format:\n")
        # gist_logfile.write('# Data calculated by the VMD volmap function\n')
        gist_logfile.write(
            'object 1 class gridpositions counts %d %d %d\n' %
            (g.grid.shape[0], g.grid.shape[1], g.grid.shape[2]))
        gist_logfile.write(
            'origin %.1f %.1f %.1f\n' %
            (g.origin[0], g.origin[1], g.origin[2]))
        gist_logfile.write('delta %.1f 0 0\n' % (g.spacing[0]))
        gist_logfile.write('delta 0 %.1f 0\n' % (g.spacing[1]))
        gist_logfile.write('delta 0 0 %.1f\n' % (g.spacing[2]))
        gist_logfile.write(
            'object 2 class gridconnections counts %d %d %d\n' %
            (g.grid.shape[0], g.grid.shape[1], g.grid.shape[2]))
        gist_logfile.write(
            'object 3 class array type double rank 0 items %d data follows\n' %
            (g.grid.shape[0] * g.grid.shape[1] * g.grid.shape[2]))
        # gist_logfile.write("#EndHeader\n")
        print "Performing energy calculations ..."
        t = time.time()
        g.calculate_grid_quantities()
        print "Took seconds.", time.time() - t
        g.writeGistData(options.prefix)
        gist_logfile.close()


def entry_point():
    main()

if __name__ == '__main__':
    entry_point()
