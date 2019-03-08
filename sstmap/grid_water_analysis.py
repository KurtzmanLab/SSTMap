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

import sys
import math
import numpy as np
import mdtraj as md

from sstmap.water_analysis import WaterAnalysis
from sstmap.utils import print_progress_bar, function_timer
import _sstmap_ext as calc
from sstmap.utils import *

GASKCAL = 0.0019872041

class GridWaterAnalysis(WaterAnalysis):
    @function_timer
    def __init__(self, topology_file, trajectory, start_frame=0, num_frames=0,
                 supporting_file=None, ligand_file=None,
                 grid_center=None, grid_dimensions=[20, 20, 20],
                 grid_resolution=[0.5, 0.5, 0.5], rho_bulk=0.0334, prefix="test"):
        """
        Initialize an instance of grid-based solvation analysis calculations.

        Parameters
        ----------
        topology_file : str
            File name of system topology.
        trajectory : str
            File name of molecular dynamics trajectory.
        start_frame : int, optional
            Index of the frame to begin calculations from. Default=0
        num_frames : int, optional
            Total number of frames to be processed. Default=0
        supporting_file : str, optional
            File name of supporting files for parsing system topology. Default=None
        ligand_file : str, optional
            PDB File containing ligand atom. Default=None
        grid_center : list, optional
            List containing x, y and z coordinates of the center of grid. Default=None
        grid_dimensions : list, optional
            List containing dimensions of grid in each direction, by default a cubic grid of
            20 points in each direction is used.
        grid_resolution : list, option
            Spacing of grid points, by default a grid spacing of 0.5 Angstrom is used.
        rho_bulk : float, optional
            Reference water density, if None, 0.0334 molecules per angstrom cubed is used.
        prefix : str, optional
            String to be used as a prefix for the output file names, if None, "test" is used.
        """

        print("Initializing ...")
        self.start_frame = start_frame
        self.num_frames = num_frames
        self.rho_bulk = float(rho_bulk)
        super(GridWaterAnalysis, self).__init__(topology_file, trajectory, supporting_file)

        self.grid_dims = np.asarray(grid_dimensions, int)
        self.resolution = grid_resolution[0]
        self.prefix = prefix
        if ligand_file is None and grid_center is None:
            sys.exit("Please provide value of the grid center as a list of x, y, z coordinates or\
                         speicify a ligand PDB file whose center would be chosen as grid center.")

        if ligand_file is not None and grid_center is None:
            # TODO: change this to a utility function
            # get ligad center
            lig = md.load_pdb(ligand_file, no_boxchk=True)
            com = np.zeros((lig.n_frames, 3))
            masses = np.ones(lig.n_atoms)
            masses /= masses.sum()
            com[0, :] = lig.xyz[0, :].astype('float64').T.dot(masses)
            grid_center = com[0, :] * 10.0
        self.voxel_vol = self.resolution ** 3.0
        # set 3D grid around the region of interest
        self.initialize_grid(grid_center, grid_resolution, grid_dimensions)
        # initialize data structures to store voxel data
        self.voxeldata, self.voxel_quarts, self.voxel_O_coords = self.initialize_voxel_data()
        # print "Reading in trajectory ..."
        # self.trj = md.load(self.trajectory, top=self.paramname)[self.start_frame: self.start_frame + self.num_frames]
        # print "Done!"

    def initialize_grid(self, center, resolution, dimensions):
        """
        Initialize grid data structure.

        Parameters
        ----------
        center : list
            List containing x, y and z coordinates of the center of grid.
        resolution : list
            Spacing of grid points in each direction.
        dimensions : list
            List containing dimensions of grid in each direction.
        """

        # set grid center, res and dimension
        print("Initializing ...")
        self.center = np.array(center, dtype=np.float_)
        self.dims = np.array(dimensions, dtype=np.int_)
        self.spacing = np.array(resolution, dtype=np.float_)
        self.gridmax = self.dims * self.spacing + 1.5        # set origin
        o = self.center - (0.5 * self.dims * self.spacing)
        self.origin = np.around(o, decimals=3)

        # set grid size (in terms of total points alog each axis)
        length = np.array(self.dims / self.spacing, dtype=np.float_)
        self.grid_size = np.ceil((length / self.spacing) + 1.0)
        self.grid_size = np.cast['uint32'](self.grid_size)

        # Finally allocate the space for the grid
        self.grid = np.zeros(self.dims, dtype=np.int_)

    def initialize_voxel_data(self):
        """
        Initializes data structures where GIST data is stored

        Returns
        -------
        voxel_array : numpy.ndarray
            A numpy array with rows equal to the number of voxels and columns corresponding
            to various properties of each voxel.
        voxel_quarts : list
            A list containing N empty list, where N is equal to the number of voxels, each empty
            list will store quaternions of each water found in the corresponding voxel during the simulation.
        voxel_coords : list
            A list containing N empty list, where N is equal to the number of voxels, each empty
            list will store coordinates of the oxygen of each water found in the corresponding voxel
             during the simulation.
        """

        v_count = 0
        voxel_array = np.zeros((self.grid.size, 35), dtype="float64")
        for index, value in np.ndenumerate(self.grid):
            _index = np.array(index, dtype=np.int32)
            point = _index * self.spacing + self.origin + 0.5 * self.spacing
            voxel_array[v_count, 1] = point[0]
            voxel_array[v_count, 2] = point[1]
            voxel_array[v_count, 3] = point[2]
            voxel_array[v_count, 0] = v_count
            v_count += 1
        voxel_quarts = [[] for i in range(voxel_array.shape[0])]
        voxel_O_coords = [[] for i in range(voxel_array.shape[0])]
        return voxel_array, voxel_quarts, voxel_O_coords

    def calculate_euler_angles(self, water, coords):
        """
        Calculates quaternion representing orientation of a water molecule.
        Adapted from: https://github.com/Amber-MD/cpptraj/blob/master/src/Action_GIST.cpp

        Parameters
        ----------
        water : tuple
            A tuple with two elements, index of the oxygen atom of current water and frame number.
        coords : numpy.ndarray
            An array of x, y, z coordinates corresponding to atom positions in current frame.
        """
        # define the lab frame of reference
        xlab = np.asarray([1.0, 0.0, 0.0])
        zlab = np.asarray([0.0, 0.0, 1.0])
        # create array for water oxygen atom coords, and append to this voxel's
        voxel_id = water[0]
        owat = coords[water[1], :]
        # create array for water's hydrogen 1 and 2
        h1wat = coords[water[1] + 1, :] - owat
        h2wat = coords[water[1] + 2, :] - owat
        # print frame_index, wat_O, owat, h1wat, h2wat
        # define water molecule's frame
        # H1 is water's x-axis, should be normalized
        h1wat /= np.linalg.norm(h1wat)
        h2wat /= np.linalg.norm(h2wat)
        ar1 = np.cross(h1wat, xlab)
        sar = np.copy(ar1)
        ar1 /= np.linalg.norm(ar1)
        dp1 = np.sum(xlab * h1wat)
        theta = np.arccos(dp1)
        sign = np.sum(sar * h1wat)
        if sign > 0:
            theta /= 2.0
        else:
            theta /= -2.0

        w1 = np.cos(theta)
        sin_theta = np.sin(theta)
        x1 = ar1[0] * sin_theta
        y1 = ar1[1] * sin_theta
        z1 = ar1[2] * sin_theta
        w2 = w1
        x2 = x1
        y2 = y1
        z2 = z1

        H_temp = np.zeros(3)
        H_temp[0] = ((w2*w2+x2*x2)-(y2*y2+z2*z2))*h1wat[0]
        H_temp[0] = (2*(x2*y2 - w2*z2)*h1wat[1]) + H_temp[0]
        H_temp[0] = (2*(x2*z2-w2*y2)*h1wat[2]) + H_temp[0]

        H_temp[1] = 2*(x2*y2 - w2*z2)* h1wat[0]
        H_temp[1] = ((w2*w2-x2*x2+y2*y2-z2*z2)*h1wat[1]) + H_temp[1]
        H_temp[1] = (2*(y2*z2+w2*x2)*h1wat[2]) +H_temp[1]

        H_temp[2] = 2*(x2*z2+w2*y2) * h1wat[0]
        H_temp[2] = (2*(y2*z2-w2*x2)*h1wat[1]) + H_temp[2]
        H_temp[2] = ((w2*w2-x2*x2-y2*y2+z2*z2)*h1wat[2]) + H_temp[2]

        H_temp2 = np.zeros(3,)
        H_temp2[0] = ((w2*w2+x2*x2)-(y2*y2+z2*z2))*h2wat[0]
        H_temp2[0] = (2*(x2*y2 + w2*z2)*h2wat[1]) + H_temp2[0]
        H_temp2[0] = (2*(x2*z2-w2*y2)+h2wat[2]) +H_temp2[0]

        H_temp2[1] = 2*(x2*y2 - w2*z2) *h2wat[0]
        H_temp2[1] = ((w2*w2-x2*x2+y2*y2-z2*z2)*h2wat[1]) +H_temp2[1]
        H_temp2[1] = (2*(y2*z2+w2*x2)*h2wat[2]) +H_temp2[1]

        H_temp2[2] = 2*(x2*z2+w2*y2)*h2wat[0]
        H_temp2[2] = (2*(y2*z2-w2*x2)*h2wat[1]) +H_temp2[2]
        H_temp2[2] = ((w2*w2-x2*x2-y2*y2+z2*z2)*h2wat[2]) + H_temp2[2]

        ar2 = np.cross(H_temp, H_temp2)
        ar2 /= np.linalg.norm(ar2)
        dp2 = np.sum(ar2 * zlab)
        theta = np.arccos(dp2)

        sar = np.cross(ar2, zlab)
        sign = np.sum(sar * H_temp)

        if sign < 0:
          theta /= 2.0
        else:
          theta /= -2.0

        w3 = np.cos(theta)
        sin_theta = np.sin(theta)
        x3 = xlab[0] * sin_theta
        y3 = xlab[1] * sin_theta
        z3 = xlab[2] * sin_theta

        w4 = w1*w3 - x1*x3 - y1*y3 - z1*z3
        x4 = w1*x3 + x1*w3 + y1*z3 - z1*y3
        y4 = w1*y3 - x1*z3 + y1*w3 + z1*x3
        z4 = w1*z3 + x1*y3 - y1*x3 + z1*w3
        self.voxel_quarts[voxel_id].extend([w4, x4, y4, z4])
        self.voxel_O_coords[voxel_id].extend(owat)

    @function_timer
    def calculate_entropy(self, num_frames=None):
        """
        Calculate solute-water translational and orientational entropy for each grid voxel using a nearest neighbors
        approach. This calculation can only be run after the voxels are populated with corresponding waters.

        Parameters
        ----------
        num_frames : int
            Number of frames processed during the analysis.

        """
        if num_frames is None:
            num_frames = self.num_frames
        calc.getNNTrEntropy(num_frames, self.voxel_vol, self.rho_bulk, 300.0, self.grid_dims, self.voxeldata, self.voxel_O_coords, self.voxel_quarts)


    def _process_frame(self, trj, energy, hbonds, entropy):
        """
        Frame wise calculation of GIST quantities.

        Parameters
        ----------
        trj :
            Molecular dynamic trajectory representing current frame.
        energy :
            If True, solute-water and water-water energies are calculated for each water in each voxel in current
            frame.
        hbonds : bool
            If True, solute-water and water-water hydrogen bonds are calculated for each water in each voxel in current
            frame.
        entropy : bool
            If True, water coordinates and quaternions are stored for each water in each voxel in current frame.
        """


        nbr_cutoff_sq = 3.5 ** 2
        trj.xyz *= 10.0
        coords = trj.xyz
        uc = trj.unitcell_vectors[0]*10.
        waters = []
        calc.assign_voxels(trj.xyz, self.dims, self.gridmax, self.origin, waters, self.wat_oxygen_atom_ids)

        distance_matrix = np.zeros((self.water_sites, self.all_atom_ids.shape[0]))

        for wat in waters:
            self.voxeldata[wat[0], 4] += 1
            ### wat[0]: Voxel index
            ### wat[1]: water molecule oxygen atom index 
            if energy or hbonds:
                e_lj_array, e_elec_array = np.copy(self.acoeff), np.copy(self.chg_product)
                ### Here it is assumed that oxygen atom is always the first atom. Is that always the case?
                ### We should probably check that somewhere during init?
                ###
                ### The self.neighbor_ids array contains atom indices of all atoms that should
                ### be considered as potential neighbors. Valid_neighbors has same shape as 
                ### self.neighbor_ids and is False at the position where index wat occures in
                ### self.neighbor_ids, otherwise it is True. neighbor_ids stores the indices of
                ### the actual neighbor candidates that will be commited to get_pairwise_distances
                ### routine and has length of self.neighbor_ids-1. wat_nbrs_shell is of length neighbor_ids
                ### and holds the shell_index of each neighbor candidate atom (0:first shell, 1: beyond first
                ### shell)
                valid_neighbors = np.ones(self.neighbor_ids.shape[0], dtype=bool)
                valid_neighbors[np.where(self.neighbor_ids==wat)] = False
                neighbor_ids   = self.neighbor_ids[valid_neighbors]
                wat_nbrs_shell = self.wat_nbrs_shell[valid_neighbors]
                calc.get_pairwise_distances(wat, self.all_atom_ids, 
                                            np.array([nbr_cutoff_sq]), neighbor_ids, wat_nbrs_shell,
                                            coords, uc, distance_matrix, 0)
                wat_nbrs = self.all_atom_ids[np.where(wat_nbrs_shell==0)]
                self.voxeldata[wat[0], 19] += wat_nbrs.shape[0]
                calc.calculate_energy(wat[1], distance_matrix, e_elec_array, e_lj_array, self.bcoeff)

                if self.prot_atom_ids.shape[0] != 0:
                    #self.voxeldata[wat[0], 13] += np.sum(e_lj_array[:, self.prot_atom_ids])
                    #self.voxeldata[wat[0], 13] += np.sum(e_elec_array[:, self.prot_atom_ids])
                    self.voxeldata[wat[0], 13] += np.sum(e_lj_array[:, self.non_water_atom_ids])
                    self.voxeldata[wat[0], 13] += np.sum(e_elec_array[:, self.non_water_atom_ids])

                self.voxeldata[wat[0], 15] += np.sum(
                    e_lj_array[:, self.wat_oxygen_atom_ids[0]:wat[1]]) + np.sum(e_lj_array[:, wat[1] + self.water_sites:])
                self.voxeldata[wat[0], 15] += np.sum(
                    e_elec_array[:, self.wat_oxygen_atom_ids[0]:wat[1]]) + np.sum(
                    e_elec_array[:, wat[1] + self.water_sites:])
                e_nbr_list = [np.sum(e_lj_array[:, wat_nbrs + i] + e_elec_array[:, wat_nbrs + i]) for i in
                              range(self.water_sites)]
                self.voxeldata[wat[0], 17] += np.sum(e_nbr_list)


                ### Might be usefull for API to have the neighbors and shell
                ### indices available.
                self.wat_nbrs_shell[valid_neighbors] = wat_nbrs_shell
                self.neighbor_ids[valid_neighbors] = neighbor_ids
                """
                ###DEBUG START###
                elj_sw = np.sum(e_lj_array[:, :self.wat_oxygen_atom_ids[0]])
                eelec_sw = np.sum(e_elec_array[:, :self.wat_oxygen_atom_ids[0]])
                elj_ww = np.sum(e_lj_array[:, self.wat_oxygen_atom_ids[0]:wat[1]]) + np.sum(e_lj_array[:, wat[1] + 1:])
                eelec_ww = np.sum(e_elec_array[:, self.wat_oxygen_atom_ids[0]:wat[1]]) + np.sum(e_elec_array[:, wat[1] + self.water_sites:])
                e_nbr_list = [np.sum(e_lj_array[:, wat_nbrs + i] + e_elec_array[:, wat_nbrs + i]) for i in xrange(self.water_sites)]
                enbr = np.sum(e_nbr_list)
                print "Calc: ", elj_sw, eelec_sw, elj_ww, eelec_ww, enbr
                distance_matrix = np.sqrt(distance_matrix)
                energy_lj, energy_elec = self.calculate_energy(distance_matrix)
                test_1 = np.sum(energy_lj[:self.wat_oxygen_atom_ids[0]:])
                test_2 = np.sum(energy_elec[:, self.non_water_atom_ids])
                test_3 = np.nansum(energy_lj[self.wat_oxygen_atom_ids[0]:])
                test_4 = np.sum(energy_elec[:, self.wat_atom_ids[0]:wat[1]]) + np.sum(energy_elec[:, wat[1] + self.water_sites:])
                test_5 = 0.0
                test_5 += np.sum(energy_lj[self.wat_oxygen_atom_ids[0]:][(wat_nbrs - self.wat_oxygen_atom_ids[0]) / self.water_sites])
                for i in range(self.water_sites):
                    test_5 += np.sum(energy_elec[:, wat_nbrs + i])
                print "Ref: ", test_1, test_2, test_3, test_4, test_5
                ###DEBUG END###
                """

                # H-bond calculations
                if hbonds:
                    prot_nbrs_all = self.prot_atom_ids[
                        np.where(distance_matrix[0, :][self.prot_atom_ids] <= nbr_cutoff_sq)]
                    prot_nbrs_hb = prot_nbrs_all[np.where(self.prot_hb_types[prot_nbrs_all] != 0)]

                    if wat_nbrs.shape[0] > 0:
                        hb_ww = self.calculate_hydrogen_bonds(trj, wat[1], wat_nbrs)
                        acc_ww = hb_ww[:, 0][np.where(hb_ww[:, 0] == wat[1])].shape[0]
                        don_ww = hb_ww.shape[0] - acc_ww
                        self.voxeldata[wat[0], 25] += hb_ww.shape[0]
                        self.voxeldata[wat[0], 31] += don_ww
                        self.voxeldata[wat[0], 33] += acc_ww
                        if wat_nbrs.shape[0] != 0 and hb_ww.shape[0] != 0:
                            self.voxeldata[wat[0], 21] += wat_nbrs.shape[0] / hb_ww.shape[0]

                    if prot_nbrs_hb.shape[0] > 0:
                        hb_sw = self.calculate_hydrogen_bonds(trj, wat[1], prot_nbrs_hb, water_water=False)
                        acc_sw = hb_sw[:, 0][np.where(hb_sw[:, 0] == wat[1])].shape[0]
                        don_sw = hb_sw.shape[0] - acc_sw
                        self.voxeldata[wat[0], 23] += hb_sw.shape[0]
                        self.voxeldata[wat[0], 27] += don_sw
                        self.voxeldata[wat[0], 29] += acc_sw

            if entropy:
                self.calculate_euler_angles(wat, coords[0, :, :])

    @function_timer
    def calculate_grid_quantities(self, energy=True, entropy=True, hbonds=True):
        """
        Performs grid-based solvation thermodynamics and structure calculations by iterating
        over frames in the trajectory.

        Parameters
        ----------
        energy : bool, optional

        entropy :

        hbonds :

        Returns
        -------

        """
        print_progress_bar(0, self.num_frames)
        if not self.topology_file.endswith(".h5"):
            topology = md.load_topology(self.topology_file)
        read_num_frames = 0
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
                    self._process_frame(trj, energy, hbonds, entropy)
                    read_num_frames += 1
            if read_num_frames < self.num_frames:
                print(("{0:d} frames found in the trajectory, resetting self.num_frames.".format(read_num_frames)))
                self.num_frames = read_num_frames

        # Normalize voxel quantities
        for voxel in range(self.voxeldata.shape[0]):
            if self.voxeldata[voxel, 4] > 1.0:
                self.voxeldata[voxel, 14] = self.voxeldata[voxel, 13] / (self.voxeldata[voxel, 4] * 2.0)
                self.voxeldata[voxel, 13] /= (self.num_frames * self.voxel_vol * 2.0)
                self.voxeldata[voxel, 16] = self.voxeldata[voxel, 15] / (self.voxeldata[voxel, 4] * 2.0)
                self.voxeldata[voxel, 15] /= (self.num_frames * self.voxel_vol * 2.0)
                if self.voxeldata[voxel, 19] > 0.0:
                    self.voxeldata[voxel, 18] = self.voxeldata[voxel, 17] / (self.voxeldata[voxel, 19] * 2.0)
                    self.voxeldata[voxel, 17] /= (self.num_frames * self.voxel_vol * self.voxeldata[voxel, 19] * 2.0)
                for i in range(19, 35, 2):
                    self.voxeldata[voxel, i + 1] = self.voxeldata[voxel, i] / self.voxeldata[voxel, 4]
                    self.voxeldata[voxel, i] /= (self.num_frames * self.voxel_vol)
            else:
                self.voxeldata[voxel, 13] *= 0.0
                self.voxeldata[voxel, 15] *= 0.0
                if self.voxeldata[voxel, 19] > 0.0:
                    self.voxeldata[voxel, 17] *= 0.0
                for i in range(19, 35, 2):
                    self.voxeldata[voxel, i] *= 0.0

        # Calculate entropies
        if entropy:
            self.calculate_entropy(num_frames=self.num_frames)

    @function_timer
    def write_data(self, prefix=None):
        """
        Write a tabular summary of GIST data for each voxel as a text file.

        Parameters
        ----------
        prefix : str, optional
            Prefix to be used for naming output dx files, if None, then prefix set at the time of initializing GIST
            object is used.

        """

        if prefix == None:
            prefix = self.prefix
        print("Writing voxel data ...")
        with open(prefix + "_gist_data.txt", "w") as f:
            gist_header = "voxel x y z nwat gO dTStr-dens dTStr-norm dTSor-dens dTSor-norm dTSsix-dens dTSsix-norm Esw-dens Esw-norm Eww-dens Eww-norm Eww-nbr-dens Eww-nbr-norm Nnbr-dens Nnbr-norm fHB-dens fHB-norm Nhbsw_dens Nhbsw_norm Nhbww_dens Nhbww_norm Ndonsw_dens Ndonsw_norm Naccsw_dens Naccsw_norm Ndonww_dens Ndonww_norm Naccww_dens Naccww_norm\n"
            f.write(gist_header)
            formatted_output_occupied_voxels = "{0[0]:.0f} {0[1]:.3f} {0[2]:.3f} {0[3]:.3f} {0[4]:.0f} {0[5]:.6f} "
            formatted_output_one_voxels = formatted_output_occupied_voxels
            formatted_output_empty_voxels = "{0[0]:.0f} {0[1]:.3f} {0[2]:.3f} {0[3]:.3f} {0[4]:.0f} {0[5]:.0f} "
            for q in range(7, 35):
                formatted_output_occupied_voxels += "{0[%d]:.6f} " % q
                formatted_output_empty_voxels += "{0[%d]:.0f} " % q
                if q in [19, 20]:
                    formatted_output_one_voxels += "{0[%d]:.3f} " % q
                elif q in [7, 8, 11, 12]:
                    formatted_output_one_voxels += "{0[%d]:.6f} " % q
                else:
                    formatted_output_one_voxels += "{0[%d]:.0f} " % q
            formatted_output_occupied_voxels += "\n"
            formatted_output_one_voxels += "\n"
            formatted_output_empty_voxels += "\n"
            for k in range(self.voxeldata.shape[0]):
                if self.voxeldata[k, 4] == 0.0:
                    f.write(formatted_output_empty_voxels.format(self.voxeldata[k, :]))
                elif self.voxeldata[k, 4] == 1.0:
                    mask_one_voxel_data = np.zeros(self.voxeldata[k, :].shape[0])
                    mask_one_voxel_data[list(range(0, 6)) + [7, 8, 11, 12] + list(range(19, 21))] = self.voxeldata[k, list(range(0, 6)) + [7, 8, 11, 12] + list(range(19, 21))]
                    f.write(formatted_output_one_voxels.format(mask_one_voxel_data))
                else:
                    f.write(formatted_output_occupied_voxels.format(self.voxeldata[k, :]))

    @function_timer
    def generate_dx_files(self, prefix=None):
        """
        Write GIST grid quantities in Data Explorer (DX) file format.

        Parameters
        ----------
        prefix : str, optional
            Prefix to be used for naming output dx files, if None, then prefix set at the time of initializing GIST
            object is used.
        """

        if prefix == None:
            prefix = self.prefix
        print("Generating dx files ...")
        gist_header = "voxel x y z nwat gO gH dTStrans-dens dTStrans-norm dTSorient-dens dTSorient-norm dTSsix-dens dTSsix-norm Esw-dens Esw-norm Eww-dens Eww-norm Eww-nbr-dens Eww-nbr-norm neighbor-dens neighbor-norm fHB-dens fHB-norm Nhbww-dens Nhbww-norm Nhbsw-dens Nhbsw-norm Ndonsw-dens Ndonsw-norm Naccsw-dens Naccsw-norm Ndonww-dens Ndonww-norm Naccww-dens Naccww-norm\n"
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
        dx_header += 'object 3 class array type double rank 0 items %d data follows\n' % (
            self.grid.shape[0] * self.grid.shape[1] * self.grid.shape[2])
        dx_file_objects = []

        data_keys = gist_header.strip("\n").split()

        for data_field, title in enumerate(data_keys):
            if data_field > 4 and data_field % 2 == 1 and title != "gH":
                f = open(prefix + "_" + title + ".dx", 'w')
                f.write(dx_header)
                dx_file_objects.append(f)
            else:
                dx_file_objects.append(None)

        for k in range(1, len(self.voxeldata) + 1):
            # print "writing data for voxel: ", k
            for column_i in range(5, len(data_keys), 2):
                dx_file_objects[column_i].write(
                    "%g " % (self.voxeldata[k - 1][column_i]))
                if k % 3 == 0:
                    dx_file_objects[column_i].write("\n")

        for f in dx_file_objects:
            if f is not None:
                f.close()

    def print_system_summary(self):
        """
        Print a summary of the molecular system.
        """

        print("System information:")
        print(("\tParameter file: %s\n" % self.topology_file))
        print(("\tTrajectory: %s\n" % self.trajectory))
        print(("\tFrames: %d, Total Atoms: %d, Waters: %d, Solute Atoms: %d\n" \
              % (self.num_frames, self.all_atom_ids.shape[0], self.wat_oxygen_atom_ids.shape[0],
                 self.non_water_atom_ids.shape[0])))
        # print "\tWater Model: %s\n"
        print("Grid information:")
        print(("\tGIST grid center: %5.3f %5.3f %5.3f\n" % (self.center[0], self.center[1], self.center[2])))
        print(("\tGIST grid dimensions: %i %i %i\n" % (self.dims[0], self.dims[1], self.dims[2])))
        print(("\tGIST grid spacing: %5.3f A^3\n" % (self.spacing[0])))

    def print_calcs_summary(self, num_frames=None):
        """
        Print a summary of GIST calculations

        Parameters
        ----------
        num_frames : int, option
            Number of frames processed during the simulation
        """

        if num_frames is None:
            num_frames = self.num_frames
        print("Summary of main calculations:")
        nwat_grid = 0.0
        Eswtot = 0.0
        Ewwtot = 0.0
        for k in self.voxeldata:
            if k[4] > 1.0:
                nwat_grid += k[4] / (num_frames * self.voxel_vol)
                # print k[11]
                Eswtot += k[13]
                Ewwtot += k[15]

        nwat_grid *= self.voxel_vol
        Eswtot *= self.voxel_vol * 2.0
        Ewwtot *= self.voxel_vol
        print(("Number of frames processed: %d" % num_frames))
        print(("\tAverage number of water molecules over the grid: %d" % nwat_grid))
        print(("\tTotal Solute-Water Energy over the grid: %.6f" % Eswtot ))
        print(("\tTotal Water-Water Energy over the grid: %.6f" % Ewwtot))
