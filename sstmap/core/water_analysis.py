"""
This module implements base class and utility classes for wateranalysistools package. Base class WaterAnalayisBase is inherited by main classess implementing solvation analysis algorithms and
provides methods for reading
"""

# import statements
import sys
import time
import os
from functools import wraps

import numpy as np
import scipy
import mdtraj as md
import parmed as pmd

DON_ACC_LIST = ["oxygen", "nitrogen", "sulfur"]

def function_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print ("Total time running %s: %2.2f seconds" %
               (function.func_name, t1-t0))
        return result
    return function_timer

class WaterAnalysis(object):
    """docstring for WaterAnalysis"""

    def __init__(self, topology_file, trajectory, start_frame=0,
                 num_frames=0, desmond_helper_file=None):
        self.topology_file = topology_file
        self.trajectory = trajectory
        first_frame = md.load_frame(self.trajectory, 0, top=self.topology_file)
        self.topology = first_frame.topology
        self.start_frame = start_frame
        self.num_frames = num_frames
        self.box_type = "Unspecified"
        orthogonal = np.allclose(md.load_frame(self.trajectory, 0, top=self.topology_file).unitcell_angles, 90)
        if orthogonal:
            self.box_type = "Orthorhombic"
        else:
            sys.exit("Only orthorhombic periodic boxes are currently supported.")
        
        if desmond_helper_file is not None:
            self.desmond_helper_file = desmond_helper_file
        self.rho_bulk = 0.0334
        self._generate_indices()
        self._generate_non_bonded_params()

    def _generate_indices(self):
        """
        Returns index arrays.
        """
        # use mdtraj to load in first frame to obtain indices
        # FIXME: add IO error handling for both topology and trajectory
        # reading.
        self.all_atom_ids = self.topology.select("all")
        self.wat_atom_ids = self.topology.select("water")
        self.wat_oxygen_atom_ids = self.topology.select("water and name O")
        # NOTE: The name is misleading, this is all solute atoms, including
        # ions in the system
        self.non_water_atom_ids = self.topology.select("not water")
        ######################################################
        # Not sure what this does
        wat_O_sites = []
        for i, at in enumerate(self.wat_oxygen_atom_ids):
            if i < len(self.wat_oxygen_atom_ids) - 1:
                # print i, self.wat_oxygen_atom_ids[i],
                # self.wat_oxygen_atom_ids[i + 1] - self.wat_oxygen_atom_ids[i]
                wat_O_sites.append(
                    self.wat_oxygen_atom_ids[i + 1] - self.wat_oxygen_atom_ids[i])
            else:
                # print i, self.wat_oxygen_atom_ids[i], self.wat_atom_ids[-1] -
                # self.wat_oxygen_atom_ids[i] + 1
                wat_O_sites.append(
                    self.wat_atom_ids[-1] - self.wat_oxygen_atom_ids[i] + 1)
        self.wat_O_sites = np.asarray(wat_O_sites)
        #######################################################
        # FIXME: is non water atom id just the protein atoms or everything else
        # apart from water
        # Obtain H-bond typing info
        acc_list = []
        don_list = []
        acc_don_list = []
        # To speed up calculations, we will pre-generate donor atom, hydrogen
        # pairs
        self.don_H_pair_dict = {}
        # obtain a list of non-water bonds
        non_water_bonds = [(bond[0].index, bond[1].index)
                           for bond in self.topology.bonds if bond[0].residue.name != "HOH"]
        # iterate over solute atom ids
        for at in self.non_water_atom_ids:
            # obtain bonds associated with donors or acceptors
            if self.topology.atom(at).element.name in DON_ACC_LIST:
                bonds_of_at = [bond for bond in non_water_bonds if at in bond]
                if self.topology.atom(at).element.name == "nitrogen":
                    # if a nitrogen atom is bonded to a hydrogn atom, save donor-H pair and added to donors
                    # print at, bonds_of_at
                    don_h_pairs = []
                    for at1, at2 in bonds_of_at:
                        if self.topology.atom(at2).element.name == "hydrogen":
                            don_h_pairs.append([at1, at2])
                        if self.topology.atom(at1).element.name == "hydrogen":
                            don_h_pairs.append([at2, at1])
                    if len(
                            don_h_pairs) != 0 and at not in self.don_H_pair_dict.keys():
                        don_list.append(at)
                        self.don_H_pair_dict[at] = don_h_pairs
                    # if no bonds with hydrogen found, add to acceptors
                    if len(don_h_pairs) == 0:
                        acc_list.append(at)
                if self.topology.atom(at).element.name in ["oxygen", "sulfur"]:
                    # if an oxygen or a sulfur atom is bonded to a hydrogen,
                    # add to the list of acceptor-donors and save donor-H pair
                    don_h_pairs = []
                    for at1, at2 in bonds_of_at:
                        if self.topology.atom(at2).element.name == "hydrogen":
                            don_h_pairs.append([at1, at2])
                        if self.topology.atom(at1).element.name == "hydrogen":
                            don_h_pairs.append([at2, at1])
                    if len(
                            don_h_pairs) != 0 and at not in self.don_H_pair_dict.keys():
                        acc_don_list.append(at)
                        self.don_H_pair_dict[at] = don_h_pairs
                    # if no bonds with hydrogen found, add to acceptors
                    if len(don_h_pairs) == 0:
                        acc_list.append(at)

        self.solute_acc_ids = np.array(acc_list, dtype=np.int)
        self.solute_acc_don_ids = np.array(acc_don_list, dtype=np.int)
        self.solute_don_ids = np.array(don_list, dtype=np.int)
        self.prot_hb_types = np.zeros(
            len(self.non_water_atom_ids), dtype=np.int_)
        for at_id in self.solute_acc_ids:
            self.prot_hb_types[at_id] = 1
        for at_id in self.solute_don_ids:
            self.prot_hb_types[at_id] = 2
        for at_id in self.solute_acc_don_ids:
            self.prot_hb_types[at_id] = 3

    def _generate_non_bonded_params(self):
        """
        Retunrs vdw params and charges for each atom in the list.
        """

        # use parmed to get parameters
        # the way parameters are provided by parmed is essentially
        parmed_topology_object = pmd.load_file(self.topology_file)
        vdw = []
        chg = []
        for at in self.all_atom_ids:
            # print at, self.param_class.atoms[at].charge*18.2223,
            # self.param_class.atoms[at].sigma,
            # self.param_class.atoms[at].epsilon
            vdw.append([parmed_topology_object.atoms[at].sigma,
                        parmed_topology_object.atoms[at].epsilon])
            chg.append(parmed_topology_object.atoms[at].charge)
        # FIXME: Is this multiplication step applicable to all topologies
        self.vdw = np.asarray(vdw)
        self.chg = np.asarray(chg) * 18.2223
        self.water_sites = self.wat_oxygen_atom_ids[1] - self.wat_oxygen_atom_ids[0]
        water_sig = self.vdw[self.wat_oxygen_atom_ids[0]][0]
        water_eps = self.vdw[self.wat_oxygen_atom_ids[0]][1]        
        self.water_water_acoeff = 4*water_eps*(water_sig**12)
        self.water_water_bcoeff = -4*water_eps*(water_sig**6)

        solute_water_sig = 0.5*(water_sig + self.vdw[self.non_water_atom_ids, 0])
        solute_water_eps = np.sqrt(water_eps*self.vdw[self.non_water_atom_ids, 1])
        self.solute_water_acoeff = 4*solute_water_eps*(solute_water_sig**12)
        self.solute_water_bcoeff = -4*solute_water_eps*(solute_water_sig**6)

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
        return energy_lj*0.5, energy_elec*0.5


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
        angles = md.utils.in_units_of(md.compute_angles(traj, angle_triplets), "radians", "degrees")
        angles_ww = angles[0, 0:water_nbrs.shape[0]*4]
        angles_sw = angles[0, water_nbrs.shape[0]*4:]
        angle_triplets_ww = angle_triplets[:water_nbrs.shape[0]*4]
        angle_triplets_sw = angle_triplets[water_nbrs.shape[0]*4:]
        hbonds_ww = angle_triplets_ww[np.where(angles_ww <= 30.0)]
        hbonds_sw = angle_triplets_sw[np.where(angles_sw <= 30.0)]
        return (hbonds_ww, hbonds_sw)




##########################################################################
# Class and methods for 'efficient' neighbor search                                                             #
##########################################################################

class NeighborSearch(object):

    def __init__(self, xyz, dist):
        """
        Class for fast queries of coordinates that are within distance <dist>
        of specified coordinate. This class must first be initialized from an
        array of all available coordinates, and a distance threshold. The
        query() method can then be used to get a list of points that are within
        the threshold distance from the specified point.
        """
        # create an array of indices around a cubic grid
        self.neighbors = []
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                for k in (-1, 0, 1):
                    self.neighbors.append((i, j, k))
        self.neighbor_array = np.array(self.neighbors, np.int)

        self.min_ = np.min(xyz, axis=0)
        self.cell_size = np.array([dist, dist, dist], np.float)
        cell = np.array((xyz - self.min_) / self.cell_size)  # , dtype=np.int)
        # create a dictionary with keys corresponding to integer representation
        # of transformed XYZ's
        self.cells = {}
        for ix, assignment in enumerate(cell):
            # convert transformed xyz coord into integer index (so coords like
            # 1.1 or 1.9 will go to 1)
            indices = assignment.astype(int)
            # create interger indices
            t = tuple(indices)

            # NOTE: a single index can have multiple coords associated with it
            # if this integer index is already present
            if t in self.cells:
                # obtain its value (which is a list, see below)
                xyz_list, trans_coords, ix_list = self.cells[t]
                # append new xyz to xyz list associated with this entry
                xyz_list.append(xyz[ix])
                # append new transformed xyz to transformed xyz list associated
                # with this entry
                trans_coords.append(assignment)
                # append new array index
                ix_list.append(ix)
            # if this integer index is encountered for the first time
            else:
                # create a dictionary key value pair,
                # key: integer index
                # value: [[list of x,y,z], [list of transformed x,y,z], [list
                # of array indices]]
                self.cells[t] = ([xyz[ix]], [assignment], [ix])

        self.dist_squared = dist * dist

    def query_nbrs_single_point(self, point):
        """
        Given a coordinate point, return all point indexes (0-indexed) that
        are within the threshold distance from it.
        """
        cell0 = np.array((point - self.min_) / self.cell_size, dtype=np.int)
        tuple0 = tuple(cell0)
        near = []
        for index_array in tuple0 + self.neighbor_array:
            t = tuple(index_array)
            if t in self.cells:
                xyz_list, trans_xyz_list, ix_list = self.cells[t]
                for (xyz, ix) in zip(xyz_list, ix_list):
                    diff = xyz - point
                    if np.dot(diff, diff) <= self.dist_squared and float(
                            np.dot(diff, diff)) > 0.0:
                        # near.append(ix)
                        # print ix, np.dot(diff, diff)
                        near.append(ix)
        return near

    def query_point_and_distance(self, point):
        """
        Given a coordinate point, return all point indexes (0-indexed) that
        are within the threshold distance from it.
        """
        cell0 = np.array((point - self.min_) / self.cell_size, dtype=np.int)
        tuple0 = tuple(cell0)
        near = []
        for index_array in tuple0 + self.neighbor_array:
            t = tuple(index_array)
            if t in self.cells:
                xyz_list, trans_xyz_list, ix_list = self.cells[t]
                for (xyz, ix) in zip(xyz_list, ix_list):
                    diff = xyz - point
                    if np.dot(diff, diff) <= self.dist_squared and float(
                            np.dot(diff, diff)) > 0.0:
                        # near.append(ix)
                        # print ix, np.dot(diff, diff)
                        near.append((ix, np.sqrt(np.dot(diff, diff))))
        return near

    def query_nbrs_multiple_points(self, points):
        """
        Given a coordinate point, return all point indexes (0-indexed) that
        are within the threshold distance from it.
        shape of points has to be (n_lig_atoms, 3)
        """
        near = []
        for point in points:
            cell0 = np.array(
                (point - self.min_) / self.cell_size,
                dtype=np.int)
            tuple0 = tuple(cell0)

            for index_array in tuple0 + self.neighbor_array:
                t = tuple(index_array)
                if t in self.cells:
                    xyz_list, trans_xyz_list, ix_list = self.cells[t]
                    for (xyz, ix) in zip(xyz_list, ix_list):
                        diff = xyz - point
                        if np.dot(diff, diff) <= self.dist_squared and float(
                                np.dot(diff, diff)) > 0.0:
                            # near.append(ix)
                            # print ix, np.dot(diff, diff)
                            if ix not in near:
                                near.append(ix)
        return near

#*********************************************************************************************#
