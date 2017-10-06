##############################################################################
# SSTMap: A Python library for the calculation of water structure and 
#         thermodynamics on solute surfaces from molecular dynamics 
#         trajectories.
# Copyright 2016-2017 Lehman College City University of New York and the Authors
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
# License along with SSTMap. If not, see <http://www.gnu.org/licenses/>.
##############################################################################
"""
This module contains implementation of a parent class for water analysis in 
molecular dynamics simulation trajectories. The data and methods in this class
are common to site-based and grid-based analysis classes.

Please reference the following if you use this code in your research:
[1] Haider K, Wickstrom L, Ramsey S, Gilson MK and Kurtzman T. Enthalpic Breakdown 
of Water Structure on Protein Active Site Surfaces. J Phys Chem B. 120:8743-8756, 
2016. http://dx.doi.org/10.1021/acs.jpcb.6b01094.
"""

##############################################################################
# Imports
##############################################################################
import sys
import os

import numpy as np
import mdtraj as md
import parmed as pmd
from parmed.charmm import CharmmParameterSet
##############################################################################
# Globals
##############################################################################

DON_ACC_LIST = ["oxygen", "nitrogen", "sulfur"]
_WATER_RESNAMES = ['H2O', 'HHO', 'OHH', 'HOH',  'OH2', 'SOL', 'WAT', 'TIP', 'TIP2', 'TIP3', 'TIP4', 'T3P', 'T4P', 'T5P']
#_WAT_NBR_AVG = ["T3P": 5.25]
ANGLE_CUTOFF_RAD = 0.523599
requirements = {   
    "prmtop": ["prmtop", "", "lorentz-bertholot"],
    "parm7": ["parm7", "", "lorentz-bertholot"],
    "psf": ["toppar", "Please provide a folder named toppar that contains charmm parameter/topology files.", "lorentz-bertholot"],
    "gro": ["top", "Please provide graomcs .top file corresponding to your system and also make sure that .itp files are present in the directory where calculations are being run.", "lorentz-bertholot"],
    "pdb": ["txt", "Please provide a text file containing non-bonded parameters for your system.", "geometric"],
    }


##############################################################################
# WaterAnalysis class definition
##############################################################################

class WaterAnalysis(object):
    """Parent class for setting up water analysis calculations in molecular 
    dynamics trajectories.
    """

    def __init__(self, topology_file, trajectory, start_frame=0,
                 num_frames=0, supporting_file=None, rho_bulk=None):
        """Initialize WaterAnalysis object for a trajectory and
        corresponding topology file.
        
        Parameters
        ----------
        topology_file : string
            Filename of the system topology file.
        trajectory : string
            Filename of the molecular dynamics trajectory.
        start_frame : int, optional
            The frame index from which the calculations will begin. Default: 0
        num_frames : int, optional
            The total number of frames or the length of simulation over which 
            calculations will be performed. Default: 0
        supporting_file : None, optional
            Filename of additional file containing non-bonded parameters for
            every particle in the system. Default: None
        rho_bulk : float
            Reference bulk water density to be used in calculations. Default: None
        """

        self.topology_file = topology_file
        self.trajectory = trajectory
        self.supporting_file = supporting_file
        self.start_frame = start_frame
        #assert num_frames >= 100, "A minimum of 100 frames are required for analysis."
        self.num_frames = num_frames
        self.check_topology_requiremnts(self.topology_file, self.supporting_file)
        first_frame = md.load_frame(self.trajectory, self.start_frame, top=self.topology_file)
        assert first_frame.unitcell_lengths is not None, "Could not detect unit cell information."
        self.topology = first_frame.topology
        self.box_type = "Unspecified"
        orthogonal = False
        try:
            orthogonal = np.allclose(md.load_frame(self.trajectory, 0, top=self.topology_file).unitcell_angles, 90)
            if orthogonal:
                self.box_type = "Orthorhombic"
        except Exception as e:
            print("WARNING: Only orthorhombic periodic boxes are currently supported.")
        self.rho_bulk = rho_bulk
        if self.rho_bulk is None:
            self.rho_bulk = 0.0334
        super_wat_select_exp = ""
        for i, wat_res in enumerate(_WATER_RESNAMES):
            if i < len(_WATER_RESNAMES) - 1:
                super_wat_select_exp += "resname %s or " % wat_res
            else:
                super_wat_select_exp += "resname %s" % wat_res        
        self.all_atom_ids = self.topology.select("all")
        self.wat_atom_ids = self.topology.select("water")
        self.prot_atom_ids = self.topology.select("protein")
        if self.wat_atom_ids.shape[0] == 0:
            self.wat_atom_ids = self.topology.select(super_wat_select_exp)
        assert (self.wat_atom_ids.shape[0] != 0), "Unable to recognize waters in the system!"
        assert (self.topology.atom(self.wat_atom_ids[0]).name == "O"), "Failed while constructing water oxygen atom indices!"
        self.wat_oxygen_atom_ids = np.asarray([atom for atom in self.wat_atom_ids if self.topology.atom(atom).name == "O"])
        self.non_water_atom_ids = np.setdiff1d(self.all_atom_ids, self.wat_atom_ids)
        assert (self.wat_atom_ids.shape[0] + self.non_water_atom_ids.shape[0] == self.all_atom_ids.shape[0]), "Failed to partition atom indices in the system correctly!"

    def check_topology_requiremnts(self, top_file, support_file):
        """
        Performs a check on supplied topology and supporting files to determine if the
        required files for this topology format are available for calculations and if
        checks are successful assigns combination rule corresponding to the format.

        Parameters
        ----------
        """

        topology_extension = top_file.split(".")[-1]

        if topology_extension not in requirements.keys():
            message = """SSTMap currently does not support %s topology file type.
            If this is a non-standard forcefiled, consider using a PDB file as a topplogy
            and provide a text file containing non-bonded parameters for each atom in your system.
            See sstmap.org for more details.
            """ % topology_extension
            sys.exit(message)

        required_support = requirements[topology_extension][0]
        if required_support == topology_extension:
            self.supporting_file = self.topology_file

        else:
            if support_file is None:
                message_1 = """SSTMap requires %s as a supporting file/data for %s parameter format.
                Please provide it as an argument to supporting_file argument or if you are running
                run_gist or run_hsa, provide it as an argument to -p flag.
                """ % (required_support, topology_extension)
                message_2 = """More specifically, %s""" % requirements[topology_extension][1]
                sys.exit(message_1 + message_2)
        self.comb_rule = requirements[topology_extension][-1]

    def assign_hb_types(self):
        """Generates index arrays for atoms of different types in the system, assign
        a hydrogen-bond type to each atom and generate a dictionary of H-bond donors
        where indices of each connected hydrogen are stored for each donor.

        Notes
        -----
        Several np.ndarray objects are generated and assigned as attributes 
        to WaterAnalysis object.
        """

        # Obtain H-bond typing info
        self.topology.create_standard_bonds()
        acc_list = []
        don_list = []
        acc_don_list = []
        # To speed up calculations, we will pre-generate donor atom, hydrogen
        # pairs
        
        self.don_H_pair_dict = {}
        #self.new_dict = {}
        # obtain a list of non-water bonds
        non_water_bonds = [(bond[0].index, bond[1].index)
                           for bond in self.topology.bonds if bond[0].residue.name not in _WATER_RESNAMES]
        dist_pairs = []
        keys_all = []

        for at in self.prot_atom_ids:
            # obtain bonds associated with donors or acceptors
            if self.topology.atom(at).element.name in DON_ACC_LIST:
                bonds_of_at = []
                for bond in non_water_bonds:
                    if at in bond and bond not in bonds_of_at:
                        bonds_of_at.append(bond)

                if self.topology.atom(at).element.name == "nitrogen":
                    # if a nitrogen atom is bonded to a hydrogn atom, save donor-H pair and added to donors
                    # print at, bonds_of_at
                    don_h_pairs = []
                    for at1, at2 in bonds_of_at:
                        if self.topology.atom(at2).element.name == "hydrogen":
                            don_h_pairs.append([at1, at2])
                        if self.topology.atom(at1).element.name == "hydrogen":
                            don_h_pairs.append([at2, at1])
                    #if len(don_h_pairs) != 0 and at not in list(self.don_H_pair_dict.keys()):
                    #    don_list.append(at)
                    #    self.don_H_pair_dict[at] = don_h_pairs
                    if len(don_h_pairs) != 0:
                        keys_all.append(at)
                        for bond in don_h_pairs:
                            dist_pairs.append(bond)
                        if at not in don_list:
                            don_list.append(at)
                    # if no bonds with hydrogen found, add to acceptors
                    else:
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
                    #if len(don_h_pairs) != 0 and at not in list(self.don_H_pair_dict.keys()):
                    #    acc_don_list.append(at)
                    #    self.don_H_pair_dict[at] = don_h_pairs
                    if len(don_h_pairs) != 0:
                        keys_all.append(at)
                        for bond in don_h_pairs:
                            dist_pairs.append(bond)
                        if at not in acc_don_list:
                            acc_don_list.append(at)
                    # if no bonds with hydrogen found, add to acceptors
                    else:
                        acc_list.append(at)

        for index, pair in enumerate(dist_pairs):
            if pair[0] not in self.don_H_pair_dict.keys():
                self.don_H_pair_dict[pair[0]] = [[pair[0], pair[1]]]
            else:
                self.don_H_pair_dict[pair[0]].append([pair[0], pair[1]])

        self.solute_acc_ids = np.array(acc_list, dtype=np.int)
        self.solute_acc_don_ids = np.array(acc_don_list, dtype=np.int)
        self.solute_don_ids = np.array(don_list, dtype=np.int)
        self.prot_hb_types = np.zeros(len(self.non_water_atom_ids), dtype=np.int_)

        for at_id in self.solute_acc_ids:
            self.prot_hb_types[at_id] = 1
        for at_id in self.solute_don_ids:
            self.prot_hb_types[at_id] = 2
        for at_id in self.solute_acc_don_ids:
            self.prot_hb_types[at_id] = 3

    def generate_nonbonded_params(self):
        """
        Obtains non-bonded parameters for each atom in the system.
        """

        # use parmed to get parameters
        self.water_sites = self.wat_oxygen_atom_ids[1] - self.wat_oxygen_atom_ids[0]

        vdw = []
        chg = []
        if self.supporting_file.endswith(".txt"):
            nb_data = np.loadtxt(self.supporting_file)
            for c in nb_data[:, 0]:
                chg.append(c)
            for v in nb_data[:, 1:]:
                vdw.append(v)
            chg = np.asarray(chg)

        elif self.topology_file.endswith(".psf"):
            parmed_topology_object = pmd.load_file(self.topology_file)
            param_dir = os.path.abspath(self.supporting_file)
            param_files = [os.path.join(param_dir, f) for f in os.listdir(param_dir) 
                if os.path.isfile(os.path.join(param_dir, f)) and f.endswith((".rtf", ".top", ".par", ".prm", ".inp", ".str"))]
            params = CharmmParameterSet(*param_files)
            try:
                parmed_topology_object.load_parameters(params)
            except Exception as e:
                print(e)
            for at in self.all_atom_ids:
                vdw.append([parmed_topology_object.atoms[at].sigma,
                            parmed_topology_object.atoms[at].epsilon])
                chg.append(parmed_topology_object.atoms[at].charge)

        # .prmtop, .gro, .
        else:
            parmed_topology_object = pmd.load_file(self.supporting_file)
            for at in self.all_atom_ids:
                vdw.append([parmed_topology_object.atoms[at].sigma,
                            parmed_topology_object.atoms[at].epsilon])
                chg.append(parmed_topology_object.atoms[at].charge)

        # User provided charges are supposed to be in correct units.
        if not self.supporting_file.endswith(".txt"):
            chg = np.asarray(chg) * 18.2223
        vdw = np.asarray(vdw)
        water_chg = chg[self.wat_atom_ids[0:self.water_sites]].reshape(self.water_sites, 1)
        self.chg_product = water_chg * np.tile(chg[self.all_atom_ids], (self.water_sites, 1))
        water_sig = vdw[self.wat_atom_ids[0:self.water_sites], 0].reshape(self.water_sites, 1)
        water_eps = vdw[self.wat_atom_ids[0:self.water_sites], 1].reshape(self.water_sites, 1)
        self.acoeff, self.bcoeff = self.apply_combination_rules(water_sig, water_eps, vdw, self.comb_rule)
        """
        # for debuging
        water_sig = vdw[self.wat_oxygen_atom_ids[0]][0]
        water_eps = vdw[self.wat_oxygen_atom_ids[0]][1]
        self.water_water_acoeff = 4 * water_eps * (water_sig ** 12)
        self.water_water_bcoeff = -4 * water_eps * (water_sig ** 6)
        solute_water_sig, solute_water_eps = self.apply_combination_rules_old(water_sig, water_eps, vdw, self.comb_rule)
        self.solute_water_acoeff = 4 * solute_water_eps * (solute_water_sig ** 12)
        self.solute_water_bcoeff = -4 * solute_water_eps * (solute_water_sig ** 6)
        """


    def apply_combination_rules(self, water_sig, water_eps, vdw, rule=None):
        """

        Args:
            water_sig: 
            water_eps: 
            vdw: 
            rule: 

        :return:
        """
        mixed_sig, mixed_eps = None, None
        acoeff, bcoeff = None, None
        if rule is None or rule == "lorentz-bertholot":
            mixed_sig = 0.5*(water_sig + vdw[self.all_atom_ids, 0])
            mixed_eps = np.sqrt(water_eps * vdw[self.all_atom_ids, 1])
        if rule == "geometric":
            mixed_sig = np.sqrt(water_sig * vdw[self.all_atom_ids, 0])
            mixed_eps = np.sqrt(water_eps * vdw[self.all_atom_ids, 1])
        acoeff = 4 * mixed_eps * (mixed_sig**12)
        bcoeff = 4 * mixed_eps * (mixed_sig**6)
        return acoeff, bcoeff

    def apply_combination_rules_old(self, water_sig, water_eps, vdw, rule=None):
        """

        Args:
            water_sig:
            water_eps:
            vdw:
            rule:

        Returns:

        """

        solute_water_sig, solute_water_eps = None, None
        if rule is None or rule == "lorentz-bertholot":
            solute_water_sig = 0.5*(water_sig + vdw[self.non_water_atom_ids, 0])
            solute_water_eps = np.sqrt(water_eps * vdw[self.non_water_atom_ids, 1])
        if rule == "geometric":
            solute_water_sig = np.sqrt(water_sig * vdw[self.non_water_atom_ids, 0])
            solute_water_eps = np.sqrt(water_eps * vdw[self.non_water_atom_ids, 1])

        return solute_water_sig, solute_water_eps

    def calculate_energy(self, distance_matrix):
        """Calculates total interaction energy of a water molecule with the rest of the
        system from the distance matrix and non-bonded parameter attributes of the WaterAnalysis object.

        
        Parameters
        ----------
        distance_matrix : np.ndarray, float, shape=(K, N)
            A matrix of inter-atomic distance, where K is the number of partciles in
            the water molecule and N is the total number of atoms in the system

        Returns
        -------
        energy_lj : np.ndarray, float, shape=(1, N_lj)
            Array of lennard-Jones interaction energies of the water oxygen against all solute
            particles and water oxygen atoms (N_lj) 
        energy_elec : np.ndarray, float, shape=(K, N)
            Array of electrostatic interaction energies of the water molecule against all atoms
            in the system.        
        """
        with np.errstate(invalid='ignore', divide='ignore'):
            wat_wat_dist_6 = distance_matrix[0, :][self.wat_oxygen_atom_ids] ** -6
            wat_solute_dist_6 = distance_matrix[0, :][self.non_water_atom_ids] ** -6
            wat_wat_dist_12 = distance_matrix[0, :][self.wat_oxygen_atom_ids] ** -12
            wat_solute_dist_12 = distance_matrix[0, :][self.non_water_atom_ids] ** -12
            #energy_sw_lj = (self.solute_water_acoeff*np.power(wat_solute_dist, -12)) + (self.solute_water_bcoeff*np.power(wat_solute_dist, -6))
            #energy_ww_lj = (self.water_water_acoeff*np.power(wat_wat_dist, -12)) + (self.water_water_bcoeff*np.power(wat_wat_dist, -6))
            energy_sw_lj = (self.solute_water_acoeff * wat_solute_dist_12) + (self.solute_water_bcoeff*wat_solute_dist_6)
            #energy_ww_lj = (self.water_water_acoeff*wat_wat_dist_12) + (self.water_water_bcoeff*wat_wat_dist_6)
            energy_ww_lj = (self.water_water_acoeff * wat_wat_dist_12) + (self.water_water_bcoeff*wat_wat_dist_6)

            water_chg = self.chg[self.wat_atom_ids[0:self.water_sites]].reshape(self.water_sites, 1)
            energy_elec = water_chg*np.tile(self.chg[self.all_atom_ids], (self.water_sites, 1))
            energy_elec *= distance_matrix ** -1
            #for i in xrange(self.wat_atom_ids.shape[0]):
            #    print "ref elec dist between current water atoms and %i: " % (self.wat_atom_ids[i])
            #    print distance_matrix[:, i]

            #print "Solute-water LJ Energy of this water: ", np.nansum(energy_sw_lj)
            #print "Solute-water Elec Energy of this water: ", np.sum(energy_elec[:, self.non_water_atom_ids])
            #print "Water-water LJ Energy of this water: ", np.nansum(energy_ww_lj)/2.0
            #print "Water-water Elec Energy of this water: ", (np.sum(energy_elec[:, self.wat_atom_ids[0]:water_id])/2.0) + (np.sum(energy_elec[:, water_id + self.water_sites:])/2.0)
            energy_lj = np.concatenate((energy_sw_lj, energy_ww_lj), axis=0)
        return energy_lj, energy_elec

    def calculate_hydrogen_bonds(self, traj, water, water_nbrs, solute_nbrs):
        """Calculates hydrogen bonds made by a water molecule with its first shell
        water and solute neighbors.
        
        Parameters
        ----------
        traj : md.trajectory
            MDTraj trajectory object for which hydrogen bonds are to be calculates. 
        water : int
            The index of water oxygen atom
        water_nbrs : np.ndarray, int, shape=(N^{ww}_nbr, )
            Indices of the water oxygen atoms in the first solvation shell of the water molecule.
        solute_nbrs : np.ndarray, int, shape=(N^{sw}_nbr, )
            Indices of thesolute atoms in the first solvation shell of the water molecule.
        
        Returns
        -------
        (hbonds_ww, hbonds_sw) : tuple
            A tuple consisting of two np.ndarray objects for water-water and solute-water
            hydrogen bonds. A hydrogen bond is represented by an array of indices
            of three atom particpating in the hydrogen bond, [Donor, H, Acceptor]
        """
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
        angles = md.compute_angles(traj, np.asarray(angle_triplets))
        angles[np.isnan(angles)] = 0.0
        angles_ww = angles[0, 0:water_nbrs.shape[0]*4]
        angles_sw = angles[0, water_nbrs.shape[0]*4:]
        angle_triplets_ww = angle_triplets[:water_nbrs.shape[0]*4]
        angle_triplets_sw = angle_triplets[water_nbrs.shape[0]*4:]
        hbonds_ww = angle_triplets_ww[np.where(angles_ww <= ANGLE_CUTOFF_RAD)]
        hbonds_sw = angle_triplets_sw[np.where(angles_sw <= ANGLE_CUTOFF_RAD)]
        return (hbonds_ww, hbonds_sw)

    def calculate_hydrogen_bonds2(self, coords, water, water_nbrs, solute_nbrs):
        """Calculates hydrogen bonds made by a water molecule with its first shell
        water and solute neighbors.
        
        Parameters
        ----------
        traj : md.trajectory
            MDTraj trajectory object for which hydrogen bonds are to be calculates. 
        water : int
            The index of water oxygen atom
        water_nbrs : np.ndarray, int, shape=(N^{ww}_nbr, )
            Indices of the water oxygen atoms in the first solvation shell of the water molecule.
        solute_nbrs : np.ndarray, int, shape=(N^{sw}_nbr, )
            Indices of thesolute atoms in the first solvation shell of the water molecule.
        
        Returns
        -------
        (hbonds_ww, hbonds_sw) : tuple
            A tuple consisting of two np.ndarray objects for water-water and solute-water
            hydrogen bonds. A hydrogen bond is represented by an array of indices
            of three atom particpating in the hydrogen bond, [Donor, H, Acceptor]
        """
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
        angle_triplets_ww = angle_triplets[:water_nbrs.shape[0]*4]
        angle_triplets_sw = angle_triplets[water_nbrs.shape[0]*4:]
        
        angles = np.asarray([self.calc_angle(coords, triplet) for triplet in angle_triplets])
        angles_ww = angles[0:water_nbrs.shape[0]*4]
        angles_sw = angles[water_nbrs.shape[0]*4:]

        hbonds_ww = angle_triplets_ww[np.where(angles_ww <= 30.0)]
        hbonds_sw = angle_triplets_sw[np.where(angles_sw <= 30.0)]
        return hbonds_ww, hbonds_sw

    def calc_angle(self, coords, triplet):
        """
        """
        dp = np.dot(coords[0, triplet[0], :] - coords[0, triplet[1], :], coords[0, triplet[2], :] - coords[0, triplet[1], :])
        norms = [np.linalg.norm(coords[0, triplet[0], :] - coords[0, triplet[1], :]), np.linalg.norm(coords[0, triplet[2], :] - coords[0, triplet[1], :])]
        cos_angle = dp / (norms[0] * norms[1])
        return np.degrees(np.arccos(cos_angle))

