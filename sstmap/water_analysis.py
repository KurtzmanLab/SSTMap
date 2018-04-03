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
This module contains implementation of a parent class for water analysis in 
molecular dynamics simulation trajectories. The data and methods in this class
are common to site-based and grid-based analysis classes.

Please reference the following if you use this code in your research:

[1] Haider K, Cruz A, Ramsey S, Gilson M. and Kurtzman T. Solvation Structure and Thermodynamic Mapping (SSTMap): An
Open-Source, Flexible Package for the Analysis of Water in Molecular Dynamics Trajectories. J. Chem. Theory Comput.
(10.1021/acs.jctc.7b00592) 2017.
[2] Crystal N. Nguyen, Michael K. Gilson, Tom Young. Structure and Thermodynamics of Molecular Hydration via Grid
Inhomogeneous Solvation Theory. eprint arXiv:1108.4876, (2011).
[3] Crystal N. Nguyen, Tom Kurtzman Young, and Michael K. Gilson. Grid inhomogeneous solvation theory: hydration
structure and thermodynamics of the miniature receptor cucurbit[7]uril. J. Chem. Phys. 137, 044101 (2012)
[4] Haider K, Wickstrom L, Ramsey S, Gilson MK and Kurtzman T. Enthalpic Breakdown of Water Structure on Protein Active
Site Surfaces. J Phys Chem B. 120:8743-8756, (2016). http://dx.doi.org/10.1021/acs.jpcb.6b01094.
"""

##############################################################################
# Imports
##############################################################################
import numpy as np
import parmed as pmd
import mdtraj as md
from parmed.charmm import CharmmParameterSet
from sstmap.utils import *

##############################################################################
# Globals
##############################################################################

DON_ACC_LIST = ["oxygen", "nitrogen", "sulfur"]
_WATER_RESNAMES = ['H2O', 'HHO', 'OHH', 'HOH',  'OH2', 'SOL', 'WAT', 'TIP', 'TIP2', 'TIP3', 'TIP4', 'T3P', 'T4P', 'T5P']
ANGLE_CUTOFF_RAD = 0.523599
requirements = {   
    "prmtop": ["prmtop", "", "lorentz-bertholot"],
    "parm7": ["parm7", "", "lorentz-bertholot"],
    "psf": ["toppar", "Please provide a folder named toppar that contains charmm parameter/topology files.",
            "lorentz-bertholot"],
    "gro": ["top", "Please provide graomcs .top file corresponding to your system and also make sure that .itp files "
                   "are present in the directory where calculations are being run. To get a list of .itp files being "
                   "used by gromacs topology file, type $grep #include ", "lorentz-bertholot"],
    "pdb": ["txt", "Please provide a text file containing non-bonded parameters for your system.", "geometric"],
    }


##############################################################################
# WaterAnalysis class definition
##############################################################################

class WaterAnalysis(object):
    """Parent class for setting up water analysis calculations in molecular 
    dynamics trajectories.
    """

    def __init__(self, topology_file, trajectory, supporting_file=None):
        """Initialize WaterAnalysis object for a trajectory and
        corresponding topology file.
        
        Parameters
        ----------
        topology_file : string
            Filename of the system topology file.
        trajectory : string
            Filename of the molecular dynamics trajectory.
        supporting_file : None, optional
            Filename of additional file containing non-bonded parameters for
            every particle in the system. Default: None
        """
        # Check sanity checks on files
        if not os.path.exists(topology_file) or not os.path.exists(trajectory):
            raise IOError("File %s or %s does not exist." % (topology_file, trajectory))
        self.topology_file = topology_file
        self.trajectory = trajectory

        # Check if correct supporting file is provided.
        self.supporting_file = supporting_file
        topology_extension = self.topology_file.split(".")[-1]
        required_support = requirements[topology_extension][0]
        self.comb_rule = None
        if required_support == topology_extension:
            self.supporting_file = self.topology_file
            self.comb_rule = requirements[topology_extension][-1]
        else:
            if topology_extension not in list(requirements.keys()):
                message = """SSTMap currently does not support %s topology file type.
                If this is a non-standard force-filed, consider using a PDB file as a topplogy
                and provide a text file containing non-bonded parameters for each atom in your system.
                See sstmap.org for more details.
                """ % topology_extension
                sys.exit(message)
            else:
                self.supporting_file = supporting_file
                self.comb_rule = requirements[topology_extension][-1]

        # Create Parmed topology object and perform sanity check on PBC's in the trajectory
        first_frame = md.load_frame(self.trajectory, 0, top=self.topology_file)
        assert first_frame.unitcell_lengths is not None, "Could not detect unit cell information."
        self.topology = first_frame.topology

        # Create index arrays for iteration over groups of atoms and perform some sanity checks on system topology
        super_wat_select_exp = ""
        for i, wat_res in enumerate(_WATER_RESNAMES):
            if i < len(_WATER_RESNAMES) - 1:
                super_wat_select_exp += "resname %s or " % wat_res
            else:
                super_wat_select_exp += "resname %s" % wat_res        
        self.all_atom_ids = self.topology.select("all")
        self.prot_atom_ids = self.topology.select("protein")
        self.wat_atom_ids = self.topology.select("water")
        if self.wat_atom_ids.shape[0] == 0:
            self.wat_atom_ids = self.topology.select(super_wat_select_exp)
        assert (self.wat_atom_ids.shape[0] != 0), \
            "Unable to recognize water residues in the system!"
        assert (self.topology.atom(self.wat_atom_ids[0]).name == "O"), \
            "Failed while constructing water oxygen atom indices!"
        self.wat_oxygen_atom_ids = np.asarray([atom for atom in self.wat_atom_ids if self.topology.atom(atom).name == "O"])
        self.water_sites = self.wat_oxygen_atom_ids[1] - self.wat_oxygen_atom_ids[0]
        for i in self.wat_oxygen_atom_ids:
            O, H1, H2 = self.topology.atom(i).name[0], self.topology.atom(i + 1).name[0], self.topology.atom(i + 2).name[0]
            if O != "O" or H1 != "H" or H2 != "H":
                sys.exit("Water molecules in the topology must be organized as Oxygen, Hydrogen, Hydrogen, Virtual-sites.")
        self.non_water_atom_ids = np.setdiff1d(self.all_atom_ids, self.wat_atom_ids)
        # ions or ligands
        self.non_prot_atom_ids = np.setdiff1d(self.non_water_atom_ids, self.prot_atom_ids)
        # if no protein, then set other solute to protein index variable for energy calculation purposes
        if self.prot_atom_ids.shape[0] == 0:
            self.prot_atom_ids = self.non_water_atom_ids
        assert (self.wat_atom_ids.shape[0] + self.non_water_atom_ids.shape[0] == self.all_atom_ids.shape[0]), \
            "Failed to partition atom indices in the system correctly!"

        # Obtain non-bonded parameters for the system
        print("Obtaining non-bonded parameters for the system ...")
        self.chg_product, self.acoeff, self.bcoeff = self.generate_nonbonded_params()
        assert self.chg_product.shape == self.acoeff.shape == self.bcoeff.shape, \
            "Mismatch in non-bonded parameter matrices, exiting."
        print("Done.")

        # Assign a hydrogen bond to atoms
        print("Assigning hydrogen bond types ...")
        self.don_H_pair_dict = {}
        self.prot_hb_types = np.zeros(len(self.all_atom_ids), dtype=np.int_)
        self.solute_acc_ids, self.solute_don_ids, self.solute_acc_don_ids = self.assign_hb_types()
        print("Done.")

    @function_timer
    def assign_hb_types(self):
        """Assigns a hydrogen-bond type to each atom and updates a dictionary of H-bond donors
        whose keys are donor atom ids and values are the indices of each connected hydrogen are stored for each donor.

        Returns
        -------
        solute_acc_ids : numpy.ndarray
            An array of indices corresponding to solute acceptor atoms
        solute_don_ids : numpy.ndarray
            An array of indices corresponding to solute donor atoms
        solute_acc_don_ids : numpy.ndarray
            An array of indices corresponding to solute atoms that are both acceptors and donors

        Notes
        -----
        The following attributes of the object are also updated.
        self.don_H_pair_dict:
            This dictionary is populated with keys that are indices of atoms in solute_don_ids and in
            solute_acc_don_ids. The value for each key is a list of tuples where each tuple is the pair
            of atom indices, first index is the donor atom and second index is the covalently-bonded hydrogen atom.
        self.prot_hb_types:
            Array size equal number of solute atoms and each value is the numeric h-bond type, using the
            following scheme; 0=non_hb, 1=acceptor, 2=donor, 3=both.
        """

        self.topology.create_standard_bonds()
        acc_list = []
        don_list = []
        acc_don_list = []
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
            if pair[0] not in list(self.don_H_pair_dict.keys()):
                self.don_H_pair_dict[pair[0]] = [[pair[0], pair[1]]]
            else:
                self.don_H_pair_dict[pair[0]].append([pair[0], pair[1]])

        solute_acc_ids = np.array(acc_list, dtype=np.int)
        solute_acc_don_ids = np.array(acc_don_list, dtype=np.int)
        solute_don_ids = np.array(don_list, dtype=np.int)

        for at_id in solute_acc_ids:
            self.prot_hb_types[at_id] = 1
        for at_id in solute_don_ids:
            self.prot_hb_types[at_id] = 2
        for at_id in solute_acc_don_ids:
            self.prot_hb_types[at_id] = 3

        return solute_acc_ids, solute_don_ids, solute_acc_don_ids

    @function_timer
    def generate_nonbonded_params(self):
        """
        Generates non-bonded parameters for energy calculations.

        Returns
        -------
        chg_product : numpy.ndarray
            An N_sites x N_particles matrix where N_sites is the number of sites in the water model and N_atoms
            is the total number of particles in the system, each entry of the matrix is the product of the charges
            q_i*q_j used for the calculation of electrostatic interactions.
        acoeff : numpy.ndarray
            An N_sites x N_particles matrix where N_sites is the number of sites in the water model and N_atoms
            is the total number of particles in the system, each entry of the matrix is the A coefficient in the
            AB form of Lennard Jones potential.
        bcoeff  : numpy.ndarray
            An N_sites x N_particles matrix where N_sites is the number of sites in the water model and N_atoms
            is the total number of particles in the system, each entry of the matrix is the B coefficient in the
            AB form of Lennard Jones potential.
        """

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
        else:
            parmed_topology_object = pmd.load_file(self.supporting_file)
            for at in self.all_atom_ids:
                vdw.append([parmed_topology_object.atoms[at].sigma,
                            parmed_topology_object.atoms[at].epsilon])
                chg.append(parmed_topology_object.atoms[at].charge)

        # User provided charges are assumed to be in correct units.
        # TODO: Make the units for charges explicit in docstring
        if not self.supporting_file.endswith(".txt"):
            chg = np.asarray(chg) * 18.2223
        vdw = np.asarray(vdw)
        water_chg = chg[self.wat_atom_ids[0:self.water_sites]].reshape(self.water_sites, 1)
        chg_product = water_chg * np.tile(chg[self.all_atom_ids], (self.water_sites, 1))

        water_sig = vdw[self.wat_atom_ids[0:self.water_sites], 0].reshape(self.water_sites, 1)
        water_eps = vdw[self.wat_atom_ids[0:self.water_sites], 1].reshape(self.water_sites, 1)
        mixed_sig, mixed_eps = None, None
        if self.comb_rule is None or self.comb_rule == "lorentz-bertholot":
            mixed_sig = 0.5 * (water_sig + vdw[self.all_atom_ids, 0])
            mixed_eps = np.sqrt(water_eps * vdw[self.all_atom_ids, 1])
        if self.comb_rule == "geometric":
            mixed_sig = np.sqrt(water_sig * vdw[self.all_atom_ids, 0])
            mixed_eps = np.sqrt(water_eps * vdw[self.all_atom_ids, 1])

        if mixed_eps is not None and mixed_sig is not None:
            acoeff = 4 * mixed_eps * (mixed_sig**12)
            bcoeff = 4 * mixed_eps * (mixed_sig**6)
        else:
            raise Exception("Couldn't assign vdw params")
        return chg_product, acoeff, bcoeff


    def calculate_hydrogen_bonds(self, traj, water, nbrs, water_water=True):
        """Calculates hydrogen bonds made by a water molecule with its first shell
        water and solute neighbors.

        Parameters
        ----------
        traj : md.trajectory
            MDTraj trajectory object for which hydrogen bonds are to be calculates.
        water : int
            The index of water oxygen atom
        nbrs : np.ndarray
            Indices of the water oxygen or solute atoms in the first solvation shell of the water molecule.
        water_water : bool
            Boolean for whether water-water or solute-water hbonds are desired

        Returns
        -------
        hbonds : np.ndarray
            Array of hydrogen bonds where each hydrogen bond is represented by an array of indices
            of three participating atoms [Donor, H, Acceptor]
        """
        hbond_data = []
        angle_triplets = []
        if water_water:
            for wat_nbr in nbrs:
                angle_triplets.extend(
                    [[water, wat_nbr, wat_nbr + 1], [water, wat_nbr, wat_nbr + 2], [wat_nbr, water, water + 1],
                     [wat_nbr, water, water + 2]])
        else:
            for solute_nbr in nbrs:
                if self.prot_hb_types[solute_nbr] == 1 or self.prot_hb_types[solute_nbr] == 3:
                    angle_triplets.extend([[solute_nbr, water, water + 1], [solute_nbr, water, water + 2]])
                if self.prot_hb_types[solute_nbr] == 2 or self.prot_hb_types[solute_nbr] == 3:
                    for don_H_pair in self.don_H_pair_dict[solute_nbr]:
                        angle_triplets.extend([[water, solute_nbr, don_H_pair[1]]])

        angle_triplets = np.asarray(angle_triplets)
        angles = md.compute_angles(traj, angle_triplets)
        angles[np.isnan(angles)] = 0.0
        hbonds = angle_triplets[np.where(angles[0, :] <= ANGLE_CUTOFF_RAD)]
        return hbonds


    def water_nbr_orientations(self, traj, water, nbrs):
        """Calculates orientations of the neighboring water molecules of a given water molecule. The orientation is
        defined as the miniumum of four possible Donor-H-Acceptor angles.

        Parameters
        ----------
        traj : md.trajectory
            MDTraj trajectory object for which hydrogen bonds are to be calculates.
        water : int
            The index of water oxygen atom
        nbrs : np.ndarray, int, shape=(N^{ww}_nbr, )
            Indices of the water oxygen or solute atoms in the first solvation shell of the water molecule.

        Returns
        -------
        wat_orientations : np.ndarray
            Array of angles between the given and each one of its neighbors
        """

        angle_triplets = []
        for wat_nbr in nbrs:
            angle_triplets.extend(
                [[water, wat_nbr, wat_nbr + 1], [water, wat_nbr, wat_nbr + 2], [wat_nbr, water, water + 1],
                 [wat_nbr, water, water + 2]])
        angle_triplets = np.asarray(angle_triplets)
        angles = md.compute_angles(traj, angle_triplets)
        angles[np.isnan(angles)] = 0.0
        wat_orientations = [np.rad2deg(np.min(angles[0, i*4:(i*4)+4])) for i in range(nbrs.shape[0])]
        return wat_orientations
