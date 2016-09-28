"""
docstring for this class
"""
import sys
import time
import os

import numpy as np
from scipy import stats
import mdtraj as md
import parmed as pmd

from water_analysis import WaterAnalysisBase, NeighborSearch
DEG_PER_RAD = 180.0 / np.pi
# root class definition


class ProteinWaterAnalysis(WaterAnalysisBase):
    """docstring for SiteWaterAnalysis"""

    def __init__(self, topology_file, trajectory, start_frame=0, num_frames=0,
                 near_ligand=True, desmond_helper_file=None, ligand_file=None, prefix="test"):

        super(
            ProteinWaterAnalysis,
            self).__init__(
            topology_file,
            trajectory,
            start_frame,
            num_frames,
            desmond_helper_file)
        self.ligand_file = ligand_file
        self.prefix = prefix
        # generate an extra list of solute atom ids that within a cutoff
        # distance from the ligand
        self.near_lig = np.zeros((len(self.non_water_atom_ids)), dtype=bool)
        # By default (and mostly) H-bonding groups in the vicinity of the
        # ligand are processed.
        if near_ligand:
            ligand = md.load_pdb(self.ligand_file)
            ligand_coords = ligand.xyz[0, :, :] * 10.0
            first_frame = md.load_frame(
                self.trajectory, 0, top=self.topology_file)
            solute_pos = first_frame.xyz[0, self.non_water_atom_ids, :] * 10.0

            search_space = NeighborSearch(solute_pos, 5.0)
            near_indices = search_space.query_nbrs_multiple_points(
                ligand_coords)
            solute_nbrs = [self.non_water_atom_ids[nbr_index]
                           for nbr_index in near_indices]
            for nbr in solute_nbrs:
                self.near_lig[nbr] = True
        # if near_ligand flag is set to False i.e., the user is interested in the whole surface
        # we set all elements of self.near_lig to True
        else:
            self.near_lig[:] = True
        # Set up dictionary to store the results of calculations
        self.prot_hbond_data = {}
        for hb_group in np.concatenate(
                (self.solute_acc_ids, self.solute_acc_don_ids, self.solute_don_ids)):
            hb_group_name = str(self.topology.atom(hb_group))
            # construct name
            self.prot_hbond_data[hb_group] = [
                hb_group_name, np.zeros(6, dtype="float64")]

    def calculate_protein_water_hbonds(self):
        # FIXME: Protein-water h-bonds are calculated with the assumption that atom order in topology file for water is
        # water O, H, H, dummy atoms. Check topology to confirm this.
        # generate concatenated arrays of don + acc_don and acc + acc_don ids
        combined_don_accdon_ids = np.concatenate(
            (self.solute_don_ids, self.solute_acc_don_ids))
        combined_acc_accdon_ids = np.concatenate(
            (self.solute_acc_ids, self.solute_acc_don_ids))
        for i in xrange(self.start_frame, self.start_frame + self.num_frames):
            print "Processing frame: ", i + 1, "..."
            # Load frame and obtain coordinates
            frame = md.load_frame(self.trajectory, i, top=self.topology)
            pos = frame.xyz[0, :, :] * 10.0
            # get water oxygen atom coordinates and construct search space
            oxygen_pos = pos[self.wat_oxygen_atom_ids]
            water_search_space = NeighborSearch(oxygen_pos, 3.5)
            # get solute atom coordinates and construct search space
            solute_pos = pos[self.non_water_atom_ids]
            solute_search_space = NeighborSearch(solute_pos, 3.5)
            # print "Processing acceptors"
            # iterate over solute acceptor atoms
            for solute_acceptor in combined_acc_accdon_ids:
                # check if current acceptor is near ligand
                if self.near_lig[solute_acceptor]:
                    # print "\tacceptor: ", solute_acceptor, self.prot_hbond_data[solute_acceptor][0]
                    # obtain water neighbors (oxygens)
                    wat_nbr_indices = water_search_space.query_nbrs_single_point(pos[
                                                                                 solute_acceptor])
                    nbr_wat_oxygens = [self.wat_oxygen_atom_ids[
                        nbr_index] for nbr_index in wat_nbr_indices]
                    # update water neighbors for this acceptor
                    self.prot_hbond_data[solute_acceptor][
                        1][0] += len(nbr_wat_oxygens)
                    # obtain neighboring solute atoms
                    solute_nbr_indices = solute_search_space.query_nbrs_single_point(pos[
                                                                                     solute_acceptor])
                    solute_nbr_atom_ids = [self.non_water_atom_ids[
                        nbr_index] for nbr_index in solute_nbr_indices]
                    # extract donors
                    nbr_solute_donors = np.intersect1d(
                        solute_nbr_atom_ids, combined_don_accdon_ids)
                    # iterate over each neighbor water oxygen for this acceptor and determine H-bonds with water
                    # print "\t\tWater neighbors:", nbr_wat_oxygens
                    for wat_O in nbr_wat_oxygens:
                        hbangle_combinations = np.asarray(
                            [[solute_acceptor, wat_O, wat_O + 1], [solute_acceptor, wat_O, wat_O + 2]])
                        theta_list = md.compute_angles(
                            frame, hbangle_combinations) * DEG_PER_RAD
                        hbangle = np.min(theta_list)
                        if hbangle <= 30.0:
                            self.prot_hbond_data[solute_acceptor][1][2] += 1
                    # print "\t\tSolute neighbors:", nbr_solute_donors
                    # iterate over each neighbor solute donor for this acceptor
                    # and determine internal H-bonds
                    for donor in nbr_solute_donors:
                        for don_H_pair in self.don_H_pair_dict[donor]:
                            hbangle = md.compute_angles(frame, np.asarray(
                                [[solute_acceptor, donor, don_H_pair[1]]])) * DEG_PER_RAD
                            if hbangle <= 30.0:
                                self.prot_hbond_data[
                                    solute_acceptor][1][4] += 1
                    # print len(nbr_wat_oxygens), nbr_solute_donors
            # iterate over solute donor atoms
            # print "Processing donors"
            for solute_donor in combined_don_accdon_ids:
                # check if current donor is near ligand
                if self.near_lig[solute_donor]:
                    # print "\tdonor: ", solute_donor, self.prot_hbond_data[solute_donor][0]
                    # print solute_donor, self.prot_hbond_data[solute_donor][0],
                    # obtain water neighbors (oxygens)
                    wat_nbr_indices = water_search_space.query_nbrs_single_point(pos[
                                                                                 solute_donor])
                    nbr_wat_oxygens = [self.wat_oxygen_atom_ids[
                        nbr_index] for nbr_index in wat_nbr_indices]
                    # update water neighbors for this donor (skip donors in acc_don_ids to avoid over-counting)
                    # if solute_donor not in self.solute_acc_don_ids:
                    self.prot_hbond_data[solute_donor][
                        1][0] += len(nbr_wat_oxygens)
                    # obtain neighboring solute atoms
                    # print "\t\tWater neighbors:", nbr_wat_oxygens
                    solute_nbr_indices = solute_search_space.query_nbrs_single_point(pos[
                                                                                     solute_donor])
                    solute_nbr_atom_ids = [self.non_water_atom_ids[
                        nbr_index] for nbr_index in solute_nbr_indices]
                    # extract acceptors
                    nbr_solute_acceptors = np.intersect1d(
                        solute_nbr_atom_ids, combined_acc_accdon_ids)
                    # print "\t\tSolute neighbors:", nbr_solute_acceptors
                    # iterate over each neighbor water oxygen for this donor
                    for solute_don_H_pair in self.don_H_pair_dict[
                            solute_donor]:
                        for wat_O in nbr_wat_oxygens:
                            hbangle = md.compute_angles(frame, np.asarray(
                                [[wat_O, solute_donor, solute_don_H_pair[1]]])) * DEG_PER_RAD
                            if hbangle <= 30.0:
                                self.prot_hbond_data[solute_donor][1][3] += 1
                                # break
                        for acceptor in nbr_solute_acceptors:
                            hbangle = md.compute_angles(frame, np.asarray(
                                [[acceptor, solute_donor, solute_don_H_pair[1]]])) * DEG_PER_RAD
                            if hbangle <= 30.0:
                                self.prot_hbond_data[solute_donor][1][5] += 1
                                # break
                        # print len(nbr_wat_oxygens), nbr_solute_acceptors

        for hb_group in self.prot_hbond_data:
            if self.near_lig[hb_group]:
                self.prot_hbond_data[hb_group][1][0] /= self.num_frames
                self.prot_hbond_data[hb_group][1][2] /= self.num_frames
                self.prot_hbond_data[hb_group][1][3] /= self.num_frames
                self.prot_hbond_data[hb_group][1][4] /= self.num_frames
                self.prot_hbond_data[hb_group][1][5] /= self.num_frames

    def writeHBsummary(self):
        f = open(self.prefix + "_prot_hb_summary.txt", "w")
        header = "atom_index atom_name wat_nbrs hbsw_tot hbss_tot hbsw_acc hbsw_don hbss_acc hbss_don\n"
        f.write(header)
        for hb_group in self.prot_hbond_data:
            if self.near_lig[hb_group]:
                d = self.prot_hbond_data[hb_group][1]
                l = "%i %s %f %f %f %f %f %f %f\n" % (hb_group, self.prot_hbond_data[hb_group][
                                                      0], d[0], d[2] + d[3], d[4] + d[5], d[2], d[3], d[4], d[5])
                f.write(l)
        f.close()

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
    if len(sys.argv[1:]) == 0:
        print "no argument given!"
        parser.print_help()
    else:
        print "Setting up calculations..."
        p = ProteinWaterAnalysis(
            options.prmtop,
            options.trjname,
            start_frame=options.start_frame,
            num_frames=options.frames,
            ligand_file=options.ligand,
            prefix=options.prefix)
        p.calculate_protein_water_hbonds()
        p.writeHBsummary()


def entry_point():
    main()

if __name__ == '__main__':
    entry_point()
