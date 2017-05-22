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
# License along with SSTMap. If not, see <http://www.gnu.org/licenses/>.
##############################################################################

import sys
import os
import time
from functools import wraps
import _sstmap_ext as calc

import numpy as np
from scipy import stats
import mdtraj as md
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import cm

#rcParams['font.family'] = 'serif'
#rcParams['font.serif'] = ['Cambria Math'] + rcParams['font.serif']

##############################################################################
# Utilities
##############################################################################

def function_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print ("Total time running %s: %2.2f seconds" %
               (function.__name__, t1-t0))
        return result
    return function_timer

def print_progress_bar (count, total):
    """
    Create and update progress bar during a loop.
    
    Parameters
    ----------
    iteration : int
        The number of current iteration, used to calculate current progress. 
    total : int
        Total number of iterations
    
    Notes
    -----
        Based on:
        http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = "=" * filled_len + ' ' * (bar_len - filled_len)

    sys.stdout.write('Progress |%s| %s%s Done\r' % (bar, percents, '%'))
    sys.stdout.flush()
    if count == total: 
        print



def plot_enbr_distribution(data_dir, site_indices=None, nbr_norm=False, ref_data=None, ref_nbrs=None):
    """
    Generate an Enbr plot for an arbitrary list of sites. First site should be the reference system.
    sites: a list of keys which represent site labels
    data: a dictionary of sites
    x_values: data points on x-axis
    nbr_norm: Normalize by number of neighbors
    outname: name of output file
    
    Parameters
    ----------
    data_dir : TYPE
        Description
    site_indices : None, optional
        Description
    nbr_norm : bool, optional
        Description
    ref_data : None, optional
        Description
    ref_nbrs : None, optional
        Description
    """
    enbr_files = []
    enbr = {}
    ref_enbr = None
    nbr_files = []
    nbr_values = []

    if not os.path.isdir(data_dir):
        sys.exit(
            "Data directory not found, please check path of the directory again.")

    if site_indices is None:
        enbr_files = [
            f for f in os.listdir(data_dir) if f.endswith("Ewwnbr.txt")]
        if nbr_norm:
            nbr_files = [
                f for f in os.listdir(data_dir) if f.endswith("Nnbrs.txt")]
    else:
        enbr_files = [f for f in os.listdir(data_dir) if f.endswith(
            "Ewwnbr.txt") and int(f[0:3]) in site_indices]
        if nbr_norm:
            nbr_files = [f for f in os.listdir(data_dir) if f.endswith(
                "Nnbrs.txt") and int(f[0:3]) in site_indices]

    for index, file in enumerate(enbr_files):
        site_i = int(file[0:3])
        enbr[site_i] = np.loadtxt(data_dir + "/" + file)
        if nbr_norm:
            nbrs = np.loadtxt(data_dir + "/" + nbr_files[index])
            nbr_values.append(np.sum(nbrs) /nbrs.shape[0])
    if ref_data is not None:
        ref_enbr = np.loadtxt(ref_data)
        if nbr_norm:
            ref_enbr *= ref_nbrs

    for index, site_i in enumerate(enbr.keys()):
        print("Generating Enbr plot for: ", site_i, enbr_files[index])
        # Get x and p_x for current site
        site_enbr = enbr[site_i]*0.5
        x_low, x_high = -5.0, 3.0
        enbr_min, enbr_max = np.min(site_enbr), np.max(site_enbr)
        if enbr_min < x_low:
            x_low = enbr_min
        if enbr_max > x_high:
            x_high = enbr_max

        x = np.linspace(x_low, x_high)
        kernel = stats.gaussian_kde(site_enbr)
        p_x = kernel.evaluate(x)
        if nbr_norm:
            site_nbrs = nbr_values[index]
            p_x *= site_nbrs
        # Get x and p_x for reference site, if available
        p_x_ref = None
        if ref_enbr is not None:
            kernel = stats.gaussian_kde(ref_enbr)
            p_x_ref = kernel.evaluate(x)
        # Set up plot
        fig, ax = plt.subplots(1)
        fig.set_size_inches(3, 3)
        plt.xlim(x_low, x_high)
        plt.ylim(0.0, np.max(p_x) + 0.1)
        start, end = ax.get_ylim()
        ax.yaxis.set_ticks(np.arange(start, end, 0.2))
        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(start, end, 2.0))
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
        x_label = r'$\mathit{E_{n} (kcal/mol)}$'
        y_label = r'$\mathit{\rho(E_{n})}$'
        if nbr_norm:
            y_label = r'$\mathit{\rho(E_{n})N^{nbr}}$'
        ax.set_xlabel(x_label, size=14)
        ax.set_ylabel(y_label, size=14)
        ax.yaxis.tick_left()
        ax.xaxis.tick_bottom()
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        plt.minorticks_on()
        plt.tick_params(which='major', width=1, length=4, direction='in')
        plt.tick_params(which='minor', width=1, length=2, direction='in')
        plt.tick_params(axis='x', labelsize=12)
        plt.tick_params(axis='y', labelsize=12)
        plt.plot(
            x, p_x, antialiased=True, linewidth=1.0, color="red", label=site_i)
        if p_x_ref is not None:
            plt.plot(x, p_x_ref, antialiased=True, linewidth=1.0,
                     color="green", label="Reference")
        fig_name = "%03d_" % site_i
        plt.legend(loc='upper right', prop={'size': 10}, frameon=False)
        plt.tight_layout()
        plt.savefig(data_dir + "/" + fig_name + "Enbr_plot.png", dpi=300)
        plt.close()


def plot_rtheta_distribution(data_dir, site_indices=None):

    """
    Parameters
    ----------
    data_dir : TYPE
        Description
    site_indices : None, optional
        Description
    
    """
    rtheta_files = []
    rtheta_data = {}

    print(data_dir)
    if not os.path.isdir(data_dir):
        sys.exit(
            "Data directory not found, please check path of the directory again.")

    if site_indices is None:
        rtheta_files = [
            f for f in os.listdir(data_dir) if f.endswith("r_theta.txt")]
    else:
        rtheta_files = [f for f in os.listdir(data_dir) if f.endswith(
            "r_theta.txt") and int(f[0:3]) in site_indices]

    for index, file in enumerate(rtheta_files):
        site_i = int(file[0:3])
        rtheta_data[site_i] = np.loadtxt(data_dir + "/" + file)

    integ_counts = 16.3624445886
    for index, site_i in enumerate(rtheta_data.keys()):
        print("Generating r_theta plot for: ", site_i, rtheta_files[index])
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        theta = rtheta_data[site_i][:, 0]
        r = rtheta_data[site_i][:, 1]
        #Nnbr = len(r)/nwat
        # print nwat, Nnbr
        # generate index matrices
        X, Y = np.mgrid[0:130:131j, 2.0:6.0:41j]
        # generate kernel density estimates
        values = np.vstack([theta, r])
        kernel = stats.gaussian_kde(values)
        positions = np.vstack([X.ravel(), Y.ravel()])
        Z = np.reshape(kernel(positions).T, X.shape)
        Z *= integ_counts*0.1
        #Z /= integ_counts
        sum_counts_kernel = 0
        # print kernel.n
        # correct Z
        for i in range(0, Y.shape[1]):
            d = Y[0, i]
            # get shell_vol
            d_low = d - 0.1
            vol = (4.0 / 3.0) * np.pi * (d**3)
            vol_low = (4.0 / 3.0) * np.pi * (d_low**3)
            shell_vol = vol - vol_low

            counts_bulk = 0.0329*shell_vol
            sum_counts_kernel += np.sum(Z[:, i])
            #Z[:,i] /= counts_bulk
            Z[:, i] = Z[:, i],counts_bulk

        print(sum_counts_kernel)
        legend_label = "%03d_" % site_i
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0.5,
                        antialiased=True, alpha=1.0, cmap=cm.coolwarm, label=legend_label)
        x_label = r"$\theta^\circ$"
        y_label = r"$r (\AA)$"
        ax.set_xlabel(x_label)
        ax.set_xlim(0, 130)
        ax.set_ylabel(y_label)
        ax.set_ylim(2.0, 6.0)
        z_label = r'$\mathrm{P(\theta, \AA)}$'
        ax.set_zlabel(z_label)
        #ax.legend(legend_label, loc='upper left', prop={'size':6})
        #ax.set_zlim(0.0, 0.15)
        plt.savefig(data_dir + "/" + legend_label + "rtheta_plot.png", dpi=300)
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

def read_gist_summary(gist_data_file):
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

def write_watpdb_from_list(traj, filename, water_id_list, full_water_res=False):
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
    
    pdb_line_format = "{0:6}{1:>5}  {2:<3}{3:<1}{4:>3} {5:1}{6:>4}{7:1}   {8[0]:>8.3f}{8[1]:>8.3f}{8[2]:>8.3f}{9:>6.2f}{10:>6.2f}{11:>12s}\n"
    ter_line_format = "{0:3}   {1:>5}      {2:>3} {3:1}{4:4} \n"
    pdb_lines = []
    # write form the list of (water, frame) tuples
    coords = md.utils.in_units_of(traj.xyz, "nanometers", "angstroms")
    # at_index, wat in enumerate(water_id_list):
    at = 1
    res = 1
    with open(filename + ".pdb", 'w') as f:
        for i in range(0,len(water_id_list)):
            wat = water_id_list[i]
            at_index = at #% 10000
            res_index = res % 10000
            #wat_coords = md.utils.in_units_of(
            #    coords[wat[0], wat[1], :], "nanometers", "angstroms")
            wat_coords = coords[wat[0], wat[1], :]
            #chain_id = possible_chains[chain_id_index]
            chain_id = "A"
            pdb_line = pdb_line_format.format(
                "ATOM", at_index, "O", " ", "WAT", chain_id, res_index, " ", wat_coords, 0.00, 0.00, "O")
            #pdb_lines.append(pdb_line)
            f.write(pdb_line)
        
            if full_water_res:
                #H1_coords = md.utils.in_units_of(
                #    coords[wat[0], wat[1] + 1, :], "nanometers", "angstroms")
                H1_coords = coords[wat[0], wat[1] + 1, :]
                pdb_line_H1 = pdb_line_format.format("ATOM", at_index + 1, "H1", " ", "WAT", chain_id, res_index, " ", H1_coords, 0.00, 0.00, "H")
                #pdb_lines.append(pdb_line_H1)
                f.write(pdb_line_H1)
                #H2_coords = md.utils.in_units_of(
                #    coords[wat[0], wat[1] + 2, :], "nanometers", "angstroms")
                H2_coords = coords[wat[0], wat[1] + 2, :]
                pdb_line_H2 = pdb_line_format.format("ATOM", at_index + 2, "H2", " ", "WAT", chain_id, res_index, " ", H2_coords, 0.00, 0.00, "H")
                #pdb_lines.append(pdb_line_H2)
                f.write(pdb_line_H2)
                at += 3
                res += 1
            else:
                at += 1
                res += 1
            if res_index == 9999:
                ter_line = ter_line_format.format(
                    "TER", at, "WAT", chain_id, res_index)
                at = 1
                #pdb_lines.append(ter_line)
    #pdb_lines.append("END")
    #np.savetxt(filename + ".pdb", np.asarray(pdb_lines), fmt="%s")



def write_watpdb_from_coords(traj, filename, wat_coords):
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
    
    pdb_line_format = "{0:6}{1:>5}  {2:<3}{3:<1}{4:>3} {5:1}{6:>4}{7:1}   {8[0]:>8.3f}{8[1]:>8.3f}{8[2]:>8.3f}{9:>6.2f}{10:>6.2f}{11:>12s}\n"
    ter_line_format = "{0:3}   {1:>5}      {2:>3} {3:1}{4:4} \n"
    pdb_lines = ["REMARK Initial number of clusters: N/A\n"]
    # write form the list of (water, frame) tuples
    for at in range(len(wat_coords)):
        wat_coord = wat_coords[at]
        at_index = at % 10000
        res_index = at % 10000
        chain_id = "A"
        pdb_line = pdb_line_format.format(
            "ATOM", at_index, "O", " ", "WAT", chain_id, res_index, " ", wat_coord, 0.00, 0.00, "O")
        pdb_lines.append(pdb_line)
        if res_index == 9999:
            ter_line = ter_line_format.format(
                "TER", at_index, "WAT", chain_id, res_index)
            pdb_lines.append(ter_line)
    
    with open(filename + ".pdb", "w") as f:
        f.write("".join(pdb_lines))

class NeighborSearch(object):
    """
    Class for relatively fast queries of coordinates within a distance
    of specified coordinate. 
    """
    def __init__(self, xyz, dist):
        """Initialize a NeighborSearch object by providing an array of
        coordinates and a distance threshold.

        Parameters
        ----------
        xyz : np.ndarray, float, shape=(N, 3)
            A multidmimensional array of three dimensional coordinates
        dist : float
            A distance cutoff to identify points within this distance of the
            query point.
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
        Given a coordinate point, return all point indexes (0-indexed) and 
        corresponding distances that are within the threshold distance from it.
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
                (point - self.min_) / self.cell_size, dtype=np.int)
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



def main():
    #plot_enbr_distribution("test_1_hsa_data", site_indices=[1, 3], nbr_norm=True)
    plot_rtheta_distribution(
        "test_1_angular_structure_data", site_indices=[1, 3])


if __name__ == '__main__':
    main()
