from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
import sys
import os

import numpy as np
from scipy import stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import cm


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
        print(os.listdir(data_dir))
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
            nbr_values.append(old_div(np.sum(nbrs),nbrs.shape[0]))
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
            vol = (old_div(4.0,3.0))*np.pi*(d**3)
            vol_low = (old_div(4.0,3.0))*np.pi*(d_low**3)
            shell_vol = vol - vol_low

            counts_bulk = 0.0329*shell_vol
            sum_counts_kernel += np.sum(Z[:, i])
            #Z[:,i] /= counts_bulk
            Z[:, i] = old_div(Z[:, i],counts_bulk)

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


def main():
    directory = "/Users/kamranhaider/sstmap_testcases/test_utils_plotting/casp3_cluster_ene_data/"
    ref = "/Users/kamranhaider/sstmap_testcases/test_utils_plotting/bulkwat_Enbr"
    plot_enbr_distribution(directory, site_indices=[0, 4], ref_data=ref, ref_nbrs=5.25)

if __name__ == '__main__':
    main()


