import sys
import os

import numpy as np
from scipy import stats

def plot_enbr_distribution(data_dir, site_indices = [], nbr_norm = False, ref_data = None, ref_nbrs = None):
    """
    generate an Enbr plot for an arbitrary list of sites. First site should be the reference system.
    sites: a list of keys which represent site labels
    data: a dictionary of sites
    x_values: data points on x-axis
    nbr_norm: Normalize by number of neighbors
    outname: name of output file
    """
    enbr_files = []
    enbr_data = {}
    ref_enbr_data = None
    nbr_files = []
    nbr_values = []

    if not os.path.isdir(data_dir):
        sys.exit("Data directory not found, please check path of the directory again.")

    if len(site_indices) == 0:
        enbr_files = [f for f in os.listdir(data_dir) if f.endswith("Ewwnbr.txt")]
        if nbr_norm:
            nbr_files = [f for f in os.listdir(data_dir) if f.endswith("Nnbrs.txt")]
    else:
        enbr_files = [f for f in os.listdir(data_dir) if f.endswith("Ewwnbr.txt") and int(f[0:3]) in site_indices]
        if nbr_norm:
            nbr_files = [f for f in os.listdir(data_dir) if f.endswith("Nnbrs.txt") and int(f[0:3]) in site_indices]

    print nbr_files
    for index, file in enumerate(enbr_files):
        site_i = int(file[0:3])
        enbr_data[site_i] = np.loadtxt(data_dir + "/" + file)
        if nbr_norm:
            nbrs = np.loadtxt(data_dir + "/" + nbr_files[index])
            nbr_values.append(np.sum(nbrs)/nbrs.shape[0])
    if ref_data is not None:
        ref_enbr_data = np.loadtxt(ref_data)


    for inex, site_i in enumerate(enbr_data.keys()):
        print "Generating Enbr plot for: ", site_i
        print nbr_values[index]
    """
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
def plot_rtheta_distribution(data_dir, site_indices):
    
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
    plot_enbr_distribution("test_1_hsa_data", site_indices=[1, 3], nbr_norm=True)


if __name__ == '__main__':
    main()
