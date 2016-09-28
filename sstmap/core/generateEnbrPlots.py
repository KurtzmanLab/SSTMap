# system imports
import os
from optparse import OptionParser

# Third party imports
import numpy as np
from scipy import stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from optparse import OptionParser
#import seaborn as sns
#sns.set_style("white")
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Cambria Math'] + rcParams['font.serif'] 


x_range = np.linspace(-4.0, 2.5)

def genPlot(data, save_dir, x_values = x_range, nbr_norm = False):
	"""
	generate an Enbr plot for an arbitrary list of sites. First site should be the reference system.
	sites: a list of keys which represent site labels
	data: a dictionary of sites
	x_values: data points on x-axis
	nbr_norm: Normalize by number of neighbors
	outname: name of output file
	"""
	colors = ["green", "blue", "red", "orange"]
	for k in data.keys():
		"generating Enbr plot for: ", k
		y_values = data[k][0]*data[k][1]
		for val in range(len(x_values)):
			if x_values[val] <= -3.26:
				 y_values[val] = 0.0

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

def getHSAData(hsa_data_file):
    '''
    Returns a dictionary with hydration site index as keys and a list of various attributes as values.
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


data_dict = {}


parser = OptionParser()
parser.add_option("-i", "--input_hsa_summary", dest="hsa_summary_file", type="string", help="Input HSA summary file")
parser.add_option("-d", "--data_directory", dest="data_directory", type="string", help="Location of Enbr data")
(options, args) = parser.parse_args()

#hsa_data = getHSAData(options.input_hsa_summary)
Enbr_files = []

for f in os.listdir(options.data_directory):
    if f.endswith("pair_ene"):
        clust_number = f[0:3]
        filename = options.data_directory + "/" + f
        data = np.loadtxt(filename)
        print "generating kernel for ", f
        kernel = stats.gaussian_kde(data)
        enbr_values = kernel.evaluate(x_range)
        data_dict[clust_number] = [enbr_values, 1.0]


# generate Figure 1 middle plot
print "generating plots"
#Fig1 = ["Pure Water", "Methane", "CB7", "Phe"]
genPlot(data_dict, options.data_directory)
