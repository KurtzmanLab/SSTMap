from __future__ import absolute_import
import sys
sys.path.append("../")

from grid_water_analysis import GridWaterAnalysis
from site_water_analysis import SiteWaterAnalysis
from utils import print_progress_bar

testsystem_dir = "../testsystems/"

def read_tip3p_desmond():
	"""Can sstmap read a pure water tip3p system?
	"""
	top = testsystem_dir + "tip3p_desmond/tip3p.pdb"
	traj = testsystem_dir + "tip3p_desmond/tip3p_desmond.nc"
	nbr_param_file = testsystem_dir + "tip3p_desmond/tip3p_cms_nb_parms.txt"
	clusters = testsystem_dir + "tip3p_desmond/clustercenterfile.pdb"

	hsa = SiteWaterAnalysis(top, traj, start_frame=0, num_frames=1000,
							clustercenter_file=clusters, desmond_helper_file=nbr_param_file,
							prefix="tip3p_desmond")
	hsa.print_system_summary()



read_tip3p_desmond()