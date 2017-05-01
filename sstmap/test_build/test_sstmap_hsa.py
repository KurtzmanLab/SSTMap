import sys
import os
#sys.path.append("/Users/kamranhaider/Dropbox/SSTMap/sstmap")

from sstmap.site_water_analysis import SiteWaterAnalysis
from readtestsystems import *

testsystem_dir = os.path.abspath("testsystems") + "/"

def main():
    # run a test arg calculation
    #clusters_arg = testsystem_dir + "arg_amber/ref_clustercenterfile.pdb"
    arg_hsa = read_arg_amber_hsa()
    arg_hsa.initialize_hydration_sites()
    arg_hsa.print_system_summary()
    arg_hsa.calculate_site_quantities(num_frames=100, entropy=False)
    arg_hsa.write_calculation_summary()
    arg_hsa.write_data()

if __name__ == '__main__':
    status = main()
