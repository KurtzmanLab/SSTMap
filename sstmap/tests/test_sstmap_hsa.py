from __future__ import absolute_import
import sys
import os
#sys.path.append("../")

from sstmap.site_water_analysis import SiteWaterAnalysis


#from readtestsystems import *
testsystem_dir = os.path.abspath("testsystems") + "/"

def main():
    # run a test arg calculation
    #clusters_arg = testsystem_dir + "arg_amber/ref_clustercenterfile.pdb"
    arg_hsa = read_arg_amber_hsa()
    arg_hsa.initialize_hydration_sites()
    arg_hsa.print_system_summary()
    arg_hsa.calculate_site_quantities(num_frames=100)
    arg_hsa.write_calculation_summary()
    arg_hsa.write_data()
    arg_hsa.calculate_angular_structure(site_indices=[1, 2, 3, 4], num_frames=100)
    arg_hsa.calculate_lonranged_ww_energy(site_indices=[1, 2, 3, 4], num_frames=100)
    # add test for long_range and r_theta
    # run a test casp calculation

if __name__ == '__main__':
    status = main()
