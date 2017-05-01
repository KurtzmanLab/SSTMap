import sys
import os

from sstmap.grid_water_analysis import GridWaterAnalysis
from readtestsystems import *

def main():
    # run a test arg calculation
    #clusters_casp = testsystem_dir + "casp3_desmond/ref_clustercenterfile.pdb"
    arg_gist = read_arg_amber_gist()
    arg_gist.print_system_summary()
    arg_gist.calculate_grid_quantities(num_frames=1000)
    arg_gist.write_data()
    arg_gist.generate_dx_files()
    

if __name__ == '__main__':
    status = main()
