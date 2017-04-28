from __future__ import absolute_import
import sys
import os
sys.path.append("../")

from grid_water_analysis import GridWaterAnalysis
from site_water_analysis import SiteWaterAnalysis
from readtestsystems import *


def main():
    # run a test bulk water simulation
    tip3p_hsa = read_tip3p()
    tip3p_hsa.initialize_hydration_sites()
    tip3p_hsa.print_system_summary()
    tip3p_hsa.calculate_site_quantities(start_frame=0, num_frames=10, entropy=False)
    tip3p_hsa.write_calculation_summary()
    tip3p_hsa.write_data()

if __name__ == '__main__':
    status = main()
