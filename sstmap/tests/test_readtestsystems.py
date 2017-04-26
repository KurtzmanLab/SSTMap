from __future__ import absolute_import
import sys
sys.path.append("../")

from grid_water_analysis import GridWaterAnalysis
from site_water_analysis import SiteWaterAnalysis
from utils import *

testsystem_dir = "../testsystems/"

def read_tip3p_desmond():
    """Can sstmap read a pure water tip3p system?
    """
    top = testsystem_dir + "tip3p_desmond/tip3p.pdb"
    traj = testsystem_dir + "tip3p_desmond/tip3p_desmond.nc"
    nb_param_file = testsystem_dir + "tip3p_desmond/tip3p_cms_nb_parms.txt"
    clusters = testsystem_dir + "tip3p_desmond/clustercenterfile.pdb"

    hsa = SiteWaterAnalysis(top, traj, start_frame=0, num_frames=1000,
                            clustercenter_file=clusters, desmond_helper_file=nb_param_file,
                            prefix="tip3p_desmond")
    return hsa

def read_casp3_desmond_with_clusters():
    """
    Can sstmap read a pure water tip3p system in amber
    """
    top = testsystem_dir + "casp3_desmond/casp3.pdb"
    traj = "/Users/kamranhaider/sstmap_testcases/casp3_desmond/simulation/casp3_pd_2/casp3_converted.nc"
    nb_param_file = testsystem_dir + "casp3_desmond/casp3_nb_parms.txt"
    ligand = testsystem_dir + "casp3_desmond/ligand.pdb"
    clusters = testsystem_dir + "casp3_desmond/ref_clustercenterfile.pdb"
    hsa = SiteWaterAnalysis(top, traj, start_frame=0, num_frames=10000, 
                            ligand_file=ligand,
                            clustercenter_file=clusters, 
                            desmond_helper_file=nb_param_file,
                            prefix="casp")
    return hsa

def read_casp3_desmond_without_clusters():
    """
    Can sstmap read a pure water tip3p system in amber
    """
    top = testsystem_dir + "casp3_desmond/casp3.pdb"
    traj = "/Users/kamranhaider/sstmap_testcases/casp3_desmond/simulation/casp3_pd_2/casp3_converted.nc"
    nb_param_file = testsystem_dir + "casp3_desmond/casp3_nb_parms.txt"
    ligand = testsystem_dir + "casp3_desmond/ligand.pdb"
    hsa = SiteWaterAnalysis(top, traj, start_frame=0, num_frames=10000, 
                            ligand_file=ligand, desmond_helper_file=nb_param_file,
                            prefix="casp")
    return hsa

def read_arg_amber_without_clusters():
    """Read Arg amber test system
    """
    top = testsystem_dir + "arg_amber/arg.prmtop"
    #traj = testsystem_dir + "arg_amber/md1ns.nc"
    traj = "/Users/kamranhaider/sstmap_testcases/arg_amber/md10ns.nc"
    ligand = testsystem_dir + "arg_amber/ligand.pdb"
    hsa = SiteWaterAnalysis(top, traj, start_frame=0, num_frames=10000, 
                            ligand_file=ligand,
                            prefix="arg")
    return hsa

def read_arg_amber_with_clusters():
    """Read Arg amber test system
    """
    top = testsystem_dir + "arg_amber/arg.prmtop"
    traj = testsystem_dir + "arg_amber/md1ns.nc"
    clusters = testsystem_dir + "arg_amber/ref_clustercenterfile.pdb"
    ligand = testsystem_dir + "arg_amber/ligand.pdb"
    hsa = SiteWaterAnalysis(top, traj, start_frame=0, num_frames=1000, 
                            clustercenter_file=clusters, ligand_file=ligand,
                            prefix="arg")
    return hsa

if __name__ == '__main__':
    #print("Generating reference system ...")
    #ref_hsa = read_casp3_desmond_with_clusters()
    #ref_hsa.initialize_hydration_sites()
    print("Generating test system ...")
    #test_hsa = read_casp3_desmond_without_clusters()
    test_hsa = read_casp3_desmond_with_clusters()
    test_hsa.initialize_hydration_sites()
    test_hsa.print_system_summary()
    test_hsa.calculate_site_quantities(start_frame=0, num_frames=100)
    test_hsa.write_calculation_summary()
    #hsa_arg.write_data()

    #hsa_arg.initialize_hydration_sites()
    #hsa_arg.calculate_site_quantities(start_frame=0, num_frames=100)
    #hsa_arg.write_calculation_summary()
    #hsa_arg.write_data()
    #print("Test 2")
    #hsa_casp = read_casp3_desmond_without_clusters()   
    #test_write_wat_pdb(hsa_casp)
    
