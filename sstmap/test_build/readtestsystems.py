import sys
import os
#sys.path.append("/Users/kamranhaider/Dropbox/SSTMap/sstmap")

from sstmap.grid_water_analysis import GridWaterAnalysis
from sstmap.site_water_analysis import SiteWaterAnalysis

def read_tip3p():
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

def read_casp3_desmond_hsa(clusters=None):
    """
    Can sstmap read a pure water tip3p system in amber
    """
    top = testsystem_dir + "casp3_desmond/casp3.pdb"
    traj = "/Users/kamranhaider/sstmap_testcases/casp3_desmond/simulation/casp3_pd_2/casp3_converted.nc"
    nb_param_file = testsystem_dir + "casp3_desmond/casp3_nb_parms.txt"
    ligand = testsystem_dir + "casp3_desmond/ligand.pdb"
    hsa = SiteWaterAnalysis(top, traj, start_frame=0, num_frames=10000, 
                            ligand_file=ligand,
                            clustercenter_file=clusters, 
                            desmond_helper_file=nb_param_file,
                            prefix="casp")
    return hsa

def read_casp3_desmond_gist():
    """
    Can sstmap read a pure water tip3p system in amber
    """
    top = testsystem_dir + "casp3_desmond/casp3.pdb"
    traj = "/Users/kamranhaider/sstmap_testcases/casp3_desmond/simulation/casp3_pd_2/casp3_converted.nc"
    nb_param_file = testsystem_dir + "casp3_desmond/casp3_nb_parms.txt"
    ligand = testsystem_dir + "casp3_desmond/ligand.pdb"
    hsa = SiteWaterAnalysis(top, traj, start_frame=0, num_frames=10000, 
                            ligand_file=ligand,
                            clustercenter_file=clusters, 
                            desmond_helper_file=nb_param_file,
                            prefix="casp")
    return hsa

def read_arg_amber_hsa(clusters=None):
    """Read Arg amber test system
    """
    top = testsystem_dir + "arg_amber/arg.prmtop"
    #traj = testsystem_dir + "arg_amber/md1ns.nc"
    traj = testsystem_dir + "arg_amber/md10ns.nc"

    ligand = testsystem_dir + "arg_amber/ligand.pdb"
    hsa = SiteWaterAnalysis(top, traj, start_frame=0, num_frames=10000, 
                            clustercenter_file=clusters, ligand_file=ligand,
                            prefix="arg")
    return hsa

def read_arg_amber_gist():
    """Read Arg amber test system
    """
    top = testsystem_dir + "arg_amber/arg.prmtop"
    traj = testsystem_dir + "arg_amber/md1ns.nc"
    ligand = testsystem_dir + "arg_amber/ligand.pdb"
    gist = GridWaterAnalysis(top, traj, start_frame=0, num_frames=1000, 
                            ligand_file=ligand,
                            grid_dimensions=[20.0, 20.0, 20.0], grid_center=[22.27, 20.46, 19.05],
                            prefix="arg")
    return gist

