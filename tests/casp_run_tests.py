from __future__ import absolute_import
from site_water_analysis import SiteWaterAnalysis

if test_1:
    hsa = SiteWaterAnalysis("/Users/kamranhaider/arg_gist_calcs/mdtraj/arg.prmtop", "/Users/kamranhaider/arg_gist_calcs/mdtraj/md10ns.ncdf",
                            start_frame=0, num_frames=1000, ligand_file="/Users/kamranhaider/arg_gist_calcs/mdtraj/ligand_arg.pdb", prefix="arg")
    # hsa.print_system_summary()
    # hsa.calculate_site_quantities()
    #hsa.calculate_angular_structure([1, 3])
    #hsa.calculate_lonranged_ww_energy([1, 3])
    # hsa.write_summary()
    # hsa.write_data()
else:
    hsa = SiteWaterAnalysis("/Users/kamranhaider/arg_gist_calcs/mdtraj/arg.prmtop", "/Users/kamranhaider/arg_gist_calcs/mdtraj/md10ns.ncdf",
                            start_frame=0, num_frames=100, cluster_center_file="/Users/kamranhaider/arg_gist_calcs/03_organizewaters/clustercenterfile.pdb", prefix="test_2")
    # hsa.print_system_summary()
    #hsa.calculate_angular_structure([1, 3])
    #hsa.calculate_lonranged_ww_energy([1, 3])
    # hsa.calculate_site_quantities()
    # hsa.write_summary()
    # hsa.write_data()
#top = "/Users/kamranhaider/simulation_data_from_microway_home/casp3/casp3_pd_2/casp3_pd_2-out_converted.gro"
top = "/Users/kamranhaider/simulation_data_from_microway_home/casp3/casp3_pd_2/casp3_pd_2-out.pdb"
#traj = "/Users/kamranhaider/simulation_data_from_microway_home/casp3/casp3_pd_2/casp3_converted.nc"
traj = "/Users/kamranhaider/water_structure_data_davinci/casp3/clustering/01_wkeeper/casp3_converted.nc"
des = "/Users/kamranhaider/Downloads/casp3_pd_2-out_cms_nb_parms.txt"
lig = "/Users/kamranhaider/water_structure_data_davinci/casp3/clustering/01_wkeeper/ligand.pdb"
clusters = "/Users/kamranhaider/water_structure_data_davinci/casp3/clustering/03_organizewaters/clustercenterfile.pdb"
#hsa = SiteWaterAnalysis(top, traj, start_frame=0, num_frames=100, cluster_center_file=clusters, desmond_helper_file=des, prefix="test_casp")
#hsa = SiteWaterAnalysis(top, traj, start_frame=0, num_frames=10000, ligand_file=lig, desmond_helper_file=des, prefix="casp")
# hsa.print_system_summary()
# hsa.calculate_site_quantities()
# hsa.write_summary()
# hsa.write_data()
