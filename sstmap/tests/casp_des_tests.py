from __future__ import absolute_import
from site_water_analysis import SiteWaterAnalysis

test_dir = "/Users/kamranhaider/sstmap_testcases/casp3_desmond/simulation/casp3_pd_2/"
top = test_dir + "casp3_pd_2-out.pdb"
traj = test_dir + "casp3_converted.nc"
lig = test_dir + "ligand.pdb"
clusters = test_dir + "test_casp_clustercenterfile.pdb"
des = test_dir + "casp3_pd_2-out_cms_nb_parms.txt"
# Test 1
# create SiteWaterAnalysis object and print summary
# includes clustering
#hsa = SiteWaterAnalysis(top, traj, 
#                        start_frame=0, num_frames=10000, 
#                        ligand_file=lig, desmond_helper_file=des,
#                        prefix=test_dir+"test_casp")
#hsa.print_system_summary()

# Test 2
# create SiteWaterAnalysis object and print summary
# skip clustering
# hsa = SiteWaterAnalysis(top, traj, 
#                        start_frame=0, num_frames=10000, 
#                        cluster_center_file=clusters, desmond_helper_file=des,
#                        prefix="test_casp")
#hsa.print_system_summary()

#hsa.calculate_site_quantities()
#hsa.write_summary()
#hsa.write_data()
#hsa.calculate_angular_structure()
#hsa.calculate_lonranged_ww_energy([1, 3])
# hsa.write_summary()
# hsa.write_data()
# Test 3
#clusters_old = test_dir + "clustercenterfile.pdb"
#hsa = SiteWaterAnalysis(top, traj, 
#                        start_frame=2000, num_frames=10000, 
#                        cluster_center_file=clusters_old, desmond_helper_file=des,
#                        prefix="old_casp")
#hsa.print_system_summary()

#hsa.calculate_site_quantities()
#hsa.write_summary()
#hsa.write_data()

# Test 3
single_clusters = test_dir + "single_clustercenterfile.pdb"
hsa = SiteWaterAnalysis(top, traj,
                        start_frame=2000, num_frames=100,
                        cluster_center_file=single_clusters,
                        desmond_helper_file=des,
                        prefix="single_casp")
#hsa.print_system_summary()

hsa.calculate_site_quantities()
hsa.write_summary()
#hsa.write_data()
