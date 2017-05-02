import _sstmap_entropy as e1
import _sstmap_probableconfig as e2
#b.bruteclust_run("/Users/kamranhaider/Dropbox/SSTMap/sstmap/tests/clustercenterfile.pdb", "/Users/kamranhaider/Dropbox/SSTMap/sstmap/tests/within5Aofligand.pdb")


e1.run_bruteclust("/Users/kamranhaider/Dropbox/SSTMap/sstmap/test_build/clustercenterfile.pdb", "/Users/kamranhaider/Dropbox/SSTMap/sstmap/test_build/within5Aofligand.pdb")
e1.run_kdhsa102("/Users/kamranhaider/Dropbox/SSTMap/sstmap/test_build/cluster.000001.pdb", "cluster.000001.pdb")
e2.run_6dimprob("/Users/kamranhaider/Dropbox/SSTMap/sstmap/test_build/cluster.000001.pdb")