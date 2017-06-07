"""
Test scripts
General purpose: Take in two files as input and run the tests
Should be able to import the module and run for hsa and gist calculations run
as part of installation.

Place it insdie the test suite, at the end of test scripts import and run tests
Use numpy testing module.

Test quantities: divide into three groups
"""


import sys
import numpy as np
import numpy.testing as npt
quantities = ["voxel", "xcoord", "ycoord", "zcoord",
                "n_wat", "g_O", "dTStrans-dens", "dTStrans-norm", "dTSorient-dens", 
                "dTSorient-norm", "Esw-dens", "Esw-norm", "Eww-dens", "Eww-norm-unref", 
                "neighbor-dens", "neighbor-norm"]
class TestGistOutput():
    
    def __init__(self, test_data, ref_data):
        self.test_data = test_data
        self.ref_data = ref_data

    def test_grid(self):
        
        passed = True
        try:
            #npt.assert_equal(self.test_data.shape, self.ref_data.shape)
            npt.assert_almost_equal(self.test_data[:, 1:4], self.ref_data[:, 1:4], decimal=3)
        except Exception as e:
            print e
            passed = False

        return passed

    def test_voxel_number(self):

        passed = True
        try:
            npt.assert_equal(self.test_data.shape, self.ref_data.shape)
        except Exception as e:
            print e
            passed = False

        return passed

    def test_quantity(self, quantity_index):

        passed = True
        try:
            npt.assert_array_almost_equal(self.test_data[:, quantity_index], self.ref_data[:, quantity_index], decimal=2)
        except Exception as e:
            print e
            passed = False

        return passed


def read_hsa_sstmap(sstmap_gist_summary):
    columns_to_read = [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 17, 18]
    sstmap_data = np.loadtxt(sstmap_gist_summary, skiprows=1, usecols=columns_to_read)
    #sstmap_data = sstmap_data[np.where(sstmap_data[:, 4] != 1.0)]
    return np.round(sstmap_data, 3)

def parse_args():
    """Parse the command line arguments and perform some validation on the
    arguments
    Returns
    -------
    args : argparse.Namespace
        The namespace containing the arguments
    """
    parser = ArgumentParser(
        description='''Run tests of HSA calculations against validated output.''')

    parser.add_argument('-t', '--test_hsa_summary', required=True, type=str,
                        help='''Summary file of HSA calculation to be tested.''')
    parser.add_argument('-r', '--ref_hsa_summary', required=True, type=str,
                        help='''A reference summary file of HSA calculation''')
    args = parser.parse_args()
    return args

def run_all_gist_tests(test_data_file, ref_data_file):

    test_result = {1: "Passed", 0: "Failed"}
    test_data = read_gist_sstmap(test_data_file)
    ref_data = read_gist_cpptraj(ref_data_file)
    diff_nwat = []
    for row in xrange(test_data.shape[0]):
        if test_data[row, 4] <= 1:
            test_data[row, 6:14] *= 0.0
        # record voxels with different water number but exclude them for tests
        else:
            if abs(int(test_data[row, 4]) - int(ref_data[row, 4])) >= 1:
                diff_nwat.append([test_data[row, :], ref_data[row, :]])
                test_data[row, 4:] *= 0.0
                ref_data[row, 4:] *= 0.0
    # Run tests
    print "Checking grid and voxel placement ...",
    testcase = TestGistOutput(test_data, ref_data)
    result = testcase.test_voxel_number()
    print test_result[bool(result)]
    result = testcase.test_grid()
    print test_result[bool(result)]
    print "Coverage for remaining tests.", 100*float(test_data.shape[0] - len(diff_nwat))/test_data.shape[0]
    print "Checking quantities ..."
    test_num = 0
    
    for quantity_index in xrange(4, test_data.shape[1]):
        result = testcase.test_quantity(quantity_index)
        print "Testing %s ... %s" % (quantities[quantity_index], test_result[bool(result)])
        

def main():
    args = parse_args()
    run_all_gist_tests(args.test_gist_summary, args.ref_gist_summary)

if __name__ == '__main__':
    main()
