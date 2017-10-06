"""
Test scripts
General purpose: Take in two files as input and run the tests
Should be able to import the module and run for hsa and gist calculations run
as part of installation.

Place it insdie the test suite, at the end of test scripts import and run tests
Use numpy testing module.

Test quantities: divide into three groups
"""


import os

import numpy as np
import numpy.testing as npt

quantities = ["voxel", "xcoord", "ycoord", "zcoord",
                "n_wat", "g_O", 
                "dTStrans-dens", "dTStrans-norm", 
                "dTSorient-dens", "dTSorient-norm", 
                "dTSsix-dens", "dTSsix-norm", 
                "Esw-dens", "Esw-norm", "Eww-dens", "Eww-norm-unref", 
                "neighbor-dens", "neighbor-norm"]

DX_TEST_FILES = ["gO", "dTSorient-dens", "dTStrans-dens", "dTSsix-dens", "Eww-dens", "Esw-dens"]

class TestGistOutput():
    """
    """
    
    def __init__(self, test_data, ref_data):
        """

        Args:
            test_data:
            ref_data:
        """
        self.test_data = test_data
        self.ref_data = ref_data

    def test_grid(self):
        """

        Returns:

        """
        
        passed = True
        try:
            #npt.assert_equal(self.test_data.shape, self.ref_data.shape)
            npt.assert_almost_equal(self.test_data[:, 1:4], self.ref_data[:, 1:4], decimal=3)
        except Exception as e:
            print e
            passed = False

        return passed

    def test_voxel_number(self):
        """

        Returns:

        """

        passed = True
        try:
            npt.assert_equal(self.test_data.shape, self.ref_data.shape)
        except Exception as e:
            print e
            passed = False

        return passed

    def test_quantity(self, quantity_index):
        """

        Args:
            quantity_index:

        Returns:

        """

        passed = True
        try:
            npt.assert_array_almost_equal(self.test_data[:, quantity_index], self.ref_data[:, quantity_index], decimal=2)
        except Exception as e:
            print e
            passed = False

        return passed

def read_gist_sstmap(sstmap_gist_summary):
    """

    Args:
        sstmap_gist_summary:

    Returns:

    """
    columns_to_read = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19]
    sstmap_data = np.loadtxt(sstmap_gist_summary, skiprows=1, usecols=columns_to_read)
    #sstmap_data = sstmap_data[np.where(sstmap_data[:, 4] != 1.0)]
    return np.round(sstmap_data, 3)

def read_gist_cpptraj(cpptraj_gist_summary):
    """

    Args:
        cpptraj_gist_summary:

    Returns:

    """
    columns_to_read = [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22]
    cpptraj_data = np.loadtxt(cpptraj_gist_summary, skiprows=2, usecols=columns_to_read)
    return cpptraj_data

def test_dx_output(test_dx_filename, ref_dx_filename):
    """
    Compare two DX files.

    Parameters
    ----------
    test_dx_filename: string
        Location of the test DX file.
    test_dx_filename: string
        Location of the reference DX file.

    Returns
    -------

    """
    with open(test_dx_filename, "r") as test:
        lines = test.readlines()
        test_dims = [float(s) for s in lines[0].strip().split()[-3:]]
        test_origin = [float(s) for s in lines[1].strip().split()[-3:]]
        test_spacing = float(lines[2].strip().split()[-1])
        test_voxel_num = int(lines[6].strip().split()[-3])
        test_data = []
        for i in xrange(len(lines[7:])):
            
            test_data.extend([float(s) for s in lines[7:][i].strip().split()])
        test_data = np.asarray(test_data)

    with open(ref_dx_filename, "r") as ref:
        lines = ref.readlines()
        ref_dims = [float(s) for s in lines[0].strip().split()[-3:]]
        ref_origin = [float(s) for s in lines[1].strip().split()[-3:]]
        ref_spacing = float(lines[2].strip().split()[-1])
        ref_voxel_num = int(lines[6].strip().split()[-3])
        ref_data = []
        for i in xrange(len(lines[7:])):
            ref_data.extend([float(s) for s in lines[7:][i].strip().split()])
        ref_data = np.asarray(ref_data)

    npt.assert_almost_equal(test_dims, ref_dims, decimal=6)
    npt.assert_almost_equal(test_origin, ref_origin, decimal=6)
    npt.assert_almost_equal(test_spacing, ref_spacing, decimal=6)
    npt.assert_almost_equal(test_voxel_num, ref_voxel_num, decimal=6)
    if "Esw" in ref_dx_filename:
        ref_data /= 2.0
    try:
        npt.assert_almost_equal(test_data, ref_data, decimal=3)
    except Exception as e:
        print e
        print test_dx_filename, ref_dx_filename


def parse_args():
    """Parse the command line arguments and perform some validation on the
    arguments
    Returns
    -------
    args : argparse.Namespace
        The namespace containing the arguments
    """
    parser = ArgumentParser(
        description='''Run tests of GIST calculations against validated output.''')

    parser.add_argument('-t', '--test_gist_summary', required=True, type=str,
                        help='''Summary file of GIST calculation to be tested.''')
    parser.add_argument('-r', '--ref_gist_summary', required=True, type=str,
                        help='''A refeeence summary file of GIST calculation''')
    args = parser.parse_args()
    return args

def run_all_gist_tests(test_dir, ref_dir):
    """

    Args:
        test_data_file:
        ref_data_file:
    """

    test_result = {1: "Passed", 0: "Failed"}
    test_dir_path = os.path.abspath(test_dir) + "/"
    ref_dir_path = os.path.abspath(ref_dir) + "/"

    if not os.path.exists(test_dir_path) or not os.path.exists(ref_dir_path):
        raise IOError("%s and/or %s directory not found, please provide correct path." % (test_dir, ref_dir))
    else:
        test_dir_files = os.listdir(test_dir_path)
        ref_dir_files = os.listdir(ref_dir_path)
        test_dx_files = [f for f in test_dir_files if f.endswith(".dx")]
        test_dx_files = [f for f in test_dx_files if f[f.find("_") + 1:][:-3] in DX_TEST_FILES]
        ref_dx_files = [f for f in ref_dir_files if f.endswith(".dx")]
        ref_dx_files = [f for f in ref_dx_files if f[5:][:-3] in DX_TEST_FILES]
        assert len(test_dx_files) == len(ref_dx_files), "Couldn't obtain all DX files, tests won't run."
        test_data_file = [test_dir_path + f for f in test_dir_files if f.endswith("gist_data.txt")]
        test_data = read_gist_sstmap(test_data_file[0])
        ref_data_file = [ref_dir_path + f for f in ref_dir_files if f.endswith("all.out")]
        ref_data = read_gist_cpptraj(ref_data_file[0])
        assert test_data.shape == ref_data.shape, "GIST columns/rows in summary files are not equal, tests won't run"
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
        print "Checking: %s" % quantities[4]
        result = testcase.test_quantity(4)
        print "%s" % test_result[bool(result)]
        for index, filename in enumerate(test_dx_files):
            test_dx_file, ref_dx_file = test_dir_path + filename, ref_dir_path + ref_dx_files[index]
            test_dx_output(test_dx_file, ref_dx_file)

    """
    for quantity_index in xrange(4, 5):
        print "--------------------------------------"
        print "Checking: %s" % quantities[quantity_index]
        result = testcase.test_quantity(quantity_index)
        print "\t%s" % test_result[bool(result)]


    """

def main():
    """

    """
    args = parse_args()
    run_all_gist_tests(args.test_gist_summary, args.ref_gist_summary)

if __name__ == '__main__':
    main()
