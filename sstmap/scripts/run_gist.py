from argparse import ArgumentParser
from sstmap.water_analysis import requirements
from sstmap.grid_water_analysis import GridWaterAnalysis
import os
import sys

#param_message = ""
#for k in requirements.keys():



def parse_args():
    """Parse the command-line arguments and check if input args are valid.

    Returns
    -------
    args : argparse.Namespace
        The namespace containing the arguments
    """
    parser = ArgumentParser(description='''Run SSTMap grid-based (GIST) calculations through command-line.''')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input_top', required=True, type=str, default=None,
                          help='''Input toplogy File.''')
    required.add_argument('-t', '--input_traj', required=True, type=str, default=None,
                          help='''Input trajectory file.''')
    required.add_argument('-l', '--ligand', required=True, type=str, default=None,
                          help='''Input ligand PDB file.''')
    required.add_argument('-g', '--grid_dim', required=True, nargs=3, type=float, default=[20.0, 20.0, 20.0],
                          help='''grid dimensions e.g., 10 10 10''')
    required.add_argument('-f', '--num_frames', required=False, type=int, default=10000,
                          help='''Total number of frames to process.''')
    parser._action_groups.append(parser._action_groups.pop(1))
    parser.add_argument('-p', '--param_file', required=False, type=str, default=None,
                          help='''Additional parameter files, specific for MD package''')
    parser.add_argument('-s', '--start_frame', required=False, type=int, default=0,
                          help='''Starting frame.''')
    parser.add_argument('-d', '--bulk_density', required=False, type=float, default=0.0334,
                        help='''Bulk density of the water model.''')
    parser.add_argument('-b', '--calc_hbonds', required=False, type=bool, default=False,
                        help='''True or False for whether to calculate h-bonds during calculations.''')
    parser.add_argument('-o', '--output_prefix', required=False, type=str,
                          help='''Prefix for all the results files.''', default="gist")
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    file_arguments = [args.input_top, args.input_traj, args.ligand]
    files_present = [os.path.isfile(f) for f in file_arguments]
    for index, present in enumerate(files_present):
        if not present:
            sys.exit("%s not found. Please make sure it exits or give the correct path." % file_arguments[index])
    if args.param_file is not None:
       if not os.path.exists(args.param_file):# or not os.path.isdir(args.param_file):
             sys.exit("%s not found. Please make sure it exits or give the correct path." % args.param_file)
    return args

def main():
    args = parse_args()
    g = GridWaterAnalysis(args.input_top, args.input_traj,
                          start_frame=args.start_frame, num_frames=args.num_frames,
                          ligand_file=args.ligand, supporting_file=args.param_file,
                          grid_dimensions=args.grid_dim,
                          rho_bulk=args.bulk_density, prefix=args.output_prefix)
    g.print_system_summary()
    g.calculate_grid_quantities(hbonds=args.calc_hbonds)
    g.print_calcs_summary()
    g.write_data()
    g.generate_dx_files()

def entry_point():
    main()

if __name__ == '__main__':
    entry_point()
