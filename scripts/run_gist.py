from argparse import ArgumentParser
from sstmap.grid_water_analysis import GridWaterAnalysis


def parse_args():
    """Parse the command line arguments and perform some validation on the
    arguments
    Returns
    -------
    args : argparse.Namespace
        The namespace containing the arguments
    """
    parser = ArgumentParser(
        description='''Run GIST calculations through command-line.''')
    required = parser.add_argument_group('required arguments')

    required.add_argument('-i', '--input_parm', required=True, type=str,
                          help='''Input toplogy File.''')
    required.add_argument('-t', '--input_traj', required=True, type=str,
                          help='''Input trajectory file.''')
    required.add_argument('-l', '--ligand', required=True, type=str,
                          help='''Input ligand PDB file.''')
    required.add_argument('-g', '--grid_dim', required=True, nargs=3, type=float,
                          help='''grid dimensions''')
    required.add_argument('-f', '--num_frames', required=False, type=int,
                          help='''Total number of frames to process.''')
    parser.add_argument('-s', '--start_frame', required=False, type=int,
                        help='''Starting frame.''')
    parser.add_argument('-o', '--output_prefix', required=False, type=str,
                        help='''Prefix for all the results files.''')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    g = GridWaterAnalysis(args.input_parm, args.input_traj,
                          start_frame=args.start_frame, num_frames=args.num_frames,
                          ligand_file=args.ligand, prefix=args.output_prefix,
                          grid_dimensions=args.grid_dim)
    g.print_system_summary()
    g.calculate_grid_quantities()
    g.print_calcs_summary()
    g.write_data()
    g.generate_dx_files()

def entry_point():
    main()

if __name__ == '__main__':
    entry_point()
