from argparse import ArgumentParser
from sstmap.core.grid_water_analysis import GridWaterAnalysis
def parse_args():
    """Parse the command-line arguments and check if input args are valid.

    Returns
    -------
    args : argparse.Namespace
        The namespace containing the arguments
    """
    parser = ArgumentParser(description='''Run grid-based SSTMap calculations through command-line.''')
    required = parser.add_argument_group('required arguments')    
    required.add_argument('-i', '--input_parm', required=True, type=str,
                          help='''Input toplogy File.''')
    required.add_argument('-t', '--input_traj', required=True, type=str,
                          help='''Input trajectory file.''')
    required.add_argument('-l', '--ligand', required=True, type=str,
                          help='''Input ligand PDB file.''')
    required.add_argument('-g', '--grid_dim', required=True, nargs = 3, type=float,
                          help='''grid dimensions''')
    parser._action_groups.append(parser._action_groups.pop(1))
    parser.add_argument('-f', '--num_frames', type=int,
                          help='''Total number of frames to process.''')
    parser.add_argument('-s', '--start_frame', type=int,
                          help='''Starting frame.''')
    parser.add_argument('-o', '--output_prefix', type=str,
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
    g.process_grid()
    g.print_calcs_summary()
    g.write_data(args.output_prefix)
    g.generate_dx_files(args.output_prefix)

def entry_point():
    main()

if __name__ == '__main__':
    entry_point()
