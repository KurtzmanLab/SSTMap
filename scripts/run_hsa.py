from argparse import ArgumentParser
from sstmap.site_water_analysis import SiteWaterAnalysis


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

    parser.add_argument('-i', '--input_parm', required=True, type=str,
                        help='''Input toplogy File.''')
    parser.add_argument('-t', '--input_traj', required=True, type=str,
                        help='''Input trajectory file.''')
    parser.add_argument('-c', '--clusters', required=True, type=str,
                        help='''PDB file containing cluster centers.''')
    parser.add_argument('-f', '--num_frames', required=False, type=int,
                        help='''Total number of frames to process.''')
    parser.add_argument('-s', '--start_frame', required=False, type=int,
                        help='''Starting frame.''')
    parser.add_argument('-o', '--output_prefix', required=False, type=str,
                        help='''Prefix for all the results files.''')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    h = SiteWaterAnalysis(args.input_parm, args.input_traj,
                          start_frame=args.start_frame, num_frames=args.num_frames,
                          cluster_center_file=args.clusters, prefix=args.output_prefix)
    h.initialize_hydration_sites()
    h.print_system_summary()
    h.calculate_site_quantities()
    h.write_calculation_summary()
    h.write_data()


def entry_point():
    main()

if __name__ == '__main__':
    entry_point()
