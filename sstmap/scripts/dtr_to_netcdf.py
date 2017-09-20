from argparse import ArgumentParser
import mdtraj as md


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

    parser.add_argument('-i', '--input_parm', required=False, type=str,
                        help='''Input toplogy File.''')
    parser.add_argument('-t', '--input_traj', required=True, type=str,
                        help='''Input trajectory file.''')
    parser.add_argument('-o', '--output_prefix', required=False, type=str,
                        help='''Prefix for all the results files.''')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    print("Reading in trajectory ...")
    traj = md.load_dtr(args.input_traj, top=args.input_parm)
    print traj
    print("Outputting NETCDF ...")
    traj.save_netcdf(args.output_prefix + "_converted.nc")
    print("Outputting PDB file of frame 1 ...")
    traj[0].save_pdb(args.output_prefix + "_converted.pdb")
    print("Done")


def entry_point():
    main()


if __name__ == '__main__':
    entry_point()
