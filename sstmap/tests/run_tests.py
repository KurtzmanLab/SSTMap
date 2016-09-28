import sys

"""
Run some test calculations to check if each class is implemented correctly.
"""


def testSiteBased(topology, trajectory, ligand,
                  n_frame, s_frame, sites, prefix):
    try:
        from wateranalysistools import SiteWaterAnalysis
    except ImportError:
        raise ImportError(
            "Could not import modules, please check installation")
    else:
        h = SiteWaterAnalysis(
            topology,
            trajectory,
            start_frame=s_frame,
            num_frames=n_frame,
            ligand_file=ligand,
            cluster_center_file=sites,
            prefix=prefix)
        # h.generate_clusters()
        h.calculate_site_quantities()
        h.write_summary()
        h.write_data()


def testProteinBased(topology, trajectory, ligand,
                     n_frame, s_frame, sites, prefix):
    try:
        from wateranalysistools import ProteinWaterAnalysis
    except ImportError:
        raise ImportError(
            "Could not import modules, please check installation")
    else:
        p = ProteinWaterAnalysis(
            topology,
            trajectory,
            start_frame=s_frame,
            num_frames=n_frame,
            ligand_file=sites,
            prefix=prefix,
            near_ligand=False)
        p.calculate_protein_water_hbonds()
        p.writeHBsummary()


def testGridBased(topology, trajectory, ligand,
                  n_frame, s_frame, sites, prefix):
    try:
        from wateranalysistools import GridWaterAnalysis
    except ImportError:
        raise ImportError(
            "Could not import modules, please check installation")
    else:
        g = GridWaterAnalysis(
            topology,
            trajectory,
            start_frame=s_frame,
            num_frames=n_frame,
            ligand_file=sites,
            prefix=prefix,
            near_ligand=False)


def main():
    topology = "testsystems/arg_amber_100ps/arg.prmtop"
    trajectory = "testsystems/arg_amber_100ps/arg_100ps.nc"
    hydration_sites = "testsystems/arg_amber_100ps/clustercenterfile.pdb"
    ligand = "testsystems/arg_amber_100ps/ligand.pdb"
    prefix = "arg_test"
    n_frame = 100
    s_frame = 0

    print "Testing SiteBased water analysis for Arg in TIP3P, using 100 ps simulation."
    #testSiteBased(topology, trajectory, ligand, n_frame, s_frame, hydration_sites, prefix)
    print "Skipping"
    print "Testing ProteinBased water analysis for Arg in TIP3P, using 100 ps simulation."
    #testProteinBased(topology, trajectory, ligand, n_frame, s_frame, hydration_sites, prefix)
    print "Skipping."
    print "Testing GridBased water analysis for Arg in TIP3P, using 100 ps simulation."
    testSiteBased(topology, trajectory, ligand, n_frame, s_frame, prefix)

if __name__ == "__main__":
    sys.exit(main())
