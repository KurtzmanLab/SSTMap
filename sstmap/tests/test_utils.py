
import mdtraj as md
from sstmap.site_water_analysis import SiteWaterAnalysis
from readtestsystems import *
from sstmap.utils import *

def test_write_wat_pdb(hsa):
    """
    Test write water PDB function
    """
    ligand = md.load_pdb(hsa.ligand)
    ligand_coords = ligand.xyz[0, :, :]
    binding_site_atom_indices = list(range(ligand_coords.shape[0]))

    trj = md.load(hsa.trajectory, top=hsa.topology)
    for i_frame in range(trj.n_frames):
        for pseudo_index in range(ligand_coords.shape[0]):
            trj.xyz[i_frame, pseudo_index,
                    :] = ligand_coords[pseudo_index, :]

    binding_site_waters = md.compute_neighbors(
        trj, 0.50, binding_site_atom_indices,
        haystack_indices=hsa.wat_oxygen_atom_ids)
    water_id_frame_list = [(i, nbr) for i in
                           range(len(binding_site_waters))
                           for nbr in binding_site_waters[i]]
    print("Writing all waters within 5 Angstrom of ligand to pdb file.")
    write_watpdb_from_list(trj, hsa.prefix + "_within5Aofligand",
                        water_id_frame_list, full_water_res=True)

def main():
    print("Test 1")
    hsa_arg = read_arg_amber_hsa(clusters="arg_amber/ref_clustercenterfile.pdb")   
    test_write_wat_pdb(hsa_arg)
    
if __name__ == '__main__':
    status = main()
    sys.exit(status)