import sys
from sstmap.utils import read_hsa_summary

hsa_quantities = {}
gist_quantities = {}
hsa_header ='index x y z nwat occupancy \
            Esw EswLJ EswElec Eww EwwLJ \
            EwwElec Etot Ewwnbr TSsw TSww \
            TStot Nnbrs Nhbww Nhbsw Nhbtot \
            f_hb_ww f_enc Acc_ww Don_ww Acc_sw \
            Don_sw solute_acceptors solute_donors'
gist_header = 'index x y z nwat gO gH \
            dTStr_dens dTStr_norm dTSor_dens dTSor_norm \
            Esw_dens Esw_norm Eww_dens Eww_norm Eww_nbr_dens \
            Eww_nbr_norm Nnnbr_dens Nnbr_norm fHB_dens \
            fHB_norm fenc_dens fenc_norm Nhbsw_dens Nhbsw_norm \
            Nhbww_dens Nhbww_norm donsw_dens donsw_norm accsw_dens \
            accsw_norm donww_dens donww_norm accww_dens accww_norm'
for i, q in enumerate(hsa_header.split()):
    hsa_quantities[q] = i
for i, q in enumerate(gist_header.split()):
    gist_quantities[q] = i




def compare_hsa(test_results, known_results):
    test_hsa = read_hsa_summary(test_results)
    known_hsa = read_hsa_summary(known_results)
    print known_hsa
    single_site_comparison(test_hsa[0], known_hsa[0], 'Ewwnbr')

def single_site_comparison(test_site_data, known_site_data, quantity):
    pass

def single_voxel_comparison(test_site_data, known_site_data, quantity):
    pass

def main():
    compare_hsa(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()