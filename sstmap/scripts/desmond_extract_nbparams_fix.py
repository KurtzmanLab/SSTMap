import sys
import numpy as np
from schrodinger.application.desmond.cms import Vdw
import schrodinger.application.desmond.ffiostructure as ffiostructure
from schrodinger import structure


def write_nonbonded_parameters(cms_file):
    vdw = []  # combined list of all vdw params from all ct's
    chg = []  # combined list of all charges from all ct's
    ct_list = [e for e in ffiostructure.CMSReader(cms_file)]
    # this means this works only on cms files with separate CT blocks
    struct_ct = ct_list[1:]
    # get total number of solute atoms
    for ct in struct_ct:
        # print ct.atom_total, len(ct.ffio.site), len(ct.ffio.pseudo)
        # else: # do the normal parsing
        ct_vdw = []  # list of vdw objects for each ct
        ct_chg = []
        ct_pseudo_vdw = []
        ct_pseudo_chg = []
        n_atomic_sites = 0
        n_pseudo_sites = 0
        vdw_type = {}  # dict of vdw types, Vdw object in this list are uninitialized
        for e in ct.ffio.vdwtype:
            vdw_type[e.name] = Vdw((e.name,), e.funct, (e.c1, e.c2,))
            #print (e.name,), e.funct, (e.c1, e.c2,)

        # for each site (i.e., an atom in most cases)
        for e in ct.ffio.site:
            # print e.type.lower()
            if e.type.lower() == 'pseudo':
                # add to vdw list for this ct
                ct_pseudo_vdw.append(vdw_type[e.vdwtype])
                ct_pseudo_chg.append(e.charge)
                n_pseudo_sites += 1
            else:
                # add to vdw list for this ct
                ct_vdw.append(vdw_type[e.vdwtype])
                ct_chg.append(e.charge)
                n_atomic_sites += 1

            # print e.index, e.charge
                # check if this site belongs to pseudoparticle, if yes raise corresponding number
        # print n_atomic_sites, n_pseudo_sites
        ct_vdw *= int(ct.atom_total / n_atomic_sites)
        ct_chg *= int(ct.atom_total / n_atomic_sites)
        vdw.extend(ct_vdw)
        chg.extend(ct_chg)
        if n_pseudo_sites != 0:
            ct_pseudo_vdw *= int(len(ct.ffio.pseudo) / n_pseudo_sites)
            ct_pseudo_chg *= int(len(ct.ffio.pseudo) / n_pseudo_sites)
            # print int(ct.atom_total / len( ct.ffio.site ))
            vdw.extend(ct_pseudo_vdw)
            chg.extend(ct_pseudo_chg)

    chg = np.asarray(chg) * 18.2223
    vdw_params = []
    for v in vdw:
        vdw_params.extend([v.c])
    vdw_params = np.asarray(vdw_params)
    with open(cms_file[0:-4] + "_cms_nb_parms.txt", "w") as f:
        for i in range(vdw_params.shape[0]):
            f.write("{0:.8f} {1:.8f} {2:.8f}\n".format(
                chg[i], vdw_params[i, 0], vdw_params[i, 1]))
    # return (chg, vdw_params)


if (__name__ == '__main__'):
    cms_file = sys.argv[1]
    write_nonbonded_parameters(cms_file)
