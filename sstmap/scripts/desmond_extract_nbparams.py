import sys
import numpy
from schrodinger.application.desmond.cms import Vdw
import schrodinger.application.desmond.ffiostructure as ffiostructure
from schrodinger import structure


def write_nonbonded_parameters(cms_file):
    # obtain LJ and Elec params
    vdw = []  # combined list of all vdw params from all ct's
    chg = []  # combined list of all charges from all ct's
    ct_list = [e for e in ffiostructure.CMSReader(cms_file)]
    # this means this works only on cms files with separate CT blocks
    struct_ct = ct_list[1:]
    # get tota number of solute atoms

    for ct in struct_ct:
        ct_chg = []
        ct_vdw = []
        # dict of vdw types, Vdw object in this list are uninitialized
        vdw_type = {}
        for e in ct.ffio.vdwtype:
            vdw_type[e.name] = Vdw((e.name,), e.funct, (e.c1, e.c2,))
            #print (e.name,), e.funct, (e.c1, e.c2,)
        for e in ct.ffio.site:  # for each site (i.e., an atom in most cases)
            ct_vdw.append(vdw_type[e.vdwtype])  # add to vdw list for this ct
            ct_chg.append(e.charge)
            # print e.index, e.charge
        ct_vdw *= int(ct.atom_total / len(ct.ffio.site))
        ct_chg *= int(ct.atom_total / len(ct.ffio.site))
        # print int(ct.atom_total / len( ct.ffio.site ))
        vdw.extend(ct_vdw)
        chg.extend(ct_chg)

    chg = numpy.asarray(chg) * 18.2223
    vdw_params = []
    # print len(chg)
    # print len(all_at_ids)
    for v in vdw:
        vdw_params.extend([v.c])
    vdw_params = numpy.asarray(vdw_params)
    # print vdw_params.shape
    # print chg.shape
    with open(cms_file[0:-4] + "_cms_nb_parms.txt", "w") as f:
        for i in range(vdw_params.shape[0]):
            f.write("{0:.8f} {1:.8f} {2:.8f}\n".format(
                chg[i], vdw_params[i, 0], vdw_params[i, 1]))
    # return (chg, vdw_params)


if (__name__ == '__main__'):
    cms_file = sys.argv[1]
    write_nonbonded_parameters(cms_file)
