/*
  C helper module for intensive calcs.
#=============================================================================================
# VERSION CONTROL INFORMATION
#=============================================================================================
__version__ = "$Revision: $ $Date: $"
# $Date: $
# $Revision: $
# $LastChangedBy: $
# $HeadURL: $
# $Id: $

#=============================================================================================

#=============================================================================================
# INSTALLATION INSTRUCTIONS
#=============================================================================================

  To compile on Linux:

  gcc -O3 -lm -fPIC -shared -I(directory with Python.h) -I(directory with numpy/arrayobject.h) -o _sstmap_ext_old.so _sstmap_ext_old.c

  For a desmond installation of python 2.5 (change path up to desmond directory, rest should be the same):
  
  gcc -O3 -lm -fPIC -shared -I /home/kamran/desmond/mmshare-v24012/lib/Linux-x86_64/include/python2.7 -I /home/kamran/desmond/mmshare-v24012/lib/Linux-x86_64/lib/python2.7/site-packages/numpy/core/include/ -o _sstmap_ext_old.so _sstmap_ext_old.c
  
  For a default installation of python 2.5:
  
  gcc -O3 -lm -fPIC -shared -I/usr/local/include/python2.5 -I/usr/local/lib/python2.5/site-packages/numpy/core/include -o _sstmap_ext_old.so _sstmap_ext_old.c

/Users/Kamran/anaconda/include/python2.7/Python.h
/Users/Kamran/anaconda/pkgs/numpy-1.9.0-py27_0/lib/python2.7/site-packages/numpy/core/include/
*/
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/arrayobject.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


double dist_mic(double x1, double x2, double x3, double y1, double y2, double y3, double b1, double b2, double b3) {
    /* Method for obtaining inter atom distance using minimum image convention
     */
    double dx, dy, dz;
    dx = x1-y1;
    dy = x2-y2;
    dz = x3-y3;
    if (dx > b1/2.0) dx -= b1; 
    else if (dx < -b1/2.0) dx += b1; 
    if (dy > b2/2.0) dy -= b2;
    else if (dy < -b2/2.0) dy += b2;
    if (dz > b3/2.0) dz -= b3; 
    else if (dz < -b3/2.0) dz += b3;

    return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
    }
    
double dist(double x1, double x2, double x3, double y1, double y2, double y3) {
    /* Method for Euclidean distance between two points
     */
    double dx, dy, dz;
    dx = x1-y1;
    dy = x2-y2;
    dz = x3-y3;
    return sqrt(pow(dx, 2)+ pow(dy, 2)+ pow(dz, 2));
    }

/*
Calculates electrostatic energy of a query water molecule against a set of target atoms

*/



PyObject *_sstmap_ext_old_processGrid(PyObject *self, PyObject *args)
{

    // declare variables received from Python
    int doEnergy, doEntropy, doHbonds, i_frame;
    PyArrayObject *watSites;
    PyArrayObject *coords, *periodic_box;
    PyArrayObject *vdw, *chg;
    PyArrayObject *grid_dim, *grid_max, *grid_orig, *voxel_wat_map, *voxel_data;
    PyArrayObject *wat_oxygen_ids, *prot_at_ids, *all_ids;
    PyArrayObject *hb_types;

    // Parse arguments received from Python
    if (!PyArg_ParseTuple(args, "iiiiO!O!O!O!O!O!O!O!O!O!O!O!O!O!:processGrid", 
                            &doEnergy, &doEntropy, &doHbonds, &i_frame,
                            &PyArray_Type,&watSites,
                            &PyArray_Type, &coords, &PyArray_Type, &periodic_box,
                            &PyArray_Type, &vdw, &PyArray_Type, &chg,
                            &PyArray_Type, &grid_dim, &PyArray_Type, &grid_max, &PyArray_Type, &grid_orig, &PyArray_Type, &voxel_wat_map, &PyArray_Type, &voxel_data,
                            &PyArray_Type, &wat_oxygen_ids, &PyArray_Type, &prot_at_ids, &PyArray_Type, &all_ids,
                            &PyArray_Type, &hb_types))
    {
        // raise argument parsing exception
        return NULL;
    }
    //FIXME: add consistency checks for python arrays
    // declare local variables
    double grid_max_x, grid_max_y, grid_max_z;
    double grid_orig_x, grid_orig_y, grid_orig_z;
    int grid_dim_x, grid_dim_y, grid_dim_z;
    double wat_translated_x, wat_translated_y, wat_translated_z;
    int n_wat, n_solute;
    int i_wat, j_wat, j_solute;

    int grid_index_x, grid_index_y, grid_index_z;


    grid_max_x = *(double *)PyArray_GETPTR1(grid_max, 0);
    grid_max_y = *(double *)PyArray_GETPTR1(grid_max, 1);
    grid_max_z = *(double *)PyArray_GETPTR1(grid_max, 2);
    grid_orig_x = *(double *)PyArray_GETPTR1(grid_orig, 0);
    grid_orig_y = *(double *)PyArray_GETPTR1(grid_orig, 1);
    grid_orig_z = *(double *)PyArray_GETPTR1(grid_orig, 2);
    grid_dim_x = *(int *)PyArray_GETPTR1(grid_dim, 0);
    grid_dim_y = *(int *)PyArray_GETPTR1(grid_dim, 1);
    grid_dim_z = *(int *)PyArray_GETPTR1(grid_dim, 2);

    // retrieve unit cell lengths for this frame
    float *b_x, *b_y, *b_z;
    b_x = (float *) PyArray_GETPTR2(periodic_box, 0, 0);
    b_y = (float *) PyArray_GETPTR2(periodic_box, 0, 1);
    b_z = (float *) PyArray_GETPTR2(periodic_box, 0, 2);
    //printf("Unit cell dimensions: %f %f %f\n", *b_x, *b_y, *b_z);

    n_wat = PyArray_DIM(wat_oxygen_ids, 0);
    n_solute = PyArray_DIM(prot_at_ids, 0);

    double *wat_eps, *wat_sig; // vdw parameters
    double wat_comb_eps, wat_comb_sig, Acoeff, Bcoeff;


    wat_sig = (double *) PyArray_GETPTR2(vdw, *(int *) PyArray_GETPTR1(wat_oxygen_ids, 0), 0);
    wat_eps = (double *) PyArray_GETPTR2(vdw, *(int *) PyArray_GETPTR1(wat_oxygen_ids, 0), 1);
    wat_comb_sig = (*wat_sig + *wat_sig)/2;
    wat_comb_eps = sqrt((*wat_eps)*(*wat_eps));
    Acoeff = 4*wat_comb_eps*pow(wat_comb_sig, 12);
    Bcoeff = -4*wat_comb_eps*pow(wat_comb_sig, 6);
    //printf("Number of waters: %i\n", n_wat);
    // iterate over frames
    // iterate over water atoms
    //printf("Iterating over frame: %i\n", i_frame + 1);
    for (i_wat = 0; i_wat < n_wat; i_wat++)
    {
        int wat_id; // id of current water
        float *wat_x, *wat_y, *wat_z; // coordinates
        double *wat_chg;
        // get water ID to and use it to get x, y, z coordinates and vdw params
        wat_id = *(int *) PyArray_GETPTR1(wat_oxygen_ids, i_wat); // obtain atom index for this atom
        wat_x = (float *) PyArray_GETPTR2(coords, wat_id, 0);
        wat_y = (float *) PyArray_GETPTR2(coords, wat_id, 1); 
        wat_z = (float *) PyArray_GETPTR2(coords, wat_id, 2);
        wat_translated_x = *wat_x - grid_orig_x;
        wat_translated_y = *wat_y - grid_orig_y;
        wat_translated_z = *wat_z - grid_orig_z;

        wat_chg = (double *) PyArray_GETPTR1(chg, wat_id);
        //printf("Processing water: %i with coordinates %f %f %f\n", wat_id, *wat_x, *wat_y, *wat_z);
        //printf("water oxygen ID %i and coordinates %f %f %f\n", wat_id, *wat_x, *wat_y, *wat_z);
        // check if the distance between wateer coordinates and grid origin is less than the max grid point

        if (wat_translated_x <= grid_max_x && wat_translated_y <= grid_max_y && wat_translated_z <= grid_max_z &&
            wat_translated_x >= -1.5 && wat_translated_y >= -1.5 && wat_translated_z >= -1.5)
        {
            if (wat_translated_x >= 0 && wat_translated_y >= 0 && wat_translated_z >= 0)
            {
                // transform water coordinates in units of grid dimensions
                grid_index_x = (int) (wat_translated_x/0.5);
                grid_index_y = (int) (wat_translated_y/0.5);
                grid_index_z = (int) (wat_translated_z/0.5);
                // check if water coords (in grid dimensions) are less than grid dimensions in each direction
                if (grid_index_x < grid_dim_x && grid_index_y < grid_dim_y && grid_index_z < grid_dim_z)
                {
                    int voxel_id;
                    int wat_H_index;
                    double dist, dist6, dist12;
                    double E_ww = 0.0;
                    double E_sw = 0.0;
                    double E_nbr = 0.0;
                    double n_nbr = 0;
                    
                    float *wat_H_x, *wat_H_y, *wat_H_z;
                    double *wat_H_chg;

                    int wat_nbr_index = 8;
                    int wat_sites = *(int *) PyArray_GETPTR1(watSites, i_wat);
                    
                    // obtain the voxel ID for this water
                    voxel_id = (grid_index_x*grid_dim_y + grid_index_y)*grid_dim_z + grid_index_z;
                    //voxel_ = (gridindex[0]*griddim_[1] + gridindex[1])*griddim_[2] + gridindex[2];
                    //printf("voxel id generated from: %i, %i %i\n", );
                    // increment water count for this voxel
                    *(double *)PyArray_GETPTR2(voxel_data, voxel_id, 4) += 1.0;
                    //
                    *(int *)PyArray_GETPTR3(voxel_wat_map, i_frame, i_wat, 0) = wat_id;
                    *(int *)PyArray_GETPTR3(voxel_wat_map, i_frame, i_wat, 1) = voxel_id;
                    *(int *)PyArray_GETPTR3(voxel_wat_map, i_frame, i_wat, 2) = 1;
                    //printf("Water atom ID %i with coordinates %f %f %f assigned to voxel %i.\n", wat_id, *wat_x, *wat_y, *wat_z, voxel_id);
                                        // iterate over all other water oxygen atoms to get Lennard-Jones component of energy                         
                    for (j_wat = 0; j_wat < n_wat; j_wat++)
                    {
                        int j_wat_id;
                        int j_wat_sites;
                        float *j_wat_x, *j_wat_y, *j_wat_z;
                        double pair_energy = 0.0;
                        int isNbr = 0;
                        int j_wat_H_index;
                        float *j_wat_H_x, *j_wat_H_y, *j_wat_H_z;
                        double *j_wat_chg, *j_wat_H_chg;
                        // get water ID to and use it to get x, y, z coordinates and vdw params
                        j_wat_id = *(int *) PyArray_GETPTR1(wat_oxygen_ids, j_wat);

                        j_wat_sites = *(int *) PyArray_GETPTR1(watSites, j_wat);
                        if (j_wat_id == wat_id) continue;
                        //printf("Curent pair %i %i\n", wat_id, j_wat_id);
                        //printf("sites %i %i\n", wat_sites, j_wat_sites);
                        j_wat_x = (float *) PyArray_GETPTR2(coords, j_wat_id, 0);
                        j_wat_y = (float *) PyArray_GETPTR2(coords, j_wat_id, 1); 
                        j_wat_z = (float *) PyArray_GETPTR2(coords, j_wat_id, 2);
                        j_wat_chg = (double *) PyArray_GETPTR1(chg, j_wat_id);
                        // calculate vdw energy
                        dist = dist_mic(*wat_x, *wat_y, *wat_z, *j_wat_x, *j_wat_y, *j_wat_z, *b_x, *b_y, *b_z);
                        if (dist <= 3.5) isNbr = 1;
                        dist6 = pow(dist, 6);
                        dist12 = dist6 * dist6;
                        pair_energy += (Acoeff/dist12)+(Bcoeff/dist6);
                        // calculate O-O elec energy
                        pair_energy += ((*wat_chg)*(*j_wat_chg))/dist;
                        //printf("O-O charges %f %f\n", *wat_chg, *j_wat_chg);
                        // calculate elec energy with rest of the j water molecule
                        // iterate over non-oxygen atom of water j
                        for (j_wat_H_index = j_wat_id + 1; j_wat_H_index < j_wat_id + j_wat_sites; j_wat_H_index++)
                        {
                            //printf("Current pair %i %i \n", j_wat_H_index, wat_id);
                            j_wat_H_x = (float *) PyArray_GETPTR2(coords, j_wat_H_index, 0);
                            j_wat_H_y = (float *) PyArray_GETPTR2(coords, j_wat_H_index, 1); 
                            j_wat_H_z = (float *) PyArray_GETPTR2(coords, j_wat_H_index, 2);
                            j_wat_H_chg = (double *) PyArray_GETPTR1(chg, j_wat_H_index);
                            dist = dist_mic(*wat_x, *wat_y, *wat_z, *j_wat_H_x, *j_wat_H_y, *j_wat_H_z, *b_x, *b_y, *b_z);
                            pair_energy += ((*wat_chg)*(*j_wat_H_chg))/dist;
                            //printf("%f\n", ((*wat_chg)*(*j_wat_H_chg))/dist);
                            // iterate over non-oxygen atom of water i
                            for (wat_H_index = wat_id + 1; wat_H_index < wat_id + wat_sites; wat_H_index++)
                            {
                                //printf("Current pair %i %i \n", j_wat_H_index, wat_H_index);
                                //printf("Water atom id: %i\n", wat_H_index);
                                wat_H_x = (float *) PyArray_GETPTR2(coords, wat_H_index, 0);
                                wat_H_y = (float *) PyArray_GETPTR2(coords, wat_H_index, 1); 
                                wat_H_z = (float *) PyArray_GETPTR2(coords, wat_H_index, 2);
                                wat_H_chg = (double *) PyArray_GETPTR1(chg, wat_H_index);
                                dist = dist_mic(*wat_H_x, *wat_H_y, *wat_H_z, *j_wat_H_x, *j_wat_H_y, *j_wat_H_z, *b_x, *b_y, *b_z);
                                pair_energy += ((*wat_H_chg)*(*j_wat_H_chg))/dist;
                                //printf("%f\n", ((*wat_H_chg)*(*j_wat_H_chg))/dist);
                            } 
                        }
                        // iterate over non-oxygen atom of water i 
                        for (wat_H_index = wat_id + 1; wat_H_index < wat_id + wat_sites; wat_H_index++)
                        {
                            //printf("Current pair %i %i \n", wat_H_index, j_wat_id);
                            //printf("Water atom id: %i\n", wat_H_index);
                            wat_H_x = (float *) PyArray_GETPTR2(coords, wat_H_index, 0);
                            wat_H_y = (float *) PyArray_GETPTR2(coords, wat_H_index, 1); 
                            wat_H_z = (float *) PyArray_GETPTR2(coords, wat_H_index, 2);
                            wat_H_chg = (double *) PyArray_GETPTR1(chg, wat_H_index);
                            dist = dist_mic(*wat_H_x, *wat_H_y, *wat_H_z, *j_wat_x, *j_wat_y, *j_wat_z, *b_x, *b_y, *b_z);
                            pair_energy += ((*wat_H_chg)*(*j_wat_chg))/dist;
                            //printf("%f\n", ((*wat_H_chg)*(*wat_chg))/dist);
                        }
                        //printf("LJ energy between water oxygen IDs %i %i = %f\n", wat_id, j_wat_id, E_sw_LJ);
                        if (isNbr == 1)
                        {
                            *(int *)PyArray_GETPTR3(voxel_wat_map, i_frame, i_wat, 6) += 1;
                            *(int *)PyArray_GETPTR3(voxel_wat_map, i_frame, i_wat, wat_nbr_index) = j_wat_id;
                            E_nbr += pair_energy;
                            *(double *)PyArray_GETPTR2(voxel_data, voxel_id, 18) += 1;
                            n_nbr += 1.0;
                            wat_nbr_index += 1;
                        }
                        E_ww += pair_energy;
                        //printf("pair energy between water oxygen IDs %i %i = %f\n", wat_id, j_wat_id, pair_energy);
                    } // end loop over water atoms
                    *(double *)PyArray_GETPTR2(voxel_data, voxel_id, 14) += E_ww;
                    if (n_nbr >= 1.0)
                    {
                        //printf("Enbr for water %i = %f\n", wat_id, E_nbr);
                        *(double *)PyArray_GETPTR2(voxel_data, voxel_id, 16) += E_nbr/n_nbr;
                    }
                    else
                    {
                        //printf("Enbr for water %i = %f\n", wat_id, E_nbr);
                        *(double *)PyArray_GETPTR2(voxel_data, voxel_id, 16) += E_nbr;
                    }

                    // iterate over all other water oxygen atoms to get Lennard-Jones component of energy                         
                    for (j_solute = 0; j_solute < n_solute; j_solute++)
                    {
                        int j_solute_id;
                        //char *hb_type;
                        float *j_solute_x, *j_solute_y, *j_solute_z;
                        double *j_solute_eps, *j_solute_sig, *j_solute_chg;
                        double comb_sig, comb_eps, aij, bij;
                        int hb_type;
                        // get water ID to and use it to get x, y, z coordinates and vdw params
                        j_solute_id = *(int *) PyArray_GETPTR1(prot_at_ids, j_solute);
                        j_solute_x = (float *) PyArray_GETPTR2(coords, j_solute_id, 0);
                        j_solute_y = (float *) PyArray_GETPTR2(coords, j_solute_id, 1); 
                        j_solute_z = (float *) PyArray_GETPTR2(coords, j_solute_id, 2);
                        j_solute_chg = (double *) PyArray_GETPTR1(chg, j_solute_id);
                        j_solute_sig = (double *) PyArray_GETPTR2(vdw, j_solute_id, 0);
                        j_solute_eps = (double *) PyArray_GETPTR2(vdw, j_solute_id, 1);

                        comb_sig = (*wat_sig + *j_solute_sig)/2.0;
                        comb_eps = sqrt((*wat_eps)*(*j_solute_eps));
                        //printf("index %i and Soluet atom ID %i and coordinates %f %f %f\n", j_solute, j_solute_id, *j_solute_x, *j_solute_y, *j_solute_z);
                        //printf("Eps and Sig for current pair: %f %f\n", comb_eps, comb_sig);
                        aij = 4*comb_eps*pow(comb_sig, 12);
                        bij = -4*comb_eps*pow(comb_sig, 6);
                        //printf("Acoeff and Bcoeff for current pair: %f %f\n", aij, bij);
                        // calculate vdw energy
                        dist = dist_mic(*wat_x, *wat_y, *wat_z, *j_solute_x, *j_solute_y, *j_solute_z, *b_x, *b_y, *b_z);
                        dist6 = pow(dist, 6);
                        dist12 = dist6 * dist6;
                        E_sw += (aij/dist12)+(bij/dist6);
                        E_sw += ((*wat_chg)*(*j_solute_chg))/dist;
                        // iterate over non-oxygen atom of water i 
                        for (wat_H_index = wat_id + 1; wat_H_index < wat_id + wat_sites; wat_H_index++)
                        {
                            //printf("Current pair %i %i \n", wat_H_index, j_wat_id);
                            //printf("Water atom id: %i\n", wat_H_index);
                            wat_H_x = (float *) PyArray_GETPTR2(coords, wat_H_index, 0);
                            wat_H_y = (float *) PyArray_GETPTR2(coords, wat_H_index, 1); 
                            wat_H_z = (float *) PyArray_GETPTR2(coords, wat_H_index, 2);
                            wat_H_chg = (double *) PyArray_GETPTR1(chg, wat_H_index);
                            dist = dist_mic(*wat_H_x, *wat_H_y, *wat_H_z, *j_solute_x, *j_solute_y, *j_solute_z, *b_x, *b_y, *b_z);
                            E_sw += ((*wat_H_chg)*(*j_solute_chg))/dist;
                            //printf("%f\n", ((*wat_H_chg)*(*wat_chg))/dist);
                        }

                        //printf("Distance between water %i and solute atom %i = %f\n", wat_id, j_solute_id, dist);
                        if (dist < 3.5)
                        {
                            hb_type = *(int *) PyArray_GETPTR1(hb_types, j_solute);
                            if (hb_type > 0)
                            {
                                //printf("Water: %i, HB type of %i is %i and distance is %f\n", wat_id, j_solute_id, hb_type, dist);
                                *(int *)PyArray_GETPTR3(voxel_wat_map, i_frame, i_wat, 7) += 1;
                                *(int *)PyArray_GETPTR3(voxel_wat_map, i_frame, i_wat, wat_nbr_index) = j_solute_id;
                                wat_nbr_index += 1;
                            }
                        }
                    } // end loop over solute atoms
                    *(double *)PyArray_GETPTR2(voxel_data, voxel_id, 12) += E_sw;
                    //printf("LJ energy between water %i and protein = %f\n", wat_id, E_sw_LJ);

                } // end  inside grid check  
            } // end
        } // end
    } // end iterating over waters
    return Py_BuildValue("i", 1);
}



PyObject *_sstmap_ext_old_getNNEntropy(PyObject *self, PyObject *args)
{

    // declare variables received from Python
    int voxel_id, n_wat;
    PyArrayObject *voxel_wats;
    PyArrayObject *voxel_wat_Eulers;
    //PyArrayObject *coords, *periodic_box;
    PyArrayObject *voxel_data;

    // Parse arguments received from Python
    if (!PyArg_ParseTuple(args, "iiO!O!O!:processGrid", 
                            &voxel_id, &n_wat,
                            &PyArray_Type, &voxel_wats,
                            &PyArray_Type, &voxel_wat_Eulers,
                            &PyArray_Type, &voxel_data))
    {
        return NULL;
    }
    int n, l;
    int wat_i, wat_j;
    int wat_i_frame, wat_j_frame;
    double NNor, dW, wat_or_ent;
    double voxel_dTSor = 0.0;
    double wat_i_theta, wat_i_phi, wat_i_psi;
    double wat_j_theta, wat_j_phi, wat_j_psi;
    double rx, ry, rz;
    double twopi = 2*M_PI;

    // for each water in the voxel
    for (n = 0; n < n_wat; n++)
    {
        NNor = 10000;
        wat_i = *(int *) PyArray_GETPTR2(voxel_wats, n, 0);
        wat_i_frame = *(int *) PyArray_GETPTR2(voxel_wats, n, 1);
        wat_i_theta = *(double *)PyArray_GETPTR3(voxel_wat_Eulers, wat_i_frame, wat_i, 0);
        wat_i_phi = *(double *)PyArray_GETPTR3(voxel_wat_Eulers, wat_i_frame, wat_i, 1);
        wat_i_psi = *(double *)PyArray_GETPTR3(voxel_wat_Eulers, wat_i_frame, wat_i, 2);
        //printf("water %i theta phi psi %f %f %f\n", wat_i, wat_i_theta, wat_i_phi, wat_i_psi);
        for (l = 0; l < n_wat; l++)
        {
            if(l == n) continue;
            wat_j = *(int *) PyArray_GETPTR2(voxel_wats, l, 0);
            wat_j_frame = *(int *) PyArray_GETPTR2(voxel_wats, l, 1);
            //printf("Calculating orientational distancce between water: %i and %i\n", wat_i, wat_j);
            wat_j_theta = *(double *)PyArray_GETPTR3(voxel_wat_Eulers, wat_j_frame, wat_j, 0);
            wat_j_phi = *(double *)PyArray_GETPTR3(voxel_wat_Eulers, wat_j_frame, wat_j, 1);
            wat_j_psi = *(double *)PyArray_GETPTR3(voxel_wat_Eulers, wat_j_frame, wat_j, 2);
            //printf("water %i theta phi psi %f %f %f\n", wat_j, wat_j_theta, wat_j_phi, wat_j_psi);

            rx = cos(wat_j_theta) - cos(wat_i_theta);
            ry = wat_j_phi - wat_i_phi;
            rz = wat_j_psi - wat_i_psi;

            if      (ry>M_PI) ry = twopi-ry;
            else if (ry<-M_PI) ry = twopi+ry;
            if      (rz>M_PI) rz = twopi-rz;
            else if (rz<-M_PI) rz = twopi+rz;

            dW = sqrt(rx*rx + ry*ry + rz*rz);
            //printf("Calculating orientational distancce between water: %i and %i = %f\n", wat_i, wat_j, dW);
            //dR = 0.0;
            // get six-D distance
            // get translational nearest neighbor            
            if (dW>0 && dW<NNor) NNor = dW;
            // get six-D nearest neighbor            
            
        }
        //calculate orientational entropy
        
        if (NNor<9999 && NNor>0) 
        {
            //printf("Nearest neighbour orientational distance: %f\n", NNor);
            //wat_tr_ent = log(nwtot*NNtr*NNtr*NNtr/(3.0*twopi));
            //voxel_dTStr_norm += wat_tr_ent;
            wat_or_ent = log(n_wat*NNor*NNor*NNor/(3.0*twopi));
            voxel_dTSor += wat_or_ent;
            
        }

    }
    // 
    //*(double *)PyArray_GETPTR1(ent, 2) += voxel_dTSor_norm;
    return Py_BuildValue("f", voxel_dTSor);

}

PyObject *_sstmap_ext_elecE( PyObject *self, PyObject *args) // we don't have a class here but still need *self argument
    /* Method for calculation electrostatic energy between a water molecule and a set of other molecules
    */
    {
    // First we declare variables to be used in this function
    npy_intp m, n;
    int i, j;
    int *wat, *other;
    float *b_x, *b_y, *b_z, *wx, *wy, *wz, *sx, *sy, *sz;
     
    double *wc, *sc;
    
    double d;
    double e_elec = 0.0;
    // These variables are declared to store Python objects, which will arrive here wrapped in args tuple
    PyArrayObject *wat_at_ids, *other_at_ids, *coords, *charges, *box;
    // Here we parse the args tuple as it contains all the variables sent from Python
    // In this case Python sent five Numpy arrays (hence five times O!, O means a python object, ! is a pointer check)
    // name of the function as called from Python is also in the string, will be used in locating errors
    // for each O!, we need to supply a pointer type and a pointer
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!:elecE",
        &PyArray_Type, &wat_at_ids,
        &PyArray_Type, &other_at_ids,
        &PyArray_Type, &coords,
        &PyArray_Type, &charges,
        &PyArray_Type, &box))
        {
            return NULL; /* raise argument parsing exception*/
        }
    
    // A consistency check for correct data type being sent from Python
    if (PyArray_NDIM(coords) != 2)
        {
            PyErr_Format(PyExc_ValueError, "coordinate array is not of correct dimensions or type");
            return NULL;
            
        }
        
    m = PyArray_DIM(wat_at_ids, 0); // m is indec over water atoms
    n = PyArray_DIM(other_at_ids, 0); // n is indec over all other atoms
    // periodic box information
    //printf("index over all water atoms: %i index over non-water atoms %i\n", m ,n);
    b_x = (float *) PyArray_GETPTR1(box, 0); 
    b_y = (float *) PyArray_GETPTR1(box, 1);
    b_z = (float *) PyArray_GETPTR1(box, 2);
    
    // loop over each water atom and over each other atom
    for (i = 0; i < m; i++) {
        wat = (int *) PyArray_GETPTR1(wat_at_ids, i); // obtain index for this atom (this is not array index, this is unique atom id)
        //printf("iteration index: %i corresponding water atom %i\n", i, *wat);
        wx = PyArray_GETPTR2(coords, *wat, 0); // use wat to get the correct x, y, z coordinates from coord array
        wy = PyArray_GETPTR2(coords, *wat, 1); 
        wz = PyArray_GETPTR2(coords, *wat, 2);
        wc = PyArray_GETPTR1(charges, *wat); // use wat to get the correct charge from charge array
        //printf("water atom ID x y z charge: %i %5.3f %5.3f %5.3f %5.3f \n", *wat, *wx, *wy, *wz, *wc);
        for (j = 0; j < n; j++) {
            other = (int *) PyArray_GETPTR1(other_at_ids, j);
            sx = PyArray_GETPTR2(coords, *other, 0); // obtain index of this atom
            sy = PyArray_GETPTR2(coords, *other, 1); // obtain x, y, z
            sz = PyArray_GETPTR2(coords, *other, 2);
            sc = PyArray_GETPTR1(charges, *other); // charge on this atom
            d = dist_mic(*wx, *wy, *wz, *sx, *sy, *sz, *b_x, *b_y, *b_z); // distance (based on minimum image convention)
            //d = dist(*wx, *wy, *wz, *sx, *sy, *sz); // distance (based on minimum image convention)
            e_elec += ((*wc)*(*sc))/d; // Coulombic interaction calculation
        }
        
    }
    return Py_BuildValue("f", e_elec);
}

PyObject *_sstmap_ext_vdwE( PyObject *self, PyObject *args)
    {
    npy_intp m, n;
    int i, j;
    int *wat, *other;
    float *b_x, *b_y, *b_z, *wx, *wy, *wz, *sx, *sy, *sz;

    double *w_sig, *w_eps, *s_sig, *s_eps;

    double comb_sig, comb_eps;
    float d, dist6, dist12; 
    double aij, bij;
    double e_vdw = 0.0;

    PyArrayObject *wat_at_ids, *other_at_ids, *coords, *vdwparms, *box;
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!:elecE",
        &PyArray_Type, &wat_at_ids,
        &PyArray_Type, &other_at_ids,
        &PyArray_Type, &coords,
        &PyArray_Type, &vdwparms,
        &PyArray_Type, &box))
        {
            return NULL; /* raise argument parsing exception*/
        }

    if (PyArray_NDIM(coords) != 2)
        {
            PyErr_Format(PyExc_ValueError, "coordinate array is not of correct dimensions or type");
            return NULL;
            
        }
    if (PyArray_NDIM(vdwparms) != 2)
        {
            PyErr_Format(PyExc_ValueError, "vdw parmeter array is not of correct dimensions or type");
            return NULL;
            
        }

    m = PyArray_DIM(wat_at_ids, 0);
    n = PyArray_DIM(other_at_ids, 0);
    
    b_x = (float *) PyArray_GETPTR1(box, 0);
    b_y = (float *) PyArray_GETPTR1(box, 1);
    b_z = (float *) PyArray_GETPTR1(box, 2);
    
    for (i = 0; i < m; i++) {
        wat = (int *) PyArray_GETPTR1(wat_at_ids, i);
        wx = PyArray_GETPTR2(coords, *wat, 0);
        wy = PyArray_GETPTR2(coords, *wat, 1);
        wz = PyArray_GETPTR2(coords, *wat, 2);
        w_sig = PyArray_GETPTR2(vdwparms, *wat, 0);
        w_eps = PyArray_GETPTR2(vdwparms, *wat, 1);
        
        for (j = 0; j < n; j++) {
            other = (int *) PyArray_GETPTR1(other_at_ids, j);
            sx = PyArray_GETPTR2(coords, *other, 0);
            sy = PyArray_GETPTR2(coords, *other, 1);
            sz = PyArray_GETPTR2(coords, *other, 2);
            s_sig = PyArray_GETPTR2(vdwparms, *other, 0);
            s_eps = PyArray_GETPTR2(vdwparms, *other, 1);
            comb_sig = (*w_sig + *s_sig)/2;
            comb_eps = sqrt((*w_eps)*(*s_eps));
            d = dist_mic(*wx, *wy, *wz, *sx, *sy, *sz, *b_x, *b_y, *b_z);
            //d = dist(*wx, *wy, *wz, *sx, *sy, *sz);
            //printf("water atom 1 x y z atom 2 x y z dist: %i %5.3f %5.3f %5.3f %i %5.3f %5.3f %5.3f %5.3f \n", *wat, *wx, *wy, *wz, *other, *sx, *sy, *sz, d);
            dist6 = pow(d, 6);
            dist12 = dist6 * dist6;
            aij = 4*comb_eps*pow(comb_sig, 12);
            bij = -4*comb_eps*pow(comb_sig, 6);
            e_vdw +=  (aij/dist12)+(bij/dist6);
        }
        
    }
    return Py_BuildValue("f", e_vdw);
}


/* Method Table
 * Registering all the functions that will be called from Python
 */

static PyMethodDef _sstmap_ext_old_methods[] = {
    {
        "processGrid",
        (PyCFunction)_sstmap_ext_old_processGrid,
        METH_VARARGS,
        "Process grid"
    },
    
    {
        "getNNEntropy",
        (PyCFunction)_sstmap_ext_old_getNNEntropy,
        METH_VARARGS,
        "get voxel entropy"
    },
    {
        "elecE",                           // name of the fucntion called from Python
        (PyCFunction)_sstmap_ext_elecE,     // corresponding C++ function
        METH_VARARGS,
        "compute electrostatic energies"   // doc string
    },
    {
        "vdwE",
        (PyCFunction)_sstmap_ext_vdwE,
        METH_VARARGS,
        "compute LJ energies"
    },



    
    {NULL, NULL}
};

/* Initialization function for this module
 */

PyMODINIT_FUNC init_sstmap_ext_old(void) // init function has the same name as module, except with init prefix
{
    // we produce name of the module, method table and a doc string
    Py_InitModule3("_sstmap_ext_old", _sstmap_ext_old_methods, "Process GIST calcs.\n");
    import_array(); // required for Numpy initialization
}
