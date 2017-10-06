/*
  sstmap C extension module.
##############################################################################
# SSTMap: A Python library for the calculation of water structure and
#         thermodynamics on solute surfaces from molecular dynamics
#         trajectories.
# Copyright 2016-2017 Lehman College City University of New York
# and the Authors
#
# Authors: Kamran Haider
# Contributors: Steven Ramsay, Anthony Cruz Balberdy
#
# SSTMap is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SSTMap. If not, see <http://www.gnu.org/licenses/>.
##############################################################################

*/
#define _USE_MATH_DEFINES
#include <math.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/arrayobject.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>


double dist_mic(double x1, double x2, double x3, double y1, double y2, double y3, double b1, double b2, double b3) {
    /* Method for obtaining inter atom distance using minimum image convention
     */
    //printf("x1: %f, x2: %f, x3: %f\n", x1, x2, x3);
    //printf("y1: %f, y2: %f, y3: %f\n", y1, y2, y3);
    double dx, dy, dz;
    dx = x1-y1;
    dy = x2-y2;
    dz = x3-y3;
    //printf("dx: %f, dy: %f, dz: %f\n", dx, dy, dz);
    //printf("bx: %f, by: %f, bz: %f\n", b1/2.0, b2/2.0, b3/2.0);
    if (dx > b1/2.0) dx -= b1;
    else if (dx < -b1/2.0) dx += b1;
    if (dy > b2/2.0) dy -= b2;
    else if (dy < -b2/2.0) dy += b2;
    if (dz > b3/2.0) dz -= b3;
    else if (dz < -b3/2.0) dz += b3;
    //printf("dist = %f", sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2)));
    return 1.0/(sqrt((dx*dx) +(dy*dy) + (dz*dz)));
    }

double dist_mic_squared(double x1, double x2, double x3, double y1, double y2, double y3, double b1, double b2, double b3) {
    /* Method for obtaining inter atom distance using minimum image convention
     */
    //printf("x1: %f, x2: %f, x3: %f\n", x1, x2, x3);
    //printf("y1: %f, y2: %f, y3: %f\n", y1, y2, y3);
    double dx, dy, dz;
    dx = x1-y1;
    dy = x2-y2;
    dz = x3-y3;
    //printf("dx: %f, dy: %f, dz: %f\n", dx, dy, dz);
    //printf("bx: %f, by: %f, bz: %f\n", b1/2.0, b2/2.0, b3/2.0);
    if (dx > b1/2.0) dx -= b1;
    else if (dx < -b1/2.0) dx += b1;
    if (dy > b2/2.0) dy -= b2;
    else if (dy < -b2/2.0) dy += b2;
    if (dz > b3/2.0) dz -= b3;
    else if (dz < -b3/2.0) dz += b3;
    //printf("dist = %f", sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2)));
    return (dx*dx) + (dy*dy) + (dz*dz);
    }

double dist(double x1, double x2, double x3, double y1, double y2, double y3) {
    /* Method for Euclidean distance between two points
     */
    double dx, dy, dz;
    dx = x1-y1;
    dy = x2-y2;
    dz = x3-y3;
    //printf("dist = %f", sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2)));
    return sqrt(pow(dx, 2)+ pow(dy, 2)+ pow(dz, 2));
    }

double dist_squared(double x1, double x2, double x3, double y1, double y2, double y3) {
    /* Method for Euclidean distance between two points
     */
    double dx, dy, dz;
    dx = x1-y1;
    dy = x2-y2;
    dz = x3-y3;
    //printf("dist = %f", sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2)));
    return pow(dx, 2)+ pow(dy, 2)+ pow(dz, 2);
    }

/*
Calculates electrostatic energy of a query water molecule against a set of target atoms

*/
PyObject *_sstmap_ext_assign_voxels(PyObject *self, PyObject *args)
{

    int n_frames, i_frame, n_wat, i_wat;
    PyArrayObject *coords, *grid_dim, *grid_max, *grid_orig, *wat_oxygen_ids;
    PyObject *frame_data;
    // declare local variables
    double grid_max_x, grid_max_y, grid_max_z;
    double grid_orig_x, grid_orig_y, grid_orig_z;
    int grid_dim_x, grid_dim_y, grid_dim_z;
    int grid_index_x, grid_index_y, grid_index_z;    
    double wat_translated_x, wat_translated_y, wat_translated_z;
    int wat_id; // id of current water
    float *wat_x, *wat_y, *wat_z; // coordinates
    int dims[1];
    PyArrayObject *wat_data;
    PyObject *curr_voxel;
    int voxel_id;

    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!",
        &PyArray_Type, &coords,
        &PyArray_Type, &grid_dim,
        &PyArray_Type, &grid_max,
        &PyArray_Type, &grid_orig,
        &PyList_Type, &frame_data,
        &PyArray_Type, &wat_oxygen_ids
        ))
    {
        return NULL;
    }
    //Set the number of frames
    n_frames = PyArray_DIM(coords, 0);
    n_wat = PyArray_DIM(wat_oxygen_ids, 0);
    //printf("The number of frames passed = %i\n", n_frames);
    //printf("The number of waters passed = %i\n", n_wat);



    grid_max_x = *(double *)PyArray_GETPTR1(grid_max, 0);
    grid_max_y = *(double *)PyArray_GETPTR1(grid_max, 1);
    grid_max_z = *(double *)PyArray_GETPTR1(grid_max, 2);
    grid_orig_x = *(double *)PyArray_GETPTR1(grid_orig, 0);
    grid_orig_y = *(double *)PyArray_GETPTR1(grid_orig, 1);
    grid_orig_z = *(double *)PyArray_GETPTR1(grid_orig, 2);
    grid_dim_x = *(int *)PyArray_GETPTR1(grid_dim, 0);
    grid_dim_y = *(int *)PyArray_GETPTR1(grid_dim, 1);
    grid_dim_z = *(int *)PyArray_GETPTR1(grid_dim, 2);

    for (i_frame = 0; i_frame < n_frames; i_frame++)
    {
        //printf("Iterating over frame: %i\n", i_frame);
        for (i_wat = 0; i_wat < n_wat; i_wat++)
        {
            dims[0] = 2;
            // get water ID to and use it to get x, y, z coordinates and vdw params
            wat_id = *(int *) PyArray_GETPTR1(wat_oxygen_ids, i_wat); // obtain atom index for this atom
            wat_x = (float *) PyArray_GETPTR3(coords, i_frame, wat_id, 0);
            wat_y = (float *) PyArray_GETPTR3(coords, i_frame, wat_id, 1); 
            wat_z = (float *) PyArray_GETPTR3(coords, i_frame, wat_id, 2);
            wat_translated_x = *wat_x - grid_orig_x;
            wat_translated_y = *wat_y - grid_orig_y;
            wat_translated_z = *wat_z - grid_orig_z;

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
                        // obtain the voxel ID for this water
                        voxel_id = (grid_index_x*grid_dim_y + grid_index_y)*grid_dim_z + grid_index_z;
                        //voxel_ = (gridindex[0]*griddim_[1] + gridindex[1])*griddim_[2] + gridindex[2];
                        //printf("Water atom ID %i with coordinates %f %f %f assigned to voxel %i.\n", wat_id, *wat_x, *wat_y, *wat_z, voxel_id);
                        wat_data = (PyArrayObject *) PyArray_FromDims(1, dims, NPY_INT);
                        *(int *)PyArray_GETPTR1(wat_data, 0) = voxel_id;
                        *(int *)PyArray_GETPTR1(wat_data, 1) = wat_id;
                        //printf("wat_data: %d %d %d\n", voxel_id, *(int *)PyArray_GETPTR1(wat_data, 0), *(int *)PyArray_GETPTR1(wat_data, 1));
                        curr_voxel = PyList_GetItem(frame_data, i_frame);
                        PyList_Append(curr_voxel, wat_data);
                        //DECREF?
                    }
                }
            }
        } // finish iterating over waters
    }
    return Py_BuildValue("i", 1);
}

PyObject *_sstmap_ext_get_pairwise_distances(PyObject *self, PyObject *args)
{
    PyArrayObject *wat, *target_at_ids, *coords, *periodic_box, *dist_array;
    float *b_x, *b_y, *b_z;
    int wat_sites, wat_atom, wat_atom_id;
    int num_target_at, target_at, target_at_id;
    int frame = 0;
    double d;
    float *wat_x, *wat_y, *wat_z;
    float *target_at_x, *target_at_y, *target_at_z;

    if (!PyArg_ParseTuple(args, "O!O!O!O!O!",
        &PyArray_Type, &wat,
        &PyArray_Type, &target_at_ids,
        &PyArray_Type, &coords,
        &PyArray_Type, &periodic_box,
        &PyArray_Type, &dist_array
        ))
    {
        return NULL;
    }
    // do distance calc here
        // retrieve unit cell lengths for this frame
    b_x = (float *) PyArray_GETPTR2(periodic_box, 0, 0);
    b_y = (float *) PyArray_GETPTR2(periodic_box, 0, 1);
    b_z = (float *) PyArray_GETPTR2(periodic_box, 0, 2);
    //printf("Unit cell dimensions: %f %f %f\n", *b_x, *b_y, *b_z);

    wat_sites = PyArray_DIM(dist_array, 0);

    num_target_at = PyArray_DIM(target_at_ids, 0);


    for (wat_atom = 0; wat_atom < wat_sites; wat_atom++)
    {
        wat_atom_id = *(int *) PyArray_GETPTR1(wat, 1) + wat_atom;
        wat_x = (float *) PyArray_GETPTR3(coords, frame, wat_atom_id, 0);
        wat_y = (float *) PyArray_GETPTR3(coords, frame, wat_atom_id, 1); 
        wat_z = (float *) PyArray_GETPTR3(coords, frame, wat_atom_id, 2);

        for (target_at = 0; target_at < num_target_at; target_at++)
        {
            target_at_id = *(int *) PyArray_GETPTR1(target_at_ids, target_at);
            target_at_x = (float *) PyArray_GETPTR3(coords, frame, target_at_id, 0);
            target_at_y = (float *) PyArray_GETPTR3(coords, frame, target_at_id, 1); 
            target_at_z = (float *) PyArray_GETPTR3(coords, frame, target_at_id, 2);
            //printf("Iterator: %d, atom id: %d\n", target_at, target_at_id);
            //printf("Water atom coords %f %f %f\n", *wat_x, *wat_y, *wat_z);
            //printf("Target atom coords %f %f %f\n", *target_at_x, *target_at_y, *target_at_z);
            d = dist_mic_squared(*wat_x, *wat_y, *wat_z, *target_at_x, *target_at_y, *target_at_z, *b_x, *b_y, *b_z);
            //printf("Distance between %d and %d = %3.2f\n", wat_atom_id, target_at_id, d);
            *(double *)PyArray_GETPTR2(dist_array, wat_atom, target_at) += d;
        }

    }
    return Py_BuildValue("i", 1);

}
PyObject *_sstmap_ext_getNNOrEntropy(PyObject *self, PyObject *args)
{
    int nwtot, n, l;
    double NNor, dW, wat_or_ent;
    double voxel_dTSor = 0.0;
    double rx, ry, rz;
    PyArrayObject *voxel_wat_Eulers; 
    double twopi = 2*M_PI;

    // Argument parsing to reterive everything sent from Python correctly    
    if (!PyArg_ParseTuple(args, "iO!",
                            &nwtot,
                            &PyArray_Type, &voxel_wat_Eulers))
        {
            return NULL; /* raise argument parsing exception*/
        }
    // for each water in the voxel
    for (n = 0; n < nwtot; n++){
        NNor = 10000;
        for (l = 0; l < nwtot; l++){
            if(l == n) continue;
            //printf("Calculating orientational distancce between water: %i and %i\n", l, n);
            rx = cos(*(double *)PyArray_GETPTR2(voxel_wat_Eulers, l, 0)) - cos(*(double *)PyArray_GETPTR2(voxel_wat_Eulers, n, 0));
            ry = *(double *)PyArray_GETPTR2(voxel_wat_Eulers, l, 1) - *(double *)PyArray_GETPTR2(voxel_wat_Eulers, n, 1);
            rz = *(double *)PyArray_GETPTR2(voxel_wat_Eulers, l, 2) - *(double *)PyArray_GETPTR2(voxel_wat_Eulers, n, 2);
            if      (ry>M_PI) ry = twopi-ry;
            else if (ry<-M_PI) ry = twopi+ry;
            if      (rz>M_PI) rz = twopi-rz;
            else if (rz<-M_PI) rz = twopi+rz;
            dW = sqrt(rx*rx + ry*ry + rz*rz);
            //dR = 0.0;
            // get six-D distance
            // get translational nearest neighbor            
            if (dW>0 && dW<NNor) NNor = dW;
            // get six-D nearest neighbor            
        }
        //calculate translational entropy
        if (NNor<9999 && NNor>0) {
            //printf("Nearest neighbour translational distance: %f\n", NNtr);
            //wat_tr_ent = log(nwtot*NNtr*NNtr*NNtr/(3.0*twopi));
            //voxel_dTStr_norm += wat_tr_ent;
            wat_or_ent = log(nwtot*NNor*NNor*NNor/(3.0*twopi));
            voxel_dTSor += wat_or_ent;            
        }        
    }
    // 
    //*(double *)PyArray_GETPTR1(ent, 2) += voxel_dTSor_norm;
    return Py_BuildValue("f", voxel_dTSor);
}

PyObject *_sstmap_ext_getNNTrEntropy(PyObject *self, PyObject *args)
{
    float ref_dens;
    float voxel_vol;
    float temp;
    int num_frames, n0, n1;
    double dTStranstot = 0.0;
    double dTSorienttot = 0;
    double dTSt = 0.0;
    double dTSs = 0.0;
    double dTSo = 0.0;
    int nwts = 0;
    int nwtt = 0;
    double pi = 3.141592653589793;
    double twopi = 6.283185307179586;
    double gas_kcal = 0.0019872041;
    double euler_masc = 0.5772156649;

    unsigned int voxel;
    PyArrayObject *voxel_data, *grid_dims;
    PyObject *voxel_O_coords, *voxel_quarts;
    // Argument parsing to reterive everything sent from Python correctly
    if (!PyArg_ParseTuple(args, "ifffO!O!O!O!",
                            &num_frames,
                            &voxel_vol,
                            &ref_dens,
                            &temp,
                            &PyArray_Type, &grid_dims,
                            &PyArray_Type, &voxel_data,
                            &PyList_Type, &voxel_O_coords,
                            &PyList_Type, &voxel_quarts))
        {
            return NULL; /* raise argument parsing exception*/
        }

    unsigned int nx = *(int *)PyArray_GETPTR1(grid_dims, 0);
    unsigned int ny = *(int *)PyArray_GETPTR1(grid_dims, 1);
    unsigned int nz = *(int *)PyArray_GETPTR1(grid_dims, 2);
    unsigned max_voxel_index = nx * ny * nz;
    unsigned int addx = ny * nz;
    unsigned int addy = nz;
    unsigned int addz = 1;
    //PyObject *curr_voxel_coords;
    //PyObject *curr_voxel_quarts;
    //printf("grid dims: %i %i %i frames\n", nx, ny, nz);
    //printf("DEBUG2: nx %d ny %d nz %d\n", nx, ny, nz);
    for (voxel = 0; voxel < max_voxel_index; voxel++)
    {
        //
        int numplane = voxel / addx;
        double nw_total = *(double *)PyArray_GETPTR2(voxel_data, voxel, 4);
        nwtt += nw_total;

        //double W_dens = 1.0 * N_waters_[voxel] / (NFRAME_*Vvox);
        //gO[voxel] = W_dens / BULK_DENS_;
        //printf("voxel: %i %i\n", voxel, nw_total);
        double voxel_dens = 1.0 * nw_total / (num_frames * voxel_vol);
        *(double *) PyArray_GETPTR2(voxel_data, voxel, 5) += voxel_dens / ref_dens;
        //printf("DEBUG2 voxel %d gO %g\n", voxel, voxel_dens/ref_dens);
        PyObject *curr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel);
        PyObject *curr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel);
        for (n0 = 0; n0 < (int) nw_total; n0++)
        {
              double NNd = 10000;
              double NNs = 10000;
              double NNr = 10000;
              int i0 = n0 * 3; // index over O coordinates
              int q0 = n0 * 4; // index over quar

              // access oxygen coordinates
              PyObject *VX0 = PyList_GetItem(curr_voxel_coords, i0);
              double vx0 = PyFloat_AsDouble(VX0);
              PyObject *VY0 = PyList_GetItem(curr_voxel_coords, i0 + 1);
              double vy0 = PyFloat_AsDouble(VY0);
              PyObject *VZ0 = PyList_GetItem(curr_voxel_coords, i0 + 2);
              double vz0 = PyFloat_AsDouble(VZ0);
              // access quaternions
              PyObject *QW0 = PyList_GetItem(curr_voxel_quarts, q0);
              double qw0 = PyFloat_AsDouble(QW0);
              PyObject *QX0 = PyList_GetItem(curr_voxel_quarts, q0 + 1);
              double qx0 = PyFloat_AsDouble(QX0);
              PyObject *QY0 = PyList_GetItem(curr_voxel_quarts, q0 + 2);
              double qy0 = PyFloat_AsDouble(QY0);
              PyObject *QZ0 = PyList_GetItem(curr_voxel_quarts, q0 + 3);
              double qz0 = PyFloat_AsDouble(QZ0);
              for (n1 = 0; n1 < (int) nw_total; n1++)
              if ( n1 != n0)
              {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(curr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(curr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(curr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(curr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(curr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(curr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(curr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }

                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                  //printf("DEBUG1: %g\n", rR);
                  if (rR > 0 && rR < NNr) NNr = rR;
              }
              if (nw_total > 1)
              {
                  if (NNr < 9999 && NNr > 0)
                  {
                      double dbl = log(NNr * NNr * NNr * nw_total / (3.0 * twopi));
                      //printf("DEBUG1: dbl %f\n", dbl);
                      *(double *) PyArray_GETPTR2(voxel_data, voxel, 10) += dbl;
                      dTSo += dbl;
                  }
              }
              bool cannotAddZ = (nz == 0 || ( voxel%nz == nz-1 ));
              bool cannotAddY = ((nz == 0 || ny-1 == 0) || ( voxel%(nz*(ny-1)+(numplane*addx)) < nz));
              bool cannotAddX = (voxel >= addx * (nx-1) && voxel < addx * nx );
              bool cannotSubZ = (nz == 0 || voxel%nz == 0);
              bool cannotSubY = ((nz == 0 || ny == 0) || (voxel%addx < nz));
              bool cannotSubX = ((nz == 0 || ny == 0) || (voxel >= 0 && voxel < addx));
              bool boundary = ( cannotAddZ || cannotAddY || cannotAddX ||
                                cannotSubZ || cannotSubY || cannotSubX );
              //printf("DEBUG2: boundary= %d\n", boundary);
              //TODO: Replace this massive code repetition with a reusable function
              if (!boundary)
              {
                PyObject *nbr_voxel_coords, *nbr_voxel_quarts;
                double n1_total;
                /* Iterate over neighbor voxel in +Z direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel + addz, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel + addz);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel + addz);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in +Y direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel + addy, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel + addy);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel + addy);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in +X direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel + addx, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel + addx);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel + addx);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in -Z direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel - addz, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel - addz);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel - addz);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in -Y direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel - addy, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel - addy);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel - addy);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in -X direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel - addx, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel - addx);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel - addx);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in +Z +Y direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel + addz + addy, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel + addz + addy);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel + addz + addy);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in +Z -Y direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel + addz - addy, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel + addz - addy);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel + addz - addy);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in +Z +Y direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel + addz + addy, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel + addz + addy);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel + addz + addy);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in +Z -Y direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel + addz - addy, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel + addz - addy);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel + addz - addy);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in -Z +Y direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel - addz + addy, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel - addz + addy);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel - addz + addy);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in -Z -Y direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel - addz - addy, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel - addz - addy);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel - addz - addy);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in +Z +X direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel + addz + addx, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel + addz + addx);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel + addz + addx);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in +Z -X direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel + addz - addx, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel + addz - addx);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel + addz - addx);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in -Z +X direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel - addz + addx, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel - addz + addx);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel - addz + addx);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in -Z -X direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel - addz - addx, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel - addz - addx);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel - addz - addx);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in +Y +X direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel + addy + addx, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel + addy + addx);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel + addy + addx);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in +Y -X direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel + addy - addx, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel + addy - addx);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel + addy - addx);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in -Y +X direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel - addy + addx, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel - addy + addx);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel - addy + addx);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }
                /* Iterate over neighbor voxel in -Y -X direction
                */
                n1_total = *(double *)PyArray_GETPTR2(voxel_data, voxel - addy - addx, 4);
                nbr_voxel_coords = PyList_GetItem(voxel_O_coords, voxel - addy - addx);
                nbr_voxel_quarts = PyList_GetItem(voxel_quarts, voxel - addy - addx);

                for (n1 = 0; n1 != (int)n1_total; n1++)
                {
                  // access oxygen coordinates
                  int i1 = n1 * 3; // index over O coordinates
                  int q1 = n1 * 4; // index over quar
                  PyObject *VX1 = PyList_GetItem(nbr_voxel_coords, i1);
                  double vx1 = PyFloat_AsDouble(VX1);
                  PyObject *VY1 = PyList_GetItem(nbr_voxel_coords, i1 + 1);
                  double vy1 = PyFloat_AsDouble(VY1);
                  PyObject *VZ1 = PyList_GetItem(nbr_voxel_coords, i1 + 2);
                  double vz1 = PyFloat_AsDouble(VZ1);
                  // access quaternions
                  PyObject *QW1 = PyList_GetItem(nbr_voxel_quarts, q1);
                  double qw1 = PyFloat_AsDouble(QW1);
                  PyObject *QX1 = PyList_GetItem(nbr_voxel_quarts, q1 + 1);
                  double qx1 = PyFloat_AsDouble(QX1);
                  PyObject *QY1 = PyList_GetItem(nbr_voxel_quarts, q1 + 2);
                  double qy1 = PyFloat_AsDouble(QY1);
                  PyObject *QZ1 = PyList_GetItem(nbr_voxel_quarts, q1 + 3);
                  double qz1 = PyFloat_AsDouble(QZ1);

                  double dd = dist_squared(vx0, vy0, vz0, vx1, vy1, vz1);
                  if (dd < NNd && dd > 0) { NNd = dd; }
                  double rR = 2 * acos(qw0 * qw1 +
                                       qx0 * qx1 +
                                       qy0 * qy1 +
                                       qz0 * qz1 );
                  double ds = rR * rR + dd;
                  if (ds < NNs && ds > 0) {NNs = ds; }
                  //printf("DEBUG2: voxel=%i water=%i\n", voxel, n0);
                  //printf("DEBUG2: self NNd=%f NNs=%f\n", NNd, NNs);
                 }

                NNd = sqrt(NNd);
                NNs = sqrt(NNs);

                if (NNd < 3 && NNd > 0)
                {
                  double dbl = log((NNd * NNd * NNd * num_frames * 4 * pi * ref_dens) / 3);
                  *(double *) PyArray_GETPTR2(voxel_data, voxel, 8) += dbl;
                  dTSt += dbl;
                  dbl = log((NNs * NNs * NNs * NNs * NNs * NNs * num_frames * pi * ref_dens) / 48);
                  *(double *) PyArray_GETPTR2(voxel_data, voxel, 12) += dbl;
                  dTSs += dbl;
                }
              }
         } // end loop over waters in this voxel

        double dTStrans_norm = *(double *)PyArray_GETPTR2(voxel_data, voxel, 8);
        double dTSorient_norm = *(double *)PyArray_GETPTR2(voxel_data, voxel, 10);
        double dTSsix_norm = *(double *)PyArray_GETPTR2(voxel_data, voxel, 12);
        if (dTSorient_norm != 0)
        {
            *(double *)PyArray_GETPTR2(voxel_data, voxel, 10) = gas_kcal * temp * ((dTSorient_norm/nw_total) + euler_masc);
            *(double *)PyArray_GETPTR2(voxel_data, voxel, 9) = *(double *)PyArray_GETPTR2(voxel_data, voxel, 10) * nw_total / (num_frames * voxel_vol);

        }

        dTSorienttot += *(double *)PyArray_GETPTR2(voxel_data, voxel, 9);

        if (dTStrans_norm != 0)
        {
          nwts += nw_total;
          *(double *) PyArray_GETPTR2(voxel_data, voxel, 8) = gas_kcal * temp * ((dTStrans_norm / nw_total) +
                                                                     euler_masc);
          *(double *) PyArray_GETPTR2(voxel_data, voxel, 12) = gas_kcal * temp * ((dTSsix_norm / nw_total) +
                                                                   euler_masc);
        }

        *(double *) PyArray_GETPTR2(voxel_data, voxel, 7) = *(double *) PyArray_GETPTR2(voxel_data, voxel, 8) * nw_total / (num_frames * voxel_vol);
        *(double *) PyArray_GETPTR2(voxel_data, voxel, 11) = *(double *) PyArray_GETPTR2(voxel_data, voxel, 12) * nw_total / (num_frames * voxel_vol);
        dTStranstot += *(double *) PyArray_GETPTR2(voxel_data, voxel, 7);
    } // end loop over all grid points
    dTStranstot *= voxel_vol;
    dTSorienttot *= voxel_vol;
    double dTSst = 0.0;
    double dTStt = 0.0;
    if (nwts > 0)
    {
        dTSst = gas_kcal * temp * ((dTSs / nwts) + euler_masc);
        dTStt = gas_kcal * temp * ((dTSt/nwts) + euler_masc);
    }
    printf("Total referenced orientational entropy of the grid:"
                    " dTSorient = %9.5f kcal/mol, Nf=%d\n", dTSorienttot, num_frames);
    double dTSot = gas_kcal * temp * ((dTSo/nwtt) + euler_masc);
    printf("watcount in vol = %d\n", nwtt);
    printf("watcount in subvol = %d\n", nwts);
    printf("Total referenced translational entropy of the grid:"
                        " dTStrans = %9.5f kcal/mol, Nf=%d\n", dTStranstot, num_frames);
    printf("Total 6d if all one vox: %9.5f kcal/mol\n", dTSst);
    printf("Total t if all one vox: %9.5f kcal/mol\n", dTStt);
    printf("Total o if all one vox: %9.5f kcal/mol\n", dTSot);




    return Py_BuildValue("i", 0);
}


PyObject *_sstmap_ext_get_dist_matrix(PyObject *self, PyObject *args)
{
    int nwtot, n, l;
    double dR, nx, ny, nz, lx, ly, lz;
    PyArrayObject *dist_matrix;
    PyArrayObject *wat_coords;
    // Argument parsing to reterive everything sent from Python correctly    
    if (!PyArg_ParseTuple(args, "iO!O!",
                            &nwtot,
                            &PyArray_Type, &dist_matrix,
                            &PyArray_Type, &wat_coords))
        {
            return NULL; /* raise argument parsing exception*/
        }
    // for each water in the voxel
    for (n = 0; n < nwtot; n++)
    {
        nx = *(double *)PyArray_GETPTR2(wat_coords, n, 0);
        ny = *(double *)PyArray_GETPTR2(wat_coords, n, 1);
        nz = *(double *)PyArray_GETPTR2(wat_coords, n, 2);

        for (l = 0; l < nwtot; l++)
        {
            if(l == n) continue;
            //printf("Calculating orientational distancce between water: %i and %i\n", l, n);
            lx = *(double *)PyArray_GETPTR2(wat_coords, l, 0);
            ly = *(double *)PyArray_GETPTR2(wat_coords, l, 1);
            lz = *(double *)PyArray_GETPTR2(wat_coords, l, 2);
            dR = dist(nx, ny, nz, lx, ly, lz);
            *(double *)PyArray_GETPTR2(dist_matrix, n, l) = dR;
        }
    }
    return Py_BuildValue("i", 1);
}

PyObject *_sstmap_ext_calculate_energy(PyObject *self, PyObject *args)
{
    PyArrayObject *dist, *chg, *acoeff, *bcoeff;
    int solvent_at_sites, n_atoms, wat, i, j, at_i;
    double *d, *a, *b, *c;
    double d_sqrt, d_inv, d6, d12;
    //PyObject* nbrs = PyList_New(100);

    if (!PyArg_ParseTuple(args, "iO!O!O!O!",
        &wat,
        &PyArray_Type, &dist,
        &PyArray_Type, &chg,
        &PyArray_Type, &acoeff,
        &PyArray_Type, &bcoeff))
    {
        return NULL;
    }

    solvent_at_sites = PyArray_DIM(dist, 0);
    n_atoms = PyArray_DIM(dist, 1);
    //printf("The number of solvent atoms = %i\n", solvent_at_sites);
    //printf("The number of target atoms  = %i\n", n_atoms);

    // for each water in the voxel
    for (i = 0; i < solvent_at_sites; i++)
    {
        at_i = wat + i;
        for (j = 0; j < n_atoms; j++)
        {
            if (at_i == j) continue;
            d = (double *) PyArray_GETPTR2(dist, i, j);
            a = (double *) PyArray_GETPTR2(acoeff, i, j);
            b = (double *) PyArray_GETPTR2(bcoeff, i, j);
            c = (double *) PyArray_GETPTR2(chg, i, j);
            //printf("Acoeff between %i and %i is %f\n", at_i, j, *a);
            //printf("Bcoeff between %i and %i is %f\n", at_i, j, *b);
            //printf("CHP between %i and %i is %f\n", i, j, *c);
            d_sqrt = sqrt(*d);
            d_inv = 1.0 / *d;
            d6 = d_inv * d_inv * d_inv;
            d12 = d6 * d6;
            //e_vdw += (*a * d12) + (*b * d6);
            //e_elec += *c * 1/d_sqrt;
            *a *= d12;
            *a -= (*b * d6);
            *c /= d_sqrt;
            //if (at_i > 4614 && j >= 4614)
            //{
            //    printf("Dist of elec energy between %i and %i is %f\n", at_i, j, d_sqrt);
            //    //printf("Pairwise elec energy between %i and %i is %f\n", at_i, j, *c);
            //}
        }
    }

    return Py_BuildValue("i", 1);

}


/* Method Table
 * Registering all the functions that will be called from Python
 */

static PyMethodDef _sstmap_ext_methods[] = {
    {
        "assign_voxels",
        (PyCFunction)_sstmap_ext_assign_voxels,
        METH_VARARGS,
        "Process grid"
    },
    
    {
        "get_pairwise_distances",
        (PyCFunction)_sstmap_ext_get_pairwise_distances,
        METH_VARARGS,
        "get distance matrix"
    },

    {
        "getNNOrEntropy",
        (PyCFunction)_sstmap_ext_getNNOrEntropy,
        METH_VARARGS,
        "get voxel entropy"
    },    

    {
        "getNNTrEntropy",
        (PyCFunction)_sstmap_ext_getNNTrEntropy,
        METH_VARARGS,
        "get voxel entropy"
    },  
    {
        "calculate_energy",
        (PyCFunction)_sstmap_ext_calculate_energy,
        METH_VARARGS,
        "get energy"
    },
    {
        "get_dist_matrix",
        (PyCFunction)_sstmap_ext_get_dist_matrix,
        METH_VARARGS,
        "get voxel entropy"
    },
    {NULL, NULL, 0, NULL}
};

/* Initialization function for this module
 */
//PyMODINIT_FUNC
DL_EXPORT(void) init_sstmap_ext(void) // init function has the same name as module, except with init prefix
{
    // we produce name of the module, method table and a doc string
    Py_InitModule3("_sstmap_ext", _sstmap_ext_methods, "Process GIST calcs.\n");
    import_array(); // required for Numpy initialization
}