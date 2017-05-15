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
    return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
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
    return pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
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

/*
Calculates electrostatic energy of a query water molecule against a set of target atoms

*/
PyObject *_sstmap_ext_assign_voxels(PyObject *self, PyObject *args)
{

    int n_frames, i_frame, n_wat, i_wat;
    PyArrayObject *coords, *grid_dim, *grid_max, *grid_orig, *wat_oxygen_ids;
    PyObject *frame_data;

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

    // declare local variables
    double grid_max_x, grid_max_y, grid_max_z;
    double grid_orig_x, grid_orig_y, grid_orig_z;
    int grid_dim_x, grid_dim_y, grid_dim_z;
    int grid_index_x, grid_index_y, grid_index_z;
    
    double wat_translated_x, wat_translated_y, wat_translated_z;


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
            int wat_id; // id of current water
            float *wat_x, *wat_y, *wat_z; // coordinates
            int dims[1];
            dims[0] = 2;
            PyArrayObject *wat_data;
            PyObject *curr_voxel;
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
                        int voxel_id;
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
                        //free(curr_voxel);
                        //free(wat_data);
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
    float *b_x, *b_y, *b_z;
    b_x = (float *) PyArray_GETPTR2(periodic_box, 0, 0);
    b_y = (float *) PyArray_GETPTR2(periodic_box, 0, 1);
    b_z = (float *) PyArray_GETPTR2(periodic_box, 0, 2);
    //printf("Unit cell dimensions: %f %f %f\n", *b_x, *b_y, *b_z);

    int wat_sites, wat_atom, wat_atom_id;
    wat_sites = PyArray_DIM(dist_array, 0);

    int num_target_at, target_at, target_at_id;
    num_target_at = PyArray_DIM(target_at_ids, 0);

    int frame = 0;
    double d;

    for (wat_atom = 0; wat_atom < wat_sites; wat_atom++)
    {
        float *wat_x, *wat_y, *wat_z;
        wat_atom_id = *(int *) PyArray_GETPTR1(wat, 1) + wat_atom;
        wat_x = (float *) PyArray_GETPTR3(coords, frame, wat_atom_id, 0);
        wat_y = (float *) PyArray_GETPTR3(coords, frame, wat_atom_id, 1); 
        wat_z = (float *) PyArray_GETPTR3(coords, frame, wat_atom_id, 2);

        for (target_at = 0; target_at < num_target_at; target_at++)
        {
            float *target_at_x, *target_at_y, *target_at_z;
            target_at_id = *(int *) PyArray_GETPTR1(target_at_ids, target_at);
            target_at_x = (float *) PyArray_GETPTR3(coords, frame, target_at_id, 0);
            target_at_y = (float *) PyArray_GETPTR3(coords, frame, target_at_id, 1); 
            target_at_z = (float *) PyArray_GETPTR3(coords, frame, target_at_id, 2);
            //printf("Iterator: %d, atom id: %d\n", target_at, target_at_id);
            //printf("Water atom coords %f %f %f\n", *wat_x, *wat_y, *wat_z);
            //printf("Target atom coords %f %f %f\n", *target_at_x, *target_at_y, *target_at_z);
            d = dist_mic(*wat_x, *wat_y, *wat_z, *target_at_x, *target_at_y, *target_at_z, *b_x, *b_y, *b_z);
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
    int nwtot, n_frames, n, l;
    double NNtr, dR, wat_tr_ent;
    double voxel_dTStr = 0.0;
    double nx, ny, nz, lx, ly, lz;
    PyArrayObject *voxel_wat_coords; 
    double twopi = 2*M_PI;
    // Argument parsing to reterive everything sent from Python correctly    
    if (!PyArg_ParseTuple(args, "iiO!",
                            &nwtot,
                            &n_frames,
                            &PyArray_Type, &voxel_wat_coords))
        {
            return NULL; /* raise argument parsing exception*/
        }
    // for each water in the voxel
    for (n = 0; n < nwtot; n++){
        NNtr = 10000;
        nx = *(double *)PyArray_GETPTR2(voxel_wat_coords, n, 0);
        ny = *(double *)PyArray_GETPTR2(voxel_wat_coords, n, 1);
        nz = *(double *)PyArray_GETPTR2(voxel_wat_coords, n, 2);

        for (l = 0; l < nwtot; l++){
            if(l == n) continue;
            //printf("Calculating orientational distancce between water: %i and %i\n", l, n);
            lx = *(double *)PyArray_GETPTR2(voxel_wat_coords, l, 0);
            ly = *(double *)PyArray_GETPTR2(voxel_wat_coords, l, 1);
            lz = *(double *)PyArray_GETPTR2(voxel_wat_coords, l, 2);
            dR = dist(nx, ny, nz, lx, ly, lz);
            // get translational nearest neighbor            
            if (dR>0 && dR<NNtr) NNtr = dR;
        }
        //calculate translational entropy
        if (NNtr<3 && NNtr>0) {
            //printf("Nearest neighbour translational distance: %f\n", NNtr);
            //wat_tr_ent = log(nwtot*NNtr*NNtr*NNtr/(3.0*twopi));
            //voxel_dTStr_norm += wat_tr_ent;
            wat_tr_ent = log((NNtr*NNtr*NNtr*n_frames*4*twopi*0.0334)/3);
            voxel_dTStr += wat_tr_ent;            
        }        
    }
    // 
    //*(double *)PyArray_GETPTR1(ent, 2) += voxel_dTSor_norm;
    return Py_BuildValue("f", voxel_dTStr);
}

static PyObject * _sstmap_ext_writepdb(PyObject * self, PyObject * args)
{
    PyObject *data, *l;
    int i, n;

    if (!PyArg_ParseTuple(args, "iO",
        &n,
        &PyList_Type, &data))
        {
            return NULL; /* raise argument parsing exception*/
        }

    FILE *fp = fopen("within5Aofligand.pdb", "ab");
    if (fp != NULL)
    {
        for (i = 0; i < n; i++)
            {
                l = PyList_GetItem(data, i);
                fputs(l, fp);
            }
    }
    fclose(fp);
    
    return Py_BuildValue("i", 1);

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
        "Process grid"
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
        "get_dist_matrix",
        (PyCFunction)_sstmap_ext_get_dist_matrix,
        METH_VARARGS,
        "get voxel entropy"
    },
    {
        "write_pdb",
        (PyCFunction)_sstmap_ext_writepdb,
        METH_VARARGS,
        "Write function for large data"

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