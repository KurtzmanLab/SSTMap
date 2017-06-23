/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   probable.h
 * Author: stevenramsey
 *
 * Created on March 16, 2016, 9:06 AM
 */

#ifndef PROBABLE_H
#define PROBABLE_H
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <string.h>
#endif /* PROBABLE_H */

static int D = 3;

struct point {
    int dim;
    double* x;
    point();
    void set_point(double* vals);
    void set_point(const point &p);
    //~point();
    void print_point();
    void zeros();
    void ones();
    //void set_dimension(int y);


};

bool same_points (point p, point q);

double dist(const point &p, const point &q);

/*struct box {
    point lo, hi; //diagonally opposite points in the box (min, max)
    //box () {} //empty normal constructor
    void set_box(const point &mylo, const point &myhi); //copy those points to be our lo and hi
};*/



struct boxnode {
    int mom, dau1, dau2, ptlo, pthi; //these are all integers which will work to point towards the specified thing in their data structure
    point lo, hi;
    boxnode();
    void set_boxnode(point mylo, point myhi, int mymom, int myd1, int myd2, int myptlo, int mypthi);

    /*
        Feed it 2 points and the necessary indices, save those indices and create a box from the points.
        In other words this is the data structure which actively creates the box data structure, but will be used
             recursively to create the entire tree
    */
};

double dist(const boxnode &b, const point &p, int d);

struct kdtree {
    static const double BIG; //this value is a placeholder for starting box size (will be absurd)
    int dim;
    int numbox, npts; //integer counts of boxes and points
    point* pts;
    boxnode *boxes;
    int* ptindx;
    int* rptindx; //point index and reverse point index
    int* nd;
    double* dn;
    double* coord;
    //int* within1;
    //double cenn[3];
    //double maxx[3];
    //double minn[3];
    kdtree(std::vector< double > &vals);
    ~kdtree();
    //utility functions for use after tree is constructed
    double disti(int jpt, int kpt);
    int locate(point pt);
    int locate(int jpt);
    //applications to use tree
    //int nearest(point pt);
    double dnearest(point pt);
    void nnearest(int jpt, int *nn, double *dn, int n);
    static void sift_down(double *heap, int *ndx, int nn);
    int locatenear(point pt, double r, int *v, int nmax);
    //double run_tree();
    double run_tree_orient();
    double run_tree_trans(std::vector<double > &cls);
    //void run_locate();
    //void print_boxes();
    //void print_tree(int y);
    //void print_box(int y);

};

