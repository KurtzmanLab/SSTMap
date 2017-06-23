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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>

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
    return 1.0/(sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2)));
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