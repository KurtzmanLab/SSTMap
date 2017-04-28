#!/bin/bash

#cp -r $RECIPE_DIR/../../sstmap $SRC_DIR
$PYTHON setup.py clean
$PYTHON setup.py install

g++ -o /usr/local/bin/bruteclust $RECIPE_DIR/../../scripts/make_clust_brute.cpp
g++ -o /usr/local/bin/kdhsa102 $RECIPE_DIR/../../scripts/kdhsa102.cpp $RECIPE_DIR/../../scripts/kdhsa102_main.cpp
g++ -o /usr/local/bin/6dimprobable $RECIPE_DIR/../../scripts/6dimprobable.cpp $RECIPE_DIR/../../scripts/6dim_main.cpp

# Add more build steps here, if  necessary.

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
