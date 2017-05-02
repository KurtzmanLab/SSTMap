#!/bin/bash

#cp -r $RECIPE_DIR/../../sstmap $SRC_DIR
$PYTHON setup.py clean
$PYTHON setup.py install


# Add more build steps here, if  necessary.
# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
