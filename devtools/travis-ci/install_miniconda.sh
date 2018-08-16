#!/bin/bash
# Credits: Adapted from https://github.com/choderalab/pymbar/blob/master/devtools/travis-ci/install.sh
# with some modifications
pushd .
cd $HOME

# Install Miniconda
MINICONDA=Miniconda2-latest-Linux-x86_64.sh
if [[ "$TRAVIS_OS_NAME" == "osx" ]];   then MINICONDA=Miniconda2-latest-MacOSX-x86_64.sh; fi

MINICONDA_HOME=$HOME/miniconda
MINICONDA_MD5=$(curl -s https://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget -q http://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    exit 1
fi
bash $MINICONDA -b -p $MINICONDA_HOME

# Configure miniconda
export PIP_ARGS="-U"
export PATH=$MINICONDA_HOME/bin:$PATH
conda update --yes conda
conda install --yes conda-build jinja2 anaconda-client
popd
