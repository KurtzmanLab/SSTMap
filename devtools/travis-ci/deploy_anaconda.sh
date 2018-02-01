#!/usr/bin/env bash
# Deploy to anaconda solvation tools channel
conda install --yes anaconda-client
pushd .
cd $HOME/travis/miniconda3/conda-bld
ls
ls */
FILES=*/${PACKAGENAME}-*.tar.bz2
for filename in $FILES; do
    anaconda -t $CONDA_UPLOAD_TOKEN remove --force ${ORGNAME}/${PACKAGENAME}-dev/${filename}
    anaconda -t $CONDA_UPLOAD_TOKEN upload --force -u ${ORGNAME} -p ${PACKAGENAME}-dev ${filename}
done
popd

#anaconda upload /home/travis/miniconda3/conda-bld/linux-64/sstmap-1.1.0-py36_0.tar.bz2
