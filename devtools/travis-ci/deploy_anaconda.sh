#!/usr/bin/env bash
# Deploy to anaconda solvation tools channel
conda install --yes anaconda-client
pushd .
cd $HOME/miniconda/conda-bld
FILES=*/${PACKAGENAME}-dev-*.tar.bz2
for filename in $FILES; do
    anaconda -t $CONDA_UPLOAD_TOKEN remove --force ${ORGNAME}/${PACKAGENAME}-dev/${filename}
    anaconda -t $CONDA_UPLOAD_TOKEN upload --force -u ${ORGNAME} -p ${PACKAGENAME}-dev ${filename}
done
popd