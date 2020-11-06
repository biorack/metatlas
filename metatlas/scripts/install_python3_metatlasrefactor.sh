#!/bin/bash

BASEDIR='/global/common/software/m2650/nersc_python_install'
ENVPATH="/global/common/software/m2650/python3-metatlas-"$NERSC_HOST
mkdir -p $ENVPATH

#cd $BASEDIR


echo "This conda will be used."
which conda
conda env create --prefix $ENVPATH -f /global/homes/b/bpb/repos/metatlas/metatlas_env.yml

conda activate $ENVPATH
which conda
which python

echo "Installing metatlas setup.py now."
cd /global/homes/b/bpb/repos/metatlas
pip install .

echo "metatlas setup created."
echo ""
echo "to use this type:"
echo "conda activate "$ENVPATH
echo "to activate"

#which conda
cd $BASEDIR
chmod -R 755 $ENVPATH


