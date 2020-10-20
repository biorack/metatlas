#!/bin/bash
# Script to create an anaconda environment with all OpenMSI 
# packages installed

BASEDIR='/global/common/software/m2650/nersc_python_install'
CONDAPATH="/global/common/software/m2650/python-"$NERSC_HOST
echo $CONDAPATH

# Activate our conda environment
source $CONDAPATH/bin/activate

nosetests --exe -v metatlas
