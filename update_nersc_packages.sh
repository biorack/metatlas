#This must say /project/projectdirs/metatlas/anaconda/bin/python or there will be unexpected consequences
#which python

#as pasteur
#Clone metatlas source
#git clone https://github.com/biorack/metatlas.git
#this goes in /project/projectdirs/metatlas/projects

#Install to anaconda/lib/python2.7/site-packages all metatlas file and all the "REQUIRES" in setup.py
cd /project/projectdirs/metatlas/projects/metatlas
git pull
pip uninstall metatlas
pip install . --upgrade

#install cobrapy
#pip install cobra
#pip install python-libsbml
# I don't know which solvers are bunded with pip install cobra, but the github clone of cobrapy will not install any solvers.
# when you do "pip install cobra" solvers are bundled with this installation.

#installing qsopt
# this solver can be installed at nersc: https://github.com/jonls/qsopt-ex
# clone the repo from there:
# cd projects/qsopt-ex/
# ./bootstrap 
# mkdir build
# cd build/
# ../configure --prefix=/project/projectdirs/metatlas/
# make -j4
# make install
# add this to the path:
# import os
# os.environ['PATH'] = os.environ['PATH']+':/project/projectdirs/metatlas/bin'
# you can choose a solver when running a command in cobra: https://cobrapy.readthedocs.io/en/stable/solvers.html

## these might be necessary
#pip install gspread

#Install qgrid from https://github.com/quantopian/qgrid
# clone to ../metatlas/projects
# pip install .
# run the pip command from the qgrid folder.

#The remaining packages are installed like this:
#/project/projectdirs/metatlas/anaconda/bin/conda install -c https://conda.anaconda.org/rdkit rdkit

#Set permissions for metatlas users
chmod -R 750 /project/projectdirs/metatlas/anaconda/lib/python2.7/site-packages/
chgrp -R metatlas /project/projectdirs/metatlas/anaconda/lib/python2.7/site-packages/
