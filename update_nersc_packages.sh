#This must say /project/projectdirs/metatlas/anaconda/bin/python or there will be unexpected consequences
#which python

#as pasteur
#Clone metatlas source
#git clone https://github.com/biorack/metatlas.git
#this goes in /project/projectdirs/metatlas/projects

#Install to anaconda/lib/python2.7/site-packages all metatlas file and all the "REQUIRES" in setup.py
cd /project/projectdirs/metatlas/projects/metatlas
git pull
pip install .

#The remaining packages are installed like this:
#pip install git+https://github.com/quantopian/qgrid
#/project/projectdirs/metatlas/anaconda/bin/conda install -c https://conda.anaconda.org/rdkit rdkit

#Set permissions for metatlas users
chmod -R 750 /project/projectdirs/metatlas/anaconda/lib/python2.7/site-packages/
chgrp -R metatlas /project/projectdirs/metatlas/anaconda/lib/python2.7/site-packages/
