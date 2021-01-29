from __future__ import absolute_import
import sys,os
import numpy as np
import pandas as pd
sys.path.append('/global/project/projectdirs/openmsi/jupyterhub_libs/anaconda/lib/python2.7/site-packages')
from rdkit import Chem
from rdkit.Chem import PandasTools
sdf_file = '/project/projectdirs/openmsi/projects/compound_data/chembl/chembl_21.sdf.gz'
df = PandasTools.LoadSDF(sdf_file)
df.to_pickle('/project/projectdirs/openmsi/projects/compound_data/chembl/chembl.pkl')
