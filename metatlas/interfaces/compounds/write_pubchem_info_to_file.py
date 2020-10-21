
# coding: utf-8

# In[1]:


import os
import pandas as pd
import requests
import json
import requests
import multiprocessing as mp
from random import randint
from time import sleep

from rdkit import Chem

df = pd.read_csv('/project/projectdirs/metatlas/projects/unique_compounds.csv.gz')
df.rename(columns = {'monoisotopoic_mw':'monoisotopic_mw'},inplace=True)


df = df.convert_objects(convert_numeric=True)
df.keys()

inchi_keys = df[(~ pd.isnull(df.inchi_key))].inchi_key


out_dir = '/project/projectdirs/metatlas/projects/pubchem_info/'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)


# # remove those that have already been done

import glob
files = glob.glob(os.path.join(out_dir,'*.json'))
done_inchi_key = []
for f in files:
    done_inchi_key.append(os.path.basename(f).split('.')[0])
inchi_keys = list(set(inchi_keys) - set(done_inchi_key))

print len(inchi_keys)

def write_pubchem_info_to_file(ik):
    suffix = '.json'
    fname = os.path.join(out_dir, ik + suffix)
    # this query will return 
    # {"Fault": {"Message": "No records found for the given CID(s)", "Code": "PUGREST.NotFound"}}
    # if 3d data is not available from pubchem.
    # leaving it out does not get you 3d. It defaults to 2d if left out.  Have to run it twice.
    #
    # missing inchikeys will return
    # {"Fault": {"Message": "No CID found", "Code": "PUGREST.NotFound", "Details": ["No CID found that matches the given InChI key"]}}
    #
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/%s/json?record_type=3d'%ik
    response = requests.get(url)
    d = response.json()
    if 'Fault' in d.keys():
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/%s/json'%ik
        response = requests.get(url)
        d = response.json()
    if not "Fault" in d.keys():
        with open(fname, 'w') as fid:
            json.dump(d, fid)

#run this after doing parralel with 4 cores.
# any more cores will cause an error
for ik in inchi_keys:
    write_pubchem_info_to_file(ik)

# run it with 4 cores, then rerun above loop doing missing one-by-one
#pool = mp.Pool(processes=4)
#pool.map(write_pubchem_info_to_file, inchi_keys)


