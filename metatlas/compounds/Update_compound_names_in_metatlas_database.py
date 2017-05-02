import sys
sys.path.insert(0,'/global/project/projectdirs/metatlas/anaconda/lib/python2.7/site-packages' )
from metatlas import metatlas_objects as metob
import pandas as pd

df = pd.read_pickle('/project/projectdirs/openmsi/projects/ben_run_pactolus/unique_compounds_updated_pubchem_info.pkl')
df = df[df.new_common_name.str.contains('[A-Za-z0-9]+')]

import time
list_of_updates = []
for i,row in df.iterrows():
    update_dict = dict(inchi_key=row.metatlas_inchi_key, 
                       name=row.new_common_name,
                       synonyms = row.new_synonyms,
                       pubchem_compound_id = '%d'%row.new_pubchem_compound_id,
                       pubchem_url = 'http://pubchem.ncbi.nlm.nih.gov/compound/%d'%row.new_pubchem_compound_id
                    )
    list_of_updates.append(update_dict)
    

#import multiprocessing as mp
db = metob.database
compounds = db['compounds']

def update_compound_table(input_dict):
    print 'updating compound'
    L = compounds.update(input_dict, ['inchi_key'])
    print 'done updating'

#pool = mp.Pool(processes=2)
t0 = time.time()
print 'updating first compounds now'
#data = pool.map(update_compound_table, list_of_updates[:4])
for i in range(100):
    update_compound_table(list_of_updates[i])
    print time.time() - t0

