import json
import requests
import os
import pandas as pd
import time

def write_classyfire_info_to_file(input_tuple,wait=True):
    outdir = input_tuple[0]
    inchi_key = input_tuple[1]
    suffix = '.json'
    fname = os.path.join(outdir, inchi_key + suffix)
    if not os.path.isfile(fname):
        txt = get_classyfire_from_inchikey(inchi_key)
        with open(fname, 'w') as fid:
            json.dump(txt, fid)
        if wait==True:
            time.sleep(2)

def get_classyfire_from_inchikey(inchi_key):    
    """
    example: 
    http://classyfire.wishartlab.com/entities/ZCPXJLJVCGFYTC-UHFFFAOYSA-N.json
    """
    url = 'http://classyfire.wishartlab.com/entities/%s.json'%inchi_key
    response = requests.get(url)
    return response.text

def read_classyfire_json(filename):
    if os.path.isfile(filename):
        with open(filename, 'r') as fid:
            return json.load(fid)
    else:
        return None


def get_classyfire_classes(inchi_keys,basepath='/global/projectb/sandbox/metatlas/projects/classyfire_annotations/'):
    cf_hits = []
    cf_keys = ['url', 'name', 'chemont_id', 'description']
    for ik in inchi_keys:
        my_file = os.path.join(basepath,'%s.json'%ik)
        temp = {'inchi_key':ik}
        if os.path.isfile(my_file):
    #         print('%s %s'%(ik,name))
            with open(my_file,'r') as fid:
                t = fid.read()
            classykeys = [ "kingdom","superclass","class", "subclass"]
            j = json.loads(json.loads(t))

            for k in classykeys:
                if (k in j) and (j[k] is not None):
                    for m in cf_keys:
                        temp['%s_%s'%(k,m)] = j[k][m]
        cf_hits.append(temp)
    cf = pd.DataFrame(cf_hits)
    cf = cf[~pd.isna(cf['kingdom_name'])]
    cf.reset_index(inplace=True,drop=True)
    return cf

def make_new_classyfire_entries(not_found):
    ### submit new molecules to classyfire
    dummy_script = """curl -is http://classyfire.wishartlab.com/queries.json -X POST -d '{"label":"asdf","query_input":"query_str", "query_type":"STRUCTURE"}' -H "Content-Type: application/json"
    """

    cf_str = []
    max_cpds = 10
    counter = 0
    all_scripts = []
    for i,row in not_found.iterrows():
    #     cf_str.append('%s\\t%s'%(row.inchi_key,Chem.MolToSmiles(Chem.MolFromInchi(row.inchi_x),isomericSmiles=True)))
        cf_str.append('%s\\t%s'%(row.inchi_key,row.inchi_x))
        counter += 1
        if counter > max_cpds:
            all_scripts.append(dummy_script.replace('query_str','\\n'.join(cf_str)))
            counter = 0
            cf_str = []
    all_scripts.append(dummy_script.replace('query_str','\\n'.join(cf_str)))

    with open('classyfire_remaining_smiles_397.sh','w') as fid:
        for script in all_scripts:
            fid.write('%s\n'%script)
            fid.write('sleep 10s\n')

            
def make_classyfire_table(classyfire_path='/global/projectb/sandbox/metatlas/projects/classyfire_annotations/',
                          output_filename='/global/homes/b/bpb/Downloads/classyfire_all_compounds.csv.gz',
                         compression=True
                         iks=None):
    files = glob.glob(os.path.join(classyfire_path,'*.json'))
    if iks is None:
        iks = [os.path.basename(f).split('.')[0] for f in files]
    cf = get_classyfire_classes(iks)
    if compression is True:
        cf.to_csv(output_filename,index=None,compression='gzip')
    else:
        cf.to_csv(output_filename,index=None)
    return cf

def make_classyfire_network(classyfire_obo_filename='ChemOnt_2_1.obo',
                            nodetable_filename='classyfire_chemont_nodetable.tab',
                            network_filename= "classyfire_ontology_network_template.graphml",
                           parentchild_filename='classyfire_chemont_network.tab',
                           save=True):
    from copy import deepcopy
    import pandas as pd
    import numpy as np
    import networkx as nx

    with open(classyfire_obo_filename,'r') as fid:
        chemont = fid.read()
    chemont = chemont.split('[Term]')
    chemont = [c.strip() for c in chemont]
    chemont.pop(0) #remove header

    print('There are %d entries'%len(chemont))

    # make an empty dict that has all possible chemont terms
    attributes = {}
    for c in chemont:
        for a in c.split('\n'):
            attributes[a.split(': ')[0]] = np.nan

    chemont_df = []
    for c in chemont:
        chemont_df.append(deepcopy(attributes))
        chemont_df[-1]['synonym'] = ''
        for a in c.split('\n'):
            split_str = a.split(': ')
            attr = split_str[0]
            value = split_str[-1]
            # there are many synonyms for each entry.  make delimited list
            if attr in ['synonym','xref']:
                chemont_df[-1][attr] = '%s; %s'%(chemont_df[-1]['synonym'],value)

            if ' ! ' in value:
                chemont_df[-1][attr] = value.split(' ! ')[0].strip()
            else:
                chemont_df[-1][attr] = value.strip()

    chemont_df = pd.DataFrame(chemont_df)
    
    if save is True:
        chemont_df[['id','is_a']].to_csv(parentchild_filename,index=None,sep='\t')


    chemont_info = chemont_df[['id','name','def','synonym','xref','alt_id','comment']]
    chemont_info = chemont_info.drop_duplicates(['id','name'])
    chemont_info.rename(columns={'name':'ontology_name'},inplace=True)
    if save is True:
        chemont_info.to_csv(nodetable_filename,index=None,sep='\t')


    G=nx.from_pandas_edgelist(chemont_df, 'is_a', 'id')
    nx.set_node_attributes(G, chemont_info.set_index('id').to_dict('index'))
    G.remove_node('CHEMONTID:0000000')
    if save is True:
        nx.write_graphml_lxml(G,network_filename)
    return G

def make_ontology_gmt_file(classyfire_file,classyfire_file='/global/homes/b/bpb/Downloads/classyfire_all_compounds.csv.gz',gmt_file = '/Users/bpb/Downloads/classyfire_gmt_file.gmt'):
    cf = pd.read_csv(classyfire_file)#,usecols=['inchi_key','class_name'])

    cf_cols = [c for c in cf.columns if '_id' in c] + ['inchi_key']
    cf = cf[cf_cols]
    cf.fillna('',inplace=True)

    
    gmt = pd.concat([cf.melt(id_vars=[col],value_vars=['inchi_key']).rename(columns={col:'term'}) for col in cf_cols[:-1]],axis=0)
    gmt.drop(columns=['variable'],inplace=True)
    gmt.dropna(subset=['term'],inplace=True)
    gmt = gmt.groupby('term')['value'].apply(list).reset_index(name='compounds')
    gmt['dummy'] = ''
    gmt = gmt[['term','dummy','compounds']]

    exploded_columns = pd.DataFrame(gmt['compounds'].values.tolist())
    gmt.drop(columns=['compounds'],inplace=True)
    gmt = pd.concat([gmt,exploded_columns],axis=1)
    # gmt.head()
    gmt.to_csv(gmt_file,sep='\t',index=None,header=False)

    with open(gmt_file,'r') as fid:
        gmt = fid.read()

    gmt = [g.strip() for g in gmt.split('\n')]
    gmt = '\n'.join(gmt)

    with open(gmt_file,'w') as fid:
        fid.write('%s'%gmt)
#     cf_str = "HFKNJQYMAGMXTR-CAPHXMBKSA-N\tCC/C=C\\C[C@H](O)/C=C/[C@H]1[C@H](O)C[C@H](O)[C@@H]1CC(=O)CCCCC(=O)O\nDSLVJFBJCIYHLK-DHGAWWCCSA-N\tCOc1cc(-c2oc3cc(O)cc(O)c3c(=O)c2O[C@@H]2O[C@H](COC(C)=O)[C@@H](O)C(O)C2O)ccc1O"

#     payload = {'label':'value',
#                'query_input':cf_str,
#                'query_type':'STRUCTURE'}
#     r = requests.post('http://classyfire.wishartlab.com/entities/queries.json', data = payload)

# #filesize less than 60 is failure
# (58, 'HFWIMJHBCIGYFH-UHFFFAOYSA-N', 63)
# {"status":"500","error":"Internal Server Error"}
# ()
# (46, 'ZCPXJLJVCGFYTC-UHFFFAOYSA-N', 4780)
# {"status":"404","error":"Not Found"}
# ()
# (4, 'CXUXMSACCLYMBI-FPLPWBNLSA-N', 2247)
# {}
# print(df[df.filesize<60].shape)
# for f in df[df.filesize < 60].filename:
#     os.remove(f)