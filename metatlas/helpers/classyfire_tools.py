import json
import requests
import os
import pandas as pd

def write_classyfire_info_to_file(input_tuple):
    outdir = input_tuple[0]
    inchi_key = input_tuple[1]
    suffix = '.json'
    fname = os.path.join(outdir, inchi_key + suffix)
    if not os.path.isfile(fname):
        txt = get_classyfire_from_inchikey(inchi_key)
        with open(fname, 'w') as fid:
            json.dump(txt, fid)

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


def get_classyfire_classes(inchi_keys):
    cf_hits = []
    cf_keys = ['url', 'name', 'chemont_id', 'description']
    for ik in inchi_keys:
        my_file = '/global/projectb/sandbox/metatlas/projects/classyfire_annotations/%s.json'%ik
        temp = {'inchi_key':ik}
        if os.path.isfile(my_file):
    #         print('%s %s'%(ik,name))
            with open(my_file,'r') as fid:
                t = fid.read()
            classykeys = [ "kingdom","superclass","class", "subclass"]
            j = json.loads(json.loads(t))

            for k in classykeys:
                if j[k] is not None:
                    for m in cf_keys:
                        temp['%s_%s'%(k,m)] = j[k]['name']
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