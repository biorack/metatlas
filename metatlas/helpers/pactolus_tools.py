from __future__ import print_function

import sys
from metatlas import metatlas_objects as metob
from metatlas import h5_query as h5q
import qgrid

from metatlas.helpers import metatlas_get_data_helper_fun as mgd

from matplotlib import pyplot as plt
import pandas as pd
import os
import tables
import pickle
import h5py
import dill
import numpy as np
import h5py
from requests import Session
import os.path
import glob as glob
import json


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
# from rdkit.Chem.rdMolDescriptors import ExactMolWt
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG,display

def import_network(network_file='network_v1p0.cyjs'):
    with open(network_file) as data_file:    
        data = json.load(data_file)
    print(data['elements']['nodes'][0]['data'].keys())
    network = {}
    network['data'] = data
    network['x'] = []
    network['y'] = []
    network['node_id'] = []
    network['node_name'] = []
    network['node_inchi_key'] = []
    for n in data['elements']['nodes']:
        network['x'].append(n['position']['x'])
        network['y'].append(n['position']['y'])
        network['node_id'].append(float(n['data']['SUID']))
        network['node_name'].append(n['data']['compound_name'])
        network['node_inchi_key'].append(n['data']['inchi_key'])
        
    network['x'] = np.asarray(network['x'])
    network['y'] = np.asarray(network['y'])
    network['node_id'] = np.asarray(network['node_id'])
    return network

def merge_sheets(my_sheets,score_cutoff=0.1,mz_cutoff=0.1):
    dfs = []
    for sheet in my_sheets:
        df = pd.read_csv(sheet,index_col=False)
        df['source_file'] = os.path.basename(sheet)
        df = filter_hits(df)
        df = df.sort_values(by='score').drop_duplicates(subset=['inchi_key'], keep='last')
        dfs.append(df)
        
    df_all_files = pd.concat(dfs)
    #print 'making key'
    #df_all_files['inchi_key'] = df_all_files.inchi.apply(lambda x: Chem.InchiToInchiKey(str(x)))
    print(df_all_files.keys())
    #df_all_files.set_index(['inchi','inchi_key','metatlas name'],inplace=True)
    df_all_files.set_index(['inchi_key','mass'],inplace=True)
    return df_all_files

def read_pacolus_results(pactolus_file,min_score=0.0):
    """
    This is a new 20161213 version of the pactolus file readers.
    The hope is to circumvent all the copies from the original hdf5 file.
    
    Input:
    pactolus_file: the full path to a conforming pactolus search result
    
    Output:
    scan_df:
    
    tree_df:
    
    """
    with h5py.File(pactolus_file,'r') as fid:
    #read score_matrix, convert all by all matrix to lists of scores
        idx = range(fid['score_matrix'].shape[0])  
        d = {'retention time':fid['scan_metadata']['peak_rt'][idx],
            'precursor intensity':fid['scan_metadata']['peak_intensity'][idx],
            'precursor mz':fid['scan_metadata']['peak_mz'][idx],
            'polarity': fid['scan_metadata']['polarity'][idx],
            'index': idx}
        scan_df = pd.DataFrame(d)
        scan_df['filename'] = pactolus_file

        m = fid['score_matrix'][:]
        hits = []
        for mm in m:
            idx = np.where(mm>min_score)[0]
            hits.append(sorted([(mm[i],i) for i in idx])[::-1])
        df = pd.DataFrame({'scores':hits})
        b_flat = pd.DataFrame([[i, x[0], x[1]] 
                           for i, y in df.scores.apply(list).iteritems() 
                           for x in y], columns=['index','score','compound']).set_index('index')
        scan_df = scan_df.merge(b_flat, how = 'outer',left_index=True, right_index=True)

        #get a list of True/False if any hits for a compound:
        f = np.any(m.T>min_score,axis=1)
        #only do this for ones that get a hit
        idx = np.where(f)[0]#range(fid['score_matrix'].shape[1])
        lookup = fid['tree_file_lookup_table'][:]
        d = {'filename': fid['tree_file_lookup_table']['filename'][idx],
             'ms1_mass': fid['tree_file_lookup_table']['ms1_mass'][idx],
             'inchi': fid['tree_file_lookup_table']['inchi'][idx],
             'permanent_charge': fid['tree_file_lookup_table']['permanent_charge'][idx],
              'index': idx}
    #     get inchikey like this:
        d['inchi_key'] = [os.path.basename(a).split('.')[0].split('_')[-1] for a in fid['tree_file_lookup_table']['filename'][idx]]

        tree_df = pd.DataFrame(d)
        # tree_df.set_index('index',drop=True,inplace=True)

    return scan_df,tree_df

def filter_hits(df,score_cutoff=0.1,mz_cutoff=0.1):
    df = df.drop( df[df['score'] < score_cutoff].index )
    mass = df['mass']
    mz = df['precursor mz']
    adduct = 1.007276
    if len(mz)>0:
        df = df.drop( df[(abs( mass + adduct - mz ) > mz_cutoff) & (df['polarity'] == 1)].index )
        df = df.drop( df[(abs( mass - adduct - mz ) > mz_cutoff) & (df['polarity'] == 0)].index )
    return df

def join_pactolus_tables(my_sheets,score_cutoff=0.001,mz_cutoff=0.02):#,use_field='precursor intensity'):
    output_df = pd.DataFrame()
    #rewrite to have a list in each cell: do the max command last
    #pull group info from metob
    for sheet in my_sheets:
        df = pd.read_excel(sheet)
        print(df.shape)
        df = df.drop( df[df['score'] < score_cutoff].index )

        mass = df['mass'][df['polarity'] == 1]
        mz = df['precursor mz'][df['polarity'] == 1]
        adduct = 1.007276
        if len(mz)>0:
            df = df.drop( df[abs( mass + adduct - mz ) > mz_cutoff].index )

        mass = df['mass'][df['polarity'] == 0]
        mz = df['precursor mz'][df['polarity'] == 0]
        adduct = 1.007276
        if len(mz)>0:
            df = df.drop( df[abs( mass - adduct - mz ) > mz_cutoff].index )
        print(df.shape)
        for index, row in input_df.iterrows():
            if type(row['name']) != float:
                try:
                    output_df.index.tolist().index(row['name']) #see if the compound has been added before
                    if pd.notnull(output_df.loc[row['name'],os.path.basename(sheet)]):
                        if row['score'] > output_df.loc[row['name'],os.path.basename(sheet)]:
                            output_df.loc[row['name'],os.path.basename(sheet)] = row['score']
                except:
                    output_df.loc[row['name'],os.path.basename(sheet)] = row['score']
    return output_df
def is_pactolus_result_file(output_file):
    my_keys = output_file.keys()
    counter = 0
    print('checking file',output_file)
    for k in my_keys:
        if 'match_matrix' not in k:
            counter = counter +1
    try:
        score_matrix = output_file['score_matrix'][:]
        num = score_matrix.shape[0]
    except:
        num = 0
    if counter > 3:# 7:
        return True
    else:
        return False

def broadcast_hits_to_expand_dataframe(df):
    df.loc[:,'score'] = df.score.apply(np.atleast_1d)
    df.loc[:,'inchi_key'] = df.inchi_key.apply(np.atleast_1d)
    df.loc[:,'mass'] = df.mass.apply(np.atleast_1d)
    all_scores = np.hstack(df.score)
    all_inchi_key = np.hstack(df.inchi_key)
    all_mass = np.hstack(df.mass)
    all_polarity = np.hstack([[n]*len(l) for n, l in df[['polarity', 'score']].values])
    all_precursor_intensity = np.hstack([[n]*len(l) for n, l in df[['precursor intensity', 'score']].values])
    all_precursor_mz = np.hstack([[n]*len(l) for n, l in df[['precursor mz', 'score']].values])
    all_retention_time = np.hstack([[n]*len(l) for n, l in df[['retention time', 'score']].values])
    df2 = pd.DataFrame({'polarity':all_polarity,'precursor intensity':all_precursor_intensity,'precursor mz':all_precursor_mz,'retention time':all_retention_time,'score':all_scores,'inchi_key':all_inchi_key, 'mass':all_mass})
    return df2
    
def make_output_tables(target_dir,score_cutoff = 0.0,to_csv=True):
#    , overwrite=True,
#                       score_cutoff=0.1,rt_min=0,rt_max=20,intensity_min = 1e4,to_excel=True):
    """
    """
#    df_lookup = pd.DataFrame({'inchi':neutral_inchi,'metatlas name':metatlas_name,'mass':neutral_mass})

    files = glob.glob(os.path.join(target_dir,'*.h5'))
    files = [f for f in files if os.path.basename(f).startswith('pactolus_results')]
    scan_dfs = []
    compound_dfs = []
    all_dfs = []
    # my_file = '/project/projectdirs/openmsi/projects/ben_run_pactolus/rccoates/pactolus_results_20151215_RCC_C18_ACN_Phz_POS_MSMS_WCS417_PCARhizo_S2.h5'
    for my_file in files:
        do_process = False #stupid thing to only operate on results of valid files without having to check twice
        outfile = os.path.join(target_dir,'%s.csv'%os.path.basename(my_file).split('.')[0])
        if (not os.path.isfile(outfile)) or (overwrite):
            with h5py.File(my_file) as output_file:
                if is_pactolus_result_file(output_file):
                    score_matrix = output_file['score_matrix'][:]
                    #if 'inchi' in output_file['compound_metadata'].keys():
                    #    inchi = output_file['compound_metadata']['inchi'][:]
                    #else:
                    #    print score_matrix.shape
                    inchi_key = np.asarray([os.path.basename(a[0]).split('.')[0].split('_')[-1] for a in output_file['tree_file_lookup_table']]) #np.asarray(range(score_matrix.shape[1]))
                    mass = np.asarray([a[1] for a in output_file['tree_file_lookup_table']]) #np.asarray(range(score_matrix.shape[1]))
                        
#                     idx = np.argwhere(score_matrix > score_cutoff)
#                     d = {'retention time': [output_file['scan_metadata']['peak_rt'][i] for i in idx[:,0]],
#                         'precursor intensity': [output_file['scan_metadata']['peak_intensity'][i] for i in idx[:,0]],
#                         'precursor mz': [output_file['scan_metadata']['peak_mz'][i] for i in idx[:,0]],
#                         'polarity': [output_file['scan_metadata']['polarity'][i] for i in idx[:,0]]}
                    d = {'retention time': output_file['scan_metadata']['peak_rt'][:],
                            'precursor intensity': output_file['scan_metadata']['peak_intensity'][:],
                            'precursor mz': output_file['scan_metadata']['peak_mz'][:],
                            'polarity': output_file['scan_metadata']['polarity'][:]}
                    d['score'] = []            
                    d['inchi_key'] = []
                    d['mass'] = []
                    for i in range(score_matrix.shape[0]):
                        idx = np.argwhere(score_matrix[i,:] > score_cutoff).flatten()
                        d['score'].append(score_matrix[i,idx])
                        d['inchi_key'].append(inchi_key[idx])
                        d['mass'].append(mass[idx])
                    
                    df = pd.DataFrame(d)
                    df = broadcast_hits_to_expand_dataframe(df)
                    #df = pd.merge(df,df_lookup,on='inchi_key',how='outer')#ignore_index=True,axis=0)
                    do_process = True
            if do_process:
                print(os.path.basename(outfile))
                if df.shape[0]>0:
                    all_dfs.append(df)
                    if to_csv:
                        df.to_csv(outfile,index=False)
    if all_dfs:
        return all_dfs
    
    
def get_neutral_inchi_and_name(use_pickle=True):
    import pickle
    if use_pickle:
        with open('metatlas_name.pickle', 'rb') as handle:
            metatlas_name = pickle.load(handle)
        with open('neutral_inchi.pickle', 'rb') as handle:
            neutral_inchi = pickle.load(handle)
        with open('neutral_mass.pickle', 'rb') as handle:
            neutral_mass = pickle.load(handle)
    else:
        c = metob.retrieve('Compound',inchi='InChI=%',username='*')
        neutral_inchi = []
        metatlas_name = []
        neutral_mass = []
        for cc in c:
            myMol = Chem.MolFromInchi(cc.inchi.encode('utf-8'))
            myMol, neutralised = NeutraliseCharges(myMol)
            neutral_mass.append(Chem.Descriptors.ExactMolWt(myMol))
            inchi = Chem.MolToInchi(myMol)
            neutral_inchi.append( inchi )
            metatlas_name.append(cc.name)

        with open('metatlas_name.pickle', 'wb') as handle:
            pickle.dump(metatlas_name,handle)
        with open('neutral_inchi.pickle', 'wb') as handle:
            pickle.dump(neutral_inchi,handle)
        with open('neutral_inchi_key.pickle', 'wb') as handle:
            pickle.dump(neutral_inchi,handle)

        with open('neutral_mass.pickle', 'wb') as handle:
            pickle.dump(neutral_mass,handle)
    return metatlas_name,neutral_inchi, neutral_mass

""" contribution from Hans de Winter """
def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

_reactions=None
def NeutraliseCharges(mol, reactions=None):
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
#     mol = Chem.MolFromSmiles(smiles)
    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            rms_smiles = Chem.MolToSmiles(rms[0])
            mol = Chem.MolFromSmiles(rms_smiles)
    if replaced:
        return (mol, True) #Chem.MolToSmiles(mol,True)
    else:
        return (mol, False)

def check_for_failed_jobs(target_dir):
    err_files = [os.path.splitext(os.path.basename(f))[0] for f in glob.glob(os.path.join(target_dir,'*.err'))]
    sbatch_files_full_path = [f for f in glob.glob(os.path.join(target_dir,'*.sbatch'))]
    sbatch_files = [os.path.splitext(os.path.basename(f))[0].split('.')[0] for f in glob.glob(os.path.join(target_dir,'*.sbatch'))]
    hdf5_files = [os.path.splitext(os.path.basename(f))[0].replace('pactolus_results_','') for f in glob.glob(os.path.join(target_dir,'*.h5'))]
    # print hdf5_files
#     print len(err_files),len(sbatch_files),len(hdf5_files)
    failed_jobs = list(set(sbatch_files) - set(hdf5_files))
    if not failed_jobs:
        print("no failed jobs exist")
    else:
        print("failed jobs:")
        for f in failed_jobs:
            print(f)
        for j in failed_jobs:
            print("sbatch",sbatch_files_full_path[index_containing_substring(sbatch_files_full_path,j)])

def index_containing_substring(the_list, substring):
    for i, s in enumerate(the_list):
        if substring in s:
              return i

def check_job_status(do_print=True,computer = 'edison'):
    my_session = Session()
    import getpass
    usr = getpass.getuser()
    pwd = getpass.getpass("enter password for user %s: " % usr)
    r = my_session.post("https://newt.nersc.gov/newt/auth", {"username": usr, "password": pwd})
    r = my_session.get("https://newt.nersc.gov/newt/queue/%s/?user=%s"%(computer,usr))
    my_jobs = r.json()
    # print my_jobs
    if do_print:
        print("You have",len(my_jobs),"jobs running or in the queue to run")
        for i,j in enumerate(my_jobs):
            print(i,'\t',j['status'], j['name'],j['memory'],j['nodes'], j['procs'], j['timeuse'])
    return my_jobs

def create_pactolus_msms_data_container(myfiles,target_directory,min_intensity,min_rt = 1,max_rt = 22,make_container=True):
    # peak_arrayindex: This is a 2D array with the shape (num_spectra, 3). 
    # The dataset contains an index that tells us: 
    # i) the x location of each spectrum [:,0], 
    # ii) the y location of each spectrum [:,1], and 
    # iii) and the index where the spectrum starts in the peak_mz and peak_value array. 
    # In item 1/2 I first fill the array with [0,i,0] values to define unique x/y locations 
    # for each spectrum and in the second line I then create the last column with start index
    # of the spectra which is just the cumulative-sum of the length of the spectra.
    # when you create the start stop locations you will need to:
    # prepend [0] to the cummulative sums (the first spectrum starts at 0 not its length).
    # remove the last entry to make sure the array has the correct length
    # That is why I did the following:
    # np.cumsum([0] + [ ri['m/z array'].shape[0] for ri in good_list ])[:-1]
    if not os.path.exists(target_directory):
        try:
            os.makedirs(target_directory)
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    for myfile in myfiles:
        finfo = h5q.get_info(myfile)
        with tables.open_file(myfile) as fid:    
            num_pos_data = finfo['ms1_pos']['nrows'] + finfo['ms2_pos']['nrows']
            num_neg_data = finfo['ms1_neg']['nrows'] + finfo['ms2_neg']['nrows']
            do_polarity = []
            if num_pos_data > 0:
                do_polarity.append(1)
            if num_neg_data > 0:
                do_polarity.append(0)
            scan_polarity = []
            for my_polarity in do_polarity:
                container_file = os.path.join(target_directory,'container_file_polarity_%d.h5'%(my_polarity))
                if not os.path.isfile(container_file):
                    make_container=True
                if make_container:
                    data = h5q.get_data(fid,ms_level=2,polarity=my_polarity,min_rt = min_rt,max_rt=max_rt,min_precursor_intensity=min_intensity)#TODO: filter by intensity,)
                    prt,pmz,pintensity = mgd.get_unique_scan_data(data)
                    for i in range(len(pintensity)):
                        scan_polarity.append(my_polarity)
                    msms_data = mgd.organize_msms_scan_data(data,prt,pmz,pintensity)
                    fpl = {}
                    # peak_mz : This is a 1D arrays with m/z values for all the spectra stored as spectrum_1, spectrum_2 etc.
                    fpl['peak_mz'] = np.concatenate(tuple( s[:,0] for s in msms_data['spectra']), axis = -1)
                    # peak_value: This is a 1D arrays with the intensity values corresponding to the m/z values stored in peak_mz
                    fpl['peak_value'] = np.concatenate(tuple( s[:,1] for s in msms_data['spectra']), axis = -1)
                    fpl['precursor_mz'] = np.asarray(msms_data['precursor_mz'])
                    fpl['peak_arrayindex'] = np.asarray([[0, i, 0] for i,rt in enumerate(prt)]) 
                    fpl['peak_arrayindex'][:,2] = np.cumsum([0] + [ s[:,0].shape[0] for s in msms_data['spectra'] ])[:-1]
                    with h5py.File(container_file,'a') as output_file:
                        group_name = os.path.basename(myfile)
                        if group_name in output_file.keys():
                            output_file.__delitem__(group_name)
                #         if group_name not in output_file.keys():
                        output_group = output_file.create_group(group_name)
                #         else:
                #             output_group = output_file[group_name]
                        for key, value in fpl.iteritems():
                            output_group[key] = value
                        experiment_group = output_group.create_group('experiment_metadata')
                        experiment_group['filename'] = group_name
                        scan_group = output_group.create_group('scan_metadata')
                        scan_group['peak_mz'] = np.asarray(msms_data['precursor_mz'])
                        scan_group['peak_rt'] = np.asarray(msms_data['precursor_rt'])
                        scan_group['peak_intensity'] = np.asarray(msms_data['precursor_intensity'])
                        scan_group['polarity'] = np.asarray(scan_polarity) # 1 for pos and 0 for neg
                write_pactolus_job_file(myfile,container_file,my_polarity)
    return container_file




def write_pactolus_job_file(myfile,
                            container_file,
                            my_polarity,
                            new_tree_file = '/project/projectdirs/metatlas/projects/clean_pactolus_trees/tree_lookup.npy',
                            base_script_name = '/project/projectdirs/openmsi/projects/ben_run_pactolus/do_not_modify_template_pactolus_script.sh'):
    #regexp the fpl_data path to create lots of jobs:
    # /project/projectdirs/openmsi/projects/ben_run_pactolus/Pactolus_NERSC_BASTet_C18_POS_Archetypes.h5:/20150510_C18_POS_MSMS_HA13-1.h5
    # regexp the outfile test_pactolus_72_2_realtime.h5
    # regexp the log files
    #SBATCH --output=job_pactolus_realtime1_out.txt
    #SBATCH --error=job_pactolus_realtime1_err.txt
    read_pat = '/project/projectdirs/openmsi/projects/ben_run_pactolus/Pactolus_NERSC_BASTet_C18_POS_Archetypes.h5:/20150510_C18_POS_MSMS_HA13-1.h5'
    save_pat = 'test_pactolus_72_2_realtime.h5'
    out_pat = 'job_pactolus_realtime1_out.txt'
    err_pat = 'job_pactolus_realtime1_err.txt'
    
    old_tree_file = '/project/projectdirs/openmsi/projects/ben_trees/metacyc_max_depth_5.npy'
    
    
    pos_neutralizations = '[-1.00727646677,-2.0151015067699998,0.00054857990946]'
    neg_neutralizations = '[1.00727646677,2.0151015067699998,-0.00054857990946]'
    job_pat = 'job_pactolus_'
    with open(base_script_name,'r') as fid:
        base_script_text = fid.read()
    # print base_script_text



    group_name = os.path.basename(myfile)
    no_extension = group_name.split('.')[0]
    
    ##### CHANGE THIS LINE ######
    new_read_pat = '"%s:%s"'%(container_file,group_name)
    #############################
    
    new_save_pat = '"%s"'%os.path.join(os.path.dirname(container_file),'pactolus_results_' + group_name)
    new_out_pat = '"%s"'%os.path.join(os.path.dirname(container_file),no_extension + '.out')
    new_err_pat = '"%s"'%os.path.join(os.path.dirname(container_file),no_extension + '.err')
    new_job_pat = '"%s"'%no_extension
    
    replace_text = [(read_pat,new_read_pat),
                    (save_pat,new_save_pat),
                    (out_pat,new_out_pat),
                    (err_pat,new_err_pat),
                    (job_pat,new_job_pat),
                    (old_tree_file,new_tree_file)]
    
    temp_text = base_script_text
    for rt in replace_text:
        temp_text = temp_text.replace(rt[0],rt[1])
        
#     temp_text = temp_text.replace('#SBATCH --time=00:15:00','#SBATCH --time=00:45:00')
    ##### CHANGE THIS LINE ######
    if my_polarity == 0:
        temp_text = temp_text.replace(pos_neutralizations,neg_neutralizations)
    #############################    
    
    #store the job name in a seperate script file so it can be submited to queue
    #each job will be called <no_extension>.sbatch
    #jobfile will be a list of squeue <no_extension.sbatch\n
    
    ##### CHANGE THIS LINE ######
    new_job_name = '%s/%s_polarity_%d.sbatch'%(os.path.dirname(container_file),os.path.basename(myfile),my_polarity)
    #############################    
    
    with open(new_job_name,'w') as fid:
        fid.write('%s'%temp_text)
    

def submit_all_jobs(target_directory,computer='edison',usr=None,pwd=None):
    import glob
    from requests import Session
    import getpass
    if not usr:
        usr = getpass.getuser()
    if not pwd:
        pwd = getpass.getpass("enter password for user %s: " % usr)
    s = Session()
    r = s.post("https://newt.nersc.gov/newt/auth", {"username": usr, "password": pwd})
    all_files = glob.glob(os.path.join(target_directory,'*.sbatch'))
    for a in all_files:
        r = s.post("https://newt.nersc.gov/newt/queue/%s/"%computer, {"jobfile": a})
    return s
#############################    
