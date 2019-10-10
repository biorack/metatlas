#!/usr/bin/env python

import pactolus.scoring as pactolus
from pactolus.mzml_loader import mzml_to_df

import numpy as np
import glob as glob
import os
import pandas as pd
import h5py
from datetime import datetime
import multiprocessing as mp
import sys
import argparse



def t_str():
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print('%s: importing'%t_str())

def create_msms_dataframe(df):
    """
    create a dataframe organized into spectra from a raw dataframe of points
    """
    #removed polarity and hdf5_file
    if 'precursor_MZ' in df.columns: #shrapnel from early days of not using consistent naming
        df=df.rename(columns = {'precursor_MZ':'precursor_mz'})
    grouped = df.groupby(['precursor_mz','rt','precursor_intensity','collision_energy']).aggregate(lambda x: tuple(x))
    grouped.mz = grouped.mz.apply(list)
    grouped.i = grouped.i.apply(list)
    grouped = grouped.reset_index()
    grouped['spectrum'] = map(lambda x,y:(x,y),grouped['mz'],grouped['i'])
    grouped['spectrum'] = grouped['spectrum'].apply(lambda x: zip(x[0],x[1]))
    grouped.drop(['mz','i'], axis=1, inplace=True)
    return grouped

def do_score(input):
    tree_filename = input[0]
    max_depth = input[1]
    indices = input[2]
    spectra = input[3]
    ms2_neutralizations = input[4]
    params = input[5]
    tree,inchikey = get_tree(tree_filename,max_depth=max_depth)
    output = []
    for i,idx in enumerate(indices):
        try:
            score,match_matrix = pactolus.calculate_MIDAS_score(spectra[i],
                                                                tree,
                                                                params['ms2_tolerance'],
                                                                neutralizations=ms2_neutralizations,
                                                                want_match_matrix=True)
            if score > 0:
                output.append((idx,score,inchikey))#,match_matrix))
        except Exception as e:
            print('spectrum ',i,'could not be scored for',inchikey,'spectrum=',spectra[i],'error=',e)
    return output

def assemble_results(r,spectra_trees):
    list_concat  =[]
    for rr in r:
        for jj in rr:
            list_concat.append((jj[0],jj[1],jj[2]))#,jj[3]))
    return pd.DataFrame(list_concat,columns=['index','score','inchikey'])#,'match_matrix'])

def get_tree(tree_file,max_depth=None):
    """
    Given a path to an hdf5 tree file, 
    return the tree at max_depth and the inchikey
    """
    with h5py.File(tree_file,'r') as fid:
        inchikey = fid.keys()[0]
        if max_depth:
            data_key = 'FragTree_at_max_depth=%d'%max_depth
        else:
            data_key = fid[inchikey].keys()[0]
        tree = fid[inchikey][data_key][:]
    return tree,inchikey

def filtered_spectrum(spectrum,precursor_mz,tolerance):
    """
    remove items in spectrum greather than precursor mz + tolerance
    """
    return [t for t in spectrum if t[0] < (precursor_mz - tolerance)]


def setup_inputs(unique_tree_files,spectra_trees,params):
    """
    setup inputs for running pactolus in parallel
    """
    mp_setup = []
    for tree_filename,indices in unique_tree_files.items():
        temp = []
        temp.append(tree_filename)
        temp.append(spectra_trees.loc[indices[0],'max_depth'])
        temp.append(indices)
        temp_spectra = []
        for idx in indices:
            mz_arr = np.asarray(spectra_trees.loc[idx,'spectrum'])
            temp_spectra.append(mz_arr)
        temp.append(temp_spectra)
        temp.append(spectra_trees.loc[idx,'ms2_neutralizations'])
        temp.append(params)
        mp_setup.append(temp)
    return mp_setup


def make_output_filename(file,outdir):
    bname = os.path.basename(file)
    bname,extention = os.path.splitext(bname)
    return os.path.join(outdir,bname)+'.pactolus.gz'

def get_spectra_from_file(file,outdir):
    if file.split('.')[-1].lower() == 'mzml':
        raw_data = mzml_to_df(file) #returns a dict of dataframes from an mzml file
    elif (file.split('.')[-1].lower() == 'h5') | (file.split('.')[-1].lower() == 'hdf5') | (file.split('.')[-1].lower() == 'hdf'):
        raw_data = mgd.df_container_from_metatlas_file(file) #This is used when input is hdf5 file
    spectra=None
    if isinstance(raw_data['ms2_pos'],pd.DataFrame) & isinstance(raw_data['ms2_neg'],pd.DataFrame): #it has both pos and neg spectra
        spectra = pd.concate([create_msms_dataframe(raw_data['ms2_pos']),create_msms_dataframe(raw_data['ms2_neg'])])
    elif isinstance(raw_data['ms2_pos'],pd.DataFrame):
        spectra = create_msms_dataframe(raw_data['ms2_pos'])
    elif isinstance(raw_data['ms2_neg'],pd.DataFrame):
        spectra = create_msms_dataframe(raw_data['ms2_neg'])
    else:
        print('File has no MSMS data.')#, file=sys.stderr)
        open(make_output_filename(file,outdir), 'a').close() #make empty file

    return spectra

def associate_spectra_with_trees(spectra,trees,args):
    #filter out values < precursor + tolerance
    keep_spectra = []
    for i,row in spectra.iterrows():
        spectra.set_value(i,'spectrum',filtered_spectrum(row.spectrum,row.precursor_mz,args['ms1_tolerance']))

    for i,row in spectra.iterrows():
        if len(row.spectrum) > 1: #spectrum has one or more ions
            keep_spectra.append(True)
        else:
            keep_spectra.append(False)
    spectra = spectra[keep_spectra] #this will drop any spectra that don't have values.

    #Get hits to trees neutralizing by adduct
    
    tree_match = []
    for i,row in spectra.iterrows():
        if 'positive' in row.polarity:
            ms1_neutralizations = args['ms1_pos_neutralizations']
        else:
            ms1_neutralizations = args['ms1_neg_neutralizations']
        for adduct in ms1_neutralizations:
            hits = np.argwhere(abs(trees['ms1_mass'] - (row.precursor_mz-adduct)) < args['ms1_tolerance']).flatten()
            if len(hits)>0:
                for h in hits:
                    tree_match.append((i,trees['filename'][h],trees['max_depth'][h],adduct,args['ms1_tolerance']))
    spectra_trees = pd.merge(spectra,
             pd.DataFrame(tree_match,columns=['hit_index','tree_filename','max_depth','ms1_neutralization','ms1_tolerance']),how='outer',
             left_index=True,
             right_on='hit_index').drop('hit_index',1)

    # spectra_trees['score'] = 0
    # spectra_trees['match_matrix'] = ''
    # spectra_trees['match_matrix'] = spectra_trees['match_matrix'].astype(object)
    spectra_trees['ms2_tolerance'] = args['ms2_tolerance']
    spectra_trees['ms2_neutralizations'] = ''
    spectra_trees['ms2_neutralizations'] = spectra_trees['ms2_neutralizations'].astype(object)
    spectra_trees['ms2_neutralizations'] = spectra_trees.apply(lambda x: args['ms2_pos_neutralizations'] if 'positive' in x['polarity'] else args['ms2_neg_neutralizations'],axis=1)
    spectra_trees.reset_index(inplace=True,drop=True)
    return spectra_trees

def score_and_save(file,mp_data,spectra_trees,args):
    p = mp.Pool(args['num_cores'])
    r = p.map(do_score,mp_data)
    p.close()
    p.terminate()
#    r = []
#    for d in mp_data:
#        r.append(do_score(d))
#        print('done scoring')

    temp_hits = assemble_results(r,spectra_trees)
    temp = pd.merge(spectra_trees,temp_hits,how='outer',left_index=True,right_on='index').drop('index',1)
    
    print('output merged')
    # temp = temp[['precursor_MZ','rt','precursor_intensity','collision_energy','spectrum','tree_filename','ms1_neutralization','score']]
    temp = temp[temp.score>0]
    temp.to_csv(make_output_filename(file,args['outdir']), index=False, compression='gzip')

def main():
    # print command line arguments
    print('%s: starting'%t_str())
    parser = argparse.ArgumentParser(description='a command line tool for searching mzml files with pactolus.  saves a gzipped csv.')
    parser.add_argument('-i','--infile', help='mzml file to search', required=True)
    parser.add_argument('-o','--outdir', help='name of output dir to store results', required=False, default=".")
    parser.add_argument('-m2t','--ms2_tolerance', help='tolerance in Daltons for ms2', type=float,default=0.01)
    parser.add_argument('-m1t','--ms1_tolerance', help='tolerance in Daltons for ms1', type=float,default=0.01)
    parser.add_argument('-m1pn','--ms1_pos_neutralizations', help='adducts to neutralize for in ms1: 1.007276,18.033823,22.989218', type=float,nargs='+',default=[1.007276,18.033823,22.989218])
    parser.add_argument('-m2pn','--ms2_pos_neutralizations', help='ionization states to neutralize for in ms2: -1.00727646677,-2.0151015067699998,0.00054857990946', type=float,nargs='+',default=[-1.00727646677,-2.0151015067699998,0.00054857990946])
    parser.add_argument('-m1nn','--ms1_neg_neutralizations', help='adducts to neutralize for in ms1: -1.007276, 59.013851', type=float,nargs='+',default=[-1.007276, 59.013851])
    parser.add_argument('-m2nn','--ms2_neg_neutralizations', help='ionization states to neutralize for in ms2: 1.00727646677,2.0151015067699998,-0.00054857990946', type=float,nargs='+',default=[1.00727646677,2.0151015067699998,-0.00054857990946])
    parser.add_argument('-t','--tree_file', help='tree file: /project/projectdirs/metatlas/projects/clean_pactolus_trees/tree_lookup.npy', default='/project/projectdirs/metatlas/projects/clean_pactolus_trees/tree_lookup.npy')
    parser.add_argument('-n','--num_cores', help='number of cores to use for multiprocessing', type=int,default=32)
    parser.add_argument('-overwrite','--overwrite', help='Overwrite pre-existing file(s): True/False', type=bool,default=False)

    args = vars(parser.parse_args())

    outdir = os.path.abspath(args['outdir'])
    
    try: 
        os.makedirs(outdir)
    except OSError:
        if not os.path.isdir(outdir):
            raise

    trees = np.load(args['tree_file'])

    
    if os.path.isdir(args['infile']): 
        # read mzml 
        files = glob.glob(os.path.join(args['infile'],'*.mzML'))
        print(len(files),'files with mzML')
        files.extend(glob.glob(os.path.join(args['infile'],'*.mzml')))
        print(len(files),'files with mzml')
        if args['overwrite'] == False:
            files = [f for f in files if not os.path.isfile(make_output_filename(f,outdir))]

        print(len(files),'files to process')
    else:
        files = [args['infile']]

    for file in files:
        print('%s: getting spectra for %s'%(t_str(),file))
        try:
            spectra = get_spectra_from_file(file,args['outdir'])
        except:
            print('ERROR: UNKNOWN ERROR WITH FILE %s'%file)
            open(make_output_filename(file,outdir), 'a').close() #make empty file
            print('EMPTY MADE: %s'%file)
            sys.exit('ERROR EXITING: UNKNOWN ERROR WITH FILE %s'%file)

        if spectra is not None:
            print('%s: linking spectra with trees for %s'%(t_str(),file))
            spectra_trees = associate_spectra_with_trees(spectra,trees,args)
            unique_tree_files = spectra_trees.groupby('tree_filename').indices
            print('%s: setting up input for %s'%(t_str(),file))
            mp_data = setup_inputs(unique_tree_files,spectra_trees,args)
            print('%s: data is setup; doing multiprocessing with %d cores'%(t_str(),args['num_cores']))
            score_and_save(file,mp_data,spectra_trees,args)
            print('%s: finished'%t_str())


if __name__ == "__main__":
    main()
