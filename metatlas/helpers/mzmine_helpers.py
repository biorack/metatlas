from __future__ import print_function
import numpy as np
import sys
import json
import pandas as pd
import re
import glob as glob
import os
from collections import defaultdict
from xml.etree import cElementTree as ET
import multiprocessing as mp

import matplotlib.pyplot as plt

import metatlas.metatlas_objects as metob
from metatlas.helpers import metatlas_get_data_helper_fun as ma_data
from metatlas.helpers import dill2plots as dp

try:
    basestring
except NameError:  # python3
    basestring = str

# setting this very high easily causes out of memory
NUM_THREADS = 8

BATCH_FILE_PATH = '/project/projectdirs/metatlas/projects/mzmine_parameters/batch_files/'
PYTHON_BINARY = '/global/common/software/m2650/python-cori/bin/python'
BINARY_PATH = '/project/projectdirs/metatlas/projects/mzmine_parameters/MZmine'

#we don't need to request haswell on genepool partition of Cori
#remove this line
#SBATCH -C haswell


SLURM_HEADER = """#!/bin/bash
#SBATCH -N 1 -n 64
#SBATCH --error="slurm.err"
#SBATCH --output="slurm.out"
#SBATCH --qos=genepool
#SBATCH -A pkscell
#SBATCH -t 24:00:00
#SBATCH -L project

export MPLBACKEND="agg"
export HDF5_USE_FILE_LOCKING=FALSE

module load java

"""

def see_if_plots_work():
    print('I am here')
    fig = plt.figure()
    print('I made a figure')
    fig.savefig('fig.pdf')
    print('I saved a figure')


def make_figures_from_filtered_data(params,all_files,my_atlas):

    pool = mp.Pool(processes=min(NUM_THREADS, len(all_files)))
    print("# Made Pool")
    print('# Getting data 2')
    metatlas_dataset = pool.map(ma_data.get_data_for_atlas_df_and_file, all_files)
    print('# Getting data 3')
    pool.close()
    pool.join()
    # pool.terminate()
    print('# Acquired Data for plotting')
    output_dir = os.path.join(params['basedir'],'%s_%s'%(params['basename'],params['polarity']))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    atlas_identifications = dp.export_atlas_to_spreadsheet(my_atlas,os.path.join(output_dir,'atlas_export.csv'))
    # atlas_identifications = dp.export_atlas_to_spreadsheet(myAtlas,'%s/sheets/%s.csv'%(plot_location_label,myAtlas.name))
    peak_height = dp.make_output_dataframe(input_fname = '',input_dataset = metatlas_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='peak_height' , output_loc=os.path.join(output_dir,'sheets'))
    peak_area = dp.make_output_dataframe(input_fname = '',input_dataset = metatlas_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='peak_area' , output_loc=os.path.join(output_dir,'sheets'))
    mz_peak = dp.make_output_dataframe(input_fname = '',input_dataset = metatlas_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='mz_peak' , output_loc=os.path.join(output_dir,'sheets'))
    rt_peak = dp.make_output_dataframe(input_fname = '', input_dataset = metatlas_dataset,include_lcmsruns = [],exclude_lcmsruns = [],fieldname='rt_peak' , output_loc=os.path.join(output_dir,'sheets'))
    mz_centroid = dp.make_output_dataframe(input_fname = '',input_dataset = metatlas_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='mz_centroid' , output_loc=os.path.join(output_dir,'sheets'))
    rt_centroid = dp.make_output_dataframe(input_fname = '',input_dataset = metatlas_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='rt_centroid' , output_loc=os.path.join(output_dir,'sheets'))
    print('$$$$ Making ID FIGURES')
    dp.make_identification_figure_v2(input_dataset = metatlas_dataset, input_fname = '', include_lcmsruns = [],exclude_lcmsruns = params['blank_str'].append('QC'), output_loc=os.path.join(output_dir,'identification'))
    print('$$$$ Done Making FIGURES')

def clean_and_filter_mzmine_output(json_filename=None,n_peaks=1000):
    """
    This is a command line function that runs after the first feature finding is complete.
    After mzmine runs, make an atlas, get data, and filter
    The principle is to cast the mzmine peaks into three different categories:
        blank
        bad peak shape
        good no msms
        good with msms
    First, Remove in blank:
        This can be done from mzmine output
        exclude peaks where neither peak height nor peak area is appreciably greater in sample vs blank
    Second, Remove Bad Peak shape:
        seems like rt_peak being near rt@max is a good thing
        seems like looking at intensity max
        you don't want to exclude peaks like leu-ileu.
        many peaks will have multiple-humps
        A simplified logic is: is max intensity near the rt_peak? if yes, is it decaying by 10x +/- 0.5 minutes. Obviously, peaks that are +/- 0.5 minutes appart will get removed. Thus, find local minimum +/- 0.5 minutes from rt_peak. Set those as rt_min and rt_max.
    Third remove no MS/MS
        exclude if no feature in non-blank samples has msms
        
    Last: only keep n_peaks by max intensity in any sample.

    """


    if not json_filename:
        json_filename = sys.argv[1]
    with open(json_filename) as data_file:    
        params = json.load(data_file)
    file_to_convert = os.path.join(params['basedir'],'intermediate_results','%s_%s.csv'%(params['basename'],params['polarity']))
    print('# Working on %s %s'%(params['basename'],params['polarity']))
    if os.path.isfile(file_to_convert):
        # take the comprehensive csv from mzmine and make a peak-height only version of it
        df,original_mzmine = metatlas_formatted_atlas_from_mzmine_output(file_to_convert,params['polarity'],
                                                                                 make_atlas=False,min_rt=0.55,
                                                                                remove_fragments=False,
                                                                                remove_adducts=False,
                                                                                remove_clusters=False,)
        df.to_csv(file_to_convert.replace('.csv','') + '_formatted.csv',index=True) #save a simplified mzmine-like csv as a backup of all features found.
        
        #Filter features found in the blank
        df = clean_up_mzmine_dataframe(df)
        df_blank_compare = df.transpose().groupby(['b' if any([s in g.lower() for s in params['blank_str']]) else 's' for g in df.columns]).max().transpose()
        if 'b' in df_blank_compare.columns:
            df_features_not_in_blank = df_blank_compare[df_blank_compare['s'] > (params['sample_to_blank'] * df_blank_compare['b'])]
            print('# %d features total'%(df.shape[0]))
            print('# %d features removed by blank\n'%(df_blank_compare.shape[0] - df_features_not_in_blank.shape[0]))
        else:
            print('# No files have "blank" or "injbl" in their names.')
            df_features_not_in_blank = df_blank_compare
        df_features_not_in_blank.reset_index(inplace=True)
        print('#There are now %d features not in blank'%df_features_not_in_blank.shape[0])
        
        
        #Make an Atlas
        cids = []
        for j,row in df_features_not_in_blank.iterrows():
            my_mz_ref = metob.MzReference(mz=row.mz,mz_tolerance=row.mz_tolerance,detected_polarity=params['polarity'],lcms_run=None)
            use_rt = row.rt_min
            if row.rt_min < params['min_rt']:
                use_rt = params['min_rt']
            my_rt_ref = metob.RtReference(rt_peak=row.rt_peak,rt_min=use_rt,rt_max=row.rt_max,lcms_run=None)
            my_id = metob.CompoundIdentification(rt_references=[my_rt_ref],mz_references=[my_mz_ref],name=row.label)
            cids.append(my_id)
        my_atlas = metob.Atlas(name='untargeted atlas',compound_identifications=cids)
        atlas_df = ma_data.make_atlas_df(my_atlas)
        
        #Make Groups
        all_files = [f.replace('Peak height','').replace('filtered','').strip() for f in df.columns if '.mzML' in f]
        metatlas_files = []      
        for f in all_files:
            f = metob.retrieve('Lcmsruns',name=f,username='*')[-1]
            if isinstance(f,type(metob.LcmsRun())):
                metatlas_files.append(f)
            else:
                print('%s NOT FOUND'%f)
                break
        groups = metob.Group(name='untargeted group',items=metatlas_files)
        
        #Get Data
        all_files = []
        for my_file in groups.items:
            all_files.append((my_file,groups,atlas_df,my_atlas))
        # print('trying it without multiprocessing')
        # metatlas_dataset = []
        # for f in all_files:
            # metatlas_dataset.append(ma_data.get_data_for_atlas_df_and_file(f))
        metatlas_dataset = []
        for f in all_files:
            metatlas_dataset.append(ma_data.get_data_for_atlas_df_and_file(f))

        # pool = mp.Pool(processes=min(NUM_THREADS, len(all_files)))
        # metatlas_dataset = pool.map(ma_data.get_data_for_atlas_df_and_file, all_files)
        # pool.close()
        # pool.join()
        # pool.terminate()
        
        # remove peaks that aren't valid.  See pk_checker for details on validity.
        num_features = len(metatlas_dataset[0])
        num_files = len(metatlas_dataset)
        
        valid_peaks = [pk_checker(metatlas_dataset,atlas_df,compound_idx,params) for compound_idx in range(num_features)]

        # add another boolean after valid_peaks of if peak is in top 1000
        valid_peaks = peak_in_top_n(metatlas_dataset,n_peaks=n_peaks,prior_boolean=valid_peaks)

        df_filtered_peaks_atlas = atlas_df.loc[valid_peaks]
        file_to_convert = os.path.join(params['basedir'],'intermediate_results','%s_%s.csv'%(params['basename'],params['polarity']))
        df_filtered_peaks_atlas[['mz','rt_peak','label']].sort_values('rt_peak').to_csv(file_to_convert.replace('.csv','') + '_formatted_peakfiltered.csv',index=False)
        print('# There are now %d filtered peaks'%df_filtered_peaks_atlas.shape[0])
        # make a filtered atlas, get data, and make_plots
        filtered_ids = []
        for i,valid in enumerate(valid_peaks):
            if valid:
                filtered_ids.append(my_atlas.compound_identifications[i])
        my_atlas.compound_identifications = filtered_ids
        atlas_df = ma_data.make_atlas_df(my_atlas)
        print("# Atlas shape: %d and Compound Identification: %d"%(atlas_df.shape[0],len(my_atlas.compound_identifications)))
        #Get Data
        all_files = []
        for my_file in groups.items:
            all_files.append((my_file,groups,atlas_df,my_atlas))
        print("# Made filelist")
        print('# Getting data 1')
        make_figures_from_filtered_data(params,all_files,my_atlas)

        
def peak_height_df(metatlas_dataset,attribute='peak_height',zero_nans=True):
    """
    Turn a list of lists in a metatlas dataset into a 
    peak height dataframe where rows are features
    and columns are samples
    
    Valid attributes are:'mz_centroid','mz_peak',
    'num_ms1_datapoints','peak_area','peak_height',
    'rt_centroid','rt_peak'
    
    infs, nans, and nulls are converted to zero by default.
    """
    d = []
    for m in metatlas_dataset: #iterate over features
        row = []
        for mm in m: #iterate over files
            try:
                row.append(mm['data']['ms1_summary'][attribute])
            except:
                row.append(0)
        d.append(row)
    df = pd.DataFrame(d).T
    if zero_nans:
        df[pd.isnull(df)] = 0
    return df

def peak_in_top_n(metatlas_dataset,n_peaks=1000,prior_boolean=None):
    """
    
    """
    df = peak_height_df(metatlas_dataset)
    if prior_boolean is not None: #make dataframe
        top_peaks = df[prior_boolean].max(axis=1).rank(method='min',ascending=False)<=n_peaks
        df.loc[top_peaks.index,'top_peaks'] = top_peaks
        df['top_peaks'].fillna(False,inplace=True)
    else:
        top_peaks = df.max(axis=1).rank(method='min',ascending=False)<=n_peaks
        df['top_peaks'] = top_peaks
    return df['top_peaks'].tolist()



def metatlas_formatted_atlas_from_mzmine_output(filename,polarity,make_atlas=True,atlas_name=None,
    do_store=False,min_rt=None,max_rt=None,min_mz=None,mz_tolerance=8,
    max_mz=None,remove_adducts=False,remove_fragments=False,remove_clusters=False):
    # 
    '''
    Turn mzmine output into conforming metatlas_atlas input
    
    Input:
    filename: csv file from mzmine output
    polarity: (positive,negative)
    atlas_name: string describing the atlas. useful incase you want to save it later.
    
    Output:
    atlas_df: dataframe of atlas content
    myAtlas: metatlas atlas object (if make_atlas=True)
    mzmine_df: dataframe of all mzmine info (if make_atlas=False)

    '''
    
    mzmine_df = pd.read_csv(filename)
    if min_rt:
        mzmine_df = mzmine_df[mzmine_df['row retention time']>min_rt]
    if max_rt:
        mzmine_df = mzmine_df[mzmine_df['row retention time']<max_rt]
    if min_mz:
        mzmine_df = mzmine_df[mzmine_df['row m/z']>min_mz]
    if max_mz:
        mzmine_df = mzmine_df[mzmine_Df['row m/z']<max_mz]
    if remove_adducts:
        mzmine_df = mzmine_df[~mzmine_df['row identity'].str.contains('Adduct',na=False)]
    if remove_fragments:
        mzmine_df = mzmine_df[~mzmine_df['row identity'].str.contains('Fragment',na=False)]
    if remove_clusters:
        mzmine_df = mzmine_df[~mzmine_df['row identity'].str.contains('Complex',na=False)]


    def clean_adducts(x):
        x = re.sub(r'\d+\.\d{2,}', lambda m: format(float(m.group(0)), '.2f'), x)
        new_x = ';'.join(
            [s.strip() for s in pd.unique(x.split(';'))]
            )
        return x 

    metatlas_atlas = pd.DataFrame()
    metatlas_atlas['label'] = mzmine_df.apply(lambda x: '%.4f@%.2f'%(x['row m/z'],x['row retention time']),axis=1)
    mzmine_df['row identity'].fillna('',inplace=True)
    metatlas_atlas['adduct_assignments'] = mzmine_df['row identity']#.apply(clean_adducts)
    metatlas_atlas['mz'] = mzmine_df.apply(lambda x: x['row m/z'],axis=1)
    metatlas_atlas['mz_tolerance'] = mz_tolerance#metatlas_atlas.mz.apply(lambda x: mz_tolerance*float(x)/1e6)
    metatlas_atlas['rt_peak'] = mzmine_df.apply(lambda x: x['row retention time'],axis=1)
    rt_min_cols = [col for col in mzmine_df.columns if 'Peak RT start' in col]
    metatlas_atlas['rt_min'] = mzmine_df[rt_min_cols].apply(lambda x: x.min(),axis=1)
    rt_min_cols = [col for col in mzmine_df.columns if 'Peak RT end' in col]
    metatlas_atlas['rt_max'] = mzmine_df[rt_min_cols].apply(lambda x: x.max(),axis=1)
    metatlas_atlas['inchi_key'] = None
    metatlas_atlas['detected_polarity'] = polarity
    
    #tuplize the 'Identification method' and 'Name' from adducts and fragments 

    # stick on the peak height columns
    pk_height = [col for col in list(mzmine_df) if 'Peak height' in col]
    metatlas_atlas = pd.concat([metatlas_atlas,mzmine_df[pk_height]],axis=1)
    metatlas_atlas['max_intensity'] = metatlas_atlas[pk_height].max(axis=1)
    metatlas_atlas.reset_index(inplace=True)
    metatlas_atlas.drop('index',axis=1,inplace=True)
    if make_atlas:
        if not atlas_name:
            atlas_name = filename
        myAtlas = dp.make_atlas_from_spreadsheet(metatlas_atlas,
                                           atlas_name,
                                           filetype='dataframe',
                                           sheetname='',
                                           polarity = polarity,
                                           store=do_store)
        return metatlas_atlas, myAtlas
    else:
        return metatlas_atlas, mzmine_df


def make_task_and_job(params):#basedir,basename,polarity,files):
    if not os.path.exists(params['basedir']):
        os.mkdir(params['basedir'])

    xml_str = get_batch_file_template()
    d = xml_to_dict(xml_str)
#     for i,k in enumerate(d['batch']['batchstep']):
#         print(i,k['@method'])

    task = metob.MZMineTask()
    task.polarity = params['polarity']
    task.lcmsruns = params['files']
    new_d = replace_files(d,params['files'])


    task.min_peak_duration = params['min_peak_duration']
    task.max_peak_duration = params['max_peak_duration']
    task.rt_tol_perfile = params['rt_tol_perfile']
    task.rt_tol_multifile = params['rt_tol_multifile']
    task.min_peak_height = params['min_peak_height']
    task.noise_floor = params['noise_floor']
    task.mz_tolerance = params['mz_tolerance']
    task.min_sn_ratio = params['min_sn_ratio']
    task.min_rt = params['min_rt']
    task.max_rt = params['max_rt']
    task.output_csv = os.path.join(params['basedir'],'intermediate_results','%s_%s.csv'%(params['basename'],task.polarity))
    task.output_workspace = os.path.join('/dev/null','%s_%s.mzmine'%(params['basename'],task.polarity))
    task.input_xml = os.path.join(params['basedir'],'logs','%s_%s.xml'%(params['basename'],task.polarity))
    task.mzmine_launcher = get_latest_mzmine_binary()
    new_d = configure_crop_filter(new_d,task.polarity,params['files'])
    new_d = configure_mass_detection(new_d,task.noise_floor,task.polarity)
    new_d = configure_chromatogram_builder(new_d,task.min_peak_duration,task.min_peak_height,task.mz_tolerance,task.polarity,task.min_rt,task.max_rt)
    new_d = configure_peak_deconvolution(new_d,task.min_peak_height,task.min_sn_ratio,task.min_peak_duration,task.max_peak_duration)
    new_d = configure_isotope_adduct_fragment_search(new_d,task.mz_tolerance,task.rt_tol_perfile,task.polarity,task.min_peak_height)
    new_d = configure_join_aligner(new_d,task.mz_tolerance,task.rt_tol_multifile)
    new_d = configure_duplicate_filter(new_d,task.mz_tolerance,task.rt_tol_perfile)
    new_d = configure_gap_filling(new_d,task.mz_tolerance)
    new_d = configure_csv_output(new_d,task.output_csv)
    new_d = configure_workspace_output(new_d,task.output_workspace)

    t = dict_to_etree(new_d)
    indent_tree(t)
    xml_batch_str = tree_to_xml(t,filename=task.input_xml)
    job_runner = '%s %s'%(task.mzmine_launcher,task.input_xml)
    return job_runner


def create_job_script(m):
    """
    
    This is the first function that runs when a user initializes a new untargeted workflow

    """

    #setup directories
    if not os.path.isdir(m['basedir']):
        os.mkdir(m['basedir'])
    dirs_to_make = ['job_scripts','logs','intermediate_results','%s_%s'%(m['basename'],m['polarity'])]
    for d in dirs_to_make:
        if not os.path.isdir(os.path.join(m['basedir'],d)):
            os.mkdir(os.path.join(m['basedir'],d))

    job_cmd = make_task_and_job(m)#['basedir'],m['basename'],m['polarity'],m['files'])


    sbatch_file_name = os.path.join(m['basedir'],'job_scripts','%s_%s.sbatch'%(m['basename'],m['polarity']))
    denovo_sbatch_file_name = os.path.join(m['basedir'],'job_scripts','%s_%s_denovo.sbatch'%(m['basename'],m['polarity']))
    err_file_name = os.path.join(m['basedir'],'logs','%s_%s.err'%(m['basename'],m['polarity']))
    out_file_name = os.path.join(m['basedir'],'logs','%s_%s.out'%(m['basename'],m['polarity']))
    job_cmd_filtered = make_targeted_mzmine_job(m['basedir'],m['basename'],m['polarity'],m['files'])
    params_filename = os.path.join(m['basedir'],'logs','%s_%s_params.json'%(m['basename'],m['polarity']))
    new_params_filename = os.path.join(m['basedir'],'logs','%s_%s_params-used.json'%(m['basename'],m['polarity']))

    copy_params_command = "cp '%s' '%s'"%(params_filename,new_params_filename)

    # build up the command to use mzmine-csv and metatlas for heavy peak filtering
    command_line_command = 'from metatlas.helpers import mzmine_helpers as mzm; mzm.clean_and_filter_mzmine_output()'
    command_line_argument = '%s'%params_filename
    python_string = '%s -c "%s" "%s"'%(PYTHON_BINARY,command_line_command,command_line_argument)

    # build up the command to convert mzmine-filtered csv into formatted csv
    # command_line_command = 'from metatlas.helpers import mzmine_helpers as mzm; mzm.command_line_make_formatted_csv_from_filtered()'
    # second_python_string = '%s -c "%s" "%s"'%(PYTHON_BINARY,command_line_command,command_line_argument)

    with open(sbatch_file_name,'w') as fid:
        fid.write('%s\n'%SLURM_HEADER.replace('slurm.err',err_file_name).replace('slurm.out',out_file_name))
        if m['metatlas_path']:
            fid.write('\n\n%s\n\n'%m['metatlas_path'])
        if not m['mzmine_done']:
            fid.write('%s\n'%copy_params_command)
            fid.write('%s\n'%job_cmd)
        fid.write('%s\n'%python_string)

        if not m['mzmine_done']:
            fid.write('%s\n'%job_cmd_filtered)
        # fid.write('%s\n'%second_python_string)

    bad_words = ['qos', '-p','-C','-L','-t','-N']
    bad_time = '#SBATCH -t 24:00:00'
    good_time = '#SBATCH -t 24:00:00\n'

    bad_node = '-N 1 -n 64'
    good_node = '#SBATCH -N 1 -c 32\n'


    with open(sbatch_file_name) as oldfile, open(denovo_sbatch_file_name, 'w') as newfile:
        for line in oldfile:
            if not any(bad_word in line for bad_word in bad_words):
                newfile.write(line)
            if bad_time in line:
                newfile.write(good_time)
            if bad_node in line:
                newfile.write(good_node)
                newfile.write('#SBATCH --mem=300G\n')


    return sbatch_file_name

#####################################################
#####################################################
########         mzmine setup scripts        ########
#####################################################
#####################################################




def remove_duplicate_files(files):
    file_names = []
    unique_files = []
    for f in files:
        if not f.name in file_names:
            unique_files.append(f.mzml_file)
            file_names.append(f.name)
    return unique_files

def get_files(groups,filename_substring,file_filters,is_group=False):
    """
    if is_group is False, gets files from the experiment/folder name and filters with file_filters

    if is_group is True, gets files from the metatlas group name and filters with file_filters

    """

    for i,g in enumerate(groups):
        if is_group == True:
        # get files as a metatlas group
            groups = dp.select_groups_for_analysis(name = g,do_print=False,
                                                   most_recent = True,
                                                   remove_empty = True,
                                                   include_list = [], exclude_list = file_filters)#['QC','Blank'])
            new_files = []
            for each_g in groups:
                for f in each_g.items:
                    new_files.append(f)
        else:
            new_files = metob.retrieve('Lcmsruns',experiment=g,name=filename_substring,username='*')
        if i == 0:
            all_files = new_files
        else:
            all_files.extend(new_files)
        if len(new_files) == 0:
            print('##### %s has ZERO files!'%g)
    if len(file_filters) > 0:
        for i,ff in enumerate(file_filters):
            if i == 0:
                files = [f for f in all_files if not ff in f.name]
            else:
                files = [f for f in files if not ff in f.name]
    else:
        files = all_files
    files = remove_duplicate_files(files)
    return files


    

def clean_up_mzmine_dataframe(df):
    """
    remove a few stray columns and set index to metatlas like attributes.
    
    this leaves a metatlas-like index and peak-height columns as the only thing remaining.
    """
    #df['rt_min'] = df['rt_peak'] - 0.2
    df = df[[c for c in df.columns if not 'max_intensity' in c]]
    df = df[[c for c in df.columns if not 'Unnamed:' in c]]
    index_columns = ['mz','rt_peak','label','mz_tolerance','rt_min','rt_max','inchi_key','detected_polarity','adduct_assignments']
    df.set_index(index_columns,inplace=True)
    df = df[sorted(df.columns)]
    
    return df

def rt_checker(met_data,atlas_df,compound_idx,params):
    """
    simply checks is actual peak is within rt_timespan of stated peak
    """
    try:
        valid = abs(met_data[compound_idx]['data']['ms1_summary']['rt_peak'] - atlas_df.loc[compound_idx,'rt_peak']) < params['rt_timespan']
    except:
        valid = False
    return valid

def min_checker(met_data,atlas_df,compound_idx,params):
    """
    looks forward and backward by rt_timespan and requires that the measured peak height be 
    greater than minima.
    """
    try:
        measured_rt_peak = met_data[compound_idx]['data']['ms1_summary']['rt_peak']
        peak_height = met_data[compound_idx]['data']['ms1_summary']['peak_height']
        if pd.isnull(measured_rt_peak) or peak_height < params['min_intensity']:
            return False
        else:
            eic = met_data[compound_idx]['data']['eic']
            condition_1 = np.asarray(eic['rt']) > measured_rt_peak
            condition_2 = np.asarray(eic['rt']) < (measured_rt_peak + params['rt_timespan'])
            condition_3 = np.asarray(eic['rt']) < measured_rt_peak
            condition_4 = np.asarray(eic['rt']) > (measured_rt_peak - params['rt_timespan'])
            intensity = np.asarray(eic['intensity'])
            
            forward_idx = (condition_1) & (condition_2)
            if any(forward_idx):
                forward_pass = peak_height > (intensity[forward_idx].min() * params['peak_to_valley_ratio'])
            else:
                forward_pass = False
            
            backward_idx = (condition_3) & (condition_4)
            if any(backward_idx):
                backward_pass = peak_height > (intensity[backward_idx].min() * params['peak_to_valley_ratio'])
            else:
                backward_pass = False
            return ((forward_pass) and (backward_pass))
    except:
        return False


def pk_checker(metatlas_dataset,atlas_df,compound_idx,params):
    # across all files for a compound, make sure that
    # the peak height is above threshold,
    # the peak is within rt-span tolerance and
    # there is a minimal value on both sides of the peak
    # there is msms between rt_min and rt_max
    # this condition only needs to occur in one file
    has_msms = []
    for met_data in metatlas_dataset:
        try:
            has_msms.append(len(met_data[compound_idx]['data']['msms']['data']['mz'])>0)
        except:
            has_msms.append(False)
    rt_valid = [rt_checker(met_data,atlas_df,compound_idx,params) for met_data in metatlas_dataset]
    minmax_valid = [min_checker(met_data,atlas_df,compound_idx,params) for met_data in metatlas_dataset]
    is_valid = any((has_msms) and (rt_valid) and (minmax_valid))
    return is_valid



# def create_targeted_job_script(m):
#     """
#     """
#     job_cmd = make_targeted_mzmine_job(m['basedir'],m['basename'],m['polarity'],m['files'])
#     sbatch_file_name = os.path.join(m['basedir'],'%s_%s_filtered.sbatch'%(m['basename'],m['polarity']))
#     err_file_name = os.path.join(m['basedir'],'%s_%s_filtered.err'%(m['basename'],m['polarity']))
#     out_file_name = os.path.join(m['basedir'],'%s_%s_filtered.out'%(m['basename'],m['polarity']))
#     with open(sbatch_file_name,'w') as fid:
#         fid.write('%s\n'%SLURM_HEADER.replace('slurm.err',err_file_name).replace('slurm.out',out_file_name))
#         fid.write('%s\n'%job_cmd)
#     return sbatch_file_name

def make_targeted_mzmine_job(basedir,basename,polarity,files):
    if not os.path.exists(basedir):
        os.mkdir(basedir)

    xml_str = get_targeted_batch_file_template()
    d = xml_to_dict(xml_str)

    task = metob.MZMineTask()
    task.polarity = polarity
    task.lcmsruns = files
    new_d = replace_files(d,files)

    task.min_peak_duration = 0.025
    task.max_peak_duration = 30.0
    task.rt_tol_perfile = 0.015
    task.rt_tol_multifile = 0.15
    task.min_peak_height = 1e6
    task.noise_floor = 3e4
    task.mz_tolerance = 10.0
    task.min_sn_ratio = 2.0
    project_name = '%s_%s'%(basename,task.polarity)
    task.output_csv = os.path.join(basedir,'intermediate_results','%s_%s_filtered.csv'%(basename,task.polarity))
    task.output_workspace = os.path.join(basedir,project_name,'%s_%s.mzmine'%(basename,task.polarity))
    task.input_xml = os.path.join(basedir,'logs','%s_%s_filtered.xml'%(basename,task.polarity))
    
    peak_list_filename = os.path.join(basedir,'intermediate_results','%s_%s_formatted_peakfiltered.csv'%(basename,polarity))
    task.mzmine_launcher = get_latest_mzmine_binary()

    new_d = configure_crop_filter(new_d,task.polarity,files)
    new_d = configure_targeted_peak_detection(new_d,peak_list_filename,intensity_tolerance=1e-4,noise_level=1e4,mz_tolerance=20,rt_tolerance=0.5)
    new_d = configure_workspace_output(new_d,task.output_workspace)

    t = dict_to_etree(new_d)
    indent_tree(t)
    xml_batch_str = tree_to_xml(t,filename=task.input_xml)
    job_runner = '%s %s'%(task.mzmine_launcher,task.input_xml)
    return job_runner

def configure_targeted_peak_detection(new_d,peak_list_filename,intensity_tolerance=1e-4,noise_level=1e4,mz_tolerance=20,rt_tolerance=0.5):
    """
    Name suffix: Suffix to be added to the peak list name.
    
    Peak list file: Path of the csv file containing the list of peaks to be detected. The csv file should have three columns. 
    The first column should contain the expected M/Z, the second column the expected RT and the third the peak name. Each peak should be in a different row.
    
    Field separator: Character(s) used to separate fields in the peak list file.
    
    Ignore first line: Check to ignore the first line of peak list file.
    
    Intensity tolerance: This value sets the maximum allowed deviation from expected shape of a peak in chromatographic direction.
    
    Noise level: The minimum intensity level for a data point to be considered part of a chromatogram. All data points below this intensity level are ignored.
    
    MZ Tolerance: Maximum allowed m/z difference to find the peak
    
    RT tolerance: Maximum allowed retention time difference to find the peak
    """
    # Set the noise floor
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'TargetedPeakDetectionModule' in d['@method']][0]
    
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Peak list file' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%s'%peak_list_filename

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Intensity tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.6f'%(intensity_tolerance)

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Noise level' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.6f'%(noise_level)

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.6f'%(rt_tolerance)

    return new_d

def configure_crop_filter(new_d,polarity,files):
    """
    
    """
    # Set the noise floor
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'CropFilterModule' in d['@method']][0]
    
    # Set the filter string
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Raw data files' in d['@name']][0]
    if any(['FPS' in f for f in files]):
        new_d['batch']['batchstep'][idx]['parameter'][idx2]['name_pattern'] = '*FPS*'
    else:
        new_d['batch']['batchstep'][idx]['parameter'][idx2]['name_pattern'] = '*'


    # Set the polarity
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Scans' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['polarity'] = polarity.upper()
    # new_d['batch']['batchstep'][idx]['parameter'][idx2]['ms_level'] = '1-2'

    return new_d

def configure_mass_detection(new_d,noise_floor,polarity):
    """
    
    """
    # Set the noise floor
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'MassDetectionModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Mass detector' in d['@name']][0]
    idx3 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module']) if 'Centroid' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']['#text'] = '%.2f'%(noise_floor)
    
    # Set the polarity
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Scans' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['polarity'] = polarity.upper()
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ms_level'] = '1'

    
    return new_d

def configure_chromatogram_builder(new_d,min_peak_duration,min_peak_height,mz_tolerance,polarity,min_rt,max_rt):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'ChromatogramBuilderModule' in d['@method']][0]
    new_d['batch']['batchstep'][idx]['parameter']

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Min time span' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(min_peak_duration)

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Min height' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(min_peak_height)

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    
   
    # Set the polarity
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Scans' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['polarity'] = polarity.upper()
    
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['retention_time']['min'] = '%.3f'%min_rt
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['retention_time']['max'] = '%.3f'%max_rt

    return new_d


def configure_peak_deconvolution(new_d,min_peak_height,min_sn_ratio,min_peak_duration,max_peak_duration):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'DeconvolutionModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Algorithm' in d['@name']][0]
    idx3 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module']) if 'Local minimum search' in d['@name']][0]
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Chromatographic threshold' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%0.1
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Search minimum in RT range (min)' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%0.05
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Minimum relative height' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%0.001
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Minimum absolute height' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%min_peak_height
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Min ratio of peak top/edge' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%min_sn_ratio
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Peak duration range (min)' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['min'] = '%.3f'%min_peak_duration
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['max'] = '%.3f'%max_peak_duration
    
    #following deconvolution, many small peaks are created.  Filter them out here
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'PeakFilterModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Height' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['min'] = '%.3f'%(min_peak_height)
    
    return new_d

def configure_isotope_adduct_fragment_search(new_d,mz_tolerance,rt_tol_perfile,polarity,min_peak_height):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'Isotope' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_perfile)

    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'Adduct' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'RT tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_perfile)
    if polarity == 'negative':
        idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Adducts' in d['@name']][0]
        #the default is setup for positive mode adducts.
        #only change them if you are in negative mode
        for i,a in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['adduct']):
            if a['@selected'] == 'true':
                new_d['batch']['batchstep'][idx]['parameter'][idx2]['adduct'][i]['@selected'] = 'false'
            else:
                new_d['batch']['batchstep'][idx]['parameter'][idx2]['adduct'][i]['@selected'] = 'true'

    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'ComplexSearchModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_perfile)
    if polarity == 'negative':
        #the default is setup for positive mode adducts.
        #only change them if you are in negative mode
        idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Ionization method' in d['@name']][0]
        new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '[M-H]-'

    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'FragmentSearchModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_perfile)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Min MS2 peak height' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(min_peak_height)
    return new_d

def configure_join_aligner(new_d,mz_tolerance,rt_tol_multifile):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'JoinAlignerModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_multifile)
    return new_d

def configure_duplicate_filter(new_d,mz_tolerance,rt_tol_perfile):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'DuplicateFilterModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'RT tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_perfile)
    return new_d

def configure_gap_filling(new_d,mz_tolerance):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'SameRangeGapFillerModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    return new_d

def configure_workspace_output(new_d,output_workspace):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'ProjectSaveModule' in d['@method']][0]
    new_d['batch']['batchstep'][idx]['parameter']['#text'] = output_workspace
    return new_d

def configure_csv_output(new_d,output_csv):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'CSVExportModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Filename' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = output_csv
    return new_d

def indent_tree(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent_tree(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def get_targeted_batch_file_template(loc='do_not_change_batch_file_targeted_peak_list.xml'):
    """
    return string of text from the template batch file
    """
    with open(os.path.join(BATCH_FILE_PATH,loc),'r') as fid:
        file_text = fid.read()
    return file_text

def get_batch_file_template(loc='do_not_change_batch_file.xml'):
    """
    return string of text from the template batch file
    """
    with open(os.path.join(BATCH_FILE_PATH,loc),'r') as fid:
        file_text = fid.read()
    return file_text

def get_latest_mzmine_binary(system='Cori',version='most_recent'):
    """
    Returns the path to the mzmine launch script.
    Default is most recent.  Alternatively specify the folder containng version you want

    for example:
        version='MZmine-2.23'
    will use the launch script in that folder

    # To setup the most recent mzmine binary, follow these steps
    cd /project/projectdirs/metatlas/projects/mzmine_parameters/MZmine
    wget $(curl -s https://api.github.com/repos/mzmine/mzmine2/releases/latest | grep 'browser_' | cut -d\" -f4) -O mzmine_latest.zip
    unzip mzmine_latest.zip
    # change directories into latest mzmine download
    # cd MZmine-XXXX
    cp ../MZmine-2.24/startMZmine_NERSC_* .
    cd /project/projectdirs/metatlas/projects/
    chgrp -R metatlas mzmine_parameters
    chmod -R 770 mzmine_parameters 
    """
    mzmine_versions = glob.glob(os.path.join(BINARY_PATH,'*'))
    if version == 'most_recent':
        most_recent = sorted([os.path.basename(m) for m in mzmine_versions if 'MZmine-' in m])[-1]
    else:
        most_recent = 'version'
    launch_script = os.path.join(os.path.join(BINARY_PATH,most_recent),'startMZmine_NERSC_Headless_%s.sh'%system)
    if os.path.isfile(launch_script):
        return launch_script
    else:
        print('See the docstring, the launch script seems to be missing.')

def replace_files(d,file_list):
    """
    Replace files for mzmine task
    
    Inputs:
    d: an xml derived dictionary of batch commands
    file_list: a list of full paths to mzML files
    
    Outputs:
    d: an xml derived dict with new files in it
    """
    for i,step in enumerate(d['batch']['batchstep']):
        if 'RawDataImportModule' in step['@method']:
            d['batch']['batchstep'][i]['parameter']['file'] = file_list
    return d





def tree_to_xml(t,filename=None):
    """

    """
    xml_str = ET.tostring(t)
    if filename:
        with open(filename,'w') as fid:
            fid.write(xml_str)
    return xml_str

def dict_to_etree(d):
    """
    Convert a python dictionary to an xml str
    http://stackoverflow.com/questions/7684333/converting-xml-to-dictionary-using-elementtree

    Example:
    from collections import defaultdict
    from xml.etree import cElementTree as ET

    try:
        basestring
    except NameError:  # python3
        basestring = str

    #d is a python dictionary
    ET.tostring(dict_to_etree(d))

    """
    def _to_etree(d, root):
        if not d:
            pass
        elif isinstance(d, basestring):
            root.text = d
        elif isinstance(d, dict):
            for k,v in d.items():
                assert isinstance(k, basestring)
                if k.startswith('#'):
                    assert k == '#text' and isinstance(v, basestring)
                    root.text = v
                elif k.startswith('@'):
                    assert isinstance(v, basestring)
                    root.set(k[1:], v)
                elif isinstance(v, list):
                    for e in v:
                        _to_etree(e, ET.SubElement(root, k))
                else:
                    _to_etree(v, ET.SubElement(root, k))
        else: assert d == 'invalid type', (type(d), d)
    assert isinstance(d, dict) and len(d) == 1
    tag, body = next(iter(d.items()))
    node = ET.Element(tag)
    _to_etree(body, node)
    return node

def xml_to_dict(xml_str):
    """
    Convert an xml file into a python dictionary.
    http://stackoverflow.com/questions/7684333/converting-xml-to-dictionary-using-elementtree

    Example:
    from xml.etree import cElementTree as ET
    filename = '/global/homes/b/bpb/batch_params/xmlfile.xml'
    with open(filename,'r') as fid:
        xml_str = fid.read()

    d = xml_to_dict(xml_str)
    """
    t = ET.XML(xml_str)
    d = etree_to_dict(t)
    return d

def etree_to_dict(t):
    """
    Convert an xml tree into a python dictionary.
    http://stackoverflow.com/questions/7684333/converting-xml-to-dictionary-using-elementtree

    """

    d = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.iteritems():
                dd[k].append(v)
        d = {t.tag: {k:v[0] if len(v) == 1 else v for k, v in dd.iteritems()}}
    if t.attrib:
        d[t.tag].update(('@' + k, v) for k, v in t.attrib.iteritems())
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
                d[t.tag]['#text'] = text
        else:
            d[t.tag] = text
    return d
