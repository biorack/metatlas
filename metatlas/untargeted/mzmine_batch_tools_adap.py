from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import sys
import json
import pandas as pd
import re
import glob as glob
import os
from collections import defaultdict
from xml.etree import cElementTree as ET

from metatlas.datastructures import metatlas_objects as metob
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.plots import dill2plots as dp

# imports for the xml to dictionary round trip
from collections import Mapping
import six
from pathlib2 import PurePath
from six.moves import map

try:
    six.string_types
except NameError:  # python3
    six.string_types = str

BATCH_FILE_PATH = '/global/common/software/m2650/mzmine_parameters/batch_files/'
BINARY_PATH = '/global/common/software/m2650/mzmine_parameters/MZmine'


# new stuff:


    


#     <batchstep method="net.sf.mzmine.modules.peaklistmethods.orderpeaklists.OrderPeakListsModule">
#         <parameter name="Peak lists" type="BATCH_LAST_PEAKLISTS"/>
#     </batchstep>




    

    
#copy files here to keep I/O off low-performance filesystems
DATA_PATH = '/global/cscratch1/sd/bpb/raw_data'

#we don't need to request haswell on genepool partition of Cori
#remove this line
#SBATCH -C haswell

# /////////////////////////////////////////////////////////////////////
# /////////////////// REALTIME QUEUE SBATCH PARAMS ////////////////////
# /////////////////////////////////////////////////////////////////////
# SLURM_HEADER = """#!/bin/bash
# #SBATCH -t 04:00:00
# #SBATCH -C haswell
# #SBATCH -N 1
# #SBATCH --error="slurm.err"
# #SBATCH --output="slurm.out"
# #SBATCH -q realtime
# #SBATCH -A m1541
# #SBATCH --exclusive
# module load java

# """


# /////////////////////////////////////////////////////////////////////
# /////////////////// CORI REGULAR SBATCH PARAMS //////////////////////
# /////////////////////////////////////////////////////////////////////
# SLURM_HEADER = """#!/bin/bash
# #SBATCH -N 1 -c 64
# #SBATCH --exclusive
# #SBATCH --error="slurm.err"
# #SBATCH --output="slurm.out"
# #SBATCH --qos=genepool
# #SBATCH -A pkscell
# #SBATCH -t 24:00:00
# #SBATCH -L project

# """


#Alicia Clum: The best nodes we have right now are ExVivo, they are 1.5 Tb nodes and very fast you can submit there by changing to --qos=jgi_shared and adding -C skylake. Prior to submitting you must type "module load esslurm" since these nodes are controlled by a different scheduler.
# Set the python to this one:
#/global/common/software/m2650/mzmine_parameters/MZmine/MZmine-2.39/startMZmine_NERSC_Headless_Cori_exvivo.sh 
# /////////////////////////////////////////////////////////////////////
# /////////////////// SKYLAKE 1.5TB QUEUE SBATCH PARAMS ///////////////
# /////////////////////////////////////////////////////////////////////
SLURM_HEADER = """#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --error="slurm.err"
#SBATCH --output="slurm.out"
#SBATCH --qos=jgi_shared
#SBATCH -A pkscell
#SBATCH -C skylake
#SBATCH -t 48:00:00
#SBATCH -L project

"""


def calc_hit_vector(n,df):
    """
    for 0,1,2 n will be 3
    for 0,1,2,3 n will be 4
    df is the count from true_positives
    this function makes it a percentation of hits will sum n or more hits in last element
    """
    m = np.zeros((n))
    s = df['count']/df['count'].sum()
    nf = df['num_features']
    for i in s.index:
        if (i>(len(m)-1)) & (len(m) >1):
            m_idx = len(m)-1
        else:
            m_idx = nf[i] 
        m[m_idx] = m[m_idx] + s[i]
    return m

def summarize_results(n,true_pos_filename,base_path,project_path,feature_table_extension,rt_column,mz_column,headerrows,sep,mz_tolerance,rt_tolerance):
    """
    """
    path  = os.path.join(base_path,project_path)
    feature_file = glob.glob(os.path.join(path,'*%s'%feature_table_extension))
    if len(feature_file)>0:
        feature_file = feature_file[-1]
    else:
        return np.zeros(n),np.zeros(n), 0

    new_path = os.path.join(path,'true_pos_results')
    if not os.path.isdir(new_path):
        os.mkdir(new_path)

    new_basename = os.path.basename(feature_file).replace(feature_table_extension,'height.xlsx')
    output_filename = os.path.join(new_path,new_basename)
    if os.path.isfile(feature_file): #has the file been made already?
        with open(feature_file,'r') as fid:
            s = fid.read()

        if len(s)>0: #does the file have anything in it?
            df_experimental = pd.read_csv(feature_file,sep=sep,skiprows=headerrows)
            if df_experimental.shape[0] > 0: #are there any rows?
                if 'hilic' in feature_file.lower():
                    sheetname = 'HILIC_POS'
                    df_true_pos = pd.read_excel(true_pos_filename,sheet_name=sheetname)
                else:
                    sheetname = 'CSH_POS'
                    df_true_pos = pd.read_excel(true_pos_filename,sheet_name=sheetname)
                istd_count,bio_count,df_grouped,df_hits,total_count = prepare_true_positive_and_export(output_filename,df_experimental,df_true_pos,rt_column=rt_column,mz_column=mz_column,mz_tolerance=mz_tolerance,rt_tolerance=rt_tolerance)
                return calc_hit_vector(n,istd_count), calc_hit_vector(n,bio_count), total_count.loc[0,'total']
    return np.zeros(n),np.zeros(n), 0

def make_count_of_knowns(df_hits,df_true_pos):
    df_grouped = df_hits[['true_pos_index','CompoundName_truepos','experimental_feature_idx']]
    df_grouped.set_index(['CompoundName_truepos'],inplace=True)

    df_grouped = df_grouped.groupby(['true_pos_index']).count()
    df_grouped = pd.merge(df_grouped,df_true_pos,left_index=True,right_index=True,how='outer')
    df_grouped.rename(columns={'experimental_feature_idx':'num_features'},inplace=True)
    return df_grouped

def map_features_to_known(df_experimental,df_true_pos,rt_column='row retention time',mz_column='row m/z',mz_tolerance=0.01,rt_tolerance=0.1):
    feature_array = df_experimental[[mz_column,rt_column]].values
    reference_array = df_true_pos[['MZ','RT']].values

    idx = np.isclose(feature_array[:,None,:], reference_array, rtol=0.0, atol=[mz_tolerance,rt_tolerance]).all(axis=2) #order is m/z, rt, polarity
    feature_idx, reference_idx = np.where(idx)

    df_hits = df_true_pos.loc[reference_idx].copy()
    df_hits['experimental_feature_idx'] = feature_idx
    df_hits = pd.merge(df_hits,df_true_pos,left_index=True,right_index=True,how='outer',suffixes=['_x','_truepos'])
    df_hits.drop(columns=['%s_x'%c for c in df_true_pos.columns],inplace=True)
    df_hits = pd.merge(df_hits,df_experimental,how='left',left_on='experimental_feature_idx',right_index=True,suffixes=['_truepos','_experimental'])
    df_hits.index.name = 'true_pos_index'
    df_hits.reset_index(inplace=True)
    return df_hits

def summarize_count_per_type(df_grouped,cpd_type='ISTD'):
    """
    cpd_type is either 'ISTD' or 'TargetCPD'
    """
    cpd_count = df_grouped[df_grouped['Type']==cpd_type][['num_features','MZ']].groupby('num_features').count()
    cpd_count.reset_index(inplace=True)
    cpd_count.rename(columns={'MZ':'count'},inplace=True)
    return cpd_count

def prepare_true_positive_and_export(output_filename,df_experimental,df_true_pos,rt_column='row retention time',mz_column='row m/z',mz_tolerance=0.01,rt_tolerance=0.1):
    df_hits = map_features_to_known(df_experimental,df_true_pos,rt_column=rt_column,mz_column=mz_column,mz_tolerance=0.01,rt_tolerance=0.1)
    df_grouped = make_count_of_knowns(df_hits,df_true_pos)
    istd_count = summarize_count_per_type(df_grouped,cpd_type='ISTD')
    bio_count = summarize_count_per_type(df_grouped,cpd_type='TargetCPD')
    params = pd.DataFrame(columns=['mz_tolerance','rt_tolerance'],data=[[mz_tolerance,rt_tolerance]])
    total_count = pd.DataFrame(columns=['total'],data=[[df_experimental.shape[0]]])

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(output_filename, engine='xlsxwriter')

    # Write each dataframe to a different worksheet.
    params.to_excel(writer, sheet_name='params')
    istd_count.to_excel(writer, sheet_name='istd_count')
    bio_count.to_excel(writer, sheet_name='bio_count')
    df_grouped.to_excel(writer, sheet_name='df_grouped')
    df_hits.to_excel(writer, sheet_name='df_hits')
    total_count.to_excel(writer, sheet_name='total')
    return istd_count,bio_count,df_grouped,df_hits,total_count
    # Close the Pandas Excel writer and output the Excel file.
#     writer.save()

def mzmine_xml_to_csv(xml_file,csv_file=None,pop_input_files=True,return_df=True):
    """
    given an xml file, turn it into a csv
    optionally return either a dict of the steps or a dataframe of the steps
    """
    with open(xml_file,'r') as fid:
        xml_str = fid.read()

    d = xml_to_dict(xml_str)

#     t = dict_to_etree(d)
#     indent_tree(t)
#     s1 = tree_to_xml(t)
    
    # pop out the files
    if pop_input_files==True:
        raw_data_import = d['batch']['batchstep'].pop(0)
        original_file_list = raw_data_import['parameter']['file']

    # This is a dict representation of all the steps
    dflat = flatten(d,enumerate_types=(list,))

    # This is a tabular representation of all the steps
    df = pd.DataFrame([(k,v) for (k,v) in dflat.items()],columns=['parameter','value']).sort_values('parameter').set_index('parameter',drop=True)
    if csv_file is not None:
        df.to_csv(csv_file)
    
    if return_df==True:
        return df #return the dataframe of the steps
    else:
        return dflat #return the dict of the steps


def make_task_and_job(params):#basedir,basename,polarity,files):
    if not os.path.exists(params['basedir']):
        os.mkdir(params['basedir'])

    xml_str = get_batch_file_template()
    d = xml_to_dict(xml_str)

    # #     Initialize the task and give it values from the user supplied form
    task = metob.MZMineTask()
    task.polarity = params['polarity']
    task.lcmsruns = params['files']
    task.min_peak_duration = params['min_peak_duration']
    task.max_peak_duration = params['max_peak_duration']
    task.rt_tol_perfile = params['rt_tol_perfile']
    task.rt_tol_multifile = params['rt_tol_multifile']
    task.min_num_scans = params['min_num_scans']
    task.smoothing_scans = params['smoothing_scans']
    task.group_intensity_threshold = params['group_intensity_threshold']
    task.min_peak_height = params['min_peak_height']
    task.ms1_noise_level = params['ms1_noise_level']
    task.ms2_noise_level = params['ms2_noise_level']
    task.mz_tolerance = params['mz_tolerance']
    task.peak_to_valley_ratio = params['peak_to_valley_ratio']
    task.min_rt = params['min_rt']
    task.max_rt = params['max_rt']
    task.representative_isotope = params['representative_isotope']
    task.remove_isotopes = params['remove_isotopes']
    task.min_peaks_in_row = params['min_peaks_in_row']
    task.peak_with_msms = params['peak_with_msms']
    task.chromatographic_threshold = params['chromatographic_threshold']
    task.search_for_minimum_rt_range = params['search_for_minimum_rt_range']
    task.minimum_relative_height = params['minimum_relative_height']
    task.mz_range_scan_pairing = params['mz_range_scan_pairing']
    task.rt_range_scan_pairing = params['rt_range_scan_pairing']
    task.gapfill_intensity_tolerance = params['gapfill_intensity_tolerance']
    task.output_csv_height = os.path.join(params['basedir'],'%s_%s_peak_height.csv'%(params['basename'],task.polarity))
    task.output_csv_area = os.path.join(params['basedir'],'%s_%s_peak_area.csv'%(params['basename'],task.polarity))
    task.output_workspace = os.path.join(params['basedir'],'%s_%s.mzmine'%(params['basename'],task.polarity))
    task.output_mgf = os.path.join(params['basedir'],'%s_%s.mgf'%(params['basename'],task.polarity))
    task.input_xml = os.path.join(params['basedir'],'logs','%s_%s.xml'%(params['basename'],task.polarity))
    task.mzmine_launcher = get_latest_mzmine_binary(version=params['mzmine_version'])

    new_d = replace_files(d,params['files'])
    new_d = configure_crop_filter(new_d,task.polarity,params['files'],min_rt=task.min_rt,max_rt=task.max_rt)

    new_d = configure_mass_detection(new_d,task.ms1_noise_level,task.ms2_noise_level)

    new_d = configure_chromatogram_builder(new_d,task.min_num_scans,task.group_intensity_threshold,task.min_peak_height,task.mz_tolerance)
    
    new_d = configure_smoothing(new_d,task.smoothing_scans)

    new_d = configure_peak_deconvolution(new_d,
                                         task.min_peak_height,
                                         task.minimum_relative_height,
                                         task.search_for_minimum_rt_range,
                                         task.chromatographic_threshold,
                                         task.peak_to_valley_ratio,
                                         task.min_peak_duration,
                                         task.max_peak_duration)

    new_d = configure_isotope_search(new_d,
                                         task.mz_tolerance,
                                         task.rt_tol_perfile,
                                         task.representative_isotope,
                                         task.remove_isotopes,
                                    task.polarity)

    new_d = configure_join_aligner(new_d,task.mz_tolerance,task.rt_tol_multifile)

    new_d = configure_gap_filling(new_d,task.mz_tolerance,task.rt_tol_multifile,task.gapfill_intensity_tolerance)

    new_d = configure_rows_filter(new_d,task.min_peaks_in_row,task.peak_with_msms)

    new_d = configure_output(new_d,
                                 task.output_csv_height,
                                task.output_csv_area,
                                 task.output_workspace,
                                 task.output_mgf)


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
#     job_cmd_filtered = make_targeted_mzmine_job(m['basedir'],m['basename'],m['polarity'],m['files'])
    params_filename = os.path.join(m['basedir'],'logs','%s_%s_params.json'%(m['basename'],m['polarity']))
    new_params_filename = os.path.join(m['basedir'],'logs','%s_%s_params-used.json'%(m['basename'],m['polarity']))

    copy_params_command = "cp '%s' '%s'"%(params_filename,new_params_filename)


    with open(sbatch_file_name,'w') as fid:
        fid.write('%s\n'%SLURM_HEADER.replace('slurm.err',err_file_name).replace('slurm.out',out_file_name))
        fid.write('%s\n'%copy_params_command)
        fid.write('%s\n'%job_cmd)

#     bad_words = ['qos', '-p','-C','-L','-t','-N']
#     bad_time = '#SBATCH -t 24:00:00'
#     good_time = '#SBATCH -t 24:00:00\n'

#     bad_node = '-N 1 -c 64'
#     good_node = '#SBATCH -N 1 -c 64\n'


#     with open(sbatch_file_name) as oldfile, open(denovo_sbatch_file_name, 'w') as newfile:
#         for line in oldfile:
#             if not any(bad_word in line for bad_word in bad_words):
#                 newfile.write(line)
#             if bad_time in line:
#                 newfile.write(good_time)
#             if bad_node in line:
#                 newfile.write(good_node)
#                 newfile.write('#SBATCH --mem=494G\n')


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

def get_files(groups,filename_substring,file_filters,keep_strings,is_group=False,return_mzml=True):
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
            
    # only keep files that don't have substrings in list
    if len(file_filters) > 0:
        for i,ff in enumerate(file_filters):
            if i == 0:
                files = [f for f in all_files if not ff in f.name]
            else:
                files = [f for f in files if not ff in f.name]
    else:
        files = all_files
        
    # kick out any files that don't match atleast one of the keep_strings
    keep_this = []
    filter_used = [] #good to keep track if a filter isn't used.  likely a typo
    if len(keep_strings) > 0:
        for i,ff in enumerate(files):
            keep_this.append(any([True if f in ff.name else False for f in keep_strings]))
        for i,ff in enumerate(keep_strings):
            filter_used.append(any([True if ff in f.name else False for f in files]))
        if not all(filter_used):
            for i,f in enumerate(filter_used):
                if f==False:
                    print('%s keep string is not used'%keep_strings[i])
                    
        files = [files[i] for i,j in enumerate(keep_this) if j==True]
    
    files = remove_duplicate_files(files)
    return files


def make_targeted_mzmine_job(basedir,basename,polarity,files):
    if not os.path.exists(basedir):
        os.mkdir(basedir)

    xml_str = get_targeted_batch_file_template()
    d = xml_to_dict(xml_str)

    task = metob.MZMineTask()
    task.polarity = polarity
    task.lcmsruns = files
    new_d = replace_files(d,files)

    project_name = '%s_%s'%(basename,task.polarity)
    task.output_workspace = os.path.join(basedir,project_name,'%s_%s.mzmine'%(basename,task.polarity))
    task.input_xml = os.path.join(basedir,'logs','%s_%s_filtered.xml'%(basename,task.polarity))
    
    task.mzmine_launcher = get_latest_mzmine_binary()

    # new_d = configure_crop_filter(new_d,task.polarity,files)
    # new_d = configure_targeted_peak_detection(new_d,peak_list_filename,intensity_tolerance=1e-4,noise_level=1e4,mz_tolerance=20,rt_tolerance=0.5)
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

def configure_crop_filter(new_d,polarity,files,min_rt=0.01,max_rt=100,fps_string='FPS'):
    """
    
    """
    # identify the element for this change
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'CropFilterModule' in d['@method']][0]
    # Set the filter string
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Raw data files' in d['@name']][0]
    if any([fps_string in f for f in files]):
        new_d['batch']['batchstep'][idx]['parameter'][idx2]['name_pattern'] = '*FPS*'
    else:
        new_d['batch']['batchstep'][idx]['parameter'][idx2]['name_pattern'] = '*'

    # Set the polarity
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Scans' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['polarity'] = polarity.upper()
    
    #set the rt min and rt max use the same idx2 as polarity
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['retention_time'] = {'max':'%.4f'%max_rt,'min':'%.4f'%min_rt}
    
    # new_d['batch']['batchstep'][idx]['parameter'][idx2]['ms_level'] = '1-2'

    return new_d

def configure_mass_detection(new_d,ms1_noise_level=1e4,ms2_noise_level=1e2):
    """
    
    """
    # Find the module
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'MassDetectionModule' in d['@method']]
    #The first idx will be for MS1 and the second will be for MS2
    
    # Set the MS1 attributes
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[0]]['parameter']) if 'Mass detector' in d['@name']][0]
    idx3 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[0]]['parameter'][idx2]['module']) if 'Centroid' in d['@name']][0]
    new_d['batch']['batchstep'][idx[0]]['parameter'][idx2]['module'][idx3]['parameter']['#text'] = '%.2f'%(ms1_noise_level)
    
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[0]]['parameter']) if 'Scans' in d['@name']][0]
    new_d['batch']['batchstep'][idx[0]]['parameter'][idx2]['ms_level'] = '1'
    
    # Set the MS2 attributes
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[1]]['parameter']) if 'Mass detector' in d['@name']][0]
    idx3 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[1]]['parameter'][idx2]['module']) if 'Centroid' in d['@name']][0]
    new_d['batch']['batchstep'][idx[1]]['parameter'][idx2]['module'][idx3]['parameter']['#text'] = '%.2f'%(ms2_noise_level)
    
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[1]]['parameter']) if 'Scans' in d['@name']][0]
    new_d['batch']['batchstep'][idx[1]]['parameter'][idx2]['ms_level'] = '2'

    
    return new_d

def configure_smoothing(new_d,smoothing_scans):
    """
#        <batchstep method="net.sf.mzmine.modules.peaklistmethods.peakpicking.smoothing.SmoothingModule">
#         <parameter name="Peak lists" type="BATCH_LAST_PEAKLISTS"/>
#         <parameter name="Filename suffix">smoothed</parameter>
#         <parameter name="Filter width">9</parameter>
#         <parameter name="Remove original peak list">false</parameter>
#     </batchstep>
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'SmoothingModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Filter width' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(smoothing_scans)

    return new_d


def configure_chromatogram_builder(new_d,min_num_scans,group_intensity_threshold,min_peak_height,mz_tolerance):
    """
    new_d = configure_chromatogram_builder(new_d,task.min_num_scans,task.group_intensity_threshold,task.min_peak_height,task.mz_tolerance)

#     <batchstep method="net.sf.mzmine.modules.masslistmethods.ADAPchromatogrambuilder.ADAPChromatogramBuilderModule">
#         <parameter name="Raw data files" type="ALL_FILES"/>
#         <parameter name="Scans">
#             <ms_level>1</ms_level>
#         </parameter>
#         <parameter name="Mass list">masses</parameter>
#         <parameter name="Min group size in # of scans">5</parameter>
#         <parameter name="Group intensity threshold">1000000.0</parameter>
#         <parameter name="Min highest intensity">80000.0</parameter>
#         <parameter name="m/z tolerance">
#             <absolutetolerance>0.002</absolutetolerance>
#             <ppmtolerance>7.0</ppmtolerance>
#         </parameter>
#         <parameter name="Suffix">chromatograms</parameter>
#     </batchstep>


    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'ADAPChromatogramBuilderModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Min group size' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(min_num_scans)

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Group intensity threshold' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(group_intensity_threshold)

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Min highest intensity' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(min_peak_height)

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
   
    return new_d

def configure_adap_peak_deconvolution(new_d,min_peak_height,minimum_relative_height,search_for_minimum_rt_range,chromatographic_threshold,min_sn_ratio,min_peak_duration,max_peak_duration):
    """
    <parameter name="Algorithm" selected="Wavelets (ADAP)">
    
    <module name="Wavelets (ADAP)">
                <parameter name="S/N threshold">3.0</parameter>
                <parameter name="S/N estimator" selected="Intensity window SN">
                    <module name="Intensity window SN"/>
                    <module name="Wavelet Coeff. SN">
                        <parameter name="Peak width mult.">1.0</parameter>
                        <parameter name="abs(wavelet coeffs.)">true</parameter>
                    </module>
                </parameter>
                <parameter name="min feature height">4500.0</parameter>
                <parameter name="coefficient/area threshold">60.0</parameter>
                <parameter name="Peak duration range">
                    <min>0.0</min>
                    <max>0.5</max>
                </parameter>
                <parameter name="RT wavelet range">
                    <min>0.0</min>
                    <max>0.1</max>
                </parameter>
            </module>
            
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'DeconvolutionModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Algorithm' in d['@name']][0]
    idx3 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module']) if 'Local minimum search' in d['@name']][0]
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Chromatographic threshold' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%chromatographic_threshold
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Search minimum in RT range (min)' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%search_for_minimum_rt_range
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Minimum relative height' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%minimum_relative_height
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Minimum absolute height' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%min_peak_height
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Min ratio of peak top/edge' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%min_sn_ratio
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Peak duration range (min)' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['min'] = '%.3f'%min_peak_duration
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['max'] = '%.3f'%max_peak_duration
    return new_d

def configure_lms_peak_deconvolution(new_d,min_peak_height,minimum_relative_height,search_for_minimum_rt_range,chromatographic_threshold,min_sn_ratio,min_peak_duration,max_peak_duration):
    """
    <parameter name="Algorithm" selected="Local minimum search">
                <module name="Local minimum search">
                <parameter name="Chromatographic threshold">0.75</parameter>
                <parameter name="Search minimum in RT range (min)">0.02</parameter>
                <parameter name="Minimum relative height">0.002</parameter>
                <parameter name="Minimum absolute height">90000.0</parameter>
                <parameter name="Min ratio of peak top/edge">1.03</parameter>
                <parameter name="Peak duration range (min)">
                    <min>0.03</min>
                    <max>1.0</max>
                </parameter>
            </module>
            
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'DeconvolutionModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Algorithm' in d['@name']][0]
    idx3 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module']) if 'Local minimum search' in d['@name']][0]
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Chromatographic threshold' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%chromatographic_threshold
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Search minimum in RT range (min)' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%search_for_minimum_rt_range
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Minimum relative height' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%minimum_relative_height
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Minimum absolute height' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%min_peak_height
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Min ratio of peak top/edge' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%min_sn_ratio
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Peak duration range (min)' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['min'] = '%.3f'%min_peak_duration
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['max'] = '%.3f'%max_peak_duration
    return new_d

def configure_isotope_search(new_d,mz_tolerance,rt_tol_perfile,representative_isotope,remove_isotopes,polarity):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'Isotope' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_perfile)
    
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Representative isotope' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%s'%(representative_isotope)

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Remove original peaklist' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%s'%(str(remove_isotopes).lower())
    
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
    
    return new_d

def configure_join_aligner(new_d,mz_tolerance,rt_tol_multifile):
    """
    # Join aligner has these scores:
#             <parameter name="Minimum absolute intensity">3000.0</parameter>
#             <parameter name="Minimum score">0.6</parameter>
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'JoinAlignerModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_multifile)
    
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Minimum absolute intensity' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = 3000#'%.3f'%(mz_tolerance)
    
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Minimum score' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = 0.6#'%.3f'%(rt_tol_multifile)
    
    return new_d

def configure_rows_filter(new_d,min_peaks_in_row,peak_with_msms):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'RowsFilterModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Minimum peaks in a row' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%d'%min_peaks_in_row
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Minimum peaks in an isotope pattern' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%d'%min_peaks_in_row
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Keep only peaks with MS2 scan (GNPS)' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%d'%peak_with_msms
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

def configure_gap_filling(new_d,mz_tolerance,gapfill_intensity_tolerance,rt_tol_multifile):
    """
    #     <batchstep method="net.sf.mzmine.modules.peaklistmethods.gapfilling.peakfinder.multithreaded.MultiThreadPeakFinderModule">
#         <parameter name="Peak lists" type="BATCH_LAST_PEAKLISTS"/>
#         <parameter name="Name suffix">gap-filled</parameter>
#         <parameter name="Intensity tolerance">0.05</parameter>
#         <parameter name="m/z tolerance">
#             <absolutetolerance>0.001</absolutetolerance>
#             <ppmtolerance>5.0</ppmtolerance>
#         </parameter>
#         <parameter name="Retention time tolerance" type="absolute">0.03</parameter>
#         <parameter name="Remove original peak list">false</parameter>
#     </batchstep>
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'gapfilling.peakfinder' in d['@method']][0]
    
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Intensity tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(gapfill_intensity_tolerance)
    
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_multifile)
    
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    
    
    return new_d

def configure_output(new_d,output_csv_height,output_csv_area,output_workspace,output_mgf):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'CSVExportModule' in d['@method']]
    #the first will be height the second will be area

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[0]]['parameter']) if 'Filename' in d['@name']][0]
    new_d['batch']['batchstep'][idx[0]]['parameter'][idx2]['#text'] = output_csv_height
    
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[1]]['parameter']) if 'Filename' in d['@name']][0]
    new_d['batch']['batchstep'][idx[1]]['parameter'][idx2]['#text'] = output_csv_area
    
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'GNPSExportModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Filename' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = output_mgf
    
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'ProjectSaveAsModule' in d['@method']][0]
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

def get_batch_file_template(loc='bootcamp_adap_template.xml'):
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

    wget $(curl -s https://api.github.com/repos/mzmine/mzmine2/releases/v2.33 | grep 'browser_' | cut -d\" -f4) -O mzmine_latest.zip


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
    mzmine_versions = glob.glob(os.path.join(BINARY_PATH,'*' + os.path.sep))
    if version == 'most_recent':
        most_recent = sorted([os.path.basename(m) for m in mzmine_versions if 'MZmine-' in m])[-1]
    else:
        most_recent = [m.split(os.path.sep)[-2] for m in mzmine_versions if version in m][-1]
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
        print(type(d),d)
        print('\n\n\n')
        if not d:
            pass

        if type(d) is {}.values().__class__:
            d = list(d.values)
        
        if isinstance(d, six.string_types):
            root.text = d
        elif isinstance(d, dict):
            for k,v in d.items():
                assert isinstance(k, six.string_types)
                if k.startswith('#'):
                    assert k == '#text' and isinstance(v, six.string_types)
                    root.text = v
                elif k.startswith('@'):
                    assert isinstance(v, six.string_types)
                    root.set(k[1:], v)
                elif isinstance(v, list):
                    for e in v:
                        _to_etree(e, ET.SubElement(root, k))
                else:
                    _to_etree(v, ET.SubElement(root, k))
#         elif isinstance(d,dict_values):
#             d = [d]
#             _to_etree(d,ET.SubElement(root, k))
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
            for k, v in six.iteritems(dc):
                dd[k].append(v)
        d = {t.tag: {k:v[0] if len(v) == 1 else v for k, v in six.iteritems(dd)}}
    if t.attrib:
        d[t.tag].update(('@' + k, v) for k, v in six.iteritems(t.attrib))
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
                d[t.tag]['#text'] = text
        else:
            d[t.tag] = text
    return d



##########################################################
#### From Here ###########################################
### https://github.com/ianlini/flatten-dict ##############
##########################################################
##########################################################



def tuple_reducer(k1, k2):
    if k1 is None:
        return (k2,)
    else:
        return k1 + (k2,)


def path_reducer(k1, k2):
    import os.path
    if k1 is None:
        return k2
    else:
        return os.path.join(k1, k2)
    


def tuple_splitter(flat_key):
    return flat_key


def path_splitter(flat_key):
    keys = PurePath(flat_key).parts
    return keys

REDUCER_DICT = {
    'tuple': tuple_reducer,
    'path': path_reducer,
}

SPLITTER_DICT = {
    'tuple': tuple_splitter,
    'path': path_splitter,
}


def flatten(d, reducer='tuple', inverse=False, enumerate_types=()):
    """Flatten `Mapping` object.

    Parameters
    ----------
    d : dict-like object
        The dict that will be flattened.
    reducer : {'tuple', 'path', Callable}
        The key joining method. If a `Callable` is given, the `Callable` will be
        used to reduce.
        'tuple': The resulting key will be tuple of the original keys.
        'path': Use `os.path.join` to join keys.
    inverse : bool
        Whether you want invert the resulting key and value.
    enumerate_types : Sequence[type]
        Flatten these types using `enumerate`.
        For example, if we set `enumerate_types` to ``(list,)``,
        `list` indices become keys: ``{'a': ['b', 'c']}`` -> ``{('a', 0): 'b', ('a', 1): 'c'}``.

    Returns
    -------
    flat_dict : dict
    """
    enumerate_types = tuple(enumerate_types)
    flattenable_types = (Mapping,) + enumerate_types
    if not isinstance(d, flattenable_types):
        raise ValueError("argument type %s is not in the flattenalbe types %s"
                         % (type(d), flattenable_types))

    if isinstance(reducer, str):
        reducer = REDUCER_DICT[reducer]
    flat_dict = {}

    def _flatten(d, parent=None):
        key_value_iterable = enumerate(d) if isinstance(d, enumerate_types) else six.viewitems(d)
        for key, value in key_value_iterable:
            flat_key = reducer(parent, key)
            if isinstance(value, flattenable_types):
                _flatten(value, flat_key)
            else:
                if inverse:
                    flat_key, value = value, flat_key
                if flat_key in flat_dict:
                    raise ValueError("duplicated key '{}'".format(flat_key))
                flat_dict[flat_key] = value

    _flatten(d)
    return flat_dict


# def nested_set_dict(d, keys, value):
#     """Set a value to a sequence of nested keys

#     Parameters
#     ----------
#     d : Mapping
#     keys : Sequence[str]
#     value : Any
#     """
#     assert keys
#     key = keys[0]
#     if len(keys) == 1:
#         if key in d:
#             raise ValueError("duplicated key '{}'".format(key))
#         d[key] = value
#         return
#     d = d.setdefault(key, {})
#     nested_set_dict(d, keys[1:], value)
    
    
def nested_set_dict(d, keys, value):
    """Set a value to a sequence of nested keys

    Parameters
    ----------
    d : Mapping
    keys : Sequence[str]
    value : Any
    """
    assert keys
    key = keys[0]
    if len(keys) == 1:
        if type(d) == list:
            d.append(value)
        else:
            d[key] = value
        return

    # the type is a string so make a dict if none exists
    if type(keys[1]) == int:
        if key in d:
            pass
        else:
            d[key] = []
        d = d[key]
    elif type(key)==int:
        if (key+1) > len(d):
            d.append({})
        d = d[key]
    else:
        d = d.setdefault(key, {})
    nested_set_dict(d, keys[1:], value)


def unflatten(d, splitter='tuple', inverse=False):
    """Unflatten dict-like object.

    Parameters
    ----------
    d : dict-like object
        The dict that will be unflattened.
    splitter : {'tuple', 'path', Callable}
        The key splitting method. If a Callable is given, the Callable will be
        used to split.
        'tuple': Use each element in the tuple key as the key of the unflattened dict.
        'path': Use `pathlib.Path.parts` to split keys.

    Tester
    d1 = {'a':{'b':[{'c1':'nested1!','d1':[{'e1':'so_nested1!!!'}]},
               {'c2':'nested2!','d2':[{'e2':'so_nested2!!!'}]},
               {'c3':'nested3!','d3':[{'e3':'so_nested3!!!'}]},
               {'c4':'nested4!','d4':[{'e4':'so_nested4a!!!'},
                                      {'e4':'so_nested4b!!!'},
                                      {'e4':'so_nested4c!!!'},
                                      {'e4':'so_nested4d!!!'},
                                      {'e4':'so_nested4e!!!'}]}]}}    

    Returns
    -------
    unflattened_dict : dict
    """
    if isinstance(splitter, str):
        splitter = SPLITTER_DICT[splitter]
    
    kv = sorted([(k,v) for (k,v) in d.items()])
    unflattened_dict = {}
    for kkvv in kv:
        key_tuple = kkvv[0]
        value = kkvv[1]
        nested_set_dict(unflattened_dict, key_tuple, value)

    return unflattened_dict

