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

from metatlas.io.update_lcmsfiles_in_lims import EXTENSIONS

# imports for the xml to dictionary round trip
from collections import Mapping
import six
from pathlib2 import PurePath

import time

BATCH_FILE_PATH = '/global/common/software/m2650/mzmine_parameters/batch_files/'
BINARY_PATH = '/global/common/software/m2650/mzmine_parameters/MZmine'


# new stuff:

import sys
import os
import pathlib
import argparse
from subprocess import call

# you need this!
# https://github.com/LabKey/labkey-api-python
#
# sys.path.insert(0,'/Users/bpb/repos/labkey-api-python/')
# sys.path.insert(0,'/global/homes/b/bpb/repos/labkey-api-python/')
# import labkey as lk
from labkey.api_wrapper import APIWrapper
import hashlib
import requests
import json
import pandas as pd
# import hashlib
import numpy as np
import re
from datetime import datetime, time as dtime
from subprocess import check_output
import time
import math
# sys.path.insert(0,'/global/homes/b/bpb/repos/metatlas')
# from metatlas.untargeted import mzmine_batch_tools_adap as mzm
from metatlas.plots import dill2plots as dp
import collections
from ast import literal_eval
from copy import deepcopy
import xmltodict
import zipfile

key_file = '/global/cfs/cdirs/metatlas/labkey_user.txt'
with open(key_file,'r') as fid:
    api_key = fid.read().strip()
labkey_server='metatlas.nersc.gov'
project_name='LIMS/'
api = APIWrapper(labkey_server, project_name, use_ssl=True,api_key=api_key)

def get_parent_folders_from_lcmsruns(get_groups=False):
#     SELECT DISTINCT parent_dir FROM lcmsrun_plus
    sql = """SELECT DISTINCT parent_dir FROM lcmsrun_plus"""
    if get_groups==True:
        sql = """SELECT 
lcmsrun.mzml_file.filename AS mzml_file,
regexp_replace(lcmsrun.name, '.*/([^/]*)/[^/]*$', '\1') AS parent_dir,

split_part(regexp_replace(lcmsrun.name, '.*/[^/]*/([^/]*)$', '\1'), '_', 13) AS file_name_field_12

FROM lcmsrun"""
    schema = 'lists'
    sql_result = api.query.execute_sql(schema, sql,max_rows=1e6)
    if sql_result is None:
        print(('execute_sql: Failed to load results from ' + schema + '.' + table))
        return None
    else:
        df = pd.DataFrame(sql_result['rows'])
        df = df[[c for c in df.columns if not c.startswith('_')]]
        return df


def get_files_from_disk(directory,extension):
    """
    Get on disk with date
    """
    get_with_date = ''.join(['find %s -iname "*%s"' % (directory,extension),' -printf "%Ts SplitThat%p\n"'])
    files = check_output(get_with_date, shell=True)
    files = files.decode('utf-8').splitlines()
    files = [f.split('SplitThat') for f in files]
    dates = [int(f[0].strip()) for f in files]
    files = [f[1].strip() for f in files]
    return dates,files

def complex_name_splitter(filename,
                          extensions=set(['raw', 'tab', 'gz', 'pactolus', 'mzML', 'd','h5']),
                         strippath='/global/project/projectdirs/metatlas/raw_data'):

    #Get the filename
    basename = os.path.basename(filename)
    #Path is everything not filename
    pathname = filename.replace(basename,'')
    #Don't store the basepath since files will likely move
    pathname = pathname.replace(strippath,'')
    pathname = pathname.strip('/')

    #remove extension, but keep any internal . separeted content
    pieces = set(basename.split('.')) - extensions
    name = '.'.join(pieces)
    name = name.replace('_spectral-hits','')

    #this will be a basename that has typically two folders
    #ths should not have an extension
    new_name = os.path.join(pathname,name)
    return new_name

def hash_bytestr_iter(bytesiter, hasher, ashexstr=False):
    for block in bytesiter:
        hasher.update(block)
    return hasher.hexdigest() if ashexstr else hasher.digest()

def make_sha256(afile, blocksize=65536):
    sha = hashlib.sha256()
    with open(afile, 'rb') as f:
        while True:
            data = f.read(blocksize)
            if not data:
                break
            sha.update(data)
    return format(sha.hexdigest())
# def make_sha256(fname):
#     return hash_bytestr_iter(file_as_blockiter(open(fname, 'rb')), hashlib.sha256())

def get_acqtime_from_mzml(mzml_file):
    startTimeStamp=None
    with open(mzml_file) as mzml:
        for line in mzml:
            if 'startTimeStamp' in line:
                startTimeStamp = line.split('startTimeStamp="')[1].split('"')[0].replace('T',' ').rstrip('Z')
                break
#     print startTimeStamp
    if not '-infinity' in startTimeStamp:
        date_object = datetime.strptime(startTimeStamp, '%Y-%m-%d %H:%M:%S')
        utc_timestamp = int(time.mktime(date_object.timetuple()))
    else:
        utc_timestamp = int(0)
    return utc_timestamp


def get_table_from_lims(table,columns=None):
    if columns is None:
        sql = """SELECT * FROM %s;"""%table
    else:
        sql = """SELECT %s FROM %s;"""%(','.join(columns),table)
    # base execute_sql
    schema = 'lists'
    sql_result = api.query.execute_sql(schema, sql,max_rows=1e6)
    if sql_result is None:
        print(('execute_sql: Failed to load results from ' + schema + '.' + table))
        return None
    else:
        df = pd.DataFrame(sql_result['rows'])
        df = df[[c for c in df.columns if not c.startswith('_')]]
        return df

def update_table_in_lims(df,table,method='update',max_size=1000):

    """
    Note: Do ~1000 rows at a time.  Any more and you get a 504 error.  Maybe increasing the timeout would help.
    In the header, timeout is a variable that gets set.  Check to see what its set to.  Maybe increasing it would let
    more rows be updated at a time

    Use it like this:
    update_table_in_lims(df_lims,'mzml_files')
    whatever is in 'name' or 'Key' will replace whatever used to be there with the other columns
    """
#     if columns is None:
#         cols = df.columns
#     else:
#         cols = pd.unique([index_column] + columns)
    # One of the cols needs to be the index column (almost always: Key or Name)

    N = math.ceil(float(df.shape[0]) / max_size)
    for sub_df in np.array_split(df, N):
        payload = sub_df.to_dict('records')
        if method=='update':
            api.query.update_rows('lists', table, payload,timeout=10000)
        elif method=='insert':
            api.query.insert_rows('lists', table, payload,timeout=10000)
        elif method=='delete':
            api.query.delete_rows('lists', table, payload,timeout=10000)
        else:
            print(('ERROR: Nothing to do.  Method %s is not programmed'%method))
        print('updated %d rows in %s'%(df.shape[0],table))

def get_union_of_all_lcms_names(tables=['mzml_file','hdf5_file','pactolus_file','raw_file','spectralhits_file']):
    # sort out the lcmsrun table
    sql = ['select name from %s'%t for t in tables]
    sql = ' union '.join(sql)

#     con = lk.utils.create_server_context(labkey_server, project_name, use_ssl=True,)
    # base execute_sql
    schema = 'lists'
    sql_result = api.query.execute_sql(schema, sql,max_rows=1e6)
    if sql_result is None:
        print(('execute_sql: Failed to load results from ' + schema + '.' + table))
    else:
        return [r['name'] for r in sql_result['rows']]


def update_lcmsrun_names(tables=['mzml_file','hdf5_file','pactolus_file','raw_file','spectralhits_file']):
    #get all the names in the various raw data tables
    names = get_union_of_all_lcms_names(tables)
    #get all the names in lcmsrun (rawdata relationship) table
    lcmsruns = get_table_from_lims('lcmsrun',columns=['name'])
    lcmsruns = lcmsruns['name'].tolist()

    # this is likeley a recently uploaded file that was just created
    missing_from_lcmsruns = list(set(names) - set(lcmsruns))

    #hopefully there aren't any of these, but always good to check
    extra_in_lcmsruns = list(set(lcmsruns) - set(names))

    #add missing ones
    if len(missing_from_lcmsruns)>0:
        temp = pd.DataFrame()
        temp['name'] = missing_from_lcmsruns
        update_table_in_lims(temp,'lcmsrun',method='insert')

    #remove extra ones
    if len(extra_in_lcmsruns)>0:
        sql = """SELECT Key FROM lcmsrun where name IN (%s);"""%','.join(['\'%s\''%e for e in extra_in_lcmsruns])
        # print(sql)
        schema = 'lists'
        sql_result = api.query.execute_sql(schema, sql,max_rows=1e6)
        if sql_result is None:
            print(('execute_sql: Failed to load results from ' + schema + '.' + table))
        #     return None
        else:
            temp = pd.DataFrame(sql_result['rows'])
            temp = temp[[c for c in temp.columns if not c.startswith('_')]]
        #     return df

        if temp.shape[0]>0:
            update_table_in_lims(temp,'lcmsrun',method='delete')

    return missing_from_lcmsruns,extra_in_lcmsruns

def update_lcmsrun_matrix(file_type):
    lcmsruns = get_table_from_lims('lcmsrun',columns=['Key','name',file_type])
    lcmsruns.fillna(-1,inplace=True) #replace None indices so absolute value below has something to work on
    lcmsruns.rename(columns={file_type:'%s_existing'%file_type},inplace=True)
    data = get_table_from_lims(file_type,columns=['Key','name'])

    df = pd.merge(lcmsruns,data,on='name',how='inner')
    df.rename(columns={'Key_x':'Key','Key_y':file_type},inplace=True)
    df = df[abs(df['%s_existing'%file_type]-df[file_type])>0]
    df.drop(columns=['name','%s_existing'%file_type],inplace=True)
    print((df.shape))

    if df.shape[0]>0:
        update_table_in_lims(df,'lcmsrun',method='update')#,index_column='Key',columns=None,labkey_server='metatlas-dev.nersc.gov',project_name='/LIMS'):
    print('done updating')

def get_lcmsrun_matrix():#labkey_server='metatlas-dev.nersc.gov',project_name='/LIMS'):
    sql = 'select '
    for f in ['mzml','hdf5','raw','spectralhits','pactolus']:
        sql = '%s %s_file.filename as %s_filename,'%(sql,f,f)
    sql = '%s from lcmsrun'%sql

    # base execute_sql
    schema = 'lists'
    sql_result = api.query.execute_sql(schema, sql,max_rows=1e8)
    if sql_result is None:
        print(('execute_sql: Failed to load results from ' + schema + '.' + table))
        return None
    else:
        lcmsruns = pd.DataFrame(sql_result['rows'])
        return lcmsruns



def update_file_conversion_tasks(task,lcmsruns=None,file_conversion_tasks=None):#,labkey_server='metatlas-dev.nersc.gov',project_name='/LIMS'):
    """
    gets current tasks and current files and determines if new tasks need to be made:

    task will be:['mzml_to_hdf5','raw_to_mzml','mzml_to_spectralhits','mzml_to_pactolus']

    """
    input_type = task.split('_')[0]
    output_type = task.split('_')[-1]

    if file_conversion_tasks is None:
        file_conversion_tasks = get_table_from_lims('file_conversion_task',columns=['Key','input_file','output_file','task','status'])
    # task_idx = file_conversion_tasks['task']==task

    if lcmsruns is None:
        lcmsruns = get_lcmsrun_matrix()

    done_input_files = lcmsruns.loc[pd.notna(lcmsruns['%s_filename'%input_type]),'%s_filename'%input_type]
    done_output_files = lcmsruns.loc[pd.notna(lcmsruns['%s_filename'%output_type]),'%s_filename'%output_type]

    task_idx = file_conversion_tasks['task']==task
    inputfile_idx = file_conversion_tasks['input_file'].isin(done_input_files)
    outputfile_idx = file_conversion_tasks['output_file'].isin(done_output_files)

    # This finds where output file exists
    done_tasks_idx = (task_idx) & (outputfile_idx)
    if sum(done_tasks_idx)>0:
        update_table_in_lims(file_conversion_tasks.loc[done_tasks_idx,['Key']],'file_conversion_task',method='delete')#,labkey_server=labkey_server,project_name=project_name)
        print(('%s: There are %d tasks where output file exist and will be removed'%(task,file_conversion_tasks[done_tasks_idx].shape[0])))
    
    
    # This finds where input file is missing
    done_tasks_idx = (task_idx) & (~inputfile_idx)
    if sum(done_tasks_idx)>0:
        update_table_in_lims(file_conversion_tasks.loc[done_tasks_idx,['Key']],'file_conversion_task',method='delete')#,labkey_server=labkey_server,project_name=project_name)
        print(('%s: There are %d tasks where input file is missing and will be removed'%(task,file_conversion_tasks[done_tasks_idx].shape[0])))

    right_now_str = datetime.now().strftime("%Y%m%d %H:%M:%S")

    idx = (pd.notna(lcmsruns['%s_filename'%input_type])) & (pd.isna(lcmsruns['%s_filename'%output_type]))

    temp = pd.DataFrame()

    temp['input_file'] = lcmsruns.loc[idx,'%s_filename'%input_type]
    temp['output_file'] = temp['input_file'].apply(lambda x: re.sub('\%s$'%EXTENSIONS[input_type],'%s'%EXTENSIONS[output_type],x))
    temp['task'] = task
    temp['status'] = STATUS['initiation']
    temp['log'] = 'detected: %s'%right_now_str
    temp.reset_index(drop=True,inplace=True)
    cols = temp.columns

    temp = pd.merge(temp,file_conversion_tasks.add_suffix('_task'),left_on=['input_file','output_file'],right_on=['input_file_task','output_file_task'],how='outer',indicator=True)
    new_tasks = temp[temp['_merge']=='left_only'].copy()
    new_tasks = new_tasks[cols]
    new_tasks.reset_index(drop=True,inplace=True)

    print(("There are %d new tasks"%new_tasks.shape[0]))

    if new_tasks.shape[0]>0:
        update_table_in_lims(new_tasks,'file_conversion_task',method='insert')

def update_file_table(file_table):
    file_type = file_table.split('_')[0]
    v = GETTER_SPEC[file_type]
    print(('Getting %s files from disk'%(file_type)))
    dates,files = get_files_from_disk(PROJECT_DIRECTORY,v['extension'])
    if len(files)>0:
        df = pd.DataFrame(data={'filename':files,'file_type':file_type,'timeepoch':dates})
        df['basename'] = df['filename'].apply(os.path.basename)
        df['name'] = df['filename'].apply(complex_name_splitter) #make a name for grouping associated content
    else:
        df = pd.DataFrame()
        df['filename'] = 'None'
        df['file_type'] = file_type
        df['timeepoch'] = 0
        df['basename'] = 'None'
        df['name'] = 'None'
    
    print(('\tThere were %d files on disk'%len(files)))

    cols = ['filename','name','Key']
    df_lims = get_table_from_lims(v['lims_table'],columns=cols)
    print(('\tThere were %d files from LIMS table %s'%(df_lims.shape[0],v['lims_table'])))

        
    diff_df = pd.merge(df, df_lims,on=['filename','name'], how='outer', indicator='Exist')
    diff_df = diff_df.loc[diff_df['Exist'] != 'both'] #(left_only, right_only, or both)
    print(('\tThere are %d different'%diff_df.shape[0]))
    print('')
#         diff_df.fillna('',inplace=True)
    diff_df['parameters'] = 1

    cols = ['file_type','filename','timeepoch','basename','name']
    temp = diff_df.loc[diff_df['Exist']=='left_only',cols]
    if temp.shape[0]>0:
        update_table_in_lims(temp,file_table,method='insert')#,index_column='Key',columns=None,labkey_server='metatlas-dev.nersc.gov',project_name='/LIMS'):

    cols = ['Key','filename']
    temp = diff_df.loc[diff_df['Exist']=='right_only',cols]
    temp['Key'] = temp['Key'].astype(int)
    if temp.shape[0]>0:
        update_table_in_lims(temp,file_table,method='delete')#,index_column='Key',columns=None,labkey_server='metatlas-dev.nersc.gov',project_name='/LIMS'):
#         df.to_csv('/global/homes/b/bpb/Downloads/%s_files.tab'%k,index=None,sep='\t')

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

    
def build_untargeted_filename(output_dir,parent_dir,polarity,file_type):
    """
        file_spec = {'peak-area-mzmine':'peak-area.csv',
                'mzmine-runner':'mzmine.sh',
                'msms-mzmine':'_MSMS.mgf',
                'peak-height-mzmine':'_peak-height.csv',
                'gnps-uuid-fbmn':'_gnps-uuid.txt',
                'fbmn-runner':'fbmn.sh',
                'fbmn-sbatch':'fbmn-sbatch.sbatch',
                'mzmine-outlog':'-mzmine.out',
                'batch-params-mzmine':'_batch-params.xml',
                'quant-fbmn':'quant.csv',
                'gnps-fbmn-network':'_gnps-fbmn-network.graphml',
                'mzmine-sbatch':'mzmine-sbatch.sbatch',
                'mzmine-errlog':'-mzmine.err',
                'metadata':'_metadata.tab',
                'fbmn-errlog':'fbmn.err',
                'fbmn-outlog':'fbmn.out',
                'gnps-download':'_gnps-download.zip'}
    """
    file_spec = {'peak-area-mzmine':'peak-area.csv',
                'mzmine-runner':'_mzmine.sh',
                'msms-mzmine':'_MSMS.mgf',
                'peak-height-mzmine':'_peak-height.csv',
                'gnps-uuid-fbmn':'_gnps-uuid.txt',
                'fbmn-runner':'fbmn.sh',
                'fbmn-sbatch':'fbmn-sbatch.sbatch',
                'mzmine-outlog':'-mzmine.out',
                'batch-params-mzmine':'_batch-params.xml',
                'quant-fbmn':'quant.csv',
                'gnps-fbmn-network':'_gnps-fbmn-network.graphml',
                'mzmine-sbatch':'_mzmine-sbatch.sbatch',
                'mzmine-errlog':'-mzmine.err',
                'metadata':'_metadata.tab',
                'fbmn-errlog':'fbmn.err',
                'fbmn-outlog':'fbmn.out',
                'gnps-download':'_gnps-download.zip'}
    pathname = os.path.join(output_dir,'%s_%s'%(parent_dir,polarity))
    filename = '%s_%s%s'%(parent_dir,polarity,file_spec[file_type])
    filename = os.path.join(pathname,filename)
    return filename

def check_gnps_status(taskid):
    url = 'https://gnps.ucsd.edu/ProteoSAFe/status_json.jsp?task=%s'%taskid
    d = requests.get(url)
    try:
        d = json.loads(d.text)
        if 'status' in d:
            return d['status'], float(d['workflow_version'].split('_')[-1])
        else:
            return 'N/A',0.0
    except:
        return 'N/A',0.0

def download_gnps_graphml(taskid,outfile):
    url = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=gnps_molecular_network_graphml/" % (taskid)
    positive_graphml = "%s.graphml"%taskid
    d = requests.get(url)
    with open(outfile, "w") as fid:
        fid.write(d.text)

def submit_mzmine_jobs(polarity='positive',polarity_short='pos'):
    """
    finds initiated mzmine tasks
    submit them
    changes to running
    """
    tasktype='mzmine'
    df = get_table_from_lims('untargeted_tasks')
    update_df = []
    for i,row in df[df['%s_%s_status'%(tasktype,polarity_short)]=='01 initiation'].iterrows():
        pathname = os.path.join(row['output_dir'],'%s_%s'%(row['parent_dir'],polarity))
        runner_filename = os.path.join(pathname,'%s_%s_mzmine.sh'%(row['parent_dir'],polarity))
        if os.path.isfile(runner_filename)==True:
            with open(runner_filename,'r') as fid:
                task = call(fid.read(),shell=True)
                print(task)
            df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '04 running'
            update_df.append(i)
    if len(update_df)>0:
        cols = ['Key',
                '%s_%s_status'%(tasktype,polarity_short)]
        update_table_in_lims(df.loc[df.index.isin(update_df),cols],'untargeted_tasks',method='update')
        
def update_mzmine_status_in_untargeted_tasks(polarity='positive',polarity_short='pos'):
    """
    finds running mzmine tasks
    checks if they have output
    changes to complete if yes
    sets fbmn to initiation
    """
    tasktype='mzmine'
    df = get_table_from_lims('untargeted_tasks')
    update_df = []
    c1 = df['%s_%s_status'%(tasktype,polarity_short)]=='04 running'
    c2 = df['%s_%s_status'%(tasktype,polarity_short)]=='01 initiation'
    for i,row in df[(c1) | (c2)].iterrows():
        pathname = os.path.join(row['output_dir'],'%s_%s'%(row['parent_dir'],polarity))
        peakheight_filename = os.path.join(pathname,'%s_%s_peak-height.csv'%(row['parent_dir'],polarity))
        if os.path.isfile(peakheight_filename)==True:
            #the job was a success
            df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '07 complete'
            df.loc[i,'%s_%s_status'%('fbmn',polarity_short)] = '01 initiation'
            update_df.append(i)
    if len(update_df)>0:
        cols = ['Key',
                '%s_%s_status'%(tasktype,polarity_short),
               '%s_%s_status'%('fbmn',polarity_short)]
        update_table_in_lims(df.loc[df.index.isin(update_df),cols],'untargeted_tasks',method='update')
        
def update_fbmn_status_in_untargeted_tasks(polarity='positive',polarity_short='pos',latest_version=28.2):
    """
    finds running fbmn tasks
    checks if they have uuid
    checks their status
    changes to complete if yes
    sets spectral hits to initiation
    """
    tasktype='fbmn'
    df = get_table_from_lims('untargeted_tasks')
    update_df = []
    c1 = df['%s_%s_status'%(tasktype,polarity_short)]=='04 running'
    c2 = df['%s_%s_status'%(tasktype,polarity_short)]=='01 initiation'
    c3 = df['%s_%s_status'%(tasktype,polarity_short)]=='08 hold'
    for i,row in df[(c1) | (c2) | (c3)].iterrows():
        pathname = os.path.join(row['output_dir'],'%s_%s'%(row['parent_dir'],polarity))
        fbmn_filename = os.path.join(pathname,'%s_%s_gnps-uuid.txt'%(row['parent_dir'],polarity))
        graphml_filename = os.path.join(pathname,'%s_%s_gnps-fbmn-network.graphml'%(row['parent_dir'],polarity))
        if os.path.isfile(fbmn_filename)==True:
            #the job was submitted
            with open(fbmn_filename,'r') as fid:
                my_text = fid.read().strip()
            taskid = my_text.split('=')[-1]
            status,version = check_gnps_status(taskid)
            print('%s %.2f for %s'%(status,version,fbmn_filename))
            if version<latest_version:
                df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '01 initiation'
#                 df.loc[i,'%s_%s_status'%('gnps_msms_hits',polarity_short)] = '01 initiation'
                download_gnps_graphml(taskid,graphml_filename)
                update_df.append(i)
            elif status=='DONE':          
                df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '07 complete'
#                 df.loc[i,'%s_%s_status'%('gnps_msms_hits',polarity_short)] = '01 initiation'
                download_gnps_graphml(taskid,graphml_filename)
                update_df.append(i)
            elif status=='FAILED':          
                df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '09 error'
#                 df.loc[i,'%s_%s_status'%('gnps_msms_hits',polarity_short)] = '08 hold'
                update_df.append(i)
            elif status=='RUNNING':
                df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '04 running'
#                 df.loc[i,'%s_%s_status'%('gnps_msms_hits',polarity_short)] = '08 hold'
                update_df.append(i)
            else:
                df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '01 initiation'
#                 df.loc[i,'%s_%s_status'%('gnps_msms_hits',polarity_short)] = '08 hold'
                update_df.append(i)
#         else:
#             print('%s is not a file'%fbmn_filename)
    if len(update_df)>0:
        cols = ['Key',
                '%s_%s_status'%(tasktype,polarity_short)] # ,               '%s_%s_status'%('gnps_msms_hits',polarity_short)
        update_table_in_lims(df.loc[df.index.isin(update_df),cols],'untargeted_tasks',method='update')        

def write_new_mzmine_params(gsheet_params,my_polarity,files,basepath,parent_dir):
    """
    takes the generic mzmine parameters
    changes it to the correct polarity
    adds in the files
    saves the xml file
    """
    params = deepcopy(gsheet_params)

            #     new_basename = os.path.join(m['basedir'],'output')
            #replace the polarity in the crop filter module
    for k,v in params.items():
        if 'polarity' in k:
            params[k] = my_polarity.upper()

    # #rename all the output files
    for k,v in params.items():
        try:
            if 'placeholder_filename' in v:
                if 'gnps-job' in v:
                    params[k] = v.replace('placeholder_filename',parent_dir)
                else:
                    params[k] = v.replace('placeholder_filename',os.path.join(basepath,parent_dir))
        except TypeError:
            pass

    # #This is a good place to make the values strings.  The xml maker needs strings later on so might as well do it here
    str_d = {}
    for k,v in params.items():
        str_d[k] = str(v)

    # # #unflatten it
    param_dict_unflat = unflatten(str_d)


    new_raw_data = {'@method': 'net.sf.mzmine.modules.rawdatamethods.rawdataimport.RawDataImportModule',
                    'parameter': {'@name': 'Raw data file names',
                                  'file': files.tolist()}}
    # new_raw_data['parameter']['file'] = mzmine_things[i]['file_list']
    param_dict_unflat['batch']['batchstep'].insert(0,new_raw_data)# str_d.keys()

    xml_string = xmltodict.unparse(param_dict_unflat)

    batch_filename = '%s_batch-params.xml'%os.path.join(basepath,parent_dir)
    with open(batch_filename,'w') as fid:
        fid.write('%s'%xml_string)
    return batch_filename


def submit_fbmn_jobs(polarity='positive',polarity_short='pos',N=15):
    """
    finds initiated mzmine tasks
    submit them
    changes to running
    """
    tasktype='fbmn'
    df = get_table_from_lims('untargeted_tasks')
#     df = df[(df['parent_dir'].str.contains('202'))] #  & (df['parent_dir'].str.contains('AK'))
    update_df = []
    count = 0
    for i,row in df[df['%s_%s_status'%(tasktype,polarity_short)]=='01 initiation'].iterrows():
        pathname = os.path.join(row['output_dir'],'%s_%s'%(row['parent_dir'],polarity))
        runner_filename = os.path.join(pathname,'%s_%s_fbmn.sh'%(row['parent_dir'],polarity))
        if os.path.isfile(runner_filename)==True:
            with open(runner_filename,'r') as fid:
                task = call(fid.read(),shell=True)
                time.sleep(5)
                print(task)
            df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '04 running'
            update_df.append(i)
            count += 1
        if count==N:
            break
    if len(update_df)>0:
        cols = ['Key',
                '%s_%s_status'%(tasktype,polarity_short)]
        update_table_in_lims(df.loc[df.index.isin(update_df),cols],'untargeted_tasks',method='update')

def get_mzmine_param_dict(gdrive_file='params20190719_v2p39_IsotopeFilter_ADAP_DeDup',param_id=2):
    gsheet_params = dp.get_google_sheet(notebook_name=gdrive_file,sheet_name='Sheet1')

    new_cols = []
    for c in gsheet_params.columns:
        if c.startswith('('):
            new_cols.append(literal_eval(c))
        else:
            new_cols.append(c)
    gsheet_params.columns=new_cols

    gsheet_params = gsheet_params[gsheet_params['param_id']==param_id]
    gsheet_params = gsheet_params.to_dict(orient='records')[-1]
    gsheet_params.pop('param_id',None)

    gsheet_params = collections.OrderedDict(gsheet_params)
    return gsheet_params

gsheet_params = get_mzmine_param_dict()
gsheet_params_idx = get_mzmine_param_dict(param_id=5)

def write_mzmine_sbatch_and_runner(basepath,batch_filename,parent_dir,num_files):
    mzmine_launcher = get_latest_mzmine_binary(version='MZmine-2.39')
    sbatch_filename = '%s_mzmine-sbatch.sbatch'%os.path.join(basepath,parent_dir)
    runner_filename = '%s_mzmine.sh'%os.path.join(basepath,parent_dir)
    s = '%s %s'%(mzmine_launcher,batch_filename)
    with open(sbatch_filename,'w') as fid:
        if num_files<51:
            fid.write('%s\n%s\n'%(SLURM_HEADER.replace('slurm','%s-%s'%(os.path.join(basepath,parent_dir),'mzmine')),s))
        else:
            print('asdfasdfasdf')
            fid.write('%s\n%s\n'%(SLURM_BIGMEM_HEADER.replace('slurm','%s-%s'%(os.path.join(basepath,parent_dir),'mzmine')),s))
    with open(runner_filename,'w') as fid:
        if num_files<51:
            fid.write('sbatch %s'%sbatch_filename)
        else:
            fid.write('module load esslurm\n')
            fid.write('sbatch %s\n'%sbatch_filename)
            fid.write('module unload esslurm\n')

def write_fbmn_sbatch_and_runner(basepath,parent_dir):
    runner_filename = '%s_fbmn.sh'%os.path.join(basepath,parent_dir)

    python_binary = '/global/common/software/m2650/python3-metatlas-cori/bin/python'
    python_file = '/global/homes/b/bpb/repos/metatlas/notebooks/workspace/mzmine/send_to_gnps.py'
    python_args = '--basedir %s --basename %s --override True'%(basepath,parent_dir)
    with open(runner_filename,'w') as fid:
        fid.write('%s %s %s\n'%(python_binary,python_file,python_args))

def update_num_features():
    import subprocess
    df_tasks = get_table_from_lims('untargeted_tasks')
    # Get num features
    found = 0
    keep_rows = []
    for i,row in df_tasks.iterrows():
        for polarity in ['positive','negative']:
            peakheight = build_untargeted_filename(row['output_dir'],row['parent_dir'],polarity,'peak-height-mzmine')
            if os.path.isfile(peakheight):
                cmd = ['wc','-l','%s'%peakheight]
                result = subprocess.run(cmd, stdout=subprocess.PIPE)
        #         n = int(subprocess.check_output().split()[0])
                num = result.stdout.split()[0]
                if len(num)>0:
                    num = int(int(num) - 1)
                    if num != row['num_%s_features'%polarity[:3]]:
                        df_tasks.loc[i,'num_%s_features'%polarity[:3]] = num
                        keep_rows.append(i)
                        print(num,row['num_%s_features'%polarity[:3]],result.stderr,result.stdout)
                        print('')

    cols = [c for c in df_tasks.columns if (c.endswith('_features')) & (c.startswith('num_'))]
    cols = cols + ['Key']
    print(cols)
    if len(keep_rows)>0:
        temp = df_tasks.loc[df_tasks.index.isin(keep_rows),cols].copy()
        temp.fillna(0,inplace=True)
        update_table_in_lims(temp,'untargeted_tasks',method='update')


from pyteomics import mgf
import numpy as np
import scipy.sparse as sp
def read_mgf(filename):
    df = []
    with mgf.MGF(filename) as reader:
        for spectrum in reader:
    #         count += 1
            d = spectrum['params']
            d['spectrum'] = np.array([spectrum['m/z array'],spectrum['intensity array']])
            d['pepmass'] = d['pepmass'][0]
            df.append(d)
    ref_df = pd.DataFrame(df)
    return ref_df


def update_num_msms():
    import subprocess
    df_tasks = get_table_from_lims('untargeted_tasks')
    # Get num features
    found = 0
    keep_rows = []
    for i,row in df_tasks.iterrows():
        for polarity in ['positive','negative']:
            filename = build_untargeted_filename(row['output_dir'],row['parent_dir'],polarity,'msms-mzmine')
            if os.path.isfile(filename):
                cmd = 'cat %s | grep FEATURE_ID| wc -l'%filename
                result = subprocess.run(cmd,stdout=subprocess.PIPE, shell=True)
                num = result.stdout.strip()
                if len(num)>0:
                    num = int(num)
                    if num != row['num_%s_msms'%polarity[:3]]:
                        df_tasks.loc[i,'num_%s_msms'%polarity[:3]] = num
                        keep_rows.append(i)
                        print(num,row['num_%s_msms'%polarity[:3]],result.stderr,result.stdout)
#                         print('')

    cols = [c for c in df_tasks.columns if (c.endswith('_msms')) & (c.startswith('num_'))]
    cols = cols + ['Key']
    print(cols)
    if len(keep_rows)>0:
        temp = df_tasks.loc[df_tasks.index.isin(keep_rows),cols].copy()
        temp.fillna(0,inplace=True)
        update_table_in_lims(temp,'untargeted_tasks',method='update')
        

def update_new_untargeted_tasks(update_lims=True):
    """
    given all directories that are in the table for potential untargeted tasks
    and all the directories in the raw data folders
    return all the directories in the raw data folders
    that are not in the untargeted tasks
    
    The strip command is because there is one folder that ends in a 
    space and labkey doesn't allow this
    """
    
    
    df = get_table_from_lims('lcmsrun_plus')

    df.drop(columns=['mzml_file_container'],inplace=True)
    df.replace('',np.nan,inplace=True)

    #Check that files have been sitting around for at least 3 hours (note time zones may vary)
    time_check = df.groupby('parent_dir')['timeepoch'].max() < (time.time()-3*60*60)
    time_check_folders = time_check[time_check==True].index.tolist()


    df_untargeted = get_table_from_lims('untargeted_tasks')
    
    all_folders = df.loc[df['polarity'].isin(['POS','NEG']),'parent_dir'].unique()

    all_folders = [a.strip() for a in all_folders]
    if df_untargeted.shape[0]>0:
        folders_in_tasks = df_untargeted['parent_dir']
        folders_in_tasks = [a.strip() for a in folders_in_tasks]
    else:
        folders_in_tasks = []
    print(len(all_folders),len(folders_in_tasks))
    new_folders = np.setdiff1d(all_folders,folders_in_tasks,)
    print(len(new_folders),len(all_folders),len(folders_in_tasks))
    new_folders = list(set(new_folders) & set(time_check_folders))
    print(len(new_folders),len(all_folders),len(folders_in_tasks))
    

    # get file counts
    pos_count = df[df['polarity']=='POS'].groupby('parent_dir')['mzml_file'].count()
    neg_count = df[df['polarity']=='NEG'].groupby('parent_dir')['mzml_file'].count()
    missing = np.setdiff1d(new_folders,pos_count.index.tolist())
    for m in missing:
        pos_count[m] = 0
    missing = np.setdiff1d(new_folders,neg_count.index.tolist())
    for m in missing:
        neg_count[m] = 0
    
    outdir = '/project/projectdirs/metatlas/projects/untargeted_tasks'
    
    
#     idx1 = ~df_untargeted['mzmine_pos_status'].str.contains('complete')
#     idx2 = ~df_untargeted['mzmine_pos_status'].str.contains('not relevant')
#     idx3 = ~df_untargeted['mzmine_neg_status'].str.contains('complete')
#     idx4 = ~df_untargeted['mzmine_neg_status'].str.contains('not relevant')
#     df_untargeted_new = df_untargeted[(idx1) | (idx2)]
#     df_untargeted_new = df_untargeted[(idx3) | (idx4)]
#     print(df.shape[0])
#     df = df[df['parent_dir'].isin(df_untargeted_new['parent_dir'])]
#     print(df.shape[0])
    print('making metadata pos')
    # make the metadata sheets
    files = df[df['polarity']=='POS'].groupby('parent_dir')
    files = [(d,g) for d,g in files]
    pos_metadata_files = {}
    pos_filelist = {}
    for block in files:
        my_polarity = 'positive'
        polarity_short = 'pos'
        parent_dir = '%s_%s'%(block[0],my_polarity)
        basepath = os.path.join(outdir,parent_dir)
        if not os.path.isdir(basepath):
            os.mkdir(basepath)
        metadata_filename = '%s_%s.tab'%(parent_dir,'metadata')
        metadata_filename = os.path.join(basepath,metadata_filename)
        temp = block[1][['mzml_file','sample_group']].copy()
        temp['mzml_file'] = temp['mzml_file'].fillna('')
        temp['sample_group'] = temp['sample_group'].fillna('')
        temp['sample_group'] = temp['sample_group'].apply(lambda x: x.lower())
        ugroups = temp['sample_group'].unique()
        ugroups = [g for g in ugroups if 'exctrl' in g]
        temp.rename(columns={'mzml_file':'filename','sample_group':'ATTRIBUTE_sampletype'},inplace=True)
        temp['filename'] = temp['filename'].apply(lambda x: os.path.basename(x))
        cols = ['CONTROL','CASE','ATTRIBUTE_media']
        for i,g in enumerate(ugroups):
            if i<(len(cols)):
                temp[cols[i]] = g
#             else:
#                 print('too many controls!!! %s'%g)
            
            
        temp.to_csv(metadata_filename,sep='\t',index=False)
        pos_metadata_files[block[0]] = metadata_filename
        pos_filelist[block[0]] = block[1]['mzml_file']
        
    print('making metadata neg')
    # make the metadata sheets
    files = df[df['polarity']=='NEG'].groupby('parent_dir')
    files = [(d,g) for d,g in files]
    neg_metadata_files = {}
    neg_filelist = {}
    for block in files:
        my_polarity = 'negative'
        polarity_short = 'neg'
        parent_dir = '%s_%s'%(block[0],my_polarity)
        basepath = os.path.join(outdir,parent_dir)
        if not os.path.isdir(basepath):
            os.mkdir(basepath)
        metadata_filename = '%s_%s.tab'%(parent_dir,'metadata')
        metadata_filename = os.path.join(basepath,metadata_filename)
        temp = block[1][['mzml_file','sample_group']].copy()
        temp['mzml_file'] = temp['mzml_file'].fillna('')
        temp['sample_group'] = temp['sample_group'].fillna('')
        temp['sample_group'] = temp['sample_group'].apply(lambda x: x.lower())
        ugroups = temp['sample_group'].unique()
        ugroups = [g for g in ugroups if 'exctrl' in g]
        temp.rename(columns={'mzml_file':'filename','sample_group':'ATTRIBUTE_sampletype'},inplace=True)
        temp['filename'] = temp['filename'].apply(lambda x: os.path.basename(x))
        cols = ['CONTROL','CASE','ATTRIBUTE_media']
        for i,g in enumerate(ugroups):
            if i<(len(cols)):
                temp[cols[i]] = g
#             else:
#                 print('too many controls!!! %s'%g)
        temp.to_csv(metadata_filename,sep='\t',index=False)
        neg_metadata_files[block[0]] = metadata_filename
        neg_filelist[block[0]] = block[1]['mzml_file']
    
    print('There are %d new_folders'%len(new_folders))
    if len(new_folders)>0:
        new_folders = pd.DataFrame(data={'parent_dir':new_folders,'num_pos_files':pos_count[new_folders],'num_neg_files':neg_count[new_folders]})
        new_folders['pos_metadata_file'] = ''
        new_folders['neg_metadata_file'] = ''
        for i,row in new_folders.iterrows():
            if row['parent_dir'] in pos_metadata_files.keys():
                new_folders.loc[i,'pos_metadata_file'] = pos_metadata_files[row['parent_dir']]
                basepath = os.path.join(outdir,'%s_%s'%(row['parent_dir'],'positive'))
                parent_dir = '%s_%s'%(row['parent_dir'],'positive')
                if '_idx_' in parent_dir.lower():
                    batch_filename = write_new_mzmine_params(gsheet_params_idx,'positive',pos_filelist[row['parent_dir']],basepath,parent_dir)
                else:
                    batch_filename = write_new_mzmine_params(gsheet_params,'positive',pos_filelist[row['parent_dir']],basepath,parent_dir)
                write_mzmine_sbatch_and_runner(basepath,batch_filename,parent_dir,pos_filelist[row['parent_dir']].shape[0])
                write_fbmn_sbatch_and_runner(basepath,parent_dir)
            if row['parent_dir'] in neg_metadata_files.keys():
                new_folders.loc[i,'neg_metadata_file'] = neg_metadata_files[row['parent_dir']]
                basepath = os.path.join(outdir,'%s_%s'%(row['parent_dir'],'negative'))
                parent_dir = '%s_%s'%(row['parent_dir'],'negative')
                if '_idx_' in parent_dir.lower():
                    batch_filename = write_new_mzmine_params(gsheet_params_idx,'negative',neg_filelist[row['parent_dir']],basepath,parent_dir)
                else:
                    batch_filename = write_new_mzmine_params(gsheet_params,'negative',neg_filelist[row['parent_dir']],basepath,parent_dir)
                write_mzmine_sbatch_and_runner(basepath,batch_filename,parent_dir,neg_filelist[row['parent_dir']].shape[0])
                write_fbmn_sbatch_and_runner(basepath,parent_dir)
                
        new_folders['file_conversion_complete'] = False
        new_folders['conforming_filenames'] = False
        new_folders['mzmine_pos_status'] = '01 initiation'
        new_folders['mzmine_neg_status'] = '01 initiation'
        new_folders['fbmn_pos_status'] = '13 waiting'
        new_folders['fbmn_neg_status'] = '13 waiting'
        new_folders['gnps_msms_hits_pos_status'] = '08 hold'
        new_folders['gnps_msms_hits_neg_status'] = '08 hold'
        new_folders['output_dir'] = outdir
        new_folders['mzmine_parameter_sheet'] = 'params20190719_v2p39_IsotopeFilter_ADAP_DeDup'
        if '_idx_' in parent_dir.lower():
            new_folders['mzmine_parameter_row'] = 5
        else:
            new_folders['mzmine_parameter_row'] = 2
        new_folders['conforming_filenames'] = True
        new_folders['file_conversion_complete'] = True
        cols = [c for c in new_folders.columns if c.endswith('_pos_status')]
        new_folders.loc[new_folders['num_pos_files']==0,cols] = '12 not relevant'
        cols = [c for c in new_folders.columns if c.endswith('_neg_status')]
        new_folders.loc[new_folders['num_neg_files']==0,cols] = '12 not relevant'
        if update_lims==True:
            update_table_in_lims(new_folders,'untargeted_tasks',method='insert',max_size=1000)
    return new_folders


def count_scans_by_class(x):
    d = {}
    for i in [0.4,0.5,0.7,0.9]:
        d['MQScore_gte_%.1f'%i] = len(x.loc[x['MQScore']>=i,'#Scan#'].unique())
    for i in [3,6,9,12,15]:
        d['SharedPeaks_gte_%d'%i] = len(x.loc[x['SharedPeaks']>=i,'#Scan#'].unique())
    return pd.Series(d, index=d.keys())

def get_gnps_hits(parent_dir,output_dir,polarity,status,override=False):
    gnps_uuid_file = build_untargeted_filename(output_dir,
                                 parent_dir,
                                     polarity,
                                     'gnps-uuid-fbmn')
    gnps_zip_output_file = os.path.join(os.path.dirname(gnps_uuid_file),'%s_%s_gnps-download.zip'%(parent_dir,polarity))
    if (os.path.isfile(gnps_uuid_file)) & ('complete' in status):
        with open(gnps_uuid_file,'r') as fid:
            url = fid.read()
        gnps_uuid = url.split('=')[-1]
        if (override==True) | (not os.path.isfile(gnps_zip_output_file)):
            with zipfile.ZipFile(gnps_zip_output_file, 'r') as archive:
                hits_file = [f for f in archive.namelist() if 'DB_result' in f]
                if len(hits_file)>0:
                    hits_file = hits_file[-1]
                    hits_fid = archive.open(hits_file)
                    df = pd.read_csv(hits_fid,sep='\t')
    
                    return df
    return None



def download_file_from_server_endpoint(server_endpoint,local_file_path):
    response=requests.post(server_endpoint)
    if response.status_code==200:
        # Write the file contents in the response to a file specified by local_file_path
        with open(local_file_path,'wb') as local_file:
            for chunk in response.iter_content(chunk_size=128):
                local_file.write(chunk)
    else:
        print('Can not do this one: %s'%os.path.basename(local_file_path))

def get_gnps_zipfile(parent_dir,output_dir,polarity,status,override=False):
    gnps_uuid_file = build_untargeted_filename(output_dir,
                                 parent_dir,
                                     polarity,
                                     'gnps-uuid-fbmn')
#     print(gnps_uuid_file)
#     print(polarity)
    gnps_zip_output_file = os.path.join(os.path.dirname(gnps_uuid_file),'%s_%s_gnps-download.zip'%(parent_dir,polarity))
    if (os.path.isfile(gnps_uuid_file)) & ('complete' in status):
        if (override==True) | (not os.path.isfile(gnps_zip_output_file)):
            with open(gnps_uuid_file,'r') as fid:
                url = fid.read()
            my_uuid = url.split('=')[-1]
            gnps_url = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task=%s&view=download_cytoscape_data&show=true"%my_uuid
            print(gnps_url)
            print(gnps_zip_output_file)
            download_file_from_server_endpoint(gnps_url, gnps_zip_output_file)
            print('')
    else:
        print('FAIL',parent_dir,status)
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
# # /////////////////////////////////////////////////////////////////////
SLURM_HEADER = """#!/bin/bash
#SBATCH -t 04:00:00
#SBATCH -C haswell
#SBATCH -N 1
#SBATCH --error="slurm.err"
#SBATCH --output="slurm.out"
#SBATCH -q realtime
#SBATCH -A m1541
#SBATCH --exclusive
module load java

"""


# /////////////////////////////////////////////////////////////////////
# # /////////////////// CORI BIGMEME QUEUE SBATCH PARAMS ////////////////////
# # /////////////////////////////////////////////////////////////////////
# SLURM_HEADER = """#!/bin/bash
# #SBATCH -N 1
# #SBATCH -A m1541
# #SBATCH -t 00:30:00

# #SBATCH --clusters=escori
# #SBATCH --qos=bigmem
# #SBATCH --job-name=my_big_job
# #SBATCH --mem=550GB
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
SLURM_BIGMEM_HEADER = """#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --error="slurm.err"
#SBATCH --output="slurm.out"
#SBATCH --qos=jgi_shared
#SBATCH -A pkscell
#SBATCH -C skylake
#SBATCH -t 8:00:00
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

        if isinstance(d, str):
            root.text = d
        elif isinstance(d, dict):
            for k,v in d.items():
                assert isinstance(k, str)
                if k.startswith('#'):
                    assert k == '#text' and isinstance(v, str)
                    root.text = v
                elif k.startswith('@'):
                    assert isinstance(v, str)
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
