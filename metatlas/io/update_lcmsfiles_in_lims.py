from __future__ import absolute_import
from __future__ import print_function
import sys
import os
import argparse

# you need this!
# https://github.com/LabKey/labkey-api-python
#
# sys.path.insert(0,'/Users/bpb/repos/labkey-api-python/')
# sys.path.insert(0,'/global/homes/b/bpb/repos/labkey-api-python/')
# import labkey as lk
from labkey.api_wrapper import APIWrapper

# import requests
# import json
import pandas as pd
# import hashlib
import numpy as np
import re
from datetime import datetime, time as dtime
from subprocess import check_output
import time
import math

EXTENSIONS = {'mzml':'.mzML',
        'hdf5':'.h5',
        'spectralhits':'_spectral-hits.tab.gz',
        'pactolus':'.pactolus.gz',
        'raw':'.raw'}

STATUS = {'initiation':'01 initiation',
      'running':'04 running',
      'custom':'06 custom',
      'complete':'07 complete',
      'hold':'08 hold',
      'error':'09 error',
      'submitted':'10 submitted',
      'corrupted':'11 corrupted'}

PROJECT_DIRECTORY = '/global/project/projectdirs/metatlas/raw_data'

GETTER_SPEC = {'raw':{'extension':'.raw',
                            'lims_table':'raw_file'},
               'mzml':{'extension':'.mzml',
                            'lims_table':'mzml_file'},
               'hdf5':{'extension':'.h5',
                            'lims_table':'hdf5_file'},
               'spectralhits':{'extension':'_spectral-hits.tab.gz',
                            'lims_table':'spectralhits_file'},
              'pactolus':{'extension':'.pactolus.gz',
                            'lims_table':'pactolus_file'}}

key_file = '/global/cfs/cdirs/metatlas/labkey_user.txt'
with open(key_file,'r') as fid:
    api_key = fid.read().strip()
labkey_server='metatlas-dev.nersc.gov'
project_name='LIMS/'
api = APIWrapper(labkey_server, project_name, use_ssl=True,api_key=api_key)


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

def file_as_blockiter(afile, blocksize=65536):
    with afile:
        block = afile.read(blocksize)
        while len(block) > 0:
            yield block
            block = afile.read(blocksize)

def make_sha256(fname):
    return hash_bytestr_iter(file_as_blockiter(open(fname, 'rb')), hashlib.sha256())

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
        print('updated')

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
    df = pd.DataFrame(data={'filename':files,'file_type':file_type,'timeepoch':dates})
    df['basename'] = df['filename'].apply(os.path.basename)
    df['name'] = df['filename'].apply(complex_name_splitter) #make a name for grouping associated content
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
        
def main():
    # print command line arguments
    parser = argparse.ArgumentParser(description='a command line tool for updating the metatlas labkey lims at nersc.')
    # parser.add_argument('-m2t','--ms2_tolerance', help='tolerance in Daltons for ms2', type=float,default=0.01)
    # parser.add_argument('-m1pn','--ms1_pos_neutralizations', help='adducts to neutralize for in ms1: 1.007276,18.033823,22.989218', type=float,nargs='+',default=[1.007276,18.033823,22.989218])
    parser.add_argument('-update_file_tables','--update_file_tables', help='Update the file tables', default=False)
    parser.add_argument('-update_lcmsrun_names','--update_lcmsrun_names', help='Update the names in the lcmsruns matrix', default=False)
    parser.add_argument('-update_lcmsrun_items','--update_lcmsrun_items', help='Update the associations in the lcmsruns matrix', default=False)
    parser.add_argument('-update_fileconversion_tasks','--update_fileconversion_tasks', help='Update the file conversions tasks', default=False)

    # parser.add_argument('-n','--num_cores', help='number of cores to use for multiprocessing', type=int,default=32)
    # parser.add_argument('-overwrite','--overwrite', help='Overwrite pre-existing file(s): True/False', type=bool,default=False)
    args = vars(parser.parse_args())
    # trees = np.load(args['tree_file'])
    print(args)
    
    if str2bool(args['update_file_tables'])==True:
        # 1. Update each table individually
        tables=['mzml_file','hdf5_file','pactolus_file','raw_file','spectralhits_file']
        for t in tables:
            update_file_table(t)

    if str2bool(args['update_lcmsrun_names'])==True:
        # 2. When complete, update the rows in lcmsrun matrix
        missing_from_lcmsruns,extra_in_lcmsruns = update_lcmsrun_names()

    if str2bool(args['update_lcmsrun_items'])==True:
        # 3. populate lcmsrun table that associates all the file types under one entry
        tables=['mzml_file','hdf5_file','pactolus_file','raw_file','spectralhits_file']
        for t in tables:
            update_lcmsrun_matrix(t)

    # 4. remove any file conversion tasks that have already occured
    # TODO!!!!!# TODO!!!!!# TODO!!!!!
    # TODO!!!!!# TODO!!!!!# TODO!!!!!
    # TODO!!!!!# TODO!!!!!# TODO!!!!!
    # TODO!!!!!# TODO!!!!!# TODO!!!!!
    if str2bool(args['update_fileconversion_tasks'])==True:
        # 5. populate any file conversion tasks that need to occur
        lcmsruns = get_lcmsrun_matrix() #this could be moved up abote to step 3 and save a few queries
        file_conversion_tasks = get_table_from_lims('file_conversion_task',columns=['Key','input_file','output_file','task','status'])
        for task in ['mzml_to_hdf5','raw_to_mzml','mzml_to_spectralhits','mzml_to_pactolus']:
            update_file_conversion_tasks(task,lcmsruns=lcmsruns,file_conversion_tasks=file_conversion_tasks)

if __name__ == "__main__":
    main()

# LOOK FOR INCOMPLETE FILES MZML FILES
# # </indexedmzML>
# df = pd.read_excel('lcmsrun_2020-07-02_22-41-09.xlsx')
# bad_files = []
# for f in df['Filename']:
#     with open(f,'r') as fid:
#         mzml = fid.read().strip().endswith('</indexedmzML>')
# #     print(mzml)
#     if mzml==False:
#         bad_files.append(f)

# for f in bad_files:
#     print(f)
#     print('')
