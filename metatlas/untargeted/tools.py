import sys
sys.path.insert(0,'/global/common/software/m2650/labkey-api-python') # https://github.com/LabKey/labkey-api-python
from labkey.api_wrapper import APIWrapper
import numpy as np
import pandas as pd
import glob as glob
from pathlib2 import PurePath
import os
import requests
from bs4 import BeautifulSoup
from datetime import datetime, time
import time
import subprocess
import math
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/metatlas')
from metatlas.tools.validate_filenames import field_exists
import subprocess
import grp

BATCH_FILE_PATH = '/global/common/software/m2650/mzmine_parameters/batch_files/'
BINARY_PATH = '/global/common/software/m2650/mzmine_parameters/MZmine'

key_file = '/global/cfs/cdirs/metatlas/labkey_user.txt'
with open(key_file,'r') as fid:
    api_key = fid.read().strip()
labkey_server='metatlas.lbl.gov'
project_name='LIMS/'
api = APIWrapper(labkey_server, project_name, use_ssl=True,api_key=api_key)

SLURM_PERLMUTTER_REALTIME_HEADER = """#!/bin/bash
#SBATCH -N 1
#SBATCH --error="slurm.err"
#SBATCH --output="slurm.out"
#SBATCH -A m1541
#SBATCH -C cpu
#SBATCH --qos=realtime
#SBATCH -t 2:00:00
"""

SLURM_PERLMUTTER_HEADER = """#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --error="slurm.err"
#SBATCH --output="slurm.out"
#SBATCH -A m2650
#SBATCH -C cpu
#SBATCH --qos=regular
#SBATCH -t 3:00:00
"""

mzine_batch_params_file = "/global/common/software/m2650/mzmine_parameters/batch_files/mzmine-3.7.2-batchparams.xml"
mzine_batch_params_file_iqx = "/global/common/software/m2650/mzmine_parameters/batch_files/IQX-mzmine-3.7.2-batchparams.xml"

def start_script(script=None):
    if script is not None:
        print('\nStarting %s on:'%script, datetime.now())

def end_script(script=None):
    if script is not None:
        print('\nEnding %s on:'%script, datetime.now(),'\n')

def check_for_polarities(output_dir: str, parent_dir: str) -> list:
    """
    Check for positive and/or negative polarity directories on disk
    to decide how to iterate through project
    """
    pos_dir = os.path.join(output_dir, "%s_%s"%(parent_dir, 'positive'))
    neg_dir = os.path.join(output_dir, "%s_%s"%(parent_dir, 'negative'))
    if os.path.exists(pos_dir) and os.path.exists(neg_dir):
        return ['positive', 'negative']
    elif os.path.exists(pos_dir) and not os.path.exists(neg_dir):
        return ['positive']
    elif not os.path.exists(pos_dir) and os.path.exists(neg_dir):
        return ['negative']
    else:
        return None

def subset_df_by_status(df: pd.DataFrame, tasktype: str, status: list, inverse=False) -> pd.DataFrame:
    """
    Subset a dataframe by status
    """
    if inverse == True:
        df_subset = df[~df['%s_pos_status'%tasktype].isin(status) | ~df['%s_neg_status'%tasktype].isin(status)]
    if inverse == False:
        df_subset = df[df['%s_pos_status'%tasktype].isin(status) | df['%s_neg_status'%tasktype].isin(status)]
    return df_subset

def filter_common_bad_project_names(df: pd.DataFrame):
    df = df[~(df['parent_dir'].str.contains(' '))]
    df = df[~(df['parent_dir'].str.contains('&'))]
    df = df[~(df['parent_dir'].str.contains('Niyogi'))]
    df = df[~(df['parent_dir'].str.contains('_partial'))]
    df = df[~(df['parent_dir'].str.contains('_old'))]
    return df

def write_fbmn_tasks_to_file(task_list: list):
    """
    Takes a list of dictionaries from submit_fbmn_jobs
    and writes the fbmn task id to a file in the 
    project directory at untarageted_tasks on perlmutter
    """
    for dataset in task_list:
        experiment = dataset['experiment']
        polarity = dataset['polarity']
        task = dataset['response']['task']
        depository = '/global/cfs/cdirs/metatlas/projects/untargeted_tasks/'
        filename = f'{depository}{experiment}_{polarity}/{experiment}_{polarity}_gnps2-fbmn-task.txt'
        if task:
            with open(filename,'w') as fid:
                fid.write(f"{experiment}_{polarity}={task}\n")
                final_filename = os.path.basename(filename)
                print(f"\t\t\tGNPS2 task file written to {final_filename}")
        else:
            print("\t\t\tError: GNPS2 FBMN task ID not found")
            return

def zip_and_upload_untargeted_results(download_folder = '/global/cfs/cdirs/metatlas/projects/untargeted_outputs', \
                                      output_dir = '/global/cfs/cdirs/metatlas/projects/untargeted_tasks', \
                                      upload=True, overwrite=False, direct_input=None,min_features_admissible=0):
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    status_list = ['07 complete']
    df = subset_df_by_status(df,'fbmn',status_list)
    df = subset_df_by_status(df,'mzmine',status_list)
    if df.empty:
        print("\tNo completed untargeted results to zip and upload!")
        return
    if not df.empty:
        if overwrite == True:
            wait = 10
            print("Notice: Overwrite flag is set to True. This will overwrite any existing zip files. %s seconds to kill the script if you want to stop."%wait)
            time.sleep(wait)
        for i,row in df.iterrows():
            output_zip_archive = os.path.join(download_folder,'%s.zip'%row['parent_dir'])
            polarity_list = check_for_polarities(output_dir,row['parent_dir'])
            if overwrite==True or not os.path.exists(output_zip_archive):
                if polarity_list is None:
                    print("\tWarning! Project %s does not have a negative or a positive polarity directory. Skipping..."%row['parent_dir'])
                    continue
                neg_mzmine_file = os.path.join(output_dir, '%s_negative'%row['parent_dir'], '%s_negative_peak-height.csv'%row['parent_dir'])
                pos_mzmine_file = os.path.join(output_dir, '%s_positive'%row['parent_dir'], '%s_positive_peak-height.csv'%row['parent_dir'])
                neg_fbmn_file = os.path.join(output_dir, '%s_negative'%row['parent_dir'], '%s_negative_gnps2-fbmn-library-results.tsv'%row['parent_dir'])
                pos_fbmn_file = os.path.join(output_dir, '%s_positive'%row['parent_dir'], '%s_positive_gnps2-fbmn-library-results.tsv'%row['parent_dir'])
                if 'negative' in polarity_list and 'positive' in polarity_list:
                    if os.path.exists(neg_mzmine_file) and os.path.exists(pos_mzmine_file) and os.path.exists(neg_fbmn_file) and os.path.exists(pos_fbmn_file):
                        neg_feature_counts = check_peak_height_table(neg_mzmine_file)
                        pos_feature_counts = check_peak_height_table(pos_mzmine_file)
                        if (neg_feature_counts + pos_feature_counts) > min_features_admissible:
                            neg_directory = os.path.join(output_dir, '%s_%s'%(row['parent_dir'], 'negative'))
                            pos_directory = os.path.join(output_dir, '%s_%s'%(row['parent_dir'], 'positive'))
                            cmd = 'zip -rjuq %s %s %s'%(output_zip_archive,neg_directory,pos_directory)
                            os.system(cmd)
                            print("\tNew untargeted results in POS and NEG modes zipped for %s"%row['parent_dir'])
                            if upload == True and os.path.exists(output_zip_archive):
                                upload_to_google_drive(output_zip_archive)
                        else:
                            print("\tWarning! Project %s has less than %s combined features in NEG and POS MZmine peak height \
                                  tables. Skipping..."%(row['parent_dir'],min_features_admissible))
                            continue
                elif not 'negative' in polarity_list and 'positive' in polarity_list:
                    if os.path.exists(pos_mzmine_file) and os.path.exists(pos_fbmn_file):
                        pos_feature_counts = check_peak_height_table(pos_mzmine_file)
                        if pos_feature_counts > min_features_admissible:
                            pos_directory = os.path.join(output_dir, '%s_%s'%(row['parent_dir'], 'positive'))
                            cmd = 'zip -rjuq %s %s'%(output_zip_archive,pos_directory)
                            os.system(cmd)
                            print("\tNew untargeted results in POS mode only zipped for %s"%row['parent_dir'])
                            if upload == True and os.path.exists(output_zip_archive):
                                upload_to_google_drive(output_zip_archive)
                        else:
                            print("\tWarning! Project %s has less than %s features in POS MZmine peak height \
                                  table. Skipping..."%(row['parent_dir'],min_features_admissible))
                            continue
                elif 'negative' in polarity_list and not 'positive' in polarity_list:
                    if os.path.exists(neg_mzmine_file) and os.path.exists(neg_fbmn_file):
                        neg_feature_counts = check_peak_height_table(neg_mzmine_file)
                        if neg_feature_counts > min_features_admissible:
                            neg_directory = os.path.join(output_dir, '%s_%s'%(row['parent_dir'], 'negative'))
                            cmd = 'zip -rjuq %s %s %s'%(output_zip_archive,neg_directory,pos_directory)
                            os.system(cmd)
                            print("\tNew untargeted results in NEG mode only zipped for %s"%row['parent_dir'])
                            if upload == True and os.path.exists(output_zip_archive):
                                upload_to_google_drive(output_zip_archive)
                        else:
                            print("\tWarning! Project %s has less than %s features in NEG MZmine peak height \
                                  table. Skipping..."%(row['parent_dir'],min_features_admissible))
                            continue

def check_peak_height_table(peak_height_file):
    cmd = 'cat %s | wc -l'%peak_height_file
    result = subprocess.run(cmd,stdout=subprocess.PIPE, shell=True)
    if result:
        feature_counts = int(result.stdout.strip())
        return feature_counts
    else:
        return 0

def upload_to_google_drive(file_path):
    project_folder = os.path.basename(file_path)
    print("\t\tUploading zip to Google Drive...")
    cmd = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone copy --size-only {file_path} ben_lbl_gdrive:/untargeted_outputs'
    subprocess.check_output(cmd, shell=True)
    check_upload = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone ls ben_lbl_gdrive:/untargeted_outputs | grep "{project_folder}"'
    check_upload_out = subprocess.check_output(check_upload, shell=True)
    if check_upload_out.decode('utf-8').strip():
        print("\t\t\tGoogle Drive upload confirmed!")
    else:
        print("\t\t\tWarning: Google Drive upload check failed. Upload may not have been successful.")

def submit_quickstart_fbmn(params="",username=""):
    ##### Pick up secret file and extract password
    with open('/global/homes/m/msdata/gnps2/gnps2_%s.txt'%username,'r') as fid:
        t = fid.read()
    t = t.split('\n')[0]
    var, pword = t.split('=')
    password = pword.replace('"', '').replace("'", '').strip()
    if not password:
        print("Password is required to submit FBMN jobs")
        return

    url = "https://gnps2.org/launchworkflow"

    # First lets login with a session
    session = requests.Session()
    login_data = {"username": username, "password": password}
    session.post("https://gnps2.org/login", data=login_data)

    # Submitting the job
    response = session.post(url, data=params)
    return response.json()

def get_untargeted_status(direct_input=None):
    status_df_list = []
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    if df.empty:
        print("\tCould not find untargeted tasks to report!")
        return pd.DataFrame()
    if not df.empty:
        print("\tPrinting status for %s projects"%df.shape[0])
        for i,row in df.iterrows():
            status_df_list.append(row[['parent_dir', 'mzmine_pos_status', 'mzmine_neg_status', 'fbmn_pos_status', 'fbmn_neg_status']])
        total_status_df = pd.concat(status_df_list, axis=0)
        return total_status_df

def recursive_chown(basepath, group_name):
    gid = grp.getgrnam(group_name).gr_gid
    for root, dirs, files in os.walk(basepath):
        for dir_name in dirs:
            os.chown(os.path.join(root, dir_name), -1, gid)
        for file_name in files:
            os.chown(os.path.join(root, file_name), -1, gid)
    os.chown(basepath, -1, gid)  # Also change the ownership of the base directory itself

def download_from_url(url, target_path):
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(target_path, 'wb') as f:
            f.write(response.raw.read())
        return True
    else:
        return False

def download_fbmn_results(output_dir='/global/cfs/cdirs/metatlas/projects/untargeted_tasks',overwrite=False,direct_input=None):
    """
    finds complete fbmn tasks
    downloads the results folder
    renames and moves results to the untargeted_tasks folder
    """
    tasktype='fbmn'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    status_list = ['07 complete','09 error']
    df = subset_df_by_status(df,tasktype,status_list)
    df = subset_df_by_status(df,'mzmine',['07 complete']) # Also want to check that mzmine is complete before downloading fbmn
    if df.empty:
        print("\tNo completed FBMN data to download!")
    if not df.empty:
        for i,row in df.iterrows():
            polarity_list = check_for_polarities(output_dir,row['parent_dir'])
            if polarity_list is None:
                print("\tWarning! Project %s does not have a negative or a positive polarity directory. Skipping..."%row['parent_dir'])
                continue
            for polarity in polarity_list:
                polarity_short = polarity[:3]
                if row['%s_%s_status'%(tasktype,polarity_short)] == '12 not relevant':
                    continue
                if row['%s_%s_status'%(tasktype,polarity_short)] == '09 error':
                    print("\tWarning! FBMN task for %s %s has error status. Not downloading files."%(row['parent_dir'],polarity))
                    continue
                pathname = os.path.join(output_dir,'%s_%s'%(row['parent_dir'],polarity))
                fbmn_filename = os.path.join(pathname,'%s_%s_gnps2-fbmn-task.txt'%(row['parent_dir'],polarity))
                if os.path.isfile(fbmn_filename)==True:
                    with open(fbmn_filename,'r') as fid:
                        taskid = fid.read().split('=')[1].strip()
                    graphml_filename = os.path.join(pathname,'%s_%s_gnps2-fbmn-network.graphml'%(row['parent_dir'],polarity))
                    results_table_filename = os.path.join(pathname,'%s_%s_gnps2-fbmn-library-results.tsv'%(row['parent_dir'],polarity))
                    gnps2_link_filename = os.path.join(pathname,'%s_%s_gnps2-page-link.txt'%(row['parent_dir'],polarity))
                    if os.path.exists(graphml_filename) and os.path.exists(results_table_filename) and os.path.exists(gnps2_link_filename):
                        continue
                    print("\tDownloading FBMN results for %s with task ID %s"%(row['parent_dir'],taskid))
                    if overwrite == True or not os.path.exists(gnps2_link_filename):
                        with open(gnps2_link_filename,'w') as fid:
                            fid.write(f"https://gnps2.org/status?task={taskid}\n")
                            print("\t\tWrote GNPS2 link")
                    if overwrite == True or not os.path.exists(graphml_filename):
                        graphml_download = download_from_url("https://gnps2.org/result?task=%s&viewname=graphml&resultdisplay_type=task"%(taskid), graphml_filename)
                        if graphml_download:
                            print("\t\tDownloaded graphml file")
                        else:
                            print("\t\tError: Failed to download graphml file")
                    if overwrite == True or not os.path.exists(results_table_filename):
                        results_table_download = download_from_url("https://gnps2.org/resultfile?task=%s&file=nf_output/library/merged_results_with_gnps.tsv"%(taskid), results_table_filename)
                        if results_table_download:
                            print("\t\tDownloaded result table file")
                        else:
                            print("\t\tError: Failed to download results table")
                    if os.path.exists(graphml_filename):
                        recursive_chown(graphml_filename, 'metatlas')
                    if os.path.exists(results_table_filename):
                        recursive_chown(results_table_filename, 'metatlas')
                    if os.path.exists(gnps2_link_filename):
                        recursive_chown(gnps2_link_filename, 'metatlas')

def get_recent_mgf_files(output_dir = '/global/cfs/cdirs/metatlas/projects/untargeted_tasks', time_back=30):
    # Get the current time
    now = time.time()
    # Calculate the time time_back days ago
    time_30_days_ago = now - time_back*24*60*60
    # Convert the time to the format used by the find command
    time_30_days_ago = time.strftime('%Y%m%d', time.gmtime(time_30_days_ago))
    # Run the find command with date modified
    command = f"find {output_dir} -name '*.mgf' -type f -newermt {time_30_days_ago}"
    recent_files = subprocess.check_output(command, shell=True).decode().split('\n')
    return recent_files

def remove_contaminant_from_mgf(f):
    if 'negative' in f:
        target = 202.07880
        diff = target * 5 / 1e6
    elif 'positive' in f:
        target = 202.07770
        diff = target * 5 / 1e6
    else:
        print(f"Error: Could not determine polarity of file {f}")
        return
    
    # Open the file for reading
    with open(f, 'r') as file:
        # Read all lines from the file
        lines = file.readlines()

    # Filter the lines that contain a number within the range
    filtered_lines = [line for line in lines if not any(abs(float(num)-target)<diff for num in line.split() if num.replace('.', '', 1).isdigit())]
    
    if len(filtered_lines)<len(lines):
        print(f"\tRemoving {len(lines)-len(filtered_lines)} lines from {f}")
        # Open the file for writing
        with open(f, 'w') as file:
            # Write the lines back
            file.writelines(filtered_lines)


def get_table_from_lims(table,columns=None,max_rows=1e6):
    if columns is None:
        sql = """SELECT * FROM %s;"""%table
    else:
        sql = """SELECT %s FROM %s;"""%(','.join(columns),table)
    # base execute_sql
    schema = 'lists'
    sql_result = api.query.execute_sql(schema, sql,max_rows=max_rows,timeout=10000)
    if sql_result is None:
        print(('Execute_sql: Failed to load results from ' + schema + '.' + table))
        return None
    else:
        df = pd.DataFrame(sql_result['rows'])
        df = df[[c for c in df.columns if not c.startswith('_')]]
        return df

def update_table_in_lims(df,table,method='update',max_size=1000,pause_time=None):

    """
    Note: Do ~1000 rows at a time.  Any more and you get a 504 error.  Maybe increasing the timeout would help.
    In the header, timeout is a variable that gets set.  Check to see what its set to.  Maybe increasing it would let
    more rows be updated at a time
    
    method can be 'update','insert',or 'delete'

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
        #print('Updated %d rows in %s'%(df.shape[0],table))
        if pause_time is not None:
            time.sleep(pause_time)    # Pause this many seconds
            #print('pausing for %d seconds'%pause_time)

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
                'gnps-download':'_gnps-download.zip',
                ''msms-mzmine3':'.mgf'}
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
                'gnps-download':'_gnps-download.zip',
                'blink-hits':'_blinkhits.csv.gz',
                'blink-network-hits':'_blinknetworkhits.csv.gz',
                'simile-hits':'_similehits.csv.gz',
                'msms-mzmine3':'.mgf'}
    pathname = os.path.join(output_dir,'%s_%s'%(parent_dir,polarity))
    filename = '%s_%s%s'%(parent_dir,polarity,file_spec[file_type])
    filename = os.path.join(pathname,filename)
    return filename

def check_gnps2_status(taskid):
    url = 'https://gnps2.org/status?task=%s'%taskid
    d = requests.get(url)
    try:
        soup = BeautifulSoup(d.content, 'html.parser')
        version = soup.find("td", text="Version").find_next_sibling("td").text
        status = soup.find("td", text="Status").find_next_sibling("td").text
        if status and version:
            return status, version
        else:
            return 'N/A','N/A'
    except:
        return 'N/A','N/A'

def submit_fbmn_jobs(direct_input=None):
    tasktype = 'fbmn'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    status_list = ['13 waiting','09 error']
    df = subset_df_by_status(df,tasktype,status_list)
    df = subset_df_by_status(df,'mzmine',['07 complete']) # Also want to check that mzmine is complete before submitting fbmn
    if df.empty:
        print("\tNo new FBMN jobs to submit!")
        return
    if df.shape[0] > 20:
        print('There are too many new projects to be submitted (%s), please check if this is accurate. Exiting script.'%df.shape[0])
        return
    if not df.empty:
        task_list = []
        index_list = []
        print("\tTotal of %s projects(s) with FBMN status %s and MZmine status 'complete' to submit to GNPS2"%(df.shape[0],status_list))
        for i,row in df.iterrows():
            project = row['parent_dir']
            polarity_list = check_for_polarities(row['output_dir'],row['parent_dir'])
            if polarity_list is None:
                print("\t\tWarning! Project %s does not have a negative or a positive polarity directory. Skipping..."%row['parent_dir'])
                continue
            for polarity in polarity_list:
                polarity_short = polarity[:3]
                if row['%s_%s_status'%(tasktype,polarity_short)] == '12 not relevant' or row['%s_%s_status'%('mzmine',polarity_short)] != '07 complete':
                    print("\t\tWarning! It looks like MZmine for %s mode for %s hasn't completed yet. Skipping..."%(project,polarity))
                    continue
                if row['%s_%s_status'%(tasktype,polarity_short)] == '09 error':
                    print("\t\tWarning! %s in %s mode has an error status. Attempting to resubmit..."%(project,polarity))
                _, validate_department, _ = field_exists(PurePath(row['parent_dir']), field_num=1)
                department = validate_department.lower()
                if department =='eb':
                    department = 'egsb'
                if not department in ['jgi','egsb']:
                    print("\t\tWarning! %s does not have a valid department name in the second field. Skipping..."%project)
                    continue
                description = '%s_%s'%(project,polarity)
                spectra_file = f'USERUPLOAD/bpbowen/untargeted_tasks/{project}_{polarity}/{project}_{polarity}.mgf'
                quant_file = f'USERUPLOAD/bpbowen/untargeted_tasks/{project}_{polarity}/{project}_{polarity}_quant.csv'
                metadata_file = f'USERUPLOAD/bpbowen/untargeted_tasks/{project}_{polarity}/{project}_{polarity}_metadata.tab'
                raw_data = f'USERUPLOAD/bpbowen/raw_data/{department}/{project}'
                mgf_filename = os.path.join(row['output_dir'],'%s_%s'%(project,polarity),'%s_%s.mgf'%(project,polarity))
                mgf_lines = count_mgf_lines(mgf_filename)
                if mgf_lines == 0:
                    print("\t\tWarning! %s in %s mode has mgf file but no mgf data. Skipping..."%(project, polarity))
                    continue
                if not os.path.exists(raw_data.replace('USERUPLOAD/bpbowen/','/global/cfs/cdirs/metatlas/')):
                    print("\tWarning! %s does not have raw data. Skipping..."%project)
                    continue
                params = set_fbmn_parameters(description, quant_file, spectra_file, metadata_file, raw_data)
                job_id = submit_quickstart_fbmn(params, "bpbowen")
                task_list.append({'experiment':project,'polarity':polarity,'response':job_id})
                print("\t\tSubmitted FBMN job for %s in %s mode"%(project,polarity))
                df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '04 running'
                if len(task_list)>0:
                    write_fbmn_tasks_to_file(task_list)
                index_list.append(i)

        if len(index_list) > 0:
            index_list = list(set(index_list))
            cols = ['Key',
                    '%s_neg_status'%(tasktype),
                    '%s_pos_status'%(tasktype)]
            df = df.loc[df.index.isin(index_list),cols]
            df.replace('NaN', 0, inplace=True)
            df.fillna(0, inplace=True)
            #print("\tUpdated untargeted tasks rows:\n" + "\n".join(["\t" + line for line in df.to_string(index=False).splitlines()]))
            update_table_in_lims(df,'untargeted_tasks',method='update')
            print("\tFBMN submission complete. Use the *_gnps2-page-link.txt file to check progress")

def count_mgf_lines(mgf_file):
    if os.path.isfile(mgf_file):
        cmd = 'cat %s | grep END | grep IONS | wc -l'%mgf_file
        result = subprocess.run(cmd,stdout=subprocess.PIPE, shell=True)
        if result:
            msms_count = int(result.stdout.strip())
            return msms_count
        else:
            return 0
    else:
        return 0


def submit_mzmine_jobs(new_projects=None,direct_input=None):
    """
    finds initiated mzmine tasks
    submit them
    changes to running
    """
    tasktype = 'mzmine'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    if new_projects is not None:
        df = df[df['parent_dir'].isin(new_projects['parent_dir'].tolist())]
    status_list = ['01 initiation']
    df = subset_df_by_status(df,tasktype,status_list)
    if df.empty:
        print("\tNo new MZmine jobs to submit!")
        return
    if df.shape[0] > 20:
        print('There are too many new projects to be submitted (%s), please check if this is accurate. Exiting script.'%df.shape[0])
        return
    if not df.empty:
        index_list = []
        print("\tTotal of %s new MZmine job(s) with status %s to submit"%(df.shape[0],status_list))
        for i,row in df.iterrows():
            polarity_list = check_for_polarities(row['output_dir'],row['parent_dir'])
            if polarity_list is None:
                print("\t\tWarning! Project %s does not have a negative or a positive polarity directory. Skipping..."%row['parent_dir'])
                continue
            for polarity in polarity_list:
                polarity_short = polarity[:3]
                if row['%s_%s_status'%(tasktype,polarity_short)] == '12 not relevant' or row['%s_%s_status'%('mzmine',polarity_short)] == '07 complete':
                    continue
                pathname = os.path.join(row['output_dir'],'%s_%s'%(row['parent_dir'],polarity))
                submission_script_filename = os.path.join(pathname,'%s_%s_mzmine.sh'%(row['parent_dir'],polarity))
                if os.path.isfile(submission_script_filename)==True:
                    with open(submission_script_filename,'r') as fid:
                        print("\t\tSubmitting %s mode mzmine job for project %s"%(polarity, row['parent_dir']))
                        cmd = fid.read()
                        sbatch_output = subprocess.run(cmd, shell=True,capture_output=True,text=True)
                        print("\t\t\t%s"%sbatch_output.stdout)
                    index_list.append(i)
                    df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '04 running'
        
        if len(index_list) > 0:
            index_list = list(set(index_list))
            cols = ['Key',
                    '%s_neg_status'%(tasktype),
                    '%s_pos_status'%(tasktype)]
            df = df.loc[df.index.isin(index_list),cols]
            df.replace('NaN', 0, inplace=True)
            df.fillna(0, inplace=True)
            #print("\tUpdated untargeted tasks rows:\n" + "\n".join(["\t" + line for line in df.to_string(index=False).splitlines()]))
            update_table_in_lims(df,'untargeted_tasks',method='update')
            print("MZmine submissions complete. Use sqs to monitor progress.")
    
def set_fbmn_parameters(description, quant_file, spectra_file, metadata_file, raw_data):
    params = {
                "description": description,
                "workflowname": "feature_based_molecular_networking_workflow",
                "featurefindingtool": "MZMINE",
                "inputfeatures": quant_file,
                "inputspectra": spectra_file,
                "metadata_filename": metadata_file,
                "input_libraries": "LIBRARYLOCATION/LC/LIBRARY",
                "input_raw_spectra": raw_data,
                "input_supplemental_edges": "",
                "library_analog_search": "0",
                "library_min_cosine": "0.7",
                "library_min_matched_peaks": "3",
                "library_topk": "20",
                "min_peak_intensity": "0.0",
                "networking_max_shift": "1999",
                "networking_min_cosine": "0.7",
                "networking_min_matched_peaks": "3",
                "normalization": "None",
                "pm_tolerance": "0.01",
                "fragment_tolerance": "0.01",
                "precursor_filter": "yes",
                "api": "no"}
    return params

def update_mzmine_status_in_untargeted_tasks(direct_input=None):
    """
    finds running mzmine tasks
    checks if they have output
    changes to complete if yes
    """
    tasktype = 'mzmine'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    status_list = ['01 initiation','04 running']
    df = subset_df_by_status(df,tasktype,status_list)
    if df.empty:
        print("\tNo MZmine jobs to update!")
    if not df.empty:
        print("\tTotal of %s project(s) with MZmine status %s to attempt to update"%(df.shape[0],status_list))
        index_list = []
        for i,row in df.iterrows():
            polarity_list = check_for_polarities(row['output_dir'],row['parent_dir'])
            if polarity_list is None:
                print("\t\tWarning! Project %s does not have a negative or a positive polarity directory. Skipping..."%row['parent_dir'])
                continue
            for polarity in polarity_list:
                polarity_short = polarity[:3]
                if row['%s_%s_status'%(tasktype,polarity_short)] == '12 not relevant' or row['%s_%s_status'%(tasktype,polarity_short)] == '07 complete':
                    continue  ## This helps when one polarity is complete but the other is not
                pathname = os.path.join(row['output_dir'],'%s_%s'%(row['parent_dir'],polarity))
                old_peakheight_filename = os.path.join(pathname,'%s_%s_MSMS_quant.csv'%(row['parent_dir'],polarity))
                peakheight_filename = os.path.join(pathname,'%s_%s_peak-height.csv'%(row['parent_dir'],polarity))
                mgf_filename = os.path.join(pathname,'%s_%s.mgf'%(row['parent_dir'],polarity))
                metadata_filename = os.path.join(pathname,'%s_%s_metadata.tab'%(row['parent_dir'],polarity))
                if (os.path.isfile(mgf_filename) and os.path.isfile(metadata_filename) and (os.path.isfile(peakheight_filename) or os.path.isfile(old_peakheight_filename))):
                    print("\t\tWorking on %s in %s mode"%(row['parent_dir'],polarity))
                    print("\t\t\tAll MZmine output files found in %s directory, continuing..."%polarity)
                    print("\t\t\tCalculating feature and background counts and updating LIMS table")
                    feature_count, exctrl_count, msms_count = calculate_counts_for_lims_table(peakheight_filename,mgf_filename)
                    df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '07 complete'
                    df.loc[i,'num_%s_features'%(polarity_short)] = int(feature_count)
                    df.loc[i,'num_%s_features_in_exctrl'%(polarity_short)] = int(exctrl_count)
                    df.loc[i,'num_%s_msms'%(polarity_short)] = int(msms_count)
                    index_list.append(i)

                    if feature_count > 0:
                        print("\t\t\t\tFiltering features by peak height ratio compared to background (control) signal and writing to file")
                        df_filtered = filter_features_by_background(peakheight_filename)
                        df_filtered_path = os.path.join(pathname, '%s_%s_peak-height-filtered.csv' % (row['parent_dir'], polarity))
                        df_filtered.to_csv(df_filtered_path, index=False)
                        recursive_chown(df_filtered_path, 'metatlas')
                    else:
                        print("\t\t\t\tWarning: Not writing filtered features to file because no features were found")
        
        if len(index_list) > 0:
            index_list = list(set(index_list))
            cols = ['Key',
                    '%s_neg_status'%(tasktype),'%s_pos_status'%(tasktype),
                    'num_pos_features', 'num_neg_features',
                    'num_neg_features_in_exctrl','num_pos_features_in_exctrl',
                    'num_neg_msms','num_pos_msms']
            df = df.loc[df.index.isin(index_list),cols]
            df.replace('NaN', 0, inplace=True)
            df.fillna(0, inplace=True)
            #print("\n\tUpdated untargeted tasks rows:\n" + "\n".join(["\t" + line for line in df.to_string(index=False).splitlines()]))
            update_table_in_lims(df,'untargeted_tasks',method='update')
    print("\tMZmine status update complete.")

def filter_features_by_background(peakheight_filename=None,background_designator="ExCtrl",background_ratio=3):
    """
    Accepts a peakheight file and filters out features that have a max peak height in exp samples that is less than
    3 times the peak height in the control samples.
    """
    data_table = pd.read_csv(peakheight_filename, sep=',')
    features_to_keep = [0] # keep the header
    for i,row in data_table.iterrows():
        if i == 0:
            continue # skip the header
        exctrl_columns = [col for col in row.index if background_designator in col and 'mzML' in col]
        exctrl_max = None
        feature_max = None
        if exctrl_columns:
            exctrl_max = row[exctrl_columns].max()
        feature_columns = [col for col in row.index if background_designator not in col and 'mzML' in col]
        if feature_columns:
            feature_max = row[feature_columns].max()
        if exctrl_max and feature_max:
            if feature_max > exctrl_max*background_ratio:
                features_to_keep.append(i)
    if len(features_to_keep) <= 1:
        print("\t\t\t\t\tWarning: No features passed the background filter. Writing empty dataframe")
        data_table_filtered = data_table.head(1)
    elif len(features_to_keep) > 1:
        data_table_filtered = data_table.loc[features_to_keep]
        difference = data_table.shape[0] - data_table_filtered.shape[0]
        print("\t\t\t\t\t%s features removed by background filter"%(difference))
    return data_table_filtered

def calculate_counts_for_lims_table(peakheight_filename=None, mgf_filename=None, background_designator="ExCtrl"):
    """
    Accepts a peakheight file and mgf file and returns the number of features
    in exp and and control samples found in the peakheightfile plus number of 
    feature ids in mgf file.
    """
    data_table = pd.read_csv(peakheight_filename, sep=',')
    exctrl_count = 0
    feature_count = 0
    msms_count = 0
    
    exctrl_columns = [col for col in data_table.columns if background_designator in col]
    if exctrl_columns:
        exctrl_rows = data_table[data_table[exctrl_columns].gt(0).any(axis=1)]
        exctrl_count = int(exctrl_rows.shape[0])
        print("\t\t\t\t%s features in control samples"%(exctrl_count))
    else:
        print("\t\t\t\tWarning: No ExCtrl samples found in peak height file")

    feature_columns = [col for col in data_table.columns if background_designator not in col and 'mzML' in col]
    if feature_columns:
        feature_rows = data_table[data_table[feature_columns].gt(0).any(axis=1)]
        feature_count = int(feature_rows.shape[0])
        print("\t\t\t\t%s features in experimental samples"%(feature_count))
    if feature_count < 10:
        print("\t\t\t\tWarning: Less than 10 features found in peak height file")

    if os.path.isfile(mgf_filename):
        cmd = 'cat %s | grep FEATURE_ID| wc -l'%mgf_filename
        result = subprocess.run(cmd,stdout=subprocess.PIPE, shell=True)
        if result:
            msms_count = int(result.stdout.strip())
            print("\t\t\t\t%s msms features in mgf file"%(msms_count))

    return feature_count, exctrl_count, msms_count

def update_fbmn_status_in_untargeted_tasks(direct_input=None):
    """
    finds running fbmn tasks
    checks if they have uuid
    checks their status
    changes to complete if yes
    sets spectral hits to initiation
    """
    tasktype='fbmn'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    status_list = ['04 running','13 waiting']
    df = subset_df_by_status(df,tasktype,status_list)
    if df.empty:
        print("\tNo FBMN jobs to update!")
    if not df.empty:
        index_list = []
        print("\tTotal of %s project(s) with FBMN status %s to attempt to update"%(df.shape[0],status_list))
        for i,row in df.iterrows():
            polarity_list = check_for_polarities(row['output_dir'],row['parent_dir'])
            if polarity_list is None:
                print("\t\tWarning! Project %s does not have a negative or a positive polarity directory. Skipping..."%row['parent_dir'])
                continue
            for polarity in polarity_list:
                polarity_short = polarity[:3]
                if row['%s_%s_status'%(tasktype,polarity_short)] == '12 not relevant' or row['%s_%s_status'%(tasktype,polarity_short)] == '07 complete':
                    continue ## This helps when one polarity is complete but the other is not
                pathname = os.path.join(row['output_dir'],'%s_%s'%(row['parent_dir'],polarity))
                fbmn_filename = os.path.join(pathname,'%s_%s_gnps2-fbmn-task.txt'%(row['parent_dir'],polarity))
                if os.path.isfile(fbmn_filename)==True:
                    with open(fbmn_filename,'r') as fid:
                        my_text = fid.read().strip()
                    taskid = my_text.split('=')[-1]
                    status,version = check_gnps2_status(taskid)
                    index_list.append(i)
                    if status=='DONE':
                        if df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] != '07 complete':
                            print("\t\tUpdating FBMN job status for %s in %s mode to complete"%(row['parent_dir'],polarity))
                            df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '07 complete'
                    elif status=='FAILED':
                        if df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] != '09 error':
                            print("\t\tUpdating FBMN job status for %s in %s mode to error"%(row['parent_dir'],polarity))
                            df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '09 error'
                    elif status in ['RUNNING','QUEUED']:
                        if df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] != '04 running':
                            print("\t\tUpdating FBMN job status for %s in %s mode to running"%(row['parent_dir'],polarity))
                            df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '04 running'
                    else:
                        if df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] != '13 waiting':
                            df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '13 waiting'
        if len(index_list) > 0:
            index_list = list(set(index_list))
            cols = ['Key',
                    '%s_neg_status'%(tasktype),
                    '%s_pos_status'%(tasktype)]
            df = df.loc[df.index.isin(index_list),cols]
            df.replace('NaN', 0, inplace=True)
            df.fillna(0, inplace=True)
            #print("\tUpdated untargeted tasks rows:\n" + "\n".join(["\t" + line for line in df.to_string(index=False).splitlines()]))
            update_table_in_lims(df,'untargeted_tasks',method='update')
    print("\tFBMN status update complete.")

def write_mzmine_sbatch_and_runner(basepath,batch_filename,parent_dir,filelist_filename):
    mzmine_launcher = '/global/common/software/m2650/mzmine_parameters/MZmine/MZmine-3.7.2/bin/MZmine -threads auto -t /pscratch/sd/b/bkieft/untargeted/tmp'

    sbatch_filename = '%s_mzmine-sbatch.sbatch'%os.path.join(basepath,parent_dir)
    runner_filename = '%s_mzmine.sh'%os.path.join(basepath,parent_dir)
    s = '%s -batch %s -input %s\necho MZMINE IS DONE\n\n\n'%(mzmine_launcher,batch_filename,filelist_filename)

    with open(sbatch_filename,'w') as fid:
        fid.write('%s\n%s\n'%(SLURM_PERLMUTTER_HEADER.replace('slurm','%s-%s'%(os.path.join(basepath,parent_dir),'mzmine')),s))
    with open(runner_filename,'w') as fid:
        fid.write('sbatch %s'%sbatch_filename)

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
    df = df[pd.notna(df['mzml_file'])]
    df['basename'] = df['mzml_file'].apply(os.path.basename)
    df.sort_values('timeepoch',ascending=False,inplace=True)
    df.drop_duplicates(['basename','parent_dir'],inplace=True)
    
    df.drop(columns=['mzml_file_container'],inplace=True)
    df.replace('',np.nan,inplace=True)
    
    print("\tComparing LIMS untargeted tasks with raw files on disk.")
    #Check that files have been sitting around for at least 3 hours (note time zones may vary)
    time_old = df.groupby('parent_dir')['timeepoch'].max() < (time.time()-3*60*60)
    time_new = df.groupby('parent_dir')['timeepoch'].max() >= (time.time()-3*60*60)
    time_old_folders = time_old[time_old==True].index.tolist()
    time_new_folders = time_new[time_new==True].index.tolist()

    #Check that files have "M2" in the name
    dirs_with_m2_files = df.groupby('parent_dir')['basename'].apply(lambda x: any('_MS2_' in filename for filename in x))
    dirs_with_m2_files = dirs_with_m2_files[dirs_with_m2_files].index.tolist()

    # Print parent_dirs with files less than 3 hours old
    if time_new_folders:
        print("\t\tWarning: There are data directories with files less than 3 hours old!")
        print("\t\tSkipping these projects:", time_new_folders)

    df_untargeted = get_table_from_lims('untargeted_tasks')
    
    all_folders = df.loc[df['polarity'].isin(['POS','NEG']),'parent_dir'].unique()

    all_folders = [a.strip() for a in all_folders]
    if df_untargeted.shape[0]>0:
        folders_in_tasks = df_untargeted['parent_dir']
        folders_in_tasks = [a.strip() for a in folders_in_tasks]
    else:
        folders_in_tasks = []

    new_folders = np.setdiff1d(all_folders,folders_in_tasks)
    print("\tFiltering new projects to only those with MS2 scans and all files older than 3 hours")
    new_folders = list(set(new_folders) & set(time_old_folders) & set(dirs_with_m2_files))
    #print("\tTotal projects: ", str(len(all_folders)))
    #print("\tProjects in untargeted tasks: ", str(len(folders_in_tasks)))
    print("\t\t%s new projects added to untargeted tasks:"%str(len(new_folders)))
    for folder in new_folders:
        print("\t\t\t", folder)

    outdir = '/global/cfs/cdirs/metatlas/projects/untargeted_tasks'
 
    if len(new_folders)>0:
        print('\tCreating metadata and file lists for new project(s)')
        # make the metadata sheets
        file_counts = {}
        metadata_files = {}
        filelist = {}
        for polarity in ['positive','negative']:

            # get file counts
            file_counts[polarity] = df[df['polarity']==polarity[:3].upper()].groupby('parent_dir')['mzml_file'].count()
            missing = np.setdiff1d(new_folders,file_counts[polarity].index.tolist())
            for m in missing:
                file_counts[polarity][m] = 0

            files = df[df['polarity']==polarity[:3].upper()].groupby('parent_dir')
            files = [(d,g) for d,g in files]
            metadata_files[polarity] = {}
            filelist[polarity] = {}
            for block in files:
                if block[0] in new_folders:
                    polarity_short = polarity[:3]
                    parent_dir = '%s_%s'%(block[0],polarity)
                    basepath = os.path.join(outdir,parent_dir)
                    if not os.path.isdir(basepath):
                        os.mkdir(basepath)
                        recursive_chown(basepath, 'metatlas')
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
                    if not os.path.isfile(metadata_filename):
                        temp.to_csv(metadata_filename,sep='\t',index=False)
                    metadata_files[polarity][block[0]] = metadata_filename
                    filelist[polarity][block[0]] = block[1]['mzml_file']
    else:
        print("No new projects need output folders or metadata! Exiting script.")
        return pd.DataFrame()

    if len(new_folders)>0:
        new_folders = pd.DataFrame(data={'parent_dir':new_folders,'num_pos_files':file_counts['positive'][new_folders],'num_neg_files':file_counts['negative'][new_folders]})
        new_folders['pos_metadata_file'] = ''
        new_folders['neg_metadata_file'] = ''
        for i,row in new_folders.iterrows():
            for polarity in ['positive','negative']:
                if row['parent_dir'] in metadata_files[polarity].keys():
                    new_folders.loc[i,"{}_metadata_file".format(polarity[:3])] = metadata_files[polarity][row['parent_dir']]
                    basepath = os.path.join(outdir,'%s_%s'%(row['parent_dir'],polarity))
                    parent_dir = '%s_%s'%(row['parent_dir'],polarity)
                    
                    _, validate_machine_name, _ = field_exists(PurePath(row['parent_dir']), field_num=6)
                    if validate_machine_name in ("IQX", "IDX"):
                        print("\t\tNotice: Project %s in %s mode used IQX or IDX and will use a modified parameters file."%(row['parent_dir'],polarity))
                        mzmine_running_parameters = mzine_batch_params_file_iqx
                    else:
                        mzmine_running_parameters = mzine_batch_params_file
                    params_filename = build_untargeted_filename(outdir,row['parent_dir'],polarity,'batch-params-mzmine')
                    with open(mzmine_running_parameters,'r') as fid:
                        orig_params = fid.read()
                    new_output = os.path.join(os.path.join(outdir,'%s_%s'%(row['parent_dir'],polarity)),'%s_%s'%(row['parent_dir'],polarity))
                    custom_params = orig_params.replace('/Users/bpb/Downloads/mzmine_outputs',new_output)
                    with open(params_filename,'w') as fid:
                        fid.write('%s'%custom_params)
                    filelist_filename = os.path.join(outdir,'%s_%s'%(row['parent_dir'],polarity))
                    filelist_filename = os.path.join(filelist_filename,'%s_%s_filelist.txt'%(row['parent_dir'],polarity))          
                    file_list = [f for f in filelist[polarity][row['parent_dir']].tolist() if f is not None]
                    if len(file_list)>0:
                        with open(filelist_filename,'w') as fid:
                            fid.write('%s'%'\n'.join(file_list))
                        write_mzmine_sbatch_and_runner(basepath,params_filename,parent_dir,filelist_filename)# put file selector here)

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
            # Remove any DNA SIP jobs from untargeted tasks in LIMS
            dna_sip_projects_conditional = (new_folders['parent_dir'].str.contains('DNA-SIP')) & (new_folders['parent_dir'].str.contains('_RM'))
            new_folders.loc[dna_sip_projects_conditional, new_folders.columns.str.contains('status')] = 'not relevant'
            update_table_in_lims(new_folders,'untargeted_tasks',method='insert',max_size=1000)
    else:
        print("No new projects need metadata!")

    return new_folders








#################################################################################
############ Old functions that are not used in the current pipeline ############
#################################################################################


# from oauth2client.service_account import ServiceAccountCredentials
# from metatlas.datastructures import metatlas_objects as metob
# import six
# from collections import defaultdict
# from xml.etree import cElementTree as ET
# from io import StringIO
# from collections import Mapping
# import json
# import zipfile
# import gspread
# import scipy.stats
# from ast import literal_eval
# from rdkit import Chem
# from rdkit.Chem import rdMolDescriptors
# from rdkit.Chem import Descriptors
# from rdkit.ML.Descriptors import MoleculeDescriptors
# from pyteomics import mgf
# import tarfile
# import shutil
# import argparse

# def update_fbmn_status_in_untargeted_tasks_OLD(polarity='positive',polarity_short='pos',latest_version=28.2):
#     """
#     finds running fbmn tasks
#     checks if they have uuid
#     checks their status
#     changes to complete if yes
#     sets spectral hits to initiation
#     """
#     tasktype='fbmn'
#     df = get_table_from_lims('untargeted_tasks')
#     update_df = []
#     c1 = df['%s_%s_status'%(tasktype,polarity_short)]=='04 running'
#     c2 = df['%s_%s_status'%(tasktype,polarity_short)]=='01 initiation'
#     c3 = df['%s_%s_status'%(tasktype,polarity_short)]=='08 hold'
#     for i,row in df[(c1) | (c2) | (c3)].iterrows():
#         pathname = os.path.join(row['output_dir'],'%s_%s'%(row['parent_dir'],polarity))
#         fbmn_filename = os.path.join(pathname,'%s_%s_gnps-uuid.txt'%(row['parent_dir'],polarity))
#         graphml_filename = os.path.join(pathname,'%s_%s_gnps-fbmn-network.graphml'%(row['parent_dir'],polarity))
#         if os.path.isfile(fbmn_filename)==True:
#             #the job was submitted
#             with open(fbmn_filename,'r') as fid:
#                 my_text = fid.read().strip()
#             taskid = my_text.split('=')[-1]
#             status,version = check_gnps_status(taskid)
#             print('%s %.2f for %s'%(status,version,fbmn_filename))
#             if version<latest_version:
#                 df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '01 initiation'
# #                 df.loc[i,'%s_%s_status'%('gnps_msms_hits',polarity_short)] = '01 initiation'
#                 download_gnps_graphml(taskid,graphml_filename)
#                 update_df.append(i)
#             elif status=='DONE':
#                 df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '07 complete'
# #                 df.loc[i,'%s_%s_status'%('gnps_msms_hits',polarity_short)] = '01 initiation'
#                 download_gnps_graphml(taskid,graphml_filename)
#                 update_df.append(i)
#             elif status=='FAILED':
#                 df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '09 error'
# #                 df.loc[i,'%s_%s_status'%('gnps_msms_hits',polarity_short)] = '08 hold'
#                 update_df.append(i)
#             elif status=='RUNNING':
#                 df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '04 running'
# #                 df.loc[i,'%s_%s_status'%('gnps_msms_hits',polarity_short)] = '08 hold'
#                 update_df.append(i)
#             else:
#                 df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '01 initiation'
# #                 df.loc[i,'%s_%s_status'%('gnps_msms_hits',polarity_short)] = '08 hold'
#                 update_df.append(i)
# #         else:
# #             print('%s is not a file'%fbmn_filename)
#     if len(update_df)>0:
#         cols = ['Key',
#                 '%s_%s_status'%(tasktype,polarity_short)] # ,               '%s_%s_status'%('gnps_msms_hits',polarity_short)
#         update_table_in_lims(df.loc[df.index.isin(update_df),cols],'untargeted_tasks',method='update')


# def remove_untargeted_job(experiment,
#                           outdir='/global/cfs/cdirs/metatlas/projects/untargeted_tasks',
#                           zipdir='/global/cfs/cdirs/metatlas/projects/untargeted_outputs'):
    
#     # remove untargeted zip file
#     zip_archive = os.path.join(zipdir,'%s.zip'%experiment)
#     if os.path.isfile(zip_archive):
#         os.remove(zip_archive)
#         print('removed %s'%zip_archive)

#     # remove untargeted output directory
#     for polarity in ['positive','negative']:
#         task_dir = os.path.join(outdir,'%s_%s'%(experiment,polarity))
#         if os.path.isdir(task_dir):
#             try:
#                 shutil.rmtree(task_dir)
#                 print('removed %s'%task_dir)
#             except OSError as e:
#                 print ("Error: %s - %s." % (e.filename, e.strerror))

#     schema = 'lists'
#     # remove untargeted_tasks
#     sql = """SELECT Key, parent_dir FROM untargeted_tasks
#     WHERE parent_dir='%s';"""%(experiment)

#     sql_result = api.query.execute_sql(schema, sql,max_rows=int(1e9))
#     if len(sql_result['rows'])>0:
#         df = pd.DataFrame(sql_result['rows'])
#         df = df[['Key']]
#         update_table_in_lims(df,'untargeted_tasks',method='delete')
#         print('removed tasks')

# def get_google_sheet(notebook_name = "Sheet name",
#                      token='/global/cfs/cdirs/metatlas/projects/google_sheets_auth/ipython to sheets demo-9140f8697062.json',
#                      sheet_name = 'Sheet1',
#                     literal_cols=None):
#     """
#     Returns a pandas data frame from the google sheet.
#     Assumes header row is first row.

#     To use the token hard coded in the token field,
#     the sheet must be shared with:
#     metatlas-ipython-nersc@ipython-to-sheets-demo.iam.gserviceaccount.com
#     Unique sheet names are a requirement of this approach.

#     """
# #     scope = ['https://spreadsheets.google.com/feeds']
# #     scope = ['https://www.googleapis.com/auth/spreadsheets']
#     scope = ['https://spreadsheets.google.com/feeds', 'https://www.googleapis.com/auth/drive']
#     #this is deprecated as of january, but we have pinned the version of oauth2.
#     #see https://github.com/google/oauth2client/issues/401
# #     json_key = json.load(open(token))
# #     credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'].encode(), scope)
#     credentials = ServiceAccountCredentials.from_json_keyfile_name(token, scope)
#     #here is the new way incase the version pin is removed
#     #credentials = ServiceAccountCredentials.from_json_keyfile_name(token, scope)

#     gc = gspread.authorize(credentials)
#     wks = gc.open(notebook_name)
#     istd_qc_data = wks.worksheet(sheet_name).get_all_values()
#     headers = istd_qc_data.pop(0)
#     df = pd.DataFrame(istd_qc_data,columns=headers)

#     # Use round trip through read_csv to infer dtypes
#     s = StringIO()
#     df.to_csv(s)
#     df2 = pd.read_csv(StringIO(s.getvalue()))
#     if 'Unnamed: 0' in df2.columns:
#         df2.drop(columns=['Unnamed: 0'],inplace=True)

#     #turn list elements into lists instead of strings
#     if literal_cols is not None:
#         for col in literal_cols:
#             df2[col] = df2[col].apply(literal_eval)
#     df2 = df2.fillna('')

#     return df2

# def get_parent_folders_from_lcmsruns(get_groups=False):
# #     SELECT DISTINCT parent_dir FROM lcmsrun_plus
#     sql = """SELECT DISTINCT parent_dir FROM lcmsrun_plus"""
#     if get_groups==True:
#         sql = """SELECT 
# lcmsrun.mzml_file.filename AS mzml_file,
# regexp_replace(lcmsrun.name, '.*/([^/]*)/[^/]*$', '\1') AS parent_dir,

# split_part(regexp_replace(lcmsrun.name, '.*/[^/]*/([^/]*)$', '\1'), '_', 13) AS file_name_field_12

# FROM lcmsrun"""
#     schema = 'lists'
#     sql_result = api.query.execute_sql(schema, sql,max_rows=1e6)
#     if sql_result is None:
#         print(('execute_sql: Failed to load results from ' + schema + '.' + table))
#         return None
#     else:
#         df = pd.DataFrame(sql_result['rows'])
#         df = df[[c for c in df.columns if not c.startswith('_')]]
#         return df

# def get_acqtime_from_mzml(mzml_file):
#     startTimeStamp=None
#     with open(mzml_file) as mzml:
#         for line in mzml:
#             if 'startTimeStamp' in line:
#                 startTimeStamp = line.split('startTimeStamp="')[1].split('"')[0].replace('T',' ').rstrip('Z')
#                 break
# #     print startTimeStamp
#     if not '-infinity' in startTimeStamp:
#         date_object = datetime.strptime(startTimeStamp, '%Y-%m-%d %H:%M:%S')
#         utc_timestamp = int(time.mktime(date_object.timetuple()))
#     else:
#         utc_timestamp = int(0)
#     return utc_timestamp

# def update_file_table(file_table):
#     file_type = file_table.split('_')[0]
#     v = GETTER_SPEC[file_type]
#     print(('Getting %s files from disk'%(file_type)))
#     dates,files = get_files_from_disk(PROJECT_DIRECTORY,v['extension'])
#     if len(files)>0:
#         df = pd.DataFrame(data={'filename':files,'file_type':file_type,'timeepoch':dates})
#         df['basename'] = df['filename'].apply(os.path.basename)
#         df['name'] = df['filename'].apply(complex_name_splitter) #make a name for grouping associated content
#     else:
#         df = pd.DataFrame()
#         df['filename'] = 'None'
#         df['file_type'] = file_type
#         df['timeepoch'] = 0
#         df['basename'] = 'None'
#         df['name'] = 'None'
    
#     print(('\tThere were %d files on disk'%len(files)))

#     cols = ['filename','name','Key']
#     df_lims = get_table_from_lims(v['lims_table'],columns=cols)
#     print(('\tThere were %d files from LIMS table %s'%(df_lims.shape[0],v['lims_table'])))

        
#     diff_df = pd.merge(df, df_lims,on=['filename','name'], how='outer', indicator='Exist')
#     diff_df = diff_df.loc[diff_df['Exist'] != 'both'] #(left_only, right_only, or both)
#     print(('\tThere are %d different'%diff_df.shape[0]))
#     print('')
# #         diff_df.fillna('',inplace=True)
#     diff_df['parameters'] = 1

#     cols = ['file_type','filename','timeepoch','basename','name']
#     temp = diff_df.loc[diff_df['Exist']=='left_only',cols]
#     if temp.shape[0]>0:
#         update_table_in_lims(temp,file_table,method='insert')#,index_column='Key',columns=None,labkey_server='metatlas-dev.nersc.gov',project_name='/LIMS'):

#     cols = ['Key','filename']
#     temp = diff_df.loc[diff_df['Exist']=='right_only',cols]
#     temp['Key'] = temp['Key'].astype(int)
#     if temp.shape[0]>0:
#         update_table_in_lims(temp,file_table,method='delete')#,index_column='Key',columns=None,labkey_server='metatlas-dev.nersc.gov',project_name='/LIMS'):
# #         df.to_csv('/global/homes/b/bpb/Downloads/%s_files.tab'%k,index=None,sep='\t')

# def str2bool(v):
#     if isinstance(v, bool):
#         return v
#     if v.lower() in ('yes', 'true', 't', 'y', '1'):
#         return True
#     elif v.lower() in ('no', 'false', 'f', 'n', '0'):
#         return False
#     else:
#         raise argparse.ArgumentTypeError('Boolean value expected.')

# def download_gnps_graphml(taskid,outfile):
#     url = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=gnps_molecular_network_graphml/" % (taskid)
#     positive_graphml = "%s.graphml"%taskid
#     d = requests.get(url)
#     with open(outfile, "w") as fid:
#         fid.write(d.text)

    

# def write_fbmn_sbatch_and_runner(basepath,parent_dir):
#     runner_filename = '%s_fbmn.sh'%os.path.join(basepath,parent_dir)

#     python_binary = '/global/common/software/m2650/python3-metatlas-cori/bin/python'
#     python_file = '/global/homes/b/bpb/repos/metatlas/metatlas/untargeted/send_to_gnps.py'
#     python_args = '--basedir %s --basename %s --override True'%(basepath,parent_dir)
#     with open(runner_filename,'w') as fid:
#         fid.write('%s %s %s\n'%(python_binary,python_file,python_args))



# def read_mgf(filename):
#     df = []
#     with mgf.MGF(filename) as reader:
#         for spectrum in reader:
#     #         count += 1
#             d = spectrum['params']
#             d['spectrum'] = np.array([spectrum['m/z array'],spectrum['intensity array']])
#             d['pepmass'] = d['pepmass'][0]
#             df.append(d)
#     ref_df = pd.DataFrame(df)
#     return ref_df


# def update_num_msms():
#     df_tasks = get_table_from_lims('untargeted_tasks')
#     # Get num features
#     found = 0
#     keep_rows = []
#     for i,row in df_tasks.iterrows():
#         for polarity in ['positive','negative']:
#             filename = build_untargeted_filename(row['output_dir'],row['parent_dir'],polarity,'msms-mzmine')
#             if os.path.isfile(filename):
#                 cmd = 'cat %s | grep FEATURE_ID| wc -l'%filename
#                 result = subprocess.run(cmd,stdout=subprocess.PIPE, shell=True)
#                 num = result.stdout.strip()
#                 if len(num)>0:
#                     num = int(num)
#                     if num != row['num_%s_msms'%polarity[:3]]:
#                         df_tasks.loc[i,'num_%s_msms'%polarity[:3]] = num
#                         keep_rows.append(i)
#                         print(num,row['num_%s_msms'%polarity[:3]],result.stderr,result.stdout)
# #                         print('')

#     cols = [c for c in df_tasks.columns if (c.endswith('_msms')) & (c.startswith('num_'))]
#     cols = cols + ['Key']
#     print(cols)
#     if len(keep_rows)>0:
#         temp = df_tasks.loc[df_tasks.index.isin(keep_rows),cols].copy()
#         temp.fillna(0,inplace=True)
#         update_table_in_lims(temp,'untargeted_tasks',method='update')


# def count_scans_by_class(x):
#     d = {}
#     for i in [0.4,0.5,0.7,0.9]:
#         d['MQScore_gte_%.1f'%i] = len(x.loc[x['MQScore']>=i,'#Scan#'].unique())
#     for i in [3,6,9,12,15]:
#         d['SharedPeaks_gte_%d'%i] = len(x.loc[x['SharedPeaks']>=i,'#Scan#'].unique())
#     return pd.Series(d, index=d.keys())

# def get_gnps_hits(parent_dir,output_dir,polarity,status,override=False):
#     gnps_uuid_file = build_untargeted_filename(output_dir,
#                                  parent_dir,
#                                      polarity,
#                                      'gnps-uuid-fbmn')
#     gnps_zip_output_file = os.path.join(os.path.dirname(gnps_uuid_file),'%s_%s_gnps-download.zip'%(parent_dir,polarity))
#     # if (os.path.isfile(gnps_uuid_file)) & ('complete' in status):
#     #     with open(gnps_uuid_file,'r') as fid:
#     #         url = fid.read()
#     #     gnps_uuid = url.split('=')[-1]
#     if os.path.isfile(gnps_zip_output_file):
#         with zipfile.ZipFile(gnps_zip_output_file, 'r') as archive:
#             hits_file = [f for f in archive.namelist() if 'DB_result' in f]
#             if len(hits_file)>0:
#                 hits_file = hits_file[-1]
#                 with archive.open(hits_file) as hits_fid:
#                     try:
#                         df = pd.read_csv(hits_fid,sep='\t')
#                     except:
#                         df = None

#                 return df
#     return None



# def download_file_from_server_endpoint(server_endpoint,local_file_path):
#     response=requests.post(server_endpoint)
#     if response.status_code==200:
#         # Write the file contents in the response to a file specified by local_file_path
#         with open(local_file_path,'wb') as local_file:
#             for chunk in response.iter_content(chunk_size=128):
#                 local_file.write(chunk)
#     else:
#         print('Can not do this one: %s'%os.path.basename(local_file_path))

# def get_gnps_zipfile(parent_dir,output_dir,polarity,status,override=False):
#     gnps_uuid_file = build_untargeted_filename(output_dir,
#                                  parent_dir,
#                                      polarity,
#                                      'gnps-uuid-fbmn')
# #     print(gnps_uuid_file)
# #     print(polarity)
#     gnps_zip_output_file = os.path.join(os.path.dirname(gnps_uuid_file),'%s_%s_gnps-download.zip'%(parent_dir,polarity))
#     if (os.path.isfile(gnps_uuid_file)) & ('complete' in status):
#         if (override==True) | (not os.path.isfile(gnps_zip_output_file)):
#             with open(gnps_uuid_file,'r') as fid:
#                 url = fid.read()
#             my_uuid = url.split('=')[-1]
#             gnps_url = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task=%s&view=download_cytoscape_data&show=true"%my_uuid
#             print(gnps_url)
#             print(gnps_zip_output_file)
#             download_file_from_server_endpoint(gnps_url, gnps_zip_output_file)
#             print('')
#     else:
#         print('FAIL',parent_dir,status)
# #     <batchstep method="net.sf.mzmine.modules.peaklistmethods.orderpeaklists.OrderPeakListsModule">
# #         <parameter name="Peak lists" type="BATCH_LAST_PEAKLISTS"/>
# #     </batchstep>



# def melt_dataframe(ph,md):
# #     df.drop(columns=['filename','sample'],inplace=True)
#     # df.fillna(0.0,inplace=True)
#     feature_cols = ['feature_id','mz','rt']
#     var_cols = [c for c in ph.columns if c.endswith('mzML')]
#     df = ph.melt(id_vars=feature_cols,value_vars=var_cols)#,var_name='value')
#     # df = df[~pd.isna(df['value'])]
#     df['value'].fillna(0.0,inplace=True)
#     df['value'] = df['value'].astype(float)
#     df = pd.merge(df,md,left_on='variable',right_on='filename')
#     df.drop(columns=['variable'],inplace=True)
#     df.reset_index(inplace=True,drop=True)
#     return df

# def calc_background(df,background='exctrl',background_ratio=3.0):
#     exctrl = df[df['sampletype'].str.contains(background.lower())].groupby(['feature_id','sampletype'])['value'].max().reset_index()
#     sample = df[~df['sampletype'].str.contains(background.lower())].groupby(['feature_id','sampletype'])['value'].max().reset_index()
#     ratio_df = pd.merge(exctrl.add_suffix('_exctrl'),sample.add_suffix('_sample'),left_on='feature_id_exctrl',right_on='feature_id_sample')
#     ratio_df['ratio'] = ratio_df['value_sample']/(1+ratio_df['value_exctrl'])
#     good_features = ratio_df.loc[ratio_df['ratio']>background_ratio,'feature_id_sample'].unique()
#     all_features = ratio_df['feature_id_sample'].unique()
#     bad_features = list(set(all_features) - set(good_features))
#     rm_count = len(all_features) - len(good_features)
#     print('Please remove %d features out of %d'%(rm_count,len(all_features)))
#     if len(bad_features)==0:
#         good_features = df['feature_id'].unique()
#     return good_features,bad_features






# def calc_hit_vector(n,df):
#     """
#     for 0,1,2 n will be 3
#     for 0,1,2,3 n will be 4
#     df is the count from true_positives
#     this function makes it a percentation of hits will sum n or more hits in last element
#     """
#     m = np.zeros((n))
#     s = df['count']/df['count'].sum()
#     nf = df['num_features']
#     for i in s.index:
#         if (i>(len(m)-1)) & (len(m) >1):
#             m_idx = len(m)-1
#         else:
#             m_idx = nf[i]
#         m[m_idx] = m[m_idx] + s[i]
#     return m

# def summarize_results(n,true_pos_filename,base_path,project_path,feature_table_extension,rt_column,mz_column,headerrows,sep,mz_tolerance,rt_tolerance):
#     """
#     """
#     path  = os.path.join(base_path,project_path)
#     feature_file = glob.glob(os.path.join(path,'*%s'%feature_table_extension))
#     if len(feature_file)>0:
#         feature_file = feature_file[-1]
#     else:
#         return np.zeros(n),np.zeros(n), 0

#     new_path = os.path.join(path,'true_pos_results')
#     if not os.path.isdir(new_path):
#         os.mkdir(new_path)

#     new_basename = os.path.basename(feature_file).replace(feature_table_extension,'height.xlsx')
#     output_filename = os.path.join(new_path,new_basename)
#     if os.path.isfile(feature_file): #has the file been made already?
#         with open(feature_file,'r') as fid:
#             s = fid.read()

#         if len(s)>0: #does the file have anything in it?
#             df_experimental = pd.read_csv(feature_file,sep=sep,skiprows=headerrows)
#             if df_experimental.shape[0] > 0: #are there any rows?
#                 if 'hilic' in feature_file.lower():
#                     sheetname = 'HILIC_POS'
#                     df_true_pos = pd.read_excel(true_pos_filename,sheet_name=sheetname)
#                 else:
#                     sheetname = 'CSH_POS'
#                     df_true_pos = pd.read_excel(true_pos_filename,sheet_name=sheetname)
#                 istd_count,bio_count,df_grouped,df_hits,total_count = prepare_true_positive_and_export(output_filename,df_experimental,df_true_pos,rt_column=rt_column,mz_column=mz_column,mz_tolerance=mz_tolerance,rt_tolerance=rt_tolerance)
#                 return calc_hit_vector(n,istd_count), calc_hit_vector(n,bio_count), total_count.loc[0,'total']
#     return np.zeros(n),np.zeros(n), 0

# def make_count_of_knowns(df_hits,df_true_pos):
#     df_grouped = df_hits[['true_pos_index','CompoundName_truepos','experimental_feature_idx']]
#     df_grouped.set_index(['CompoundName_truepos'],inplace=True)

#     df_grouped = df_grouped.groupby(['true_pos_index']).count()
#     df_grouped = pd.merge(df_grouped,df_true_pos,left_index=True,right_index=True,how='outer')
#     df_grouped.rename(columns={'experimental_feature_idx':'num_features'},inplace=True)
#     return df_grouped

# def map_features_to_known(df_experimental,df_true_pos,rt_column='row retention time',mz_column='row m/z',mz_tolerance=0.01,rt_tolerance=0.1):
#     feature_array = df_experimental[[mz_column,rt_column]].values
#     reference_array = df_true_pos[['MZ','RT']].values

#     idx = np.isclose(feature_array[:,None,:], reference_array, rtol=0.0, atol=[mz_tolerance,rt_tolerance]).all(axis=2) #order is m/z, rt, polarity
#     feature_idx, reference_idx = np.where(idx)

#     df_hits = df_true_pos.loc[reference_idx].copy()
#     df_hits['experimental_feature_idx'] = feature_idx
#     df_hits = pd.merge(df_hits,df_true_pos,left_index=True,right_index=True,how='outer',suffixes=['_x','_truepos'])
#     df_hits.drop(columns=['%s_x'%c for c in df_true_pos.columns],inplace=True)
#     df_hits = pd.merge(df_hits,df_experimental,how='left',left_on='experimental_feature_idx',right_index=True,suffixes=['_truepos','_experimental'])
#     df_hits.index.name = 'true_pos_index'
#     df_hits.reset_index(inplace=True)
#     return df_hits

# def summarize_count_per_type(df_grouped,cpd_type='ISTD'):
#     """
#     cpd_type is either 'ISTD' or 'TargetCPD'
#     """
#     cpd_count = df_grouped[df_grouped['Type']==cpd_type][['num_features','MZ']].groupby('num_features').count()
#     cpd_count.reset_index(inplace=True)
#     cpd_count.rename(columns={'MZ':'count'},inplace=True)
#     return cpd_count

# def prepare_true_positive_and_export(output_filename,df_experimental,df_true_pos,rt_column='row retention time',mz_column='row m/z',mz_tolerance=0.01,rt_tolerance=0.1):
#     df_hits = map_features_to_known(df_experimental,df_true_pos,rt_column=rt_column,mz_column=mz_column,mz_tolerance=0.01,rt_tolerance=0.1)
#     df_grouped = make_count_of_knowns(df_hits,df_true_pos)
#     istd_count = summarize_count_per_type(df_grouped,cpd_type='ISTD')
#     bio_count = summarize_count_per_type(df_grouped,cpd_type='TargetCPD')
#     params = pd.DataFrame(columns=['mz_tolerance','rt_tolerance'],data=[[mz_tolerance,rt_tolerance]])
#     total_count = pd.DataFrame(columns=['total'],data=[[df_experimental.shape[0]]])

#     # Create a Pandas Excel writer using XlsxWriter as the engine.
#     writer = pd.ExcelWriter(output_filename, engine='xlsxwriter')

#     # Write each dataframe to a different worksheet.
#     params.to_excel(writer, sheet_name='params')
#     istd_count.to_excel(writer, sheet_name='istd_count')
#     bio_count.to_excel(writer, sheet_name='bio_count')
#     df_grouped.to_excel(writer, sheet_name='df_grouped')
#     df_hits.to_excel(writer, sheet_name='df_hits')
#     total_count.to_excel(writer, sheet_name='total')
#     return istd_count,bio_count,df_grouped,df_hits,total_count
#     # Close the Pandas Excel writer and output the Excel file.
# #     writer.save()



# def filter_features(experiment,rt_cutoff=1,polarity='positive',output_dir='/global/cfs/cdirs/metatlas/projects/untargeted_tasks/'):
#     mgf_file = build_untargeted_filename(output_dir,experiment,polarity,'msms-mzmine')
#     ph_file =  build_untargeted_filename(output_dir,experiment,polarity,'peak-height-mzmine')
#     md_file =  build_untargeted_filename(output_dir,experiment,polarity,'metadata')

#     ph = pd.read_csv(ph_file)
#     print(ph.shape)
#     cols = ph.columns
#     cols = [c.replace(' Peak height','') for c in cols]
#     ph.columns = cols
#     cols = [c for c in cols if not 'Unnamed:' in c]
#     ph = ph[cols]
#     rename_cols = {'row ID':'feature_id','row m/z':'mz','row retention time':'rt'}
#     ph.rename(columns=rename_cols,inplace=True)


#     md = pd.read_csv(md_file,sep='\t')
#     md.rename(columns={'ATTRIBUTE_sampletype':'sampletype'},inplace=True)
#     md.fillna('',inplace=True)
#     md['sampletype'] = md['sampletype'].astype(str)
#     md['sampletype'] = md['sampletype'].apply(lambda x: x.lower())
#     if any([e for e in  md['sampletype'].tolist() if 'exctrl' in e]):
#         md['filename'] = md['filename'].apply(lambda x: os.path.basename(x))

#     df = melt_dataframe(ph,md)
#     df = df[~df['sampletype'].str.startswith('qc')]
#     good_features,bad_features = calc_background(df,background_ratio=50)#,background='media-control')
#     print(len(good_features),len(bad_features))
#     df = df[df['feature_id'].isin(good_features)]
#     idx = df['rt']>rt_cutoff
#     print('there are %d in the solvent front that will be removed'%sum(~idx))
#     df = df[idx]
#     df['polarity'] = polarity
#     # df_range = df.groupby(['feature_id','mz','rt'])['value'].agg({'max':lambda x: x.max(),'range':lambda x: (x.max()+1) / (x.min()+1)})
#     df_avg_temp = df.groupby(['feature_id','polarity','sampletype']).mean()
#     df_avg_temp.reset_index(inplace=True)
#     cols = ['feature_id','polarity','mz','rt']
#     if 'mz' in df_avg_temp.columns:
#         df_range = df_avg_temp.groupby(cols).agg(max=('value', np.max),range=('value', lambda x: (x.max()+1)/(x.min()+1)))

#         df_range = df.groupby(['feature_id','polarity','mz','rt']).agg(max=('value', np.max),range=('value', lambda x: (x.max()+1)/(x.min()+1)))
#         df_range.reset_index(inplace=True,drop=False)

#         df_range = df_range[(df_range['max']>1e6)]
#         # df_range = df_range[(df_range['range']>10) & (df_range['max']>1e6)]
#         df_range.reset_index(inplace=True,drop=True)
#         if df_range.shape[0]>0:
#             df = pd.merge(df,df_range,how='inner',on=cols)
#             df.drop(columns=['max','range'],inplace=True)
#         else:
#             df = None
#     else:
#         df = None
#     return df


# def update_feature_table():
#     df_tasks = get_table_from_lims('untargeted_tasks',columns=['parent_dir','output_dir','mzmine_pos_status','mzmine_neg_status'])
#     df_tasks = df_tasks.melt(id_vars=['parent_dir','output_dir'],value_vars=['mzmine_pos_status','mzmine_neg_status'])
#     df_tasks = df_tasks[df_tasks['value']=='07 complete']
#     df_tasks['polarity'] = df_tasks['variable'].apply(lambda x: 'positive' if x.split('_')[1]=='pos' else 'negative')

#     sql = "SELECT DISTINCT experiment,polarity FROM untargeted_features;"
#     print(sql)
#     schema = 'lists'
#     sql_result = api.query.execute_sql(schema, sql,max_rows=int(1e9))
#     done_experiments = pd.DataFrame(sql_result['rows'])

#     if done_experiments.shape[0]>0:
#         print(df_tasks.shape)
#         df_tasks = pd.merge(df_tasks,done_experiments,left_on=['parent_dir','polarity'],right_on=['experiment','polarity'],how='left',indicator=True)
#         df_tasks = df_tasks[df_tasks['_merge']=='left_only']
#         df_tasks.drop(columns=['experiment','_merge'],inplace=True)
#         print(df_tasks.shape)


#     for i,row in df_tasks.iterrows():
#         experiment = row['parent_dir']
#         polarity = row['polarity']
#         print(experiment)
#         temp = filter_features(experiment,polarity=polarity)
#         if (temp is not None):
#             if ('CONTROL' in temp.columns):
#                 temp.drop(columns=['CONTROL'],inplace=True)
#             temp['experiment'] = experiment
#             temp['polarity'] = polarity
#             # temp.drop(columns=['filename'],inplace=True)
#             # cols = temp.columns
#             # cols = [c for c in cols if c != 'value']
#             # temp = temp.groupby(cols).mean()
#             temp.reset_index(inplace=True)
#             cols = ['feature_id','experiment','polarity','mz','rt']
#             quant_cols = ['value','filename','sampletype']
#             temp = temp[cols + quant_cols]
#             temp.sort_values('value',ascending=False,inplace=True)
#             temp.drop_duplicates(subset=cols,inplace=True)
#         else:
#             print(experiment)
#             temp = pd.DataFrame()
#             temp['experiment'] = [experiment]
#             temp['polarity'] = [polarity]
#         update_table_in_lims(temp,'untargeted_features',method='insert')
        


        
# def update_gnps_hits_table(output_dir='/global/cfs/cdirs/metatlas/projects/untargeted_tasks/'):
#     print('Updating GNPS hits')
#     df_hits_db = get_table_from_lims('gnps_hits',columns=['feature_key'])
#     df_features = get_table_from_lims('untargeted_features',columns=['Key','experiment','polarity','feature_id'],max_rows=1e9)
#     cols = ['experiment','polarity']
#     all_experiments = df_features[cols].drop_duplicates(subset=cols)

#     if df_hits_db.shape[0]>0:
#         done_experiments = df_features.loc[df_features['Key'].isin(df_hits_db['feature_key']),cols].drop_duplicates()
#         todo_experiments = pd.merge(all_experiments,done_experiments,how='left',indicator=True)
#         todo_experiments = todo_experiments[todo_experiments['_merge']=='left_only']
#         todo_experiments.drop(columns=['_merge'],inplace=True)
#     else:
#         todo_experiments = all_experiments.copy()
#     print('Getting GNPS Hits for',todo_experiments.shape[0],'experiments')

#     for i,row in todo_experiments.iterrows():
#         print(row)
#         experiment = row['experiment']
#         polarity = row['polarity']

#         hits = get_gnps_hits(experiment,output_dir,
#                                  polarity,'complete',override=True)

#         if hits is not None:
#             hits = hits[pd.notna(hits['Smiles'])]
#             hits = hits[pd.notna(hits['InChIKey'])]
#             no_hits=False # set this to true if there are not any hits to update a dummy row to the table

#             if hits.shape[0]>0:
#                 hits.rename(columns={'#Scan#':'feature_id','MQScore':'score'},inplace=True)
#                 hits.columns = [c.lower() for c in hits.columns]
#                 hits.sort_values('score',ascending=False,inplace=True)
#                 hits.drop_duplicates(subset=['feature_id','inchikey'],inplace=True)

#                 ##### Add any new compounds to the compounds table #####
#                 compounds = get_table_from_lims('gnps_compounds',columns=['inchikey'])
#                 cols_compound_db = ['inchikey','compound_name','smiles']
#                 df_compound = hits[cols_compound_db].copy()
#                 df_compound.drop_duplicates(subset=['inchikey'],inplace=True)
#                 df_compound = pd.merge(df_compound,compounds,on='inchikey',how='left',indicator=True)
#                 df_compound = df_compound[df_compound['_merge']=='left_only']
#                 if df_compound.shape[0]>0:
#                     print('updating compounds')
#                     update_table_in_lims(df_compound,'gnps_compounds',method='insert')

#                 ##### Add hits to the hits table #####
#                 hits_cols = ['feature_id','score','adduct','sharedpeaks', 'massdiff','inchikey']
#                 hits_slice = df_features[(df_features['experiment']==experiment) & (df_features['polarity']==polarity)].copy()
#                 hits_slice = pd.merge(hits_slice,hits[hits_cols],on='feature_id',how='inner')
#                 if hits_slice.shape[0]>0:
#                     hits_slice.rename(columns={'Key':'feature_key'},inplace=True)
#                     hits_slice.drop(columns=['feature_id','experiment','polarity'],inplace=True)
#                     update_table_in_lims(hits_slice,'gnps_hits',method='insert')
#                 else:
#                     no_hits = True
#             else:
#                 no_hits = True
#         else:
#             no_hits = True

#         if no_hits==True:
#             # This is a placeholder entry in the table so that jobs with no hits will be removed from future queue
#             # Add one row that has a feature
#             print('no hits:',experiment,polarity)
#             hits_slice = df_features[(df_features['experiment']==experiment) & (df_features['polarity']==polarity)].copy()
#             min_feature_key = hits_slice['Key'].min()
#             hits_slice = hits_slice[hits_slice['Key']==min_feature_key]
#             hits_slice.rename(columns={'Key':'feature_key'},inplace=True)
#             hits_slice.drop(columns=['feature_id','experiment','polarity'],inplace=True)
#             update_table_in_lims(hits_slice,'gnps_hits',method='insert')

# def update_gnps_compound_properties():
#     sql = "SELECT DISTINCT inchikey FROM gnps_compound_properties;"

#     print('NEW BATCH of compound properties to calculate')
#     schema = 'lists'
#     sql_result = api.query.execute_sql(schema, sql,max_rows=int(1e9))
#     done_compounds = pd.DataFrame(sql_result['rows'])

#     # Updating Compounds That need properties
#     df = get_table_from_lims('gnps_compounds',columns=['inchikey','smiles'])
#     print(df.shape)
#     if done_compounds.shape[0]>0:
#         df = df[~df['inchikey'].isin(done_compounds['inchikey'])]
#     print(df.shape)


#     idx = (pd.notna(df['smiles'])) & (pd.notna(df['inchikey'])) & (df['smiles']!='')
#     df = df[idx]
#     df.loc[idx,'smiles'] = df.loc[idx,'smiles'].apply(lambda x: x.strip())
#     df.loc[idx,'mol'] = df.loc[idx,'smiles'].apply(lambda x: Chem.MolFromSmiles(x))
#     idx = pd.notna(df['mol'])
#     df = df[idx]

#     print(df.shape)
#     if df.shape[0]>0:
#         descriptor_names = list(rdMolDescriptors.Properties.GetAvailableProperties())
#         get_descriptors = rdMolDescriptors.Properties(descriptor_names)

#         def mol_to_descriptors(mol):
#             descriptors = get_descriptors.ComputeProperties(mol)
#             return descriptors


#         des_list = [x[0] for x in Descriptors._descList]
#         calculator = MoleculeDescriptors.MolecularDescriptorCalculator(des_list)


#         descriptors = []

#         for i,row in df.iterrows():
#             # print(row['inchikey'],row['smiles'])
#             out = {'inchikey':row['inchikey']}
#             d = list(mol_to_descriptors(row['mol']))
#             for j,dd in enumerate(d):
#                 out['property: %s'%descriptor_names[j]] = dd
#             d = list(calculator.CalcDescriptors(row['mol']))
#             for j,dd in enumerate(d):
#                 out['descriptor: %s'%des_list[j]] = dd
#             descriptors.append(out)
#             # print('done\n')

#         descriptors = pd.DataFrame(descriptors)

#         cols = [c for c in descriptors.columns if not 'inchikey' in c]
#         descriptors = descriptors.melt(id_vars='inchikey',value_vars=cols)
#         descriptors['type'] = descriptors['variable'].apply(lambda x: 'property' if 'property' in x else 'descriptor')
#         descriptors['variable'] = descriptors['variable'].apply(lambda x: x.split(':')[-1].strip())
#         update_table_in_lims(descriptors,'gnps_compound_properties',method='insert',max_size=5000)

# def mean_confidence_interval(data, confidence=0.95):
#     data = data['value'].values
#     a = 1.0 * np.array(data)
#     n = len(a)
#     m = np.mean(a)
#     if n>1:
#         se = scipy.stats.sem(a)
#         h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
#         return pd.Series({'mean':m,'lower_bound': m-h,'upper_bound':m+h,'num_replicates':n,'standard_error':se})
#     else:
#         return pd.Series({'mean':m,'lower_bound': 0,'upper_bound':0,'num_replicates':n,'standard_error':0})


# def update_untargeted_treatments():
#     print('getting master list of features')
#     sql = """SELECT Key as Key,feature_id,experiment,polarity FROM untargeted_features"""

#     schema = 'lists'
#     sql_result = api.query.execute_sql(schema, sql,max_rows=int(1e9))
#     df = pd.DataFrame(sql_result['rows'])
#     sql = "SELECT DISTINCT feature FROM untargeted_treatments;"
#     df = df[pd.notna(df['feature_id'])]

#     print('there are %d features'%df.shape[0])

#     print('getting treatment groups')
#     sql_result = api.query.execute_sql(schema, sql,max_rows=int(1e9))
#     done_features = pd.DataFrame(sql_result['rows'])
#     print('there are %d features by treatments'%done_features.shape[0])

#     if done_features.shape[0]>0:
#         df = df[~df['Key'].isin(done_features['feature'])]
#         print('After removing done features, there are %d features'%df.shape[0])

#     if df.shape[0]>0:
#         todo_experiments = df.copy()
#         todo_experiments.drop_duplicates(subset=['experiment','polarity'],inplace=True)
#         for i,row in todo_experiments.iterrows():
#             print('working on %s %s mode'%(row['experiment'],row['polarity']))
#             new_features = filter_features(row['experiment'],rt_cutoff=1,polarity=row['polarity'])
#             if (new_features is not None) and (new_features.shape[0]>0):
#                 new_features['polarity'] = row['polarity']
#                 new_features['experiment'] = row['experiment']
#                 new_features = pd.merge(new_features,df,on=['feature_id','experiment','polarity'],how='inner')
#                 new_features.rename(columns={'Key':'feature'},inplace=True)
#                 new_features.drop(columns=['feature_id','experiment','polarity'],inplace=True)
#                 print('it has %d features by treatments'%new_features.shape[0])

#                 g = new_features.groupby(['feature','sampletype'])[['value','mz']].apply(lambda x: mean_confidence_interval(x))
#                 g.reset_index(inplace=True,drop=False)
#                 update_table_in_lims(g,'untargeted_treatments',method='insert',max_size=500)


            
    
            
# #####################################################
# #####################################################
# ########         mzmine setup scripts        ########
# #####################################################
# #####################################################




# def remove_duplicate_files(files):
#     file_names = []
#     unique_files = []
#     for f in files:
#         if not f.name in file_names:
#             unique_files.append(f.mzml_file)
#             file_names.append(f.name)
#     return unique_files

# def get_files(groups,filename_substring,file_filters,keep_strings,is_group=False,return_mzml=True):
#     """
#     if is_group is False, gets files from the experiment/folder name and filters with file_filters

#     if is_group is True, gets files from the metatlas group name and filters with file_filters

#     """

#     for i,g in enumerate(groups):
#         if is_group == True:
#         # get files as a metatlas group
#             groups = dp.select_groups_for_analysis(name = g,do_print=False,
#                                                    most_recent = True,
#                                                    remove_empty = True,
#                                                    include_list = [], exclude_list = file_filters)#['QC','Blank'])
#             new_files = []
#             for each_g in groups:
#                 for f in each_g.items:
#                     new_files.append(f)
#         else:
#             new_files = metob.retrieve('Lcmsruns',experiment=g,name=filename_substring,username='*')
#         if i == 0:
#             all_files = new_files
#         else:
#             all_files.extend(new_files)
#         if len(new_files) == 0:
#             print('##### %s has ZERO files!'%g)

#     # only keep files that don't have substrings in list
#     if len(file_filters) > 0:
#         for i,ff in enumerate(file_filters):
#             if i == 0:
#                 files = [f for f in all_files if not ff in f.name]
#             else:
#                 files = [f for f in files if not ff in f.name]
#     else:
#         files = all_files

#     # kick out any files that don't match atleast one of the keep_strings
#     keep_this = []
#     filter_used = [] #good to keep track if a filter isn't used.  likely a typo
#     if len(keep_strings) > 0:
#         for i,ff in enumerate(files):
#             keep_this.append(any([True if f in ff.name else False for f in keep_strings]))
#         for i,ff in enumerate(keep_strings):
#             filter_used.append(any([True if ff in f.name else False for f in files]))
#         if not all(filter_used):
#             for i,f in enumerate(filter_used):
#                 if f==False:
#                     print('%s keep string is not used'%keep_strings[i])

#         files = [files[i] for i,j in enumerate(keep_this) if j==True]

#     files = remove_duplicate_files(files)
#     return files


# def make_targeted_mzmine_job(basedir,basename,polarity,files):
#     if not os.path.exists(basedir):
#         os.mkdir(basedir)

#     xml_str = get_targeted_batch_file_template()
#     d = xml_to_dict(xml_str)

#     task = metob.MZMineTask()
#     task.polarity = polarity
#     task.lcmsruns = files
#     new_d = replace_files(d,files)

#     project_name = '%s_%s'%(basename,task.polarity)
#     task.output_workspace = os.path.join(basedir,project_name,'%s_%s.mzmine'%(basename,task.polarity))
#     task.input_xml = os.path.join(basedir,'logs','%s_%s_filtered.xml'%(basename,task.polarity))

#     task.mzmine_launcher = get_latest_mzmine_binary()

#     # new_d = configure_crop_filter(new_d,task.polarity,files)
#     # new_d = configure_targeted_peak_detection(new_d,peak_list_filename,intensity_tolerance=1e-4,noise_level=1e4,mz_tolerance=20,rt_tolerance=0.5)
#     new_d = configure_workspace_output(new_d,task.output_workspace)

#     t = dict_to_etree(new_d)
#     indent_tree(t)
#     xml_batch_str = tree_to_xml(t,filename=task.input_xml)
#     job_runner = '%s %s'%(task.mzmine_launcher,task.input_xml)
#     return job_runner

# def configure_targeted_peak_detection(new_d,peak_list_filename,intensity_tolerance=1e-4,noise_level=1e4,mz_tolerance=20,rt_tolerance=0.5):
#     """
#     Name suffix: Suffix to be added to the peak list name.

#     Peak list file: Path of the csv file containing the list of peaks to be detected. The csv file should have three columns.
#     The first column should contain the expected M/Z, the second column the expected RT and the third the peak name. Each peak should be in a different row.

#     Field separator: Character(s) used to separate fields in the peak list file.

#     Ignore first line: Check to ignore the first line of peak list file.

#     Intensity tolerance: This value sets the maximum allowed deviation from expected shape of a peak in chromatographic direction.

#     Noise level: The minimum intensity level for a data point to be considered part of a chromatogram. All data points below this intensity level are ignored.

#     MZ Tolerance: Maximum allowed m/z difference to find the peak

#     RT tolerance: Maximum allowed retention time difference to find the peak
#     """
#     # Set the noise floor
#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'TargetedPeakDetectionModule' in d['@method']][0]

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Peak list file' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%s'%peak_list_filename

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Intensity tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.6f'%(intensity_tolerance)

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Noise level' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.6f'%(noise_level)

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.6f'%(rt_tolerance)

#     return new_d

# def configure_crop_filter(new_d,polarity,files,min_rt=0.01,max_rt=100,fps_string='FPS'):
#     """

#     """
#     # identify the element for this change
#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'CropFilterModule' in d['@method']][0]
#     # Set the filter string
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Raw data files' in d['@name']][0]
#     if any([fps_string in f for f in files]):
#         new_d['batch']['batchstep'][idx]['parameter'][idx2]['name_pattern'] = '*FPS*'
#     else:
#         new_d['batch']['batchstep'][idx]['parameter'][idx2]['name_pattern'] = '*'

#     # Set the polarity
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Scans' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['polarity'] = polarity.upper()

#     #set the rt min and rt max use the same idx2 as polarity
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['retention_time'] = {'max':'%.4f'%max_rt,'min':'%.4f'%min_rt}

#     # new_d['batch']['batchstep'][idx]['parameter'][idx2]['ms_level'] = '1-2'

#     return new_d

# def configure_mass_detection(new_d,ms1_noise_level=1e4,ms2_noise_level=1e2):
#     """

#     """
#     # Find the module
#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'MassDetectionModule' in d['@method']]
#     #The first idx will be for MS1 and the second will be for MS2

#     # Set the MS1 attributes
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[0]]['parameter']) if 'Mass detector' in d['@name']][0]
#     idx3 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[0]]['parameter'][idx2]['module']) if 'Centroid' in d['@name']][0]
#     new_d['batch']['batchstep'][idx[0]]['parameter'][idx2]['module'][idx3]['parameter']['#text'] = '%.2f'%(ms1_noise_level)

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[0]]['parameter']) if 'Scans' in d['@name']][0]
#     new_d['batch']['batchstep'][idx[0]]['parameter'][idx2]['ms_level'] = '1'

#     # Set the MS2 attributes
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[1]]['parameter']) if 'Mass detector' in d['@name']][0]
#     idx3 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[1]]['parameter'][idx2]['module']) if 'Centroid' in d['@name']][0]
#     new_d['batch']['batchstep'][idx[1]]['parameter'][idx2]['module'][idx3]['parameter']['#text'] = '%.2f'%(ms2_noise_level)

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[1]]['parameter']) if 'Scans' in d['@name']][0]
#     new_d['batch']['batchstep'][idx[1]]['parameter'][idx2]['ms_level'] = '2'


#     return new_d

# def configure_smoothing(new_d,smoothing_scans):
#     """
# #        <batchstep method="net.sf.mzmine.modules.peaklistmethods.peakpicking.smoothing.SmoothingModule">
# #         <parameter name="Peak lists" type="BATCH_LAST_PEAKLISTS"/>
# #         <parameter name="Filename suffix">smoothed</parameter>
# #         <parameter name="Filter width">9</parameter>
# #         <parameter name="Remove original peak list">false</parameter>
# #     </batchstep>
#     """
#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'SmoothingModule' in d['@method']][0]
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Filter width' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(smoothing_scans)

#     return new_d


# def configure_chromatogram_builder(new_d,min_num_scans,group_intensity_threshold,min_peak_height,mz_tolerance):
#     """
#     new_d = configure_chromatogram_builder(new_d,task.min_num_scans,task.group_intensity_threshold,task.min_peak_height,task.mz_tolerance)

# #     <batchstep method="net.sf.mzmine.modules.masslistmethods.ADAPchromatogrambuilder.ADAPChromatogramBuilderModule">
# #         <parameter name="Raw data files" type="ALL_FILES"/>
# #         <parameter name="Scans">
# #             <ms_level>1</ms_level>
# #         </parameter>
# #         <parameter name="Mass list">masses</parameter>
# #         <parameter name="Min group size in # of scans">5</parameter>
# #         <parameter name="Group intensity threshold">1000000.0</parameter>
# #         <parameter name="Min highest intensity">80000.0</parameter>
# #         <parameter name="m/z tolerance">
# #             <absolutetolerance>0.002</absolutetolerance>
# #             <ppmtolerance>7.0</ppmtolerance>
# #         </parameter>
# #         <parameter name="Suffix">chromatograms</parameter>
# #     </batchstep>


#     """
#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'ADAPChromatogramBuilderModule' in d['@method']][0]
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Min group size' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(min_num_scans)

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Group intensity threshold' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(group_intensity_threshold)

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Min highest intensity' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(min_peak_height)

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)

#     return new_d

# def configure_adap_peak_deconvolution(new_d,min_peak_height,minimum_relative_height,search_for_minimum_rt_range,chromatographic_threshold,min_sn_ratio,min_peak_duration,max_peak_duration):
#     """
#     <parameter name="Algorithm" selected="Wavelets (ADAP)">

#     <module name="Wavelets (ADAP)">
#                 <parameter name="S/N threshold">3.0</parameter>
#                 <parameter name="S/N estimator" selected="Intensity window SN">
#                     <module name="Intensity window SN"/>
#                     <module name="Wavelet Coeff. SN">
#                         <parameter name="Peak width mult.">1.0</parameter>
#                         <parameter name="abs(wavelet coeffs.)">true</parameter>
#                     </module>
#                 </parameter>
#                 <parameter name="min feature height">4500.0</parameter>
#                 <parameter name="coefficient/area threshold">60.0</parameter>
#                 <parameter name="Peak duration range">
#                     <min>0.0</min>
#                     <max>0.5</max>
#                 </parameter>
#                 <parameter name="RT wavelet range">
#                     <min>0.0</min>
#                     <max>0.1</max>
#                 </parameter>
#             </module>

#     """
#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'DeconvolutionModule' in d['@method']][0]
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Algorithm' in d['@name']][0]
#     idx3 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module']) if 'Local minimum search' in d['@name']][0]
#     idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Chromatographic threshold' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%chromatographic_threshold
#     idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Search minimum in RT range (min)' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%search_for_minimum_rt_range
#     idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Minimum relative height' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%minimum_relative_height
#     idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Minimum absolute height' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%min_peak_height
#     idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Min ratio of peak top/edge' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%min_sn_ratio
#     idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Peak duration range (min)' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['min'] = '%.3f'%min_peak_duration
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['max'] = '%.3f'%max_peak_duration
#     return new_d

# def configure_lms_peak_deconvolution(new_d,min_peak_height,minimum_relative_height,search_for_minimum_rt_range,chromatographic_threshold,min_sn_ratio,min_peak_duration,max_peak_duration):
#     """
#     <parameter name="Algorithm" selected="Local minimum search">
#                 <module name="Local minimum search">
#                 <parameter name="Chromatographic threshold">0.75</parameter>
#                 <parameter name="Search minimum in RT range (min)">0.02</parameter>
#                 <parameter name="Minimum relative height">0.002</parameter>
#                 <parameter name="Minimum absolute height">90000.0</parameter>
#                 <parameter name="Min ratio of peak top/edge">1.03</parameter>
#                 <parameter name="Peak duration range (min)">
#                     <min>0.03</min>
#                     <max>1.0</max>
#                 </parameter>
#             </module>

#     """
#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'DeconvolutionModule' in d['@method']][0]
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Algorithm' in d['@name']][0]
#     idx3 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module']) if 'Local minimum search' in d['@name']][0]
#     idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Chromatographic threshold' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%chromatographic_threshold
#     idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Search minimum in RT range (min)' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%search_for_minimum_rt_range
#     idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Minimum relative height' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%minimum_relative_height
#     idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Minimum absolute height' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%min_peak_height
#     idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Min ratio of peak top/edge' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%min_sn_ratio
#     idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Peak duration range (min)' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['min'] = '%.3f'%min_peak_duration
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['max'] = '%.3f'%max_peak_duration
#     return new_d

# def configure_isotope_search(new_d,mz_tolerance,rt_tol_perfile,representative_isotope,remove_isotopes,polarity):
#     """

#     """
#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'Isotope' in d['@method']][0]
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_perfile)

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Representative isotope' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%s'%(representative_isotope)

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Remove original peaklist' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%s'%(str(remove_isotopes).lower())

#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'Adduct' in d['@method']][0]
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'RT tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_perfile)
#     if polarity == 'negative':
#         idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Adducts' in d['@name']][0]
#         #the default is setup for positive mode adducts.
#         #only change them if you are in negative mode
#         for i,a in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['adduct']):
#             if a['@selected'] == 'true':
#                 new_d['batch']['batchstep'][idx]['parameter'][idx2]['adduct'][i]['@selected'] = 'false'
#             else:
#                 new_d['batch']['batchstep'][idx]['parameter'][idx2]['adduct'][i]['@selected'] = 'true'

#     return new_d

# def configure_join_aligner(new_d,mz_tolerance,rt_tol_multifile):
#     """
#     # Join aligner has these scores:
# #             <parameter name="Minimum absolute intensity">3000.0</parameter>
# #             <parameter name="Minimum score">0.6</parameter>

#     """
#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'JoinAlignerModule' in d['@method']][0]
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_multifile)

# #     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Minimum absolute intensity' in d['@name']][0]
# #     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = 3000#'%.3f'%(mz_tolerance)

# #     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Minimum score' in d['@name']][0]
# #     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = 0.6#'%.3f'%(rt_tol_multifile)

#     return new_d

# def configure_rows_filter(new_d,min_peaks_in_row,peak_with_msms):
#     """

#     """
#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'RowsFilterModule' in d['@method']][0]
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Minimum peaks in a row' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%d'%min_peaks_in_row
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Minimum peaks in an isotope pattern' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%d'%min_peaks_in_row
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Keep only peaks with MS2 scan (GNPS)' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%d'%peak_with_msms
#     return new_d

# def configure_duplicate_filter(new_d,mz_tolerance,rt_tol_perfile):
#     """

#     """
#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'DuplicateFilterModule' in d['@method']][0]
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'RT tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_perfile)
#     return new_d

# def configure_gap_filling(new_d,mz_tolerance,gapfill_intensity_tolerance,rt_tol_multifile):
#     """
#     #     <batchstep method="net.sf.mzmine.modules.peaklistmethods.gapfilling.peakfinder.multithreaded.MultiThreadPeakFinderModule">
# #         <parameter name="Peak lists" type="BATCH_LAST_PEAKLISTS"/>
# #         <parameter name="Name suffix">gap-filled</parameter>
# #         <parameter name="Intensity tolerance">0.05</parameter>
# #         <parameter name="m/z tolerance">
# #             <absolutetolerance>0.001</absolutetolerance>
# #             <ppmtolerance>5.0</ppmtolerance>
# #         </parameter>
# #         <parameter name="Retention time tolerance" type="absolute">0.03</parameter>
# #         <parameter name="Remove original peak list">false</parameter>
# #     </batchstep>

#     """
#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'gapfilling.peakfinder' in d['@method']][0]

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Intensity tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(gapfill_intensity_tolerance)

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_multifile)

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)


#     return new_d

# def configure_output(new_d,output_csv_height,output_csv_area,output_workspace,output_mgf):
#     """

#     """
#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'CSVExportModule' in d['@method']]
#     #the first will be height the second will be area

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[0]]['parameter']) if 'Filename' in d['@name']][0]
#     new_d['batch']['batchstep'][idx[0]]['parameter'][idx2]['#text'] = output_csv_height

#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx[1]]['parameter']) if 'Filename' in d['@name']][0]
#     new_d['batch']['batchstep'][idx[1]]['parameter'][idx2]['#text'] = output_csv_area

#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'GNPSExportModule' in d['@method']][0]
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Filename' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = output_mgf

#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'ProjectSaveAsModule' in d['@method']][0]
#     new_d['batch']['batchstep'][idx]['parameter']['#text'] = output_workspace
#     return new_d

# def configure_csv_output(new_d,output_csv):
#     """

#     """
#     idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'CSVExportModule' in d['@method']][0]
#     idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Filename' in d['@name']][0]
#     new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = output_csv
#     return new_d

# def indent_tree(elem, level=0):
#     i = "\n" + level*"  "
#     if len(elem):
#         if not elem.text or not elem.text.strip():
#             elem.text = i + "  "
#         if not elem.tail or not elem.tail.strip():
#             elem.tail = i
#         for elem in elem:
#             indent_tree(elem, level+1)
#         if not elem.tail or not elem.tail.strip():
#             elem.tail = i
#     else:
#         if level and (not elem.tail or not elem.tail.strip()):
#             elem.tail = i

# def get_targeted_batch_file_template(loc='do_not_change_batch_file_targeted_peak_list.xml'):
#     """
#     return string of text from the template batch file
#     """
#     with open(os.path.join(BATCH_FILE_PATH,loc),'r') as fid:
#         file_text = fid.read()
#     return file_text

# def get_batch_file_template(loc='bootcamp_adap_template.xml'):
#     """
#     return string of text from the template batch file
#     """
#     with open(os.path.join(BATCH_FILE_PATH,loc),'r') as fid:
#         file_text = fid.read()
#     return file_text

# def get_latest_mzmine_binary(system='Cori',version='most_recent'):
#     """
#     Returns the path to the mzmine launch script.
#     Default is most recent.  Alternatively specify the folder containng version you want

#     for example:
#         version='MZmine-2.23'
#     will use the launch script in that folder

#     wget $(curl -s https://api.github.com/repos/mzmine/mzmine2/releases/v2.33 | grep 'browser_' | cut -d\" -f4) -O mzmine_latest.zip


#     # To setup the most recent mzmine binary, follow these steps
#     cd /project/projectdirs/metatlas/projects/mzmine_parameters/MZmine
#     wget $(curl -s https://api.github.com/repos/mzmine/mzmine2/releases/latest | grep 'browser_' | cut -d\" -f4) -O mzmine_latest.zip
#     unzip mzmine_latest.zip
#     # change directories into latest mzmine download
#     # cd MZmine-XXXX
#     cp ../MZmine-2.24/startMZmine_NERSC_* .
#     cd /project/projectdirs/metatlas/projects/
#     chgrp -R metatlas mzmine_parameters
#     chmod -R 770 mzmine_parameters
#     """
#     mzmine_versions = glob.glob(os.path.join(BINARY_PATH,'*' + os.path.sep))
#     if version == 'most_recent':
#         most_recent = sorted([os.path.basename(m) for m in mzmine_versions if 'MZmine-' in m])[-1]
#     else:
#         most_recent = [m.split(os.path.sep)[-2] for m in mzmine_versions if version in m][-1]
#     launch_script = os.path.join(os.path.join(BINARY_PATH,most_recent),'startMZmine_NERSC_Headless_%s.sh'%system)
#     if os.path.isfile(launch_script):
#         return launch_script
#     else:
#         print('See the docstring, the launch script seems to be missing.')

# def replace_files(d,file_list):
#     """
#     Replace files for mzmine task

#     Inputs:
#     d: an xml derived dictionary of batch commands
#     file_list: a list of full paths to mzML files

#     Outputs:
#     d: an xml derived dict with new files in it
#     """
#     for i,step in enumerate(d['batch']['batchstep']):
#         if 'RawDataImportModule' in step['@method']:
#             d['batch']['batchstep'][i]['parameter']['file'] = file_list
#     return d





# def tree_to_xml(t,filename=None):
#     """

#     """
#     xml_str = ET.tostring(t)
#     if filename:
#         with open(filename,'w') as fid:
#             fid.write(xml_str)
#     return xml_str



# def dict_to_etree(d):
#     """
#     Convert a python dictionary to an xml str
#     http://stackoverflow.com/questions/7684333/converting-xml-to-dictionary-using-elementtree

#     Example:
#     from collections import defaultdict
#     from xml.etree import cElementTree as ET

#     try:
#         basestring
#     except NameError:  # python3
#         basestring = str

#     #d is a python dictionary
#     ET.tostring(dict_to_etree(d))

#     """
#     def _to_etree(d, root):
#         print(type(d),d)
#         print('\n\n\n')
#         if not d:
#             pass

#         if type(d) is {}.values().__class__:
#             d = list(d.values)

#         if isinstance(d, str):
#             root.text = d
#         elif isinstance(d, dict):
#             for k,v in d.items():
#                 assert isinstance(k, str)
#                 if k.startswith('#'):
#                     assert k == '#text' and isinstance(v, str)
#                     root.text = v
#                 elif k.startswith('@'):
#                     assert isinstance(v, str)
#                     root.set(k[1:], v)
#                 elif isinstance(v, list):
#                     for e in v:
#                         _to_etree(e, ET.SubElement(root, k))
#                 else:
#                     _to_etree(v, ET.SubElement(root, k))
# #         elif isinstance(d,dict_values):
# #             d = [d]
# #             _to_etree(d,ET.SubElement(root, k))
#         else: assert d == 'invalid type', (type(d), d)
#     assert isinstance(d, dict) and len(d) == 1
#     tag, body = next(iter(d.items()))
#     node = ET.Element(tag)
#     _to_etree(body, node)
#     return node

# def xml_to_dict(xml_str):
#     """
#     Convert an xml file into a python dictionary.
#     http://stackoverflow.com/questions/7684333/converting-xml-to-dictionary-using-elementtree

#     Example:
#     from xml.etree import cElementTree as ET
#     filename = '/global/homes/b/bpb/batch_params/xmlfile.xml'
#     with open(filename,'r') as fid:
#         xml_str = fid.read()

#     d = xml_to_dict(xml_str)
#     """
#     t = ET.XML(xml_str)
#     d = etree_to_dict(t)
#     return d

# def etree_to_dict(t):
#     """
#     Convert an xml tree into a python dictionary.
#     http://stackoverflow.com/questions/7684333/converting-xml-to-dictionary-using-elementtree

#     """

#     d = {t.tag: {} if t.attrib else None}
#     children = list(t)
#     if children:
#         dd = defaultdict(list)
#         for dc in map(etree_to_dict, children):
#             for k, v in six.iteritems(dc):
#                 dd[k].append(v)
#         d = {t.tag: {k:v[0] if len(v) == 1 else v for k, v in six.iteritems(dd)}}
#     if t.attrib:
#         d[t.tag].update(('@' + k, v) for k, v in six.iteritems(t.attrib))
#     if t.text:
#         text = t.text.strip()
#         if children or t.attrib:
#             if text:
#                 d[t.tag]['#text'] = text
#         else:
#             d[t.tag] = text
#     return d



# ##########################################################
# #### From Here ###########################################
# ### https://github.com/ianlini/flatten-dict ##############
# ##########################################################
# ##########################################################



# def tuple_reducer(k1, k2):
#     if k1 is None:
#         return (k2,)
#     else:
#         return k1 + (k2,)


# def path_reducer(k1, k2):
#     import os.path
#     if k1 is None:
#         return k2
#     else:
#         return os.path.join(k1, k2)



# def tuple_splitter(flat_key):
#     return flat_key


# def path_splitter(flat_key):
#     keys = PurePath(flat_key).parts
#     return keys

# REDUCER_DICT = {
#     'tuple': tuple_reducer,
#     'path': path_reducer,
# }

# SPLITTER_DICT = {
#     'tuple': tuple_splitter,
#     'path': path_splitter,
# }


# def flatten(d, reducer='tuple', inverse=False, enumerate_types=()):
#     """Flatten `Mapping` object.

#     Parameters
#     ----------
#     d : dict-like object
#         The dict that will be flattened.
#     reducer : {'tuple', 'path', Callable}
#         The key joining method. If a `Callable` is given, the `Callable` will be
#         used to reduce.
#         'tuple': The resulting key will be tuple of the original keys.
#         'path': Use `os.path.join` to join keys.
#     inverse : bool
#         Whether you want invert the resulting key and value.
#     enumerate_types : Sequence[type]
#         Flatten these types using `enumerate`.
#         For example, if we set `enumerate_types` to ``(list,)``,
#         `list` indices become keys: ``{'a': ['b', 'c']}`` -> ``{('a', 0): 'b', ('a', 1): 'c'}``.

#     Returns
#     -------
#     flat_dict : dict
#     """
#     enumerate_types = tuple(enumerate_types)
#     flattenable_types = (Mapping,) + enumerate_types
#     if not isinstance(d, flattenable_types):
#         raise ValueError("argument type %s is not in the flattenalbe types %s"
#                          % (type(d), flattenable_types))

#     if isinstance(reducer, str):
#         reducer = REDUCER_DICT[reducer]
#     flat_dict = {}

#     def _flatten(d, parent=None):
#         key_value_iterable = enumerate(d) if isinstance(d, enumerate_types) else six.viewitems(d)
#         for key, value in key_value_iterable:
#             flat_key = reducer(parent, key)
#             if isinstance(value, flattenable_types):
#                 _flatten(value, flat_key)
#             else:
#                 if inverse:
#                     flat_key, value = value, flat_key
#                 if flat_key in flat_dict:
#                     raise ValueError("duplicated key '{}'".format(flat_key))
#                 flat_dict[flat_key] = value

#     _flatten(d)
#     return flat_dict


# # def nested_set_dict(d, keys, value):
# #     """Set a value to a sequence of nested keys

# #     Parameters
# #     ----------
# #     d : Mapping
# #     keys : Sequence[str]
# #     value : Any
# #     """
# #     assert keys
# #     key = keys[0]
# #     if len(keys) == 1:
# #         if key in d:
# #             raise ValueError("duplicated key '{}'".format(key))
# #         d[key] = value
# #         return
# #     d = d.setdefault(key, {})
# #     nested_set_dict(d, keys[1:], value)


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
#         if type(d) == list:
#             d.append(value)
#         else:
#             d[key] = value
#         return

#     # the type is a string so make a dict if none exists
#     if type(keys[1]) == int:
#         if key in d:
#             pass
#         else:
#             d[key] = []
#         d = d[key]
#     elif type(key)==int:
#         if (key+1) > len(d):
#             d.append({})
#         d = d[key]
#     else:
#         d = d.setdefault(key, {})
#     nested_set_dict(d, keys[1:], value)


# def unflatten(d, splitter='tuple', inverse=False):
#     """Unflatten dict-like object.

#     Parameters
#     ----------
#     d : dict-like object
#         The dict that will be unflattened.
#     splitter : {'tuple', 'path', Callable}
#         The key splitting method. If a Callable is given, the Callable will be
#         used to split.
#         'tuple': Use each element in the tuple key as the key of the unflattened dict.
#         'path': Use `pathlib.Path.parts` to split keys.

#     Tester
#     d1 = {'a':{'b':[{'c1':'nested1!','d1':[{'e1':'so_nested1!!!'}]},
#                {'c2':'nested2!','d2':[{'e2':'so_nested2!!!'}]},
#                {'c3':'nested3!','d3':[{'e3':'so_nested3!!!'}]},
#                {'c4':'nested4!','d4':[{'e4':'so_nested4a!!!'},
#                                       {'e4':'so_nested4b!!!'},
#                                       {'e4':'so_nested4c!!!'},
#                                       {'e4':'so_nested4d!!!'},
#                                       {'e4':'so_nested4e!!!'}]}]}}

#     Returns
#     -------
#     unflattened_dict : dict
#     """
#     if isinstance(splitter, str):
#         splitter = SPLITTER_DICT[splitter]

#     kv = sorted([(k,v) for (k,v) in d.items()])
#     unflattened_dict = {}
#     for kkvv in kv:
#         key_tuple = kkvv[0]
#         value = kkvv[1]
#         nested_set_dict(unflattened_dict, key_tuple, value)

#     return unflattened_dict
