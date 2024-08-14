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
import metatlas.tools.validate_filenames as vfn
import subprocess
import grp
import logging

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

def call_logger(log_filename: str, log_level: str, log_format: str):
    logging.basicConfig(filename=log_filename, level=log_level, format=log_format, filemode='a')

def start_script(script=None):
    """
    Kick off script with a timestamp for log
    """
    if script is not None:
        logging.info(tab_print('Successfully started %s on %s:'%(script,datetime.now()), 0))

def end_script(script=None):
    """
    End script with a timestamp for log
    """
    if script is not None:
        logging.info(tab_print('Successfully ended %s on %s:'%(script,datetime.now()), 0))

def tab_print(message="",indent_level=0):
    """
    Print a message with a specified number of tabs
    """
    return ("\t"*indent_level + message)

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
    Subset a dataframe by status in LIMS
    """
    if inverse == True:
        df_subset = df[~df['%s_pos_status'%tasktype].isin(status) | ~df['%s_neg_status'%tasktype].isin(status)]
    if inverse == False:
        df_subset = df[df['%s_pos_status'%tasktype].isin(status) | df['%s_neg_status'%tasktype].isin(status)]
    return df_subset

def filter_common_bad_project_names(df: pd.DataFrame):
    """
    Takes a df object (usually from exporting the LIMS untargeted tasks table) and removes
    projects with common strings that are in untargeted tasks which should not be run
    """
    df = df[~(df['parent_dir'].str.contains(' '))]
    df = df[~(df['parent_dir'].str.contains('&'))]
    #df = df[~(df['parent_dir'].str.contains('Niyogi'))]
    df = df[~(df['parent_dir'].str.contains('_partial'))]
    df = df[~(df['parent_dir'].str.contains('_old'))]
    return df

def write_fbmn_tasks_to_file(task_list: dict, output_dir='/global/cfs/cdirs/metatlas/projects/untargeted_tasks'):
    """
    Takes a list of dictionaries from submit_fbmn_jobs
    and writes the fbmn task id to a file in the 
    project directory at untarageted_tasks on perlmutter
    """
    experiment = task_list['experiment']
    polarity = task_list['polarity']
    task = task_list['response']['task']
    filename = os.path.join(output_dir, '%s_%s'%(experiment, polarity), '%s_%s_gnps2-fbmn-task.txt'%(experiment, polarity))
    if task:
        with open(filename,'w') as fid:
            fid.write("%s_%s=%s\n"%(experiment,polarity,task))
            final_filename = os.path.basename(filename)
            logging.info(tab_print("GNPS2 task file for %s mode written to %s"%(polarity,final_filename), 4))
    else:
        logging.warning(tab_print("Warning! GNPS2 FBMN task ID not found. File gnps2-fbmn-task.txt not written.", 4))

def zip_and_upload_untargeted_results(download_folder = '/global/cfs/cdirs/metatlas/projects/untargeted_outputs', \
                                      output_dir = '/global/cfs/cdirs/metatlas/projects/untargeted_tasks', \
                                      upload=True, overwrite_zip=False, overwrite_drive=False, direct_input=None, \
                                      min_features_admissible=0, add_documentation=True, \
                                      doc_name=None):
    """
    This function is called by export_untargeted_results.py
    
    Specify the download folder (where the zip files will be saved) and
    the output folder (where task directories are located) in case the info
    in the LIMS table is not accurate
    
    Upload is True (default) when the zip folder will be uploaded to Google Drive

    Overwrite_zip is False (default) when the zip folder will not be created if it already exists in download_folder

    Overwrite_drive is False (default) when the zipped folder will not be uploaded to Google Drive if it exists there already

    Direct_input is None (default) when all available projects will undergo zip and upload. Set direct_input
    to a csv list of project names if you only want to run this function on specific untargeted_tasks

    Min_features_admissible is 0 (default) if you want the function to zip and upload only when there are more than 0 features
    """
    if add_documentation == True and doc_name is None:
        logging.warning(tab_print("Warning! 'Add_documentation' flag is True but no documentation file name provided. Exiting script...", 1))
        sys.exit(1)
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
        overwrite_zip = True
        overwrite_drive = True
    status_list = ['07 complete']
    df = subset_df_by_status(df,'fbmn',status_list)
    df = subset_df_by_status(df,'mzmine',status_list)
    df = subset_df_by_status(df,'fbmn',['09 error'],inverse=True)
    # Don't want to subset by 'not relevant' here because it initially looks for one or the other polarity to have the given status
    if df.empty:
        logging.info(tab_print("No completed untargeted results to zip and upload!", 1))
        return
    if not df.empty:
        logging.info(tab_print("%s project(s) found with complete mzmine and fbmn status."%(df.shape[0]), 1))
        for i,row in df.iterrows():
            output_zip_archive = os.path.join(download_folder,'%s.zip'%row['parent_dir'])
            polarity_list = check_for_polarities(output_dir,row['parent_dir'])
            if overwrite_zip==False and os.path.exists(output_zip_archive):
                #logging.warning(tab_print("Warning! Zip archive for %s exists. Set overwrite_zip to True if you want to replace. Skipping zip and upload..."%(row['parent_dir']), 1))
                continue
            if overwrite_zip==True or not os.path.exists(output_zip_archive):
                if polarity_list is None:
                    logging.warning(tab_print("Warning! Project %s does not have a negative or a positive polarity directory. Skipping..."%(row['parent_dir']), 1))
                    continue
                neg_mzmine_file = os.path.join(output_dir, '%s_negative'%row['parent_dir'], '%s_negative_peak-height.csv'%row['parent_dir'])
                pos_mzmine_file = os.path.join(output_dir, '%s_positive'%row['parent_dir'], '%s_positive_peak-height.csv'%row['parent_dir'])
                neg_fbmn_file = os.path.join(output_dir, '%s_negative'%row['parent_dir'], '%s_negative_gnps2-fbmn-library-results.tsv'%row['parent_dir'])
                pos_fbmn_file = os.path.join(output_dir, '%s_positive'%row['parent_dir'], '%s_positive_gnps2-fbmn-library-results.tsv'%row['parent_dir'])
                if 'negative' in polarity_list and 'positive' in polarity_list:
                    # Check that mzmine and fbmn "marker" files exist, they've probably finished sucessfully (double-check since status need to all be 'complete')
                    if os.path.exists(neg_mzmine_file) and os.path.exists(pos_mzmine_file) and os.path.exists(neg_fbmn_file) and os.path.exists(pos_fbmn_file):
                        neg_feature_counts = check_peak_height_table(neg_mzmine_file)
                        pos_feature_counts = check_peak_height_table(pos_mzmine_file)
                        if (neg_feature_counts + pos_feature_counts) > min_features_admissible:
                            # if overwrite_zip==True and os.path.exists(output_zip_archive):
                            #     wait = 10
                            #     logging.warning(tab_print("Warning! Overwrite zip is True and %s exists. Giving %s seconds to end process before zipping over existing archive..."%(os.path.basename(output_zip_archive),wait), 1))
                            #     time.sleep(wait)
                            neg_directory = os.path.join(output_dir, '%s_%s'%(row['parent_dir'], 'negative'))
                            pos_directory = os.path.join(output_dir, '%s_%s'%(row['parent_dir'], 'positive'))
                            if add_documentation == True:
                                logging.info(tab_print("Downloading latest GNPS2 user guide documentation to add to zip...", 1))
                                doc_present = add_gnps2_documentation(download_folder=download_folder,doc_name=doc_name)
                                if doc_present:
                                    cmd = 'zip -rjq - %s %s %s >%s'%(neg_directory,pos_directory,os.path.join(download_folder,doc_name),output_zip_archive)
                                else:
                                    logging.warning(tab_print("Warning! Add documentation flag is True but GNPS2 user guide documentation not available. Not adding to zip...", 2))
                                    cmd = 'zip -rjq - %s %s >%s'%(neg_directory,pos_directory,output_zip_archive)
                            else:
                                cmd = 'zip -rjq - %s %s >%s'%(neg_directory,pos_directory,output_zip_archive)
                            os.system(cmd)
                            logging.info(tab_print("New untargeted results in %s mode(s) zipped for %s"%(polarity_list,row['parent_dir']), 1))
                            if upload == True and os.path.exists(output_zip_archive):
                                upload_to_google_drive(output_zip_archive,overwrite_drive)
                        else:
                            logging.warning(tab_print("Warning! Project %s has less than %s features in the %s peak height table(s). Skipping zip and upload..."%(row['parent_dir'],min_features_admissible,polarity_list), 1))
                            continue
                elif not 'negative' in polarity_list and 'positive' in polarity_list:
                    if os.path.exists(pos_mzmine_file) and os.path.exists(pos_fbmn_file):
                        pos_feature_counts = check_peak_height_table(pos_mzmine_file)
                        if pos_feature_counts > min_features_admissible:
                            # if overwrite_zip==True and os.path.exists(output_zip_archive):
                            #     wait = 10
                            #     logging.warning(tab_print("Warning! Overwrite zip is True and %s exists. Giving %s seconds to end process before zipping over existing archive..."%(os.path.basename(output_zip_archive),wait), 1))
                            #     time.sleep(wait)
                            pos_directory = os.path.join(output_dir, '%s_%s'%(row['parent_dir'], 'positive'))
                            if add_documentation == True:
                                logging.info(tab_print("Downloading latest GNPS2 user guide documentation to add to zip...", 1))
                                doc_present = add_gnps2_documentation(download_folder=download_folder,doc_name=doc_name)
                                if doc_present:
                                    cmd = 'zip -rjq - %s %s >%s'%(pos_directory,os.path.join(download_folder,doc_name),output_zip_archive)
                                else:
                                    logging.warning(tab_print("Warning! Add documentation flag is True but GNPS2 user guide documentation not available. Not adding to zip...", 2))
                                    cmd = 'zip -rjq - %s >%s'%(pos_directory,output_zip_archive)
                            else:
                                cmd = 'zip -rjq - %s >%s'%(pos_directory,output_zip_archive)
                            os.system(cmd)
                            logging.info(tab_print("New untargeted results in %s mode(s) zipped for %s"%(polarity_list,row['parent_dir']), 1))
                            if upload == True and os.path.exists(output_zip_archive):
                                upload_to_google_drive(output_zip_archive,overwrite_drive)
                        else:
                            logging.warning(tab_print("Warning! Project %s has less than %s features in the %s peak height table(s). Skipping zip and upload..."%(row['parent_dir'],min_features_admissible,polarity_list), 1))
                            continue
                elif 'negative' in polarity_list and not 'positive' in polarity_list:
                    if os.path.exists(neg_mzmine_file) and os.path.exists(neg_fbmn_file):
                        neg_feature_counts = check_peak_height_table(neg_mzmine_file)
                        if neg_feature_counts > min_features_admissible:
                            # if overwrite_zip==True and os.path.exists(output_zip_archive):
                            #     wait = 10
                            #     logging.warning(tab_print("Warning! Overwrite zip is True and %s exists. Giving %s seconds to end process before zipping over existing archive..."%(os.path.basename(output_zip_archive),wait), 1))
                            #     time.sleep(wait)
                            neg_directory = os.path.join(output_dir, '%s_%s'%(row['parent_dir'], 'negative'))
                            if add_documentation == True:
                                logging.info(tab_print("Downloading latest GNPS2 user guide documentation to add to zip...", 1))
                                doc_present = add_gnps2_documentation(download_folder=download_folder,doc_name=doc_name)
                                if doc_present:
                                    cmd = 'zip -rjq - %s %s >%s'%(neg_directory,os.path.join(download_folder,doc_name),output_zip_archive)
                                else:
                                    logging.warning(tab_print("Warning! Add documentation flag is True but GNPS2 user guide documentation not available. Not adding to zip...", 2))
                                    cmd = 'zip -rjq - %s >%s'%(neg_directory,output_zip_archive)
                            else:
                                cmd = 'zip -rjq - %s >%s'%(neg_directory,output_zip_archive)
                            os.system(cmd)
                            logging.info(tab_print("New untargeted results in %s mode(s) zipped for %s"%(polarity_list,row['parent_dir']), 1))
                            if upload == True and os.path.exists(output_zip_archive):
                                upload_to_google_drive(output_zip_archive,overwrite_drive)
                        else:
                            logging.warning(tab_print("Warning! Project %s has less than %s features in the %s peak height table(s). Skipping zip and upload..."%(row['parent_dir'],min_features_admissible,polarity_list), 1))
                            continue

def check_peak_height_table(peak_height_file):
    """
    Open the peak_height.csv file created by MZmine
    and count the number of lines (minus the header)
    """
    cmd = 'cat %s | wc -l'%peak_height_file
    result = subprocess.run(cmd,stdout=subprocess.PIPE, shell=True)
    if result:
        feature_counts = (int(result.stdout.strip())-1) # Subtract the header line
        return feature_counts
    else:
        return 0

def add_gnps2_documentation(download_folder: str, doc_name="Untargeted_metabolomics_GNPS2_Guide.docx"):
    cmd = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone copy --update ben_lbl_gdrive:/untargeted_outputs/{doc_name} {download_folder}'
    subprocess.check_output(cmd, shell=True)
    check_download = f'ls -l {download_folder} | grep "{doc_name}"'
    check_download_out = subprocess.check_output(check_download, shell=True)
    if check_download_out.decode('utf-8').strip():
        logging.info(tab_print("GNPS2 user guide documentation (%s) downloaded."%(doc_name), 2))
        return True
    else:
        logging.warning(tab_print("Warning! GNPS2 user guide documentation not downloaded!", 2))
        return False
    
def upload_to_google_drive(file_path: str, overwrite=False):
    """
    Upload a file (usually a zip) to the metabolomics untargeted outputs google
    drive folder using rclone 
    Check credentials in the rclone config file using /global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone config file
    or /global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone listremotes
    """
    project_folder = os.path.basename(file_path)
    if overwrite == True:
        logging.info(tab_print("Uploading zip to Google Drive...", 2))
        upload_command = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone copy --update {file_path} ben_lbl_gdrive:/untargeted_outputs'
        try:
            subprocess.check_output(upload_command, shell=True)
        except:
            logging.warning(tab_print("Warning! Google Drive upload failed with overwrite=%s with exception on %s"%(overwrite, upload_command), 3))
            return
        check_upload_command = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone ls ben_lbl_gdrive:/untargeted_outputs | grep "{project_folder}"'
        try:
            check_upload_out = subprocess.check_output(check_upload_command, shell=True)
            if check_upload_out.decode('utf-8').strip():
                logging.info(tab_print("Google Drive upload confirmed!", 3))
            else:
                logging.warning(tab_print("Warning! Google Drive upload check failed. Upload may not have been successful.", 3))
        except:
            logging.warning(tab_print("Warning! Google Drive upload failed with overwrite=%s with exception on %s"%(overwrite,check_upload_command), 3))
            return
    if overwrite == False:
        check_upload_command = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone ls ben_lbl_gdrive:/untargeted_outputs | grep "{project_folder}"'
        try:
            check_upload_out = subprocess.check_output(check_upload_command, shell=True)
            if check_upload_out.decode('utf-8').strip():
                logging.info(tab_print("Overwrite is False and Google Drive folder for %s already exists. Skipping upload..."%(project_folder), 2))
                return
            else:
                logging.info(tab_print("Uploading zip to Google Drive...", 2))
                upload_command = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone copy --ignore-existing {file_path} ben_lbl_gdrive:/untargeted_outputs'
                try:
                    subprocess.check_output(upload_command, shell=True)
                    check_upload = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone ls ben_lbl_gdrive:/untargeted_outputs | grep "{project_folder}"'
                    check_upload_out = subprocess.check_output(check_upload, shell=True)
                    if check_upload_out.decode('utf-8').strip():
                        logging.info(tab_print("Google Drive upload confirmed!", 3))
                    else:
                        logging.warning(tab_print("Warning! Google Drive upload check failed. Upload may not have been successful.", 3))
                except:
                    logging.warning(tab_print("Warning! Google Drive upload failed with overwrite=%s with exception on %s"%(overwrite,upload_command), 3))
                    return
        except:
            logging.warning(tab_print("Warning! Google Drive upload failed with overwrite=%s with exception on %s"%(overwrite,check_upload_command), 3))
            return

def submit_quickstart_fbmn(params="",username=""):
    """
    Submit FBMN jobs by passing parameters
    """
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
    """
    This function is called by check_untargeted_status.py

    Print the status of a user-defined list of projects
    by calling this function with a csv list of project names
    """
    status_df_list = []
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    if df.empty:
        logging.info(tab_print("Could not find untargeted tasks to report! Check that the input argument is a single comma-separated list", 1))
        return
    if not df.empty:
        logging.info(tab_print("Generating status output table for %s projects:"%(df.shape[0]), 1))
        for i,row in df.iterrows():
            try:
                archive = row['parent_dir'] + ".zip"
                check_upload = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone ls ben_lbl_gdrive:/untargeted_outputs | grep "{archive}"'
                check_upload_out = subprocess.check_output(check_upload, shell=True)
                if check_upload_out.decode('utf-8').strip():
                    row['gdrive_status'] = '08 uploaded'
                else:
                    row['gdrive_status'] = '14 not uploaded'
            except:
                row['gdrive_status'] = '14 not uploaded'

            status_df_list.append(row[['parent_dir', 'mzmine_pos_status', 'mzmine_neg_status', 'fbmn_pos_status', 'fbmn_neg_status', 'gdrive_status']])
        
        total_status_df = pd.DataFrame(status_df_list)
        
        print("\n",total_status_df.to_string())

def recursive_chown(basepath, group_name):
    """
    Given a path (to a file or directory) this function will recursively (if directory)
    change the group ownership to the specified group name
    """
    gid = grp.getgrnam(group_name).gr_gid
    for root, dirs, files in os.walk(basepath):
        for dir_name in dirs:
            os.chown(os.path.join(root, dir_name), -1, gid)
        for file_name in files:
            os.chown(os.path.join(root, file_name), -1, gid)
    os.chown(basepath, -1, gid)  # Also change the ownership of the base directory itself

def download_from_url(url, target_path):
    """
    This function takes a URL and a target path (on disk) and downloads the contents
    of the URL with requests.get() and writes the contents to the target path
    """
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(target_path, 'wb') as f:
            f.write(response.raw.read())
        return True
    else:
        return False

def download_fbmn_results(output_dir='/global/cfs/cdirs/metatlas/projects/untargeted_tasks',overwrite_fbmn=False,direct_input=None):
    """
    This function is called by download_fbmn_results.py

    finds complete fbmn tasks (which also have complete mzmine status)
    downloads the graphml and results table files
    renames and moves results to the untargeted_tasks folder

    Overwrite is False (default) when existing GNPS2 files in the untargeted_tasks folder will not be replaced

    Direct_input is None (default) when files from GNPS2 for all available projects will be downloaded. Set direct_input
    to a csv list of project names if you only want to run this function on specific untargeted_tasks
    """
    tasktype='fbmn'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
        overwrite_fbmn = True
    status_list = ['07 complete','09 error']
    df = subset_df_by_status(df,tasktype,status_list)
    df = subset_df_by_status(df,'mzmine',['07 complete']) # Also want to check that mzmine is complete before downloading fbmn
    if df.empty:
        logging.info(tab_print("No completed FBMN data to download!", 1))
    if not df.empty:
        for i,row in df.iterrows():
            polarity_list = check_for_polarities(output_dir,row['parent_dir'])
            if polarity_list is None:
                logging.warning(tab_print("Warning! Project %s does not have a negative or a positive polarity directory. Skipping..."%(row['parent_dir']), 1))
                continue
            for polarity in polarity_list:
                polarity_short = polarity[:3]
                if row['%s_%s_status'%(tasktype,polarity_short)] == '12 not relevant':
                    continue
                if row['%s_%s_status'%(tasktype,polarity_short)] == '09 error':
                    logging.warning(tab_print("Warning! FBMN task for %s %s has error status. Not downloading files."%(row['parent_dir'],polarity), 1))
                    continue
                pathname = os.path.join(output_dir,'%s_%s'%(row['parent_dir'],polarity))
                fbmn_filename = os.path.join(pathname,'%s_%s_gnps2-fbmn-task.txt'%(row['parent_dir'],polarity))
                if os.path.isfile(fbmn_filename)==True:
                    with open(fbmn_filename,'r') as fid:
                        taskid = fid.read().split('=')[1].strip()
                    graphml_filename = os.path.join(pathname,'%s_%s_gnps2-fbmn-network.graphml'%(row['parent_dir'],polarity))
                    results_table_filename = os.path.join(pathname,'%s_%s_gnps2-fbmn-library-results.tsv'%(row['parent_dir'],polarity))
                    gnps2_link_filename = os.path.join(pathname,'%s_%s_gnps2-page-link.txt'%(row['parent_dir'],polarity))
                    # if overwrite_fbmn==True and os.path.exists(graphml_filename) and os.path.exists(results_table_filename) and os.path.exists(gnps2_link_filename):
                    #     wait = 10
                    #     logging.warning(tab_print("Warning! Overwrite is True and all FBMN files exist on disk. Giving %s seconds to end process before downloading..."%(wait), 1))
                    #     time.sleep(wait)
                    if overwrite_fbmn==False and os.path.exists(graphml_filename) and os.path.exists(results_table_filename) and os.path.exists(gnps2_link_filename):
                        continue
                    logging.info(tab_print("Downloading FBMN results for %s with task ID %s"%(row['parent_dir'],taskid), 1))
                    if overwrite_fbmn == True or not os.path.exists(gnps2_link_filename):
                        with open(gnps2_link_filename,'w') as fid:
                            fid.write(f"https://gnps2.org/status?task={taskid}\n")
                            logging.info(tab_print("Wrote GNPS2 link", 2))
                    if overwrite_fbmn == True or not os.path.exists(graphml_filename):
                        graphml_download = download_from_url("https://gnps2.org/result?task=%s&viewname=graphml&resultdisplay_type=task"%(taskid), graphml_filename)
                        if graphml_download:
                            logging.info(tab_print("Downloaded graphml file", 2))
                        else:
                            logging.info(tab_print("Error: Failed to download graphml file", 2))
                    if overwrite_fbmn == True or not os.path.exists(results_table_filename):
                        results_table_download = download_from_url("https://gnps2.org/resultfile?task=%s&file=nf_output/library/merged_results_with_gnps.tsv"%(taskid), results_table_filename)
                        if results_table_download:
                            logging.info(tab_print("Downloaded result table file", 2))
                        else:
                            logging.info(tab_print("Error: Failed to download results table", 2))
                    if os.path.exists(graphml_filename):
                        recursive_chown(graphml_filename, 'metatlas')
                    if os.path.exists(results_table_filename):
                        recursive_chown(results_table_filename, 'metatlas')
                    if os.path.exists(gnps2_link_filename):
                        recursive_chown(gnps2_link_filename, 'metatlas')
        logging.info(tab_print("All new FBMN results downloaded.", 1))

def get_recent_mgf_files(output_dir = '/global/cfs/cdirs/metatlas/projects/untargeted_tasks', time_back=30):
    """
    This finds all recent mgf files in all folders in output_dir
    and passes them to the remove_contaimant_from_mgf function
    """
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
    """
    This accepts a file path to an mgf file and removes the contaminant
    peak(s) from the file
    """
    if 'negative' in f:
        target = 202.07880
        diff = target * 5 / 1e6
    elif 'positive' in f:
        target = 202.07770
        diff = target * 5 / 1e6
    else:
        logging.info(tab_print("Error: Could not determine polarity of file %s"%(f), 1))
        return
    
    # Open the file for reading
    with open(f, 'r') as file:
        # Read all lines from the file
        lines = file.readlines()

    # Filter the lines that contain a number within the range
    filtered_lines = [line for line in lines if not any(abs(float(num)-target)<diff for num in line.split() if num.replace('.', '', 1).isdigit())]
    
    if len(filtered_lines)<len(lines):
        removed_lines = len(lines)-len(filtered_lines)
        file_path = os.path.basename(f)
        logging.info(tab_print("Deleting %s lines from %s"%(removed_lines,file_path), 1))
        # Open the file for writing
        with open(f, 'w') as file:
            # Write the lines back
            file.writelines(filtered_lines)

def remove_mgf_contaminants(output_dir = '/global/cfs/cdirs/metatlas/projects/untargeted_tasks', time_back=30):
    """
    This script is called in run_fbmn.py before submitting to GNPS2

    This wrapper script will call the two relevant functions
    """
    recent_files = get_recent_mgf_files(output_dir=output_dir,time_back=time_back)
    for f in recent_files:
        if f:
            remove_contaminant_from_mgf(f)

def get_table_from_lims(table,columns=None,max_rows=1e6):
    """
    Downloads table from LIMS and returns a pandas dataframe.
    Table should be the name of a valid "list" in LIMS
    """
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
    """
    This function takes a alphanumeric task ID from a GNPS2 run
    and returns the status and version of the run from a static
    HTML page (accessible to anyone)
    """
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

def submit_fbmn_jobs(direct_input=None,overwrite_fbmn=False,validate_names=True,output_dir='/global/cfs/cdirs/metatlas/projects/untargeted_tasks'):
    """
    This function is called by run_fbmn.py

    finds waiting or errored fbmn tasks (which also have complete mzmine status)
    submits the tasks to GNPS2
    updates the LIMS table with status running if successful

    Direct_input is None (default) when all waiting or errorer tasks will be considered for submission. Set direct_input
    to a csv list of project names if you only want to run this function on specific untargeted_tasks
    """
    tasktype = 'fbmn'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
        overwrite_fbmn = True
    status_list = ['13 waiting','09 error']
    df = subset_df_by_status(df,tasktype,status_list)
    df = subset_df_by_status(df,'mzmine',['07 complete']) # Also want to check that mzmine is complete before submitting fbmn
    if df.empty:
        logging.info(tab_print("No new FBMN jobs to submit!", 1))
        return
    if df.shape[0] > 20:
        logging.info(tab_print('There are too many new projects to be submitted (%s), please check if this is accurate. Exiting script.'%(df.shape[0]), 1))
        return
    if not df.empty:
        index_list = []
        logging.info(tab_print("Total of %s projects(s) with FBMN status %s and MZmine status 'complete' to submit to GNPS2"%(df.shape[0],status_list), 1))
        for i,row in df.iterrows():
            project = row['parent_dir']
            logging.info(tab_print("Working on project %s:"%(project), 2))
            polarity_list = check_for_polarities(row['output_dir'],row['parent_dir'])
            if polarity_list is None:
                logging.warning(tab_print("Warning! Project %s does not have a negative or a positive polarity directory. Skipping..."%(row['parent_dir']), 3))
                continue
            for polarity in polarity_list:
                polarity_short = polarity[:3]
                if row['%s_%s_status'%(tasktype,polarity_short)] == '12 not relevant' or row['%s_%s_status'%('mzmine',polarity_short)] != '07 complete':
                    logging.warning(tab_print("Warning! It looks like MZmine for %s mode for %s hasn't completed yet. Skipping..."%(project,polarity), 3))
                    continue
                if row['%s_%s_status'%(tasktype,polarity_short)] == '09 error':
                    logging.warning(tab_print("Warning! %s in %s mode has an error status. Attempting to resubmit..."%(project,polarity), 3))

                pathname = os.path.join(row['output_dir'],'%s_%s'%(row['parent_dir'],polarity))
                fbmn_filename = os.path.join(pathname,'%s_%s_gnps2-fbmn-task.txt'%(row['parent_dir'],polarity))
                if os.path.isfile(fbmn_filename)==True and overwrite_fbmn==False:
                    continue
                # if os.path.isfile(fbmn_filename)==True and overwrite_fbmn==True:
                #     wait = 10
                #     logging.warning(tab_print("Warning! Overwrite is True and %s exists. Giving %s seconds to end process before submitting..."%(os.path.basename(fbmn_filename),wait), 2))
                #     time.sleep(wait)
                if validate_names is True:
                    _, validate_department, _ = vfn.field_exists(PurePath(row['parent_dir']), field_num=1)
                if validate_names is False:
                    validate_department = row['parent_dir'].split('_')[1]
                department = validate_department.lower()
                if department =='eb':
                    department = 'egsb'
                if not department in ['jgi','egsb']:
                    logging.warning(tab_print("Warning! %s does not have a valid department name in the second field. Skipping..."%(project), 3))
                    continue
                description = '%s_%s'%(project,polarity)
                spectra_file = f'USERUPLOAD/bpbowen/untargeted_tasks/{project}_{polarity}/{project}_{polarity}.mgf'
                quant_file = f'USERUPLOAD/bpbowen/untargeted_tasks/{project}_{polarity}/{project}_{polarity}_quant.csv'
                metadata_file = f'USERUPLOAD/bpbowen/untargeted_tasks/{project}_{polarity}/{project}_{polarity}_metadata.tab'
                raw_data = f'USERUPLOAD/bpbowen/raw_data/{department}/{project}'
                mgf_filename = os.path.join(row['output_dir'],'%s_%s'%(project,polarity),'%s_%s.mgf'%(project,polarity))
                mgf_lines = count_mgf_lines(mgf_filename)
                if mgf_lines == 0:
                    logging.warning(tab_print("Warning! %s in %s mode has mgf file but no mgf data. Skipping..."%(project, polarity), 3))
                    continue
                if not os.path.exists(raw_data.replace('USERUPLOAD/bpbowen/','/global/cfs/cdirs/metatlas/')):
                    logging.warning(tab_print("Warning! %s does not have raw data. Skipping..."%(project), 3))
                    continue
                params = set_fbmn_parameters(description, quant_file, spectra_file, metadata_file, raw_data)
                job_id = submit_quickstart_fbmn(params, "bpbowen")
                task_list = {'experiment':project,'polarity':polarity,'response':job_id}
                logging.info(tab_print("Submitted FBMN job.", 3))
                df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '04 running'
                write_fbmn_tasks_to_file(task_list,output_dir)
                index_list.append(i)            

        if len(index_list) > 0:
            index_list = list(set(index_list))
            cols = ['Key',
                    '%s_neg_status'%(tasktype),
                    '%s_pos_status'%(tasktype)]
            df = df.loc[df.index.isin(index_list),cols]
            df.replace('NaN', 0, inplace=True)
            df.fillna(0, inplace=True)
            update_table_in_lims(df,'untargeted_tasks',method='update')
            logging.info(tab_print("FBMN submission(s) complete. Use GNPS2 link (gnps2-page-link.txt) to monitor progress.", 1))

def count_mgf_lines(mgf_file):
    """
    Before submitting a job to GNPS2, check if the mgf file from mzmine has any data
    by providing its path on disk
    """
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
    This function is called by run_mzmine.py

    finds initiated mzmine tasks
    submits the tasks as mzmine jobs on perlmutter
    updates the LIMS table with status running if successful

    Direct_input is None (default) when all initated tasks will be considered for submission. Set direct_input
    to a csv list of project names if you only want to run this function on specific untargeted_tasks
    """
    tasktype = 'mzmine'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    else:
        if new_projects:
            df = df[df['parent_dir'].isin(new_projects['parent_dir'])]
        else:
            logging.info(tab_print("No new MZmine jobs to submit!", 1))
            return
    status_list = ['01 initiation']
    df = subset_df_by_status(df,tasktype,status_list)
    if df.empty:
        logging.info(tab_print("No new MZmine jobs to submit!", 1))
        return
    if df.shape[0] > 20:
        logging.info(tab_print('There are too many new projects to be submitted (%s), please check if this is accurate. Exiting script.'%(df.shape[0]), 1))
        return
    if not df.empty:
        index_list = []
        logging.info(tab_print("Total of %s new MZmine job(s) with status %s to submit"%(df.shape[0],status_list), 1))
        for i,row in df.iterrows():
            polarity_list = check_for_polarities(row['output_dir'],row['parent_dir'])
            if polarity_list is None:
                logging.warning(tab_print("Warning! Project %s does not have a negative or a positive polarity directory. Skipping..."%(row['parent_dir']), 2))
                continue
            for polarity in polarity_list:
                polarity_short = polarity[:3]
                if row['%s_%s_status'%(tasktype,polarity_short)] == '12 not relevant' or row['%s_%s_status'%('mzmine',polarity_short)] == '07 complete':
                    continue
                pathname = os.path.join(row['output_dir'],'%s_%s'%(row['parent_dir'],polarity))
                submission_script_filename = os.path.join(pathname,'%s_%s_mzmine.sh'%(row['parent_dir'],polarity))
                if os.path.isfile(submission_script_filename)==True:
                    with open(submission_script_filename,'r') as fid:
                        logging.info(tab_print("Submitting %s mode mzmine job for project %s"%(polarity, row['parent_dir']), 2))
                        cmd = fid.read()
                        sbatch_output = subprocess.run(cmd, shell=True,capture_output=True,text=True)
                        logging.info(tab_print("%s"%(sbatch_output.stdout), 3))
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
            update_table_in_lims(df,'untargeted_tasks',method='update')
            logging.info(tab_print("MZmine submission(s) complete. Use sqs to monitor progress.", 1))
    
def set_fbmn_parameters(description, quant_file, spectra_file, metadata_file, raw_data):
    """
    Hard coded parameters and user-defined parameters are formatted by passing
    the arguments for file location
    """
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

def update_mzmine_status_in_untargeted_tasks(direct_input=None,skip_update=False):
    """
    This function is called by run_mzmine.py, run_fbmn.py, download_fbmn_results.py, export_untargeted_results.py and check_untargeted_status.py
    
    finds initiated or running mzmine tasks
    checks if they have valid output
    updates the LIMS table with status complete if successful

    Direct_input is None (default) when all initated or running tasks haven't produced output yet. Set direct_input
    to a csv list of project names if you only want to run this function on specific untargeted_tasks
    """
    if skip_update:
        return
    tasktype = 'mzmine'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    status_list = ['01 initiation','04 running']
    df = subset_df_by_status(df,tasktype,status_list)
    if df.empty:
        logging.info(tab_print("No MZmine jobs to update!", 1))
    if not df.empty:
        #logging.info(tab_print("Total of %s project(s) with MZmine status %s to attempt to update"%(df.shape[0],status_list), 1))
        index_list = []
        for i,row in df.iterrows():
            polarity_list = check_for_polarities(row['output_dir'],row['parent_dir'])
            if polarity_list is None:
                logging.warning(tab_print("Warning! Project %s does not have a negative or a positive polarity directory. Skipping..."%(row['parent_dir']), 2))
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
                    logging.info(tab_print("Working on %s in %s mode"%(row['parent_dir'],polarity), 2))
                    logging.info(tab_print("All MZmine output files found in %s directory, continuing..."%(polarity), 3))
                    logging.info(tab_print("Calculating feature and background counts and updating LIMS table", 3))
                    feature_count, exctrl_count, msms_count = calculate_counts_for_lims_table(peakheight_filename,mgf_filename)
                    df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '07 complete'
                    df.loc[i,'num_%s_features'%(polarity_short)] = int(feature_count)
                    df.loc[i,'num_%s_features_in_exctrl'%(polarity_short)] = int(exctrl_count)
                    df.loc[i,'num_%s_msms'%(polarity_short)] = int(msms_count)
                    index_list.append(i)

                    if feature_count > 0:
                        logging.info(tab_print("Filtering features by peak height ratio compared to background (control) signal and writing to file", 4))
                        df_filtered = filter_features_by_background(peakheight_filename)
                        df_filtered_path = os.path.join(pathname, '%s_%s_peak-height-filtered.csv' % (row['parent_dir'], polarity))
                        df_filtered.to_csv(df_filtered_path, index=False)
                        recursive_chown(df_filtered_path, 'metatlas')
                    else:
                        logging.warning(tab_print("Warning! Not writing filtered features to file because no features were found", 4))
        
        if len(index_list) > 0:
            logging.info(tab_print("Updating statuses in LIMS table for %s projects..."%(len(index_list)), 1))
            index_list = list(set(index_list))
            cols = ['Key',
                    '%s_neg_status'%(tasktype),'%s_pos_status'%(tasktype),
                    'num_pos_features', 'num_neg_features',
                    'num_neg_features_in_exctrl','num_pos_features_in_exctrl',
                    'num_neg_msms','num_pos_msms']
            df = df.loc[df.index.isin(index_list),cols]
            df.replace('NaN', 0, inplace=True)
            df.fillna(0, inplace=True)
            update_table_in_lims(df,'untargeted_tasks',method='update')
    logging.info(tab_print("MZmine status update complete.", 1))

def filter_features_by_background(peakheight_filename=None,background_designator=["TxCtrl,ExCtrl"],background_ratio=3):
    """
    Accepts a peakheight file and filters out features that have a max peak height in exp samples that is less than
    3 times the peak height in the control samples.
    """
    data_table = pd.read_csv(peakheight_filename, sep=',')
    features_to_keep = [0] # keep the header
    for i,row in data_table.iterrows():
        if i == 0:
            continue # skip the header
        exctrl_columns = [col for col in row.index if any(designator in col for designator in background_designator) and 'mzML' in col]
        exctrl_max = None
        feature_max = None
        if exctrl_columns:
            exctrl_max = row[exctrl_columns].max()
        else:
            logging.warning(tab_print("Warning! No background samples with designation %s could be found. Not filtering and returning empty dataframe."%(background_designator), 4))
            return data_table.head(1)
        feature_columns = [col for col in row.index if any(designator not in col for designator in background_designator) and 'mzML' in col]
        if feature_columns:
            feature_max = row[feature_columns].max()
        else:
            logging.warning(tab_print("Warning! No samples with features could be found. Not filtering and returning empty dataframe.", 4))
            return data_table.head(1)
        if exctrl_max and feature_max:
            if feature_max > exctrl_max*background_ratio:
                features_to_keep.append(i)
        else:
            logging.warning(tab_print("Warning! Could not calculate per-feature max intensities for peak height table. Not filtering and returning empty dataframe.", 4))
            return data_table.head(1)
    if len(features_to_keep) <= 1: # only the header row was retained
        logging.warning(tab_print("Warning! No features passed the background filter. Returning empty dataframe.", 5))
        return data_table.head(1)
    elif len(features_to_keep) > 1:
        data_table_filtered = data_table.loc[features_to_keep]
        difference = data_table.shape[0] - data_table_filtered.shape[0]
        logging.info(tab_print("%s features removed by background filter"%(difference), 5))
        return data_table_filtered
    else:
        logging.warning(tab_print("Warning! Something went wrong with the background filter. Returning empty dataframe.", 5))
        return data_table.head(1)

def calculate_counts_for_lims_table(peakheight_filename=None, mgf_filename=None, background_designator=["TxCtrl,ExCtrl"]):
    """
    Accepts a peakheight file and mgf file and returns the number of features
    in exp and and control samples found in the peakheightfile plus number of 
    feature ids in mgf file.
    """
    data_table = pd.read_csv(peakheight_filename, sep=',')
    exctrl_count = 0
    feature_count = 0
    msms_count = 0
    
    exctrl_columns = [col for col in row.index if any(designator in col for designator in background_designator) and 'mzML' in col]
    if exctrl_columns:
        exctrl_rows = data_table[data_table[exctrl_columns].gt(0).any(axis=1)]
        exctrl_count = int(exctrl_rows.shape[0])
        logging.info(tab_print("%s features in control samples"%(exctrl_count), 4))
    else:
        logging.warning(tab_print("Warning! No ExCtrl samples found in peak height file", 4))

    feature_columns = [col for col in row.index if any(designator not in col for designator in background_designator) and 'mzML' in col]
    if feature_columns:
        feature_rows = data_table[data_table[feature_columns].gt(0).any(axis=1)]
        feature_count = int(feature_rows.shape[0])
        logging.info(tab_print("%s features in experimental samples"%(feature_count), 4))
    if feature_count < 10:
        logging.warning(tab_print("Warning! Less than 10 features found in peak height file", 4))

    if os.path.isfile(mgf_filename):
        cmd = 'cat %s | grep FEATURE_ID| wc -l'%mgf_filename
        result = subprocess.run(cmd,stdout=subprocess.PIPE, shell=True)
        if result:
            msms_count = int(result.stdout.strip())
            logging.info(tab_print("%s msms features in mgf file"%(msms_count), 4))

    return feature_count, exctrl_count, msms_count

def update_fbmn_status_in_untargeted_tasks(direct_input=None,skip_update=False):
    """
    This function is called by run_mzmine.py, run_fbmn.py, download_fbmn_results.py, export_untargeted_results.py and check_untargeted_status.py
    
    finds running or waiting fbmn tasks
    checks if they have valid output
    updates the LIMS table with status complete if successful

    Direct_input is None (default) when all running or waiting tasks haven't produced output yet. Set direct_input
    to a csv list of project names if you only want to run this function on specific untargeted_tasks
    """
    if skip_update:
        return
    tasktype='fbmn'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    status_list = ['04 running','13 waiting']
    df = subset_df_by_status(df,tasktype,status_list)
    if df.empty:
        logging.info(tab_print("No FBMN jobs to update!", 1))
    if not df.empty:
        index_list = []
        #logging.info(tab_print("Total of %s project(s) with FBMN status %s to attempt to update"%(df.shape[0],status_list), 1))
        for i,row in df.iterrows():
            polarity_list = check_for_polarities(row['output_dir'],row['parent_dir'])
            if polarity_list is None:
                logging.warning(tab_print("Warning! Project %s does not have a negative or a positive polarity directory. Skipping..."%(row['parent_dir']), 2))
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
                            logging.info(tab_print("Updating FBMN job status for %s in %s mode to complete"%(row['parent_dir'],polarity), 2))
                            df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '07 complete'
                    elif status=='FAILED':
                        if df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] != '09 error':
                            logging.info(tab_print("Updating FBMN job status for %s in %s mode to error"%(row['parent_dir'],polarity), 2))
                            df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '09 error'
                    elif status in ['RUNNING','QUEUED']:
                        if df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] != '04 running':
                            logging.info(tab_print("Updating FBMN job status for %s in %s mode to running"%(row['parent_dir'],polarity), 2))
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
            update_table_in_lims(df,'untargeted_tasks',method='update')
    logging.info(tab_print("FBMN status update complete.", 1))

def write_mzmine_sbatch_and_runner(basepath,batch_filename,parent_dir,filelist_filename):
    """
    Write the sbatch and runner files for mzmine submission via slurm
    """
    mzmine_launcher = '/global/common/software/m2650/mzmine_parameters/MZmine/MZmine-3.7.2/bin/MZmine -threads auto -t /pscratch/sd/b/bkieft/untargeted/tmp'

    sbatch_filename = '%s_mzmine-sbatch.sbatch'%os.path.join(basepath,parent_dir)
    runner_filename = '%s_mzmine.sh'%os.path.join(basepath,parent_dir)
    s = '%s -batch %s -input %s\necho MZMINE IS DONE\n\n\n'%(mzmine_launcher,batch_filename,filelist_filename)

    with open(sbatch_filename,'w') as fid:
        fid.write('%s\n%s\n'%(SLURM_PERLMUTTER_HEADER.replace('slurm','%s-%s'%(os.path.join(basepath,parent_dir),'mzmine')),s))
    with open(runner_filename,'w') as fid:
        fid.write('sbatch %s'%sbatch_filename)

def write_metadata_per_new_project(df: pd.DataFrame,control=["TxCtrl","ExCtrl"],raw_data_path="/",validate_names=True) -> list:
    """
    Takes a LIMS table (usually raw data from lcmsruns_plus) and creates
    mzmine metadata for writing a new untargeted_tasks directory
    """
    new_project_list = []
    # Iterate through each unique parent_dir
    for parent_dir in df['parent_dir'].unique():
        # Filter the DataFrame for the current parent_dir
        df_filtered = df[df['parent_dir'] == parent_dir]

        # Finish raw_data path
        if validate_names is True:
            _, validate_department, _ = vfn.field_exists(PurePath(parent_dir), field_num=1)
            department = validate_department.lower()
        if validate_names is False:
            department = parent_dir.split('_')[1].lower()
        if department =='eb':
            department = 'egsb'
        full_mzml_path = os.path.join(raw_data_path,department,parent_dir)

        # Initialize info dict
        project_dict = {'parent_dir': parent_dir}
        project_dict['positive'] = {}
        project_dict['negative'] = {}

        # Determine which polarities need metadata
        positive_file_subset = df_filtered['basename'][df_filtered['basename'].apply(
            lambda x: 'POS' in x.split('_')[9] and 
            'QC' not in x.split('_')[12] and
            'InjBl' not in x.split('_')[12] and
            'ISTD' not in x.split('_')[12])].to_list()
        if positive_file_subset:
            positive_file_list = [os.path.join(full_mzml_path, file) for file in positive_file_subset]
            positive_file_count = len(positive_file_list)
            if validate_names is True:
                file_name_validation = [vfn.validate_file_name(PurePath(file),minimal=True,print_logger=False) for file in positive_file_list]
                if all(file_name_valid == True for file_name_valid in file_name_validation) is False:
                    logging.warning(tab_print("Warning! Project %s has one or more mzML files that do not have valid format. Not adding project to untargeted tasks..."%(parent_dir), 2))
                    continue
                project_name_validation = [vfn.parent_dir_num_fields(PurePath(file)) for file in positive_file_list]
                if all(len(project_name_valid) == 0 for project_name_valid in project_name_validation) is False: # Prints empty list if project name valid
                    logging.warning(tab_print("Warning! Parent dir %s does not have valid format. Not adding project to untargeted tasks..."%(parent_dir), 2))
                    continue
        else:
            positive_file_count = 0

        negative_file_subset = df_filtered['basename'][df_filtered['basename'].apply(
            lambda x: 'NEG' in x.split('_')[9] and 
            'QC' not in x.split('_')[12] and
            'InjBl' not in x.split('_')[12] and
            'ISTD' not in x.split('_')[12])].to_list()
        if negative_file_subset:
            negative_file_list = [os.path.join(full_mzml_path, file) for file in negative_file_subset]
            negative_file_count = len(negative_file_list)
            if validate_names is True:
                file_name_validation = [vfn.validate_file_name(PurePath(file),minimal=True,print_logger=False) for file in negative_file_list]
                if all(file_name_valid == True for file_name_valid in file_name_validation) is False:
                    logging.warning(tab_print("Warning! Project %s has one or more mzML files that do not have valid format. Not adding project to untargeted tasks..."%(parent_dir), 2))
                    continue
                project_name_validation = [vfn.parent_dir_num_fields(PurePath(file)) for file in negative_file_list]
                if all(len(project_name_valid) == 0 for project_name_valid in project_name_validation) is False: # Prints empty list if project name valid
                    logging.warning(tab_print("Warning! Parent dir %s does not have valid format. Not adding project to untargeted tasks..."%(parent_dir), 2))
                    continue
        else:
            negative_file_count = 0
        
        # Counts and mzml list
        if positive_file_count > 0 and negative_file_count > 0:
            project_dict['polarities'] = ['positive', 'negative']
            logging.info(tab_print("Project %s has both polarities with %s positive and %s negative files"%(parent_dir,positive_file_count,negative_file_count), 2))
            project_dict['negative']['file_count'] = negative_file_count
            project_dict['positive']['file_count'] = positive_file_count
            project_dict['negative']['file_list'] = negative_file_list
            project_dict['positive']['file_list'] = positive_file_list
        elif positive_file_count > 0:
            logging.info(tab_print("Project %s has positive mode only with %s files"%(parent_dir,positive_file_count), 2))
            project_dict['polarities'] = ['positive']
            project_dict['positive']['file_count'] = positive_file_count
            project_dict['positive']['file_list'] = positive_file_list
        elif negative_file_count > 0:
            logging.info(tab_print("Project %s has negative mode only with %s files"%(parent_dir,negative_file_count), 2))
            project_dict['polarities'] = ['negative']        
            project_dict['negative']['file_count'] = negative_file_count
            project_dict['negative']['file_list'] = negative_file_list
        else:
            project_dict['polarities'] = [None]
            continue

        for polarity in project_dict['polarities']:
            # Create metadata dataframe
            polarity_short = polarity[:3].upper()
            metadata_df = df_filtered[['basename', 'sample_group']].copy()
            metadata_df['basename'] = metadata_df['basename'].fillna('')
            metadata_df['sample_group'] = metadata_df['sample_group'].fillna('')
            metadata_df['sample_group'] = metadata_df['sample_group'].apply(lambda x: x.lower())
            ugroups = metadata_df['sample_group'].unique()
            ugroups = [g for g in ugroups if any(c.lower() in g for c in control)]
            metadata_df.rename(columns={'basename': 'filename', 'sample_group': 'ATTRIBUTE_sampletype'}, inplace=True)
            metadata_df['filename'] = metadata_df['filename'].apply(lambda x: os.path.basename(x))
            cols = ['CONTROL', 'CASE', 'ATTRIBUTE_media']
            for i, g in enumerate(ugroups):
                if i < len(cols):
                    metadata_df.loc[:, cols[i]] = g
            metadata_df = metadata_df[metadata_df['filename'].str.split('_').str[9].str.contains(polarity_short, na=False)]
            project_dict[polarity]['metadata_df'] = metadata_df

        new_project_list.append(project_dict)

    return new_project_list

def update_new_untargeted_tasks(update_lims=True,skip_update=False,background_designator=["TxCtrl","ExCtrl"], \
                                output_dir='/global/cfs/cdirs/metatlas/projects/untargeted_tasks',
                                raw_data_dir='/global/cfs/cdirs/metatlas/raw_data', validate_names=True):
    """
    This script is called by run_mzmine.py before the untargeted pipeline kicks off

    given all directories that are in the table for potential untargeted tasks
    and all the directories in the raw data folders
    return all the directories in the raw data folders
    that are not in the untargeted tasks
    
    The strip command is because there is one folder that ends in a 
    space and labkey doesn't allow this
    """
    if skip_update:
        return
    
    # Get lcmsrun table and subset
    df = get_table_from_lims('lcmsrun_plus')
    df = df[pd.notna(df['mzml_file'])]
    df['basename'] = df['mzml_file'].apply(os.path.basename)
    df.sort_values('timeepoch',ascending=False,inplace=True)
    df.drop_duplicates(['basename','parent_dir'],inplace=True)
    df.drop(columns=['mzml_file_container'],inplace=True)
    df.replace('',np.nan,inplace=True)

    # Get untargeted tasks table
    df_untargeted = get_table_from_lims('untargeted_tasks')
    
    logging.info(tab_print("Comparing LIMS untargeted tasks with raw files on disk...", 1))

    #Check that files in lcmsruns table have been sitting around for at least 3 hours (note time zones may vary)
    time_old = df.groupby('parent_dir')['timeepoch'].max() < (time.time()-3*60*60)
    time_new = df.groupby('parent_dir')['timeepoch'].max() >= (time.time()-3*60*60)
    time_old_folders = time_old[time_old==True].index.tolist()
    time_new_folders = time_new[time_new==True].index.tolist()
    #Check that files have "M2" in the name
    dirs_with_m2_files = df.groupby('parent_dir')['basename'].apply(lambda x: any('_MS2_' in filename for filename in x))
    dirs_with_m2_files = dirs_with_m2_files[dirs_with_m2_files].index.tolist()
    # Print parent_dirs with files less than 3 hours old
    if time_new_folders:
        logging.warning(tab_print("Warning! There are data directories with files less than 3 hours old!", 2))
        logging.info(tab_print("Skipping these projects:", 3))
        for folder in time_new_folders:
            logging.info(tab_print(folder, 4))

    # Get all folders in raw data table
    all_folders = df.loc[df['polarity'].isin(['POS','NEG']),'parent_dir'].unique()
    all_folders = [a.strip() for a in all_folders]

    # Get all folders in untargeted tasks table
    if df_untargeted.shape[0]>0:
        folders_in_tasks = df_untargeted['parent_dir']
        folders_in_tasks = [a.strip() for a in folders_in_tasks]
    else:
        folders_in_tasks = []

    # Find difference between lcmsrun table and untargeted table to find projects to initiate
    logging.info(tab_print("Finding new projects to initate...", 1))
    new_folders = np.setdiff1d(all_folders,folders_in_tasks)
    new_folders = list(set(new_folders) & set(time_old_folders) & set(dirs_with_m2_files))
    if len(new_folders) == 0:
        logging.info(tab_print("No new projects to add to untargeted tasks!", 2))
        return None
    
    # Check for polarities by looking for positive and negative mzml files
    df_new = df[df['parent_dir'].isin(new_folders)]
    logging.info(tab_print("Checking for polarities in new projects and validating mzml file names...", 1))
    new_project_info_list = write_metadata_per_new_project(df=df_new,control=background_designator,raw_data_path=raw_data_dir,validate_names=validate_names)
    new_project_info_list_subset = [d for d in new_project_info_list if d.get('polarities') is not None]

    # Create metadata for new projects with relevant polarities
    if len(new_project_info_list_subset)>0:
        lims_untargeted_list = []
        for new_project_dict in new_project_info_list_subset:
            project_name = new_project_dict['parent_dir']
            polarity_list = new_project_dict['polarities']
            logging.info(tab_print("Working on %s..."%(project_name), 1))

            # Write some basic info to LIMS
            lims_untargeted_table_updater = {'parent_dir': project_name}
            lims_untargeted_table_updater['conforming_filenames'] = True
            lims_untargeted_table_updater['file_conversion_complete'] = True
            lims_untargeted_table_updater['output_dir'] = output_dir
            if validate_names is True:
                _, validate_machine_name, _ = vfn.field_exists(PurePath(project_name), field_num=6)
            if validate_names is False:
                validate_machine_name = project_name.split('_')[6]
            if validate_machine_name.lower() in ("iqx", "idx"):
                mzmine_running_parameters = mzine_batch_params_file_iqx
                mzmine_parameter = 5
            else:
                mzmine_running_parameters = mzine_batch_params_file
                mzmine_parameter = 2
            lims_untargeted_table_updater['mzmine_parameter_sheet'] = mzmine_running_parameters
            lims_untargeted_table_updater['mzmine_parameter_row'] = mzmine_parameter

            for polarity in ['positive','negative']: # Don't initiate mzmine jobs on polarities that don't have sample mzmls
                polarity_short = polarity[:3]
                parent_dir = '%s_%s'%(project_name,polarity)
                basepath = os.path.join(output_dir,parent_dir)
                metadata_header = f"{polarity_short}_metadata_file"
                mzmine_status_header = f"mzmine_{polarity_short}_status"
                fbmn_status_header = f"fbmn_{polarity_short}_status"
                file_count_header = f"num_{polarity_short}_files"
                if polarity not in polarity_list: # Write blank columns to LIMS and skip writing directory/files to disk
                    lims_untargeted_table_updater[metadata_header] = ""
                    lims_untargeted_table_updater[file_count_header] = 0
                    lims_untargeted_table_updater[mzmine_status_header] = '12 not relevant'
                    lims_untargeted_table_updater[fbmn_status_header] = '12 not relevant'
                else: # Write directory/files to disk for present polarity and fill in LIMS columns
                    logging.info(tab_print("Writing MZmine submission input files...",2))
                    logging.info(tab_print("%s metadata file (*_metadata.tab)"%(polarity), 3))
                    os.mkdir(basepath)
                    recursive_chown(basepath, 'metatlas')
                    metadata_df = new_project_dict[polarity]['metadata_df']
                    metadata_filename = os.path.join(basepath,'%s_metadata.tab'%(parent_dir))
                    metadata_df.to_csv(metadata_filename, sep='\t', index=False)
                    
                    logging.info(tab_print("%s MZmine parameter file (*_batch-params.xml)"%(polarity), 3))
                    params_filename = build_untargeted_filename(output_dir,project_name,polarity,'batch-params-mzmine')
                    with open(mzmine_running_parameters,'r') as fid:
                        orig_params = fid.read()
                    new_param_path = os.path.join(basepath,parent_dir)
                    custom_params = orig_params.replace('/Users/bpb/Downloads/mzmine_outputs',new_param_path)
                    with open(params_filename,'w') as fid:
                        fid.write('%s'%custom_params)

                    logging.info(tab_print("%s mzML path list file (*_filelist.txt)"%(polarity), 3))
                    file_list = new_project_dict[polarity]['file_list']
                    filelist_filename = os.path.join(basepath,'%s_filelist.txt'%(parent_dir))
                    with open(filelist_filename,'w') as fid:
                        for mzml in file_list:
                            fid.write('%s\n'%mzml)

                    logging.info(tab_print("%s sbatch submission script and runner files (*.sh, *.sbatch)"%(polarity), 3))
                    write_mzmine_sbatch_and_runner(basepath,params_filename,parent_dir,filelist_filename)

                    lims_untargeted_table_updater[metadata_header] = metadata_filename
                    lims_untargeted_table_updater[file_count_header] = new_project_dict[polarity]['file_count']
                    lims_untargeted_table_updater[mzmine_status_header] = '01 initiation'
                    lims_untargeted_table_updater[fbmn_status_header] = '13 waiting'

            lims_untargeted_list.append(lims_untargeted_table_updater)

        lims_untargeted_df = pd.DataFrame(lims_untargeted_list)

        if update_lims is True:
            if lims_untargeted_df.shape[0] > 0:
                logging.info(tab_print("Updating LIMS table with %s new projects..."%(lims_untargeted_df.shape[0]), 1))
                update_table_in_lims(lims_untargeted_df,'untargeted_tasks',method='insert',max_size=1000)
                logging.info(tab_print("LIMS untargeted tasks table update complete.", 2))
            else:
                logging.info(tab_print("LIMS does not need updating!", 1))
    
    logging.info(tab_print("Exported new project info for MZmine submission.", 1))
    return lims_untargeted_df