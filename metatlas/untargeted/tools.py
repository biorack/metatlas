import sys
sys.path.insert(0,'/global/common/software/m2650/labkey-api-python') # https://github.com/LabKey/labkey-api-python
from labkey.api_wrapper import APIWrapper
import numpy as np
import pandas as pd
from pathlib2 import PurePath, Path
import os
import requests
from bs4 import BeautifulSoup
from datetime import datetime, time
import time
import glob
import subprocess
import math
sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/')
import metatlas.tools.validate_filenames as vfn
import subprocess
import grp
import logging
import paramiko
import zipfile
import shutil
from typing import List, Dict, Union, Optional

key_file = '/global/cfs/cdirs/metatlas/labkey_user.txt'
with open(key_file,'r') as fid:
    api_key = fid.read().strip()
labkey_server='metatlas.jgi.doe.gov'
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

BATCH_FILE_PATH = '/global/common/software/m2650/mzmine_parameters/batch_files'
BINARY_PATH = '/global/common/software/m2650/mzmine_parameters/MZmine'
mzine_batch_params_file = f"{BATCH_FILE_PATH}/mzmine-3.7.2-batchparams.xml"
mzine_batch_params_file_iqx = f"{BATCH_FILE_PATH}/IQX-mzmine-3.7.2-batchparams.xml"
mzine_batch_params_file_pos = f"{BATCH_FILE_PATH}/POS-mzmine-3.7.2-batchparams.xml"
mzine_batch_params_file_neg = f"{BATCH_FILE_PATH}/NEG-mzmine-3.7.2-batchparams.xml"

def call_logger(log_filename: str, log_level: str, log_format: str):
    logging.basicConfig(filename=log_filename, level=log_level, format=log_format, filemode='a')

def call_logger(log_filename: str, log_level: str, log_format: str, log_to_stdout: bool):
    if log_to_stdout:
        logging.basicConfig(stream=sys.stdout, level=log_level, format=log_format)
    else:
        logging.basicConfig(filename=log_filename, level=log_level, format=log_format, filemode='a')

def start_script(script:str=None) -> str:
    """
    Kick off script with a timestamp for log
    """
    if script is not None:
        return tab_print('Successfully started %s on %s:'%(script,datetime.now()), 0)
    else:
        sys.exit("No script name provided to start_script")

def end_script(script:str=None) -> str:
    """
    End script with a timestamp for log
    """
    if script is not None:
        return tab_print('Successfully completed %s on %s:'%(script,datetime.now()), 0)
    else:
        sys.exit("No script name provided to end_script")

def tab_print(message:str="",indent_level:int=0) -> str:
    """
    Print a message with a specified number of tabs
    """
    return ("\t"*indent_level + message)

def check_for_polarities(output_dir:str=None, parent_dir:str=None) -> list:
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

def subset_df_by_status(df:pd.DataFrame=None, tasktype:str=None, status:list=[], inverse:bool=False) -> pd.DataFrame:
    """
    Subset a dataframe by status in LIMS
    """
    if inverse == True:
        df_subset = df[~df['%s_pos_status'%tasktype].isin(status) | ~df['%s_neg_status'%tasktype].isin(status)]
    if inverse == False:
        df_subset = df[df['%s_pos_status'%tasktype].isin(status) | df['%s_neg_status'%tasktype].isin(status)]
    return df_subset

def filter_common_bad_project_names(df:pd.DataFrame=None) -> pd.DataFrame:
    """
    Takes a df object (usually from exporting the LIMS untargeted tasks table) and removes
    projects with common strings that are in untargeted tasks which should not be run
    """
    df = df[~(df['parent_dir'].str.contains(' '))]
    #df = df[~(df['parent_dir'].str.contains('&'))]
    #df = df[~(df['parent_dir'].str.contains('Niyogi'))]
    #df = df[~(df['parent_dir'].str.contains('_partial'))]
    #df = df[~(df['parent_dir'].str.contains('_old'))]
    return df

def write_fbmn_tasks_to_file(
    task_list: Dict,
    output_dir: str
) -> None:
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
            logging.info(tab_print("GNPS2 task file for %s mode written to %s"%(polarity,final_filename), 3))
    else:
        logging.warning(tab_print("Warning! GNPS2 FBMN task ID not found. File gnps2-fbmn-task.txt not written.", 3))

def zip_and_upload_untargeted_results(
    download_folder: str,
    output_dir: str,
    doc_name: str,
    add_documentation: bool,
    skip_zip_upload: bool,
    abridged_filenames: bool,
    upload: bool,
    overwrite_zip: bool,
    overwrite_drive: bool,
    direct_input: Optional[str] = None,
    min_features_admissible: int = 0
) -> None:
    """    
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
    if skip_zip_upload:
        logging.info("Skipping Step 7/7: Zipping up and (optionally) uploading output folders to gdrive...")
        return
    logging.info("Step 7/7: Zipping up and (optionally) uploading output folders to gdrive...")
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    status_list = ['07 complete']
    if direct_input is None:
        df = subset_df_by_status(df,'fbmn',status_list)
        df = subset_df_by_status(df,'mzmine',status_list)
        df = subset_df_by_status(df,'mzmine',['09 error'],inverse=True) # Make sure one of the polarities doesn't have an error
        df = subset_df_by_status(df,'fbmn',['09 error'],inverse=True) # Make sure one of the polarities doesn't have an error
    if df.empty:
        logging.info(tab_print("No completed untargeted results to zip and upload!", 1))
        return
    if not df.empty:
        logging.info(tab_print("%s project(s) total with complete mzmine and fbmn status."%(df.shape[0]), 1))
        zip_count = 0
        upload_count = 0
        for i,row in df.iterrows():
            project_name = row['parent_dir']
            output_zip_archive = os.path.join(download_folder,'%s.zip'%project_name)
            polarity_list = check_for_polarities(output_dir,project_name)
            if overwrite_zip==False and os.path.exists(output_zip_archive):
                #logging.warning(tab_print("Warning! Zip archive for %s exists. Set overwrite_zip to True if you want to replace. Skipping zip and upload..."%(project_name), 1))
                continue
            if overwrite_zip==True or not os.path.exists(output_zip_archive):
                if polarity_list is None:
                    logging.warning(tab_print("Warning! Project %s does not have a negative or a positive polarity directory. Skipping..."%(project_name), 1))
                    continue
                # Create variables for possible polarity directories and files
                neg_directory = os.path.join(output_dir, '%s_%s'%(project_name, 'negative'))
                pos_directory = os.path.join(output_dir, '%s_%s'%(project_name, 'positive'))
                neg_mzmine_file = os.path.join(neg_directory, '%s_negative_peak-height.csv'%project_name)
                pos_mzmine_file = os.path.join(pos_directory, '%s_positive_peak-height.csv'%project_name)
                neg_fbmn_file = os.path.join(neg_directory, '%s_negative_gnps2-fbmn-library-results.tsv'%project_name)
                pos_fbmn_file = os.path.join(pos_directory, '%s_positive_gnps2-fbmn-library-results.tsv'%project_name)
                if 'negative' in polarity_list and 'positive' in polarity_list:
                    # Check that mzmine and fbmn "marker" files exist, they've probably finished sucessfully (double-check since status need to all be 'complete' above)
                    if os.path.exists(neg_mzmine_file) and os.path.exists(pos_mzmine_file) and os.path.exists(neg_fbmn_file) and os.path.exists(pos_fbmn_file):
                        try:
                            recursive_chown(neg_directory, 'metatlas')
                            recursive_chown(pos_directory, 'metatlas')
                        except:
                            logging.info(tab_print("Note: Could not change group ownership of %s or %s."%(neg_directory,pos_directory), 1))
                        neg_feature_counts = check_peak_height_table(neg_mzmine_file)
                        pos_feature_counts = check_peak_height_table(pos_mzmine_file)
                        if (neg_feature_counts + pos_feature_counts) > min_features_admissible:
                            neg_directory = os.path.join(output_dir, '%s_%s'%(project_name, 'negative'))
                            pos_directory = os.path.join(output_dir, '%s_%s'%(project_name, 'positive'))
                            # Get rid of the mzmine job id files
                            neg_mzmine_job_id_filename = os.path.join(neg_directory,'%s_%s_mzmine-job-id.txt'%(project_name, 'negative'))
                            pos_mzmine_job_id_filename = os.path.join(pos_directory,'%s_%s_mzmine-job-id.txt'%(project_name, 'positive'))
                            if os.path.exists(neg_mzmine_job_id_filename):
                                os.remove(neg_mzmine_job_id_filename)
                            if os.path.exists(pos_mzmine_job_id_filename):
                                os.remove(pos_mzmine_job_id_filename)

                            zip_untargeted_results(target_dirs=[neg_directory,pos_directory], abridged_filenames=abridged_filenames, \
                                                   add_documentation=add_documentation, download_folder=download_folder, doc_name=doc_name, output_zip_archive=output_zip_archive)
                            zip_count += 1

                            if upload == True and os.path.exists(output_zip_archive):
                                upload_success = upload_to_google_drive(output_zip_archive,overwrite_drive)
                                if upload_success:
                                    upload_count += 1
                        else:
                            logging.warning(tab_print("Warning! Project %s has less than %s features in the %s peak height table(s). Skipping zip and upload..."%(project_name,min_features_admissible,polarity_list), 1))
                            continue
                elif not 'negative' in polarity_list and 'positive' in polarity_list:
                    if os.path.exists(pos_mzmine_file) and os.path.exists(pos_fbmn_file):
                        try:
                            recursive_chown(pos_directory, 'metatlas')
                        except:
                            logging.info(tab_print("Note: Could not change group ownership of %s."%(pos_directory), 1))
                        pos_feature_counts = check_peak_height_table(pos_mzmine_file)
                        if pos_feature_counts > min_features_admissible:
                            pos_directory = os.path.join(output_dir, '%s_%s'%(project_name, 'positive'))
                            # Get rid of the mzmine job id files
                            pos_mzmine_job_id_filename = os.path.join(output_dir,'%s_%s'%(project_name, 'positive'),'%s_%s_mzmine-job-id.txt'%(project_name, 'positive'))
                            if os.path.exists(pos_mzmine_job_id_filename):
                                os.remove(pos_mzmine_job_id_filename)

                            zip_untargeted_results(target_dirs=[pos_directory], abridged_filenames=abridged_filenames, \
                                                   add_documentation=add_documentation, download_folder=download_folder, doc_name=doc_name, output_zip_archive=output_zip_archive)
                            zip_count += 1

                            if upload == True and os.path.exists(output_zip_archive):
                                upload_success = upload_to_google_drive(output_zip_archive,overwrite_drive)
                                if upload_success:
                                    upload_count += 1
                        else:
                            logging.warning(tab_print("Warning! Project %s has less than %s features in the %s peak height table(s). Skipping zip and upload..."%(project_name,min_features_admissible,polarity_list), 1))
                            continue
                elif 'negative' in polarity_list and not 'positive' in polarity_list:
                    if os.path.exists(neg_mzmine_file) and os.path.exists(neg_fbmn_file):
                        try:
                            recursive_chown(neg_directory, 'metatlas')
                        except:
                            logging.info(tab_print("Note: Could not change group ownership of %s."%(neg_directory), 1))
                        neg_feature_counts = check_peak_height_table(neg_mzmine_file)
                        if neg_feature_counts > min_features_admissible:
                            neg_directory = os.path.join(output_dir, '%s_%s'%(project_name, 'negative'))
                            # Get rid of the mzmine job id files
                            neg_mzmine_job_id_filename = os.path.join(output_dir,'%s_%s'%(project_name, 'negative'),'%s_%s_mzmine-job-id.txt'%(project_name, 'negative'))
                            if os.path.exists(neg_mzmine_job_id_filename):
                                os.remove(neg_mzmine_job_id_filename)

                            zip_untargeted_results(target_dirs=[neg_directory], abridged_filenames=abridged_filenames, \
                                                   add_documentation=add_documentation, download_folder=download_folder, doc_name=doc_name, output_zip_archive=output_zip_archive)
                            zip_count += 1

                            if upload == True and os.path.exists(output_zip_archive):
                                upload_success = upload_to_google_drive(output_zip_archive,overwrite_drive)
                                if upload_success:
                                    upload_count += 1
                        else:
                            logging.warning(tab_print("Warning! Project %s has less than %s features in the %s peak height table(s). Skipping zip and upload..."%(project_name,min_features_admissible,polarity_list), 1))
                            continue
        
        if zip_count == 0:
            logging.info(tab_print("No newly completed untargeted projects to be zipped.", 2))
        
        logging.info(tab_print("%s new complete untargeted projects uploaded."%(upload_count), 1))


def zip_untargeted_results(
    abridged_filenames: bool,
    doc_name: str,
    add_documentation: bool,
    target_dirs: Optional[List[str]] = None,
    download_folder: Optional[str] = None,
    output_zip_archive: Optional[str] = None
) -> None:
    if target_dirs is None:
        logging.warning(tab_print("Warning! No target directory provided for renaming untargeted results files, but rename function is set to True.", 1))
        return
    if output_zip_archive is None:
        logging.warning(tab_print("Warning! No output zip archive provided for zipping untargeted results files, but zip function is set to True.", 1))
        return
    if download_folder is None:
        logging.warning(tab_print("Warning! No download location is provided for untargeted results.", 1))
        return

    # Zip the renamed files
    if add_documentation == True:
        logging.info(tab_print("Downloading latest GNPS2 user guide documentation to add to zip...", 1))
        doc_present = add_gnps2_documentation(download_folder=download_folder,doc_name=doc_name)
        if doc_present:
            if len(target_dirs) == 2:
                cmd = 'zip -rjq - %s %s %s >%s'%(target_dirs[0],target_dirs[1],os.path.join(download_folder,doc_name),output_zip_archive)
            elif len(target_dirs) == 1:
                cmd = 'zip -rjq - %s %s >%s'%(target_dirs[0],os.path.join(download_folder,doc_name),output_zip_archive)
        else:
            logging.warning(tab_print("Warning! Add documentation flag is True but GNPS2 user guide documentation not available. Not adding to zip...", 2))
            if len(target_dirs) == 2:
                cmd = 'zip -rjq - %s %s >%s'%(target_dirs[0],target_dirs[1],output_zip_archive)
            elif len(target_dirs) == 1:
                cmd = 'zip -rjq - %s >%s'%(target_dirs[0],output_zip_archive)
    else:
        if len(target_dirs) == 2:
            cmd = 'zip -rjq - %s %s >%s'%(target_dirs[0],target_dirs[1],output_zip_archive)
        elif len(target_dirs) == 1:
            cmd = 'zip -rjq - %s >%s'%(target_dirs[0],output_zip_archive)
    os.system(cmd)
    logging.info(tab_print("New untargeted results in %s zipped"%(target_dirs), 1))
    
    if abridged_filenames is True:
        rename_untargeted_files_in_archive(output_zip_archive=output_zip_archive)
    
    # Change permissions of resulting zip
    try:
        recursive_chown(output_zip_archive, 'metatlas')
    except:
        logging.info(tab_print("Note: Could not change group ownership of %s."%(output_zip_archive), 2))


def rename_untargeted_files_in_archive(
    output_zip_archive: Optional[str] = None
) -> None:
    if output_zip_archive is None:
        logging.warning(tab_print("Warning! No output zip archive provided for renaming untargeted results files, but rename function is set to True.", 1))
        return

    project_name = os.path.splitext(os.path.basename(output_zip_archive))[0]

    if len(project_name.split('_')) >= 9:
        date = project_name.split('_')[0]
        department = project_name.split('_')[1]
        submitter = project_name.split('_')[2]
        pid = project_name.split('_')[3]
        chromatography = project_name.split('_')[7]
    else:
        logging.info(tab_print("Note! Project name %s has fewer than expected fields. Skipping renaming files before zip..."%(project_name), 1))
        return

    # Check if project name follows the standard naming convention
    if not any(substring.lower() in chromatography.lower() for substring in ['C18', 'LIPID', 'HILIC']) or \
       not date.isdigit() or len(date) != 8:
            logging.warning(tab_print("Warning! Project name %s does not follow the standard naming convention. Skipping renaming files before zip..."%(project_name), 1))
            logging.warning(tab_print("Here is what could be extracted: Date: %s, Department: %s, Submitter: %s, PID: %s, Chromatography: %s"%(date, department, submitter, pid, chromatography), 2))
            return
    else:
        new_project_name = f"{date}_{department}_{submitter}_{pid}_{chromatography}"

    # Unzip the archive and rename all files
    temp_dir = f"/tmp/{project_name}"
    os.makedirs(temp_dir, exist_ok=True)
    try:
        with zipfile.ZipFile(output_zip_archive, 'r') as zip_ref:
            zip_ref.extractall(temp_dir)

        for root, _, files in os.walk(temp_dir):
            for file in files:
                if project_name in file:
                    new_file_name = file.replace(project_name, new_project_name)
                    os.rename(os.path.join(root, file), os.path.join(root, new_file_name))

        new_zip_path = output_zip_archive + ".new"
        with zipfile.ZipFile(new_zip_path, 'w') as zip_ref:
            for root, _, files in os.walk(temp_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    arcname = os.path.relpath(file_path, temp_dir)
                    zip_ref.write(file_path, arcname)

        # Overwrite the existing archive with the new zip file
        shutil.move(new_zip_path, output_zip_archive)
    except Exception as e:
        logging.error(tab_print(f"Error processing the zip archive: {e}", 2))
    finally:
        shutil.rmtree(temp_dir)

def check_peak_height_table(peak_height_file: str) -> None:
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

def add_gnps2_documentation(
    download_folder: str,
    doc_name: str = None
) -> bool:
    if doc_name is None:
        logging.warning(tab_print("Warning! GNPS2 user guide documentation file name not provided. Not downloaded!", 2))
        return False
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
    
def upload_to_google_drive(
    file_path: str, 
    overwrite: bool = False
) -> bool:
    """
    Upload a file (usually a zip) to the metabolomics untargeted outputs google
    drive folder using rclone 
    Check credentials in the rclone config file using /global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone config file
    or /global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone listremotes
    """
    project_folder = os.path.basename(file_path)
    logging.info(tab_print("Uploading zip to Google Drive...", 2))
    if overwrite == True: # Use --update flag in rclone to overwrite remote zip if it's older than the local zip
        # Do the upload
        upload_command = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone copy --update {file_path} ben_lbl_gdrive:/untargeted_outputs'
        try:
            subprocess.check_output(upload_command, shell=True)
        except:
            logging.critical(tab_print("Warning! Google Drive upload failed on upload with overwrite=%s with exception on command %s"%(overwrite, upload_command), 3))
            return False
        # Check that upload worked
        check_upload_command = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone ls ben_lbl_gdrive:/untargeted_outputs --max-depth 1 | grep "{project_folder}"'
        try:
            check_upload_out = subprocess.check_output(check_upload_command, shell=True)
            if check_upload_out.decode('utf-8').strip():
                logging.info(tab_print("Google Drive upload confirmed!", 3))
                return True
            else:
                logging.warning(tab_print("Warning! Google Drive upload check failed. Upload may not have been successful.", 3))
                return False
        except:
            logging.warning(tab_print("Warning! Google Drive upload failed on upload check with overwrite=%s with exception on %s"%(overwrite,check_upload_command), 3))
            return False
    if overwrite == False:
        # Overwrite is False, so first check if the project folder already exists on Google Drive
        check_upload_command = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone ls ben_lbl_gdrive:/untargeted_outputs --max-depth 1 | grep "{project_folder}"'
        try:
            check_upload_out = subprocess.check_output(check_upload_command, shell=True)
            if check_upload_out.decode('utf-8').strip():
                logging.info(tab_print("Note: Overwrite is False and Google Drive folder for %s already exists. Skipping upload..."%(project_folder), 2))
                return False
        except: # Nothing was returned from the check
            upload_command = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone copy --ignore-existing {file_path} ben_lbl_gdrive:/untargeted_outputs'
            try:
                subprocess.check_output(upload_command, shell=True)
                try:
                    check_upload = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone ls ben_lbl_gdrive:/untargeted_outputs --max-depth 1 | grep "{project_folder}"'
                    check_upload_out = subprocess.check_output(check_upload, shell=True)
                    if check_upload_out.decode('utf-8').strip():
                        logging.info(tab_print("Google Drive upload confirmed!", 3))
                        return True
                    else:
                        logging.warning(tab_print("Warning! Google Drive upload check failed. Upload may not have been successful.", 3))
                        return False
                except:
                    logging.warning(tab_print("Warning! Google Drive upload check failed with overwrite=%s with exception on %s"%(overwrite,check_upload_command), 3))
                    return False
            except:
                logging.critical(tab_print("Warning! Google Drive upload command failed with overwrite=%s with exception on %s"%(overwrite,upload_command), 3))
                return False

def submit_quickstart_fbmn(
    params: str = "",
    username: str = ""
) -> dict:
    """
    Submit FBMN jobs by passing parameters
    """
    if not params or not username:
        print("Params and username are required to submit FBMN jobs")
        return
    ##### Pick up secret file and extract password
    with open('/global/cfs/cdirs/metatlas/gnps2/gnps2_%s.txt'%username,'r') as fid:
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

def get_untargeted_status(
    direct_input: str = None,
    print_recent: str = None
) -> None:
    """
    Print the status of a user-defined list of projects
    by calling this function with a csv list of project names
    """
    status_df_list = []
    tab_print("Getting full table from LIMS...", 1)
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    if df.empty:
        tab_print("Could not find untargeted tasks to report! Check that the input argument is a single comma-separated list", 1)
        return
    if not df.empty:
        # Get all projects in gdrive
        all_gdrive_projects_cmd = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone lsl ben_lbl_gdrive:/untargeted_outputs --max-depth 1'
        try:
            all_gdrive_projects = subprocess.check_output(all_gdrive_projects_cmd, shell=True)
            all_gdrive_projects_list = all_gdrive_projects.decode('utf-8').strip().split('\n')
            all_gdrive_projects_df = [row.split() for row in all_gdrive_projects_list]
            all_gdrive_projects_df = pd.DataFrame(all_gdrive_projects_df)
        except:
            tab_print("Warning! Could not get list of projects from Google Drive. Skipping status report...", 1)
            return
        for i,row in df.iterrows():
            project_name = row['parent_dir']
            archive = project_name + ".zip"
            if archive in all_gdrive_projects_df[3].values:
                row['gdrive_status'] = '08 uploaded'
                day = all_gdrive_projects_df[all_gdrive_projects_df[3] == archive][1].values[0]
                time = all_gdrive_projects_df[all_gdrive_projects_df[3] == archive][2].values[0]
                row['upload_date'] = f"{str(day)} {str(time)}"
            else:
                row['gdrive_status'] = '14 not uploaded'
                row['upload_date'] = 'NA'

            status_df_list.append(row[['Modified', 'parent_dir', 'mzmine_pos_status', 'mzmine_neg_status', 'fbmn_pos_status', 'fbmn_neg_status', 'gdrive_status','upload_date']])
    
    total_status_df = pd.DataFrame(status_df_list)
    total_status_df = total_status_df.sort_values(by='Modified', ascending=True)

    if isinstance(print_recent, int):
        print(total_status_df.tail(print_recent).to_csv(sep='\t', index=False))
    else:
        print(total_status_df.to_csv(sep='\t', index=False))

def recursive_chown(
    basepath: Union[str, bytes, os.PathLike],
    group_name: str
) -> None:
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

def download_from_url(
    url: str,
    target_path: Union[str, bytes, os.PathLike]
) -> None:
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

def download_fbmn_results(
    output_dir: str,
    overwrite_fbmn: bool,
    skip_fbmn_download: bool,
    direct_input: Optional[str] = None
) -> None:
    """
    finds complete fbmn tasks (which also have complete mzmine status)
    downloads the graphml and results table files
    renames and moves results to the untargeted_tasks folder

    Overwrite is False (default) when existing GNPS2 files in the untargeted_tasks folder will not be replaced

    Direct_input is None (default) when files from GNPS2 for all available projects will be downloaded. Set direct_input
    to a csv list of project names if you only want to run this function on specific untargeted_tasks
    """
    if skip_fbmn_download:
        logging.info('Skipping Step 6/7: Checking for completed FBMN jobs and downloading results...')
        return
    logging.info('Step 6/7: Checking for completed FBMN jobs and downloading results...')
    tasktype='fbmn'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
        df = subset_df_by_status(df,'mzmine',['07 complete']) # Also want to check that mzmine is complete before downloading fbmn
    status_list = ['07 complete','09 error']
    if direct_input is None:
        df = subset_df_by_status(df,tasktype,status_list)
        df = subset_df_by_status(df,'mzmine',['07 complete']) # Also want to check that mzmine is complete before downloading fbmn
    if df.empty:
        logging.info(tab_print("No completed FBMN data to download!", 1))
    if not df.empty:
        count = 0
        for i,row in df.iterrows():
            project_name = row['parent_dir']
            polarity_list = check_for_polarities(output_dir,project_name)
            if polarity_list is None:
                logging.warning(tab_print("Warning! Project %s does not have a negative or a positive polarity directory. Skipping..."%(project_name), 1))
                continue
            for polarity in polarity_list:
                polarity_short = polarity[:3]
                pathname = os.path.join(output_dir,'%s_%s'%(project_name,polarity))
                fbmn_filename = os.path.join(pathname,'%s_%s_gnps2-fbmn-task.txt'%(project_name,polarity))
                if os.path.isfile(fbmn_filename)==True:
                    with open(fbmn_filename,'r') as fid:
                        taskid = fid.read().split('=')[1].strip()
                    graphml_filename = os.path.join(pathname,'%s_%s_gnps2-fbmn-network.graphml'%(project_name,polarity))
                    results_table_filename = os.path.join(pathname,'%s_%s_gnps2-fbmn-library-results.tsv'%(project_name,polarity))
                    gnps2_link_filename = os.path.join(pathname,'%s_%s_gnps2-page-link.txt'%(project_name,polarity))
                    if overwrite_fbmn==False and os.path.exists(graphml_filename) and os.path.exists(results_table_filename) and os.path.exists(gnps2_link_filename):
                        continue # Skip finished projects unless overwrite is forced
                    # Bail out conditions
                    logging.info(tab_print("Working on %s mode for project %s:"%(polarity,project_name), 1))
                    if row['%s_%s_status'%(tasktype,polarity_short)] == '12 not relevant':
                        logging.info(tab_print("Bailed out because FBMN status for %s mode is '12 not relevant'. Skipping download."%(polarity), 2))
                        continue
                    if row['%s_%s_status'%(tasktype,polarity_short)] != '07 complete':
                        logging.info(tab_print("Bailed out because FBMN status for %s mode is not '07 complete'. Skipping download."%(polarity), 2))
                        continue # skip this polarity even if the other one is finished
                    if row['%s_%s_status'%("mzmine",polarity_short)] != '07 complete':
                        logging.info(tab_print("Bailed out because MZmine status for %s mode is not '07 complete'. Skipping download."%(polarity), 2))
                        continue # skip this polarity even if the other one is finished
                    if row['%s_%s_status'%(tasktype,polarity_short)] == '09 error':
                        logging.info(tab_print("Bailed out because FBMN status for %s mode is '09 error'. Not downloading files."%(polarity), 2))
                        continue

                    logging.info(tab_print("Downloading %s mode FBMN results with task ID %s"%(polarity,taskid), 2))
                    if overwrite_fbmn == True or not os.path.exists(gnps2_link_filename):
                        with open(gnps2_link_filename,'w') as fid:
                            fid.write(f"https://gnps2.org/status?task={taskid}\n")
                            logging.info(tab_print("Wrote GNPS2 link", 3))
                    if overwrite_fbmn == True or not os.path.exists(graphml_filename):
                        graphml_download = download_from_url("https://gnps2.org/result?task=%s&viewname=graphml&resultdisplay_type=task"%(taskid), graphml_filename)
                        if graphml_download:
                            logging.info(tab_print("Downloaded graphml file", 3))
                            count += 1
                        else:
                            logging.critical(tab_print("Error: Failed to download graphml file", 3))
                    if overwrite_fbmn == True or not os.path.exists(results_table_filename):
                        results_table_download = download_from_url("https://gnps2.org/resultfile?task=%s&file=nf_output/library/merged_results_with_gnps.tsv"%(taskid), results_table_filename)
                        if results_table_download:
                            logging.info(tab_print("Downloaded result table file", 3))
                            count += 1
                        else:
                            logging.critical(tab_print("Error: Failed to download results table", 3))
                    try:
                        recursive_chown(pathname, 'metatlas')
                    except:
                        logging.info(tab_print("Note: Could not change group ownership of %s"%(pathname), 3))
        if count > 0:
            logging.info(tab_print("All new FBMN results downloaded.", 1))
        else:
            logging.info(tab_print("No new FBMN results to download!", 1))

def get_recent_mgf_files(
    output_dir:str = None,
    time_back:float = 30
) -> List[str]:
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

def remove_contaminant_from_mgf(
    f: str = None
) -> None:
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
        logging.critical(tab_print("Error: Could not determine polarity of file %s"%(f), 1))
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

def remove_mgf_contaminants(
    output_dir:str = None,
    time_back:float = 30
) -> None:
    """
    This script is called in run_fbmn.py before submitting to GNPS2

    This wrapper script will call the two relevant functions
    """
    recent_files = get_recent_mgf_files(output_dir=output_dir,time_back=time_back)
    for f in recent_files:
        if f:
            remove_contaminant_from_mgf(f)

def get_table_from_lims(
    table: str,
    columns: Optional[List[str]] = None,
    max_rows: Union[int, float] = 1e6
) -> pd.DataFrame:
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

def update_table_in_lims(
    df: pd.DataFrame,
    table: str,
    method: str = 'update',
    max_size: int = 1000,
    pause_time: Optional[float] = None
) -> None:
    """
    Note: Do ~1000 rows at a time.  Any more and you get a 504 error.  Maybe increasing the timeout would help.
    In the header, timeout is a variable that gets set.  Check to see what its set to.  Maybe increasing it would let
    more rows be updated at a time
    
    method can be 'update','insert',or 'delete'

    Use it like this:
    update_table_in_lims(df_lims,'mzml_files')
    whatever is in 'name' or 'Key' will replace whatever used to be there with the other columns
    """
    N = math.ceil(float(df.shape[0]) / max_size)
    for i in range(N):
        start_idx = i * max_size
        end_idx = min((i + 1) * max_size, df.shape[0])
        sub_df = df.iloc[start_idx:end_idx]
        payload = sub_df.to_dict('records')
        if method == 'update':
            api.query.update_rows('lists', table, payload, timeout=10000)
        elif method == 'insert':
            api.query.insert_rows('lists', table, payload, timeout=10000)
        elif method == 'delete':
            api.query.delete_rows('lists', table, payload, timeout=10000)
        else:
            logging.critical(tab_print('ERROR: Nothing to do.  Method %s is not programmed' % (method), 2))
        if pause_time is not None:
            time.sleep(pause_time)

def build_untargeted_filename(
    output_dir: str,
    parent_dir: str,
    polarity: str,
    file_type: Optional[str] = None
) -> str:
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

def check_gnps2_status(taskid:str = None):
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

def mirror_mzmine_results_to_gnps2(
    project: str,
    polarity: str,
    output_dir: str,
    username: str = "bpbowen"
) -> None:
    """
    Mirrors MZmine results to GNPS2.

    This function loads a password from a file and uses it to mirror
    MZmine results to GNPS2 for a given project and polarity.

    Parameters:
    - project (str): The name of the project.
    - polarity (str): The polarity of the results (e.g., 'positive' or 'negative').
    - username (str): The username for GNPS2. Default is 'bpbowen'.
    """
    # Load password from file
    with open('/global/cfs/cdirs/metatlas/gnps2/gnps2_bpbowen.txt', 'r') as f:
        for line in f:
            if line.startswith('MYVARIABLE='):
                password = line.strip().split('=')[1].replace('"', '')
                break

    if not password:
        logging.error(tab_print("Password is required to mirror data to GNPS2. Exiting", 3))
        return "Failed"
    
    logging.info(tab_print("Mirroring MZmine results for %s to GNPS2..."%(project), 3))

    # Suppress all paramiko logs
    paramiko_logger = logging.getLogger("paramiko")
    paramiko_logger.setLevel(logging.CRITICAL)
    paramiko_logger.addHandler(logging.NullHandler())
    paramiko_logger.propagate = False

    project_directory = f"{project}_{polarity}"
    local_directory = os.path.join(output_dir, project_directory)
    remote_directory = f"/untargeted_tasks/{project_directory}"
    remote_host = "sftp.gnps2.org"
    remote_port = 443
    remote_user = username

    try:
        transport = paramiko.Transport((remote_host, remote_port))
        transport.connect(username=remote_user, password=password)
        sftp = paramiko.SFTPClient.from_transport(transport)
    except paramiko.SSHException as e:
        logging.info(tab_print(f"Failed to connect to GNPS2: {e}", 4))
        wait_time = 30
        logging.info(tab_print("Waiting %s seconds and attempting to connect again..."%(wait_time), 5))
        time.sleep(wait_time)
        try:
            transport = paramiko.Transport((remote_host, remote_port))
            transport.connect(username=remote_user, password=password)
            sftp = paramiko.SFTPClient.from_transport(transport)
        except paramiko.SSHException as e:
            logging.error(tab_print(f"Failed to connect to GNPS2 again. Skipping mirror with error: {e}", 4))
            return "Failed"
    
    try:
        sftp.mkdir(remote_directory)
    except Exception as e:
        if "mkdir file exists" not in str(e).lower():
            logging.warning(tab_print(f"Notice! Did not create {remote_directory} at GNPS2 for reason: {e}", 3))

    try:
        local_directory = Path(local_directory)
        for file_path in local_directory.rglob('*'):
            if file_path.is_file() and file_path.suffix in ('.mgf', '.csv', '.tab'):
                #logging.info("Uploading %s to GNPS2..." % file_path.name)
                local_path = str(file_path)
                remote_path = f"{remote_directory}/{file_path.name}"
                sftp.put(local_path, remote_path)
                logging.info(tab_print(f"Uploaded {file_path.name} to GNPS2...", 4))
        sftp.close()
        transport.close()
        logging.info(tab_print(f"Completed MZmine results mirror to GNPS2 for {project}...", 3))
        return "Passed"
    except:
        logging.error(tab_print(f"Failed to mirror MZmine results for {project} to GNPS2", 3))
        return "Failed"

def mirror_raw_data_to_gnps2(
    project: str,
    username: str,
    raw_data_dir: str,
    overwrite_existing_dir: bool = False,
    polarity: Optional[str] = None,
    raw_data_subdir: Optional[str] = None,
    use_polarity_subdir: Optional[bool] = False
) -> None:
    """
    Mirrors raw data (mzML files) to GNPS2.

    This function loads a password from a file and uses it to mirror
    raw data (mzML files) to GNPS2 for a given project and polarity.

    Parameters:
    - project: The name of the project.
    - polarity: The polarity of the data (e.g., 'positive' or 'negative').
    - username: The username for GNPS2. Usually is 'bpbowen'.
    - raw_data_dir: The directory containing the raw data.
    - raw_data_subdir : The subdirectory within the raw data directory. Default is None.
    """
    # Load password from file
    password = None
    try:
        with open('/global/cfs/cdirs/metatlas/gnps2/gnps2_bpbowen.txt', 'r') as f:
            for line in f:
                if line.startswith('MYVARIABLE='):
                    password = line.strip().split('=')[1].replace('"', '')
                    break
    except FileNotFoundError:
        logging.error("Password file not found. Exiting.")
        return "Failed"

    if not password:
        logging.error("Password is required to mirror data to GNPS2. Exiting.")
        return "Failed"
    
    logging.info(tab_print(f"Mirroring raw data (mzML files) for {project} to GNPS2...", 2))

    # Suppress all paramiko logs
    paramiko_logger = logging.getLogger("paramiko")
    paramiko_logger.setLevel(logging.CRITICAL)
    paramiko_logger.addHandler(logging.NullHandler())
    paramiko_logger.propagate = False

    if raw_data_subdir is None: # This means we'll have to try to infer the locations of the mzML files from the project name
        _, validate_department, _ = vfn.field_exists(PurePath(project), field_num=1)
        try:
            if validate_department is None:
                subdir = 'jgi' # Assume a raw data location if project name is not parseable
            else:
                subdir = validate_department.lower()
            if subdir == 'eb':
                subdir = 'egsb'
            check_project_raw_dir = os.path.join(raw_data_dir,subdir,project)
            if os.path.exists(check_project_raw_dir):
                local_directory = check_project_raw_dir
            else:
                possible_subdirs = ["jgi", "egsb", "akuftin", "agolini", "kblouie", "lzhan", "smkosina", "jsjordan", "mtvu", "tharwood", "bpb"]
                for subdir in possible_subdirs:
                    check_project_raw_dir = os.path.join(raw_data_dir, subdir, project)
                    if os.path.exists(check_project_raw_dir):
                        local_directory = check_project_raw_dir
                        break
                else:
                    logging.error(f"Raw data directory for {project} could not be found after trying several possible alternatives. Not mirroring to GNPS2.")
                    return "Failed"
        except Exception as e:
            logging.error(tab_print(f"Error when trying to locate mzML files on disk: {e}", 2))
            return "Failed"
    else:
        local_directory = os.path.join(raw_data_dir,raw_data_subdir,project)
        if not os.path.exists(local_directory):
            logging.error(tab_print(f"Raw data directory for {project} could not be found at user-supplied location of {raw_data_dir}/{raw_data_subdir}. Not mirroring to GNPS2.", 2))
            return "Failed"

    local_directory = Path(local_directory)
    remote_subdir = local_directory.parent.name # Use the same subdir as on perlmutter instead of inferring from the project name directly
    remote_directory = f"/raw_data/{remote_subdir}/{project}"
    if use_polarity_subdir is True:
        if polarity is None:
            logging.error(tab_print("Polarity is required when use_polarity_subdir is True. Exiting.", 2))
            return "Failed"
        else:
            polarity_directory = f"{remote_directory}/{polarity}"
    remote_host = "sftp.gnps2.org"
    remote_port = 443
    remote_user = username

    try:
        transport = paramiko.Transport((remote_host, remote_port))
        transport.connect(username=remote_user, password=password)
        sftp = paramiko.SFTPClient.from_transport(transport)
    except paramiko.SSHException as e:
        logging.info(tab_print(f"Failed to connect to GNPS2: {e}", 4))
        wait_time = 30
        logging.info(tab_print("Waiting %s seconds and attempting to connect again..."%(wait_time), 5))
        time.sleep(wait_time)
        try:
            transport = paramiko.Transport((remote_host, remote_port))
            transport.connect(username=remote_user, password=password)
            sftp = paramiko.SFTPClient.from_transport(transport)
        except paramiko.SSHException as e:
            logging.error(tab_print(f"Failed to connect to GNPS2 again. Skipping mirror with error: {e}", 4))
            return "Failed"

    try:
        sftp.mkdir(remote_directory)
    except Exception as e:
        if "mkdir file exists" in str(e).lower():
            if overwrite_existing_dir is False:
                logging.info(tab_print(f"Found existing directory {remote_directory} at GNPS2 and overwrite_existing_dir is False. Exiting.", 3))
                return "Failed"
            elif overwrite_existing_dir is True:
                logging.info(tab_print(f"Found existing directory {remote_directory} at GNPS2 and overwrite_existing_dir is True. Continuing.", 3))
        else:
            logging.warning(tab_print(f"Notice! Did not create {remote_directory} at GNPS2 for reason: {e}", 3))
            return "Failed"
        
    if use_polarity_subdir is True and polarity is not None:
        polarity_short = f"_{polarity[:3].upper()}_"
        try:
            sftp.mkdir(polarity_directory)
        except Exception as e:
            if "mkdir file exists" in str(e).lower():
                if overwrite_existing_dir is False:
                    logging.info(tab_print(f"Found existing directory {polarity_directory} at GNPS2 and overwrite_existing_dir is False. Exiting.", 3))
                    return "Failed"
            else:
                logging.warning(tab_print(f"Notice! Did not create {polarity_directory} at GNPS2 for reason: {e}", 3))
                return "Failed"
    try:
        for file_path in local_directory.rglob('*'):
            if use_polarity_subdir is True:
                polarity_short = f"_{polarity[:3].upper()}_"
                if file_path.is_file() and file_path.suffix == '.mzML' and polarity_short in file_path.name:
                    #logging.info("Uploading %s to GNPS2..." % file_path.name)
                    local_path = str(file_path)
                    remote_path = f"{polarity_directory}/{file_path.name}"
                    sftp.put(local_path, remote_path)
                    logging.info(tab_print(f"Uploaded {file_path.name} to {remote_path} at GNPS2...", 3))
            if use_polarity_subdir is False:
                if file_path.is_file() and file_path.suffix == '.mzML' in file_path.name:
                    #logging.info("Uploading %s to GNPS2..." % file_path.name)
                    local_path = str(file_path)
                    remote_path = f"{remote_directory}/{file_path.name}"
                    sftp.put(local_path, remote_path)
                    logging.info(tab_print(f"Uploaded {file_path.name} to {remote_path} at GNPS2...", 3))     

        sftp.close()
        transport.close()
        logging.info(tab_print(f"Completed raw data mirror to GNPS2 for {project}...", 2))
        return "Passed"
    except Exception as e:
        logging.error(f"Failed to mirror raw data for {project} to GNPS2: {e}")
        return "Failed"

# def DEPRACATED_check_for_mzmine_files_at_gnps2(project: str, polarity: str, username="bpbowen"):
    
#     ##### Pick up secret file and extract password
#     with open('/global/cfs/cdirs/metatlas/gnps2/gnps2_%s.txt'%username,'r') as fid:
#         t = fid.read()
#     t = t.split('\n')[0]
#     var, pword = t.split('=')
#     password = pword.replace('"', '').replace("'", '').strip()
#     if not password:
#         print("Password is required to submit FBMN jobs")
#         return

#     remote_directory="/untargeted_tasks"
#     remote_host="ftp.gnps2.org"
#     remote_port=6541
#     data_dir = project + "_" + polarity

#     ftp = ftplib.FTP()
#     ftp.connect(host=remote_host,port=remote_port,timeout=120)
#     ftp.login(user=username,passwd=password)
#     ftp.cwd(remote_directory)

#     remote_directories = ftp.nlst()

#     if data_dir in remote_directories:
#         ftp.cwd(data_dir)
#         files = ftp.nlst()
#         required_files = [f"{project}_{polarity}_peak-height.csv",
#                           f"{project}_{polarity}_metadata.tab",
#                           f"{project}_{polarity}.mgf"]
#         all_files_exist = all(file in files for file in required_files)
#         if all_files_exist:
#             return True
#         else:
#             return False
#     else:
#         return False

# def DEPRACATED_sync_mzmine_results_to_gnps2():
#     sync_cmd = '/global/common/software/m2650/infrastructure_automation/gnps2_mirror/sync_untargeted_mzmine_results_to_gnps2.sh'
#     if os.path.exists(sync_cmd):
#         logging.info(tab_print("Syncing MZmine results (.tab, .csv, .mgf) at NERSC with GNPS2 before submitting FBMN job...", 1))
#         try:
#             result = subprocess.run(sync_cmd, shell=True, capture_output=True, text=True)
#             stdout = result.stdout.splitlines()
#             if stdout:
#                 for line in stdout:
#                     logging.info(tab_print(line, 2))
#             stderr = result.stderr.splitlines()
#             if stderr:
#                 for line in stderr:
#                     if "Permission denied" not in line: # Niyogi project has permission denied error
#                         logging.critical(tab_print("Sync script returned errors:", 2))
#                         logging.critical(tab_print(line, 3))
#                         return False
#             logging.info(tab_print("All MZmine results files at NERSC synced with GNPS2.", 2))
#             return True
#         except:
#             logging.critical(tab_print("Warning! Sync command failed.", 2))
#             return False
#     else:
#         logging.critical(tab_print("Warning! Could not find the sync script.", 2))
#         return False

# def DEPRACATED_sync_raw_data_to_gnps2():
#     sync_cmd = '/global/common/software/m2650/infrastructure_automation/gnps2_mirror/sync_untargeted_raw_data_to_gnps2.sh'
#     if os.path.exists(sync_cmd):
#         logging.info(tab_print("Syncing raw data (.mzML) at NERSC with GNPS2 before submitting FBMN job...", 1))
#         try:
#             result = subprocess.run(sync_cmd, shell=True, capture_output=True, text=True)
#             stdout = result.stdout.splitlines()
#             if stdout:
#                 for line in stdout:
#                     logging.info(tab_print(line, 2))
#             stderr = result.stderr.splitlines()
#             if stderr:
#                 for line in stderr:
#                     if "Permission denied" not in line: # Niyogi project has permission denied error
#                         logging.critical(tab_print("Sync script returned errors:", 2))
#                         logging.critical(tab_print(line, 3))
#                         return False
#             logging.info(tab_print("All raw data files at NERSC synced with GNPS2.", 2))
#             return True
#         except:
#             logging.critical(tab_print("Warning! Sync command failed.", 2))
#             return False
#     else:
#         logging.critical(tab_print("Warning! Could not find the sync script.", 2))
#         return False

def submit_fbmn_jobs(
    overwrite_fbmn: bool,
    output_dir: str,
    skip_fbmn_submit: bool,
    skip_mirror_raw_data: bool,
    skip_mirror_mzmine_results: bool,
    direct_input: Optional[str] = None,
    raw_data_dir: Optional[str] = None,
    raw_data_subdir: Optional[str] = None
) -> None:
    """
    finds waiting or errored fbmn tasks (which also have complete mzmine status)
    submits the tasks to GNPS2
    updates the LIMS table with status running if successful

    Direct_input is None (default) when all waiting or errorer tasks will be considered for submission. Set direct_input
    to a csv list of project names if you only want to run this function on specific untargeted_tasks
    """
    if skip_fbmn_submit:
        logging.info('Skipping Step 5/7: Submitting new FBMN jobs to GNPS2...')
        return
    logging.info('Step 5/7: Submitting new FBMN jobs to GNPS2...')
    tasktype = 'fbmn'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    status_list = ['13 waiting','09 error']
    if direct_input is None:
        df = subset_df_by_status(df,tasktype,status_list)
        df = subset_df_by_status(df,'mzmine',['07 complete']) # Also want to check that mzmine is complete before submitting fbmn
        df = subset_df_by_status(df,'mzmine',['09 error'], inverse=True) # Do not submit if mzmine has any error statuses
    if df.empty:
        logging.info(tab_print("No new FBMN jobs to submit!", 1))
        return
    # if df.shape[0] > 20:
    #     logging.info(tab_print('There are too many new projects to be submitted (%s), please check if this is accurate. Exiting script.'%(df.shape[0]), 1))
    #     return
    if not df.empty:
        #logging.info(tab_print("Total of %s projects(s) with FBMN status %s and MZmine status ['07 complete'] to submit to GNPS2:"%(df.shape[0],status_list), 1))
        index_list = []
        for i,row in df.iterrows():
            project_name = row['parent_dir']
            polarity_list = check_for_polarities(row['output_dir'],project_name)
            if polarity_list is None:
                logging.warning(tab_print("Warning! Project %s does not have a negative or a positive polarity directory. Skipping..."%(project_name), 2))
                continue
            for polarity in polarity_list:
                polarity_short = polarity[:3]
                pathname = os.path.join(row['output_dir'],'%s_%s'%(project_name,polarity))
                fbmn_filename = os.path.join(pathname,'%s_%s_gnps2-fbmn-task.txt'%(project_name,polarity))
                
                # Bail out conditions
                if row['%s_%s_status'%(tasktype,polarity_short)] == '07 complete' and row['%s_%s_status'%('mzmine',polarity_short)] == '07 complete' and overwrite_fbmn is False:
                    continue
                logging.info(tab_print("Working on %s mode for project %s:"%(polarity,project_name), 1))
                if row['%s_%s_status'%(tasktype,polarity_short)] == '09 error':
                    if row['num_%s_files'%polarity_short] == 0:
                        logging.info(tab_print("Notice! No raw data files for %s mode. Setting status to '12 not relevant'"%(polarity), 2))
                        df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '12 not relevant'
                        continue
                    if row['num_%s_msms'%polarity_short] < 2:
                        logging.info(tab_print("Notice! Insufficient MSMS detected for %s mode. Setting status to '12 not relevant'"%(polarity), 2))
                        df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '12 not relevant'
                        continue
                    logging.warning(tab_print("Warning! Project %s mode has an error status but appears to have raw data and Mzmine results. Attempting to resubmit..."%(polarity), 2))
                    if os.path.isfile(fbmn_filename):
                        os.remove(fbmn_filename) # Remove failed task ID file in order to submit again
                if row['%s_%s_status'%(tasktype,polarity_short)] == '12 not relevant':
                    logging.info(tab_print("Bailed out because FBMN status is '12 not relevant' for %s mode"%(polarity), 2))
                    continue
                if row['%s_%s_status'%('mzmine',polarity_short)] != '07 complete':
                    logging.info(tab_print("Bailed out because MZmine status is not '07 complete' for %s mode"%(polarity), 2))
                    continue
                if os.path.isfile(fbmn_filename)==True and overwrite_fbmn==False:
                    logging.info(tab_print("Bailed out because FBMN task file already exists for %s mode and overwrite is False"%(polarity), 2))
                    continue

                if raw_data_subdir is None:
                    _, validate_department, _ = vfn.field_exists(PurePath(project_name), field_num=1)
                    try:
                        if validate_department is None:
                            gnps2_subdir = 'jgi' # Assume raw data location if project name is not paresable
                        else:
                            gnps2_subdir = validate_department.lower()
                        if gnps2_subdir == 'eb':
                            gnps2_subdir = 'egsb'
                    except:
                        logging.warning(tab_print("Warning! Could not infer department/raw data location for %s. Defaulting to 'other'. Use --raw_data_subdir to provide a custom subdirectory for the raw data."%(project_name), 2))
                        gnps2_subdir = "other"
                else:
                    gnps2_subdir = raw_data_subdir

                # Get mzmine results files and raw data to GNPS2 before starting FBMN job
                if skip_mirror_mzmine_results is False:
                    logging.info(tab_print("Ensuring MZmine results are at GNPS2 before submitting FBMN job...", 2))
                    mirror = mirror_mzmine_results_to_gnps2(project=project_name,polarity=polarity,output_dir=output_dir,username="bpbowen")
                    if mirror == "Failed":
                        logging.warning(tab_print("Skipping FBMN submission for %s mode because of MZmine results upload failure"%(polarity), 3))
                        continue
                else:
                    logging.info(tab_print("Skipping MZmine results mirroring to GNPS2...", 2))
                if skip_mirror_raw_data is False:
                    logging.info(tab_print("Ensuring raw mzML data are at GNPS2 before submitting FBMN job...", 2))
                    mirror = mirror_raw_data_to_gnps2(project=project_name,polarity=polarity,username="bpbowen",raw_data_dir=raw_data_dir,raw_data_subdir=raw_data_subdir,use_polarity_subdir=False)
                    if mirror == "Failed":
                        logging.info(tab_print("Notice! Proceeding with FBMN submission for %s mode even though raw data mirror failed"%(polarity), 3))
                else:
                    logging.info(tab_print("Skipping raw data mirroring to GNPS2 for %s mode..."%(polarity), 2))

                description = '%s_%s'%(project_name,polarity)
                spectra_file = f'USERUPLOAD/bpbowen/untargeted_tasks/{project_name}_{polarity}/{project_name}_{polarity}.mgf'
                quant_file = f'USERUPLOAD/bpbowen/untargeted_tasks/{project_name}_{polarity}/{project_name}_{polarity}_quant.csv'
                metadata_file = f'USERUPLOAD/bpbowen/untargeted_tasks/{project_name}_{polarity}/{project_name}_{polarity}_metadata.tab'
                raw_data = f'USERUPLOAD/bpbowen/raw_data/{gnps2_subdir}/{project_name}'
                mgf_filename = os.path.join(row['output_dir'],'%s_%s'%(project_name,polarity),'%s_%s.mgf'%(project_name,polarity))
                mgf_lines = count_mgf_lines(mgf_filename)
                if mgf_lines == 0:
                    logging.info(tab_print("Note! MGF file in %s mode has no MSMS data. Updating FBMN status to '12 not relevant'"%(polarity), 2))
                    df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '12 not relevant'
                    index_list.append(i)
                    continue
                params = set_fbmn_parameters(description, quant_file, spectra_file, metadata_file, raw_data)
                job_id = submit_quickstart_fbmn(params, "bpbowen")
                task_list = {'experiment':project_name,'polarity':polarity,'response':job_id}
                logging.info(tab_print("Submitted FBMN job for %s mode and set LIMS status to ['04 running']."%(polarity), 2))
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
            logging.info(tab_print("FBMN submission(s) complete. Use GNPS2 *gnps2-page-link.txt to monitor progress.", 1))

def count_mgf_lines(mgf_file:str = None) -> int:
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

def submit_mzmine_jobs(
    skip_mzmine_submit: bool,
    overwrite_mzmine: bool,
    new_projects: Optional[List[str]] = None,
    direct_input: Optional[str] = None
) -> None:
    """
    finds initiated mzmine tasks
    submits the tasks as mzmine jobs on perlmutter
    updates the LIMS table with status running if successful

    Direct_input is None (default) when all initated tasks will be considered for submission. Set direct_input
    to a csv list of project names if you only want to run this function on specific untargeted_tasks
    """
    if skip_mzmine_submit:
        logging.info('Skipping Step 3/7: Submitting new MZmine jobs...')
        return
    logging.info('Step 3/7: Submitting new MZmine jobs...')
    tasktype = 'mzmine'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    else:
        if new_projects is not None and not new_projects.empty: # this means step 1 returned a dataframe of new projects to run
            df = df[df['parent_dir'].isin(new_projects['parent_dir'])]
        else:
            logging.info(tab_print("No new MZmine jobs initialized to submit!", 1))
            return
    status_list = ['01 initiation','09 error']
    if direct_input is None:
        df = subset_df_by_status(df,tasktype,status_list)
    if df.empty:
        logging.info(tab_print("No new MZmine jobs to submit!", 1))
        return
    # if df.shape[0] > 20:
    #     logging.info(tab_print('There are too many new projects to be submitted (%s), please check if this is accurate. Exiting script.'%(df.shape[0]), 1))
    #     return
    if not df.empty:
        index_list = []
        logging.info(tab_print("Total of %s new MZmine job(s) with status %s to submit"%(df.shape[0],status_list), 1))
        for i,row in df.iterrows():
            project_name = row['parent_dir']
            polarity_list = check_for_polarities(row['output_dir'],project_name)
            if polarity_list is None:
                logging.warning(tab_print("Warning! Project %s does not have a negative or a positive polarity directory. Skipping..."%(project_name), 2))
                continue
            for polarity in polarity_list:
                polarity_short = polarity[:3]
                parent_dir = '%s_%s'%(project_name,polarity)
                pathname = os.path.join(row['output_dir'],parent_dir)
                submission_script_filename = os.path.join(pathname,'%s_mzmine.sh'%(parent_dir))
                mzmine_peak_height_file = os.path.join(row['output_dir'],'%s_%s'%(project_name,polarity),'%s_%s_peak-height.csv'%(project_name,polarity))
                mzmine_mgf_file = os.path.join(row['output_dir'],'%s_%s'%(project_name,polarity),'%s_%s.mgf'%(project_name,polarity))
                
                if row['%s_%s_status'%('mzmine',polarity_short)] == '07 complete' and overwrite_mzmine is False:
                    continue
                if os.path.exists(mzmine_peak_height_file) and os.path.exists(mzmine_mgf_file) and overwrite_mzmine is False:
                    logging.info(tab_print("Note: Overwrite is False and MZmine files for %s in %s mode already exist. Skipping submission..."%(project_name,polarity), 2))
                    continue
                if row['%s_%s_status'%(tasktype,polarity_short)] == '12 not relevant':
                    continue # skip the completed polarity if the other polarity is initiated or errored and needs to be (re)submitted
                if row['%s_%s_status'%(tasktype,polarity_short)] == '09 error':
                    if row['num_%s_files'%polarity_short] == 0:
                        logging.warning(tab_print("Warning! No raw data files for %s mode. Setting status to '12 not relevant'"%(polarity), 2))
                        df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '12 not relevant'
                        df.loc[i,'%s_%s_status'%("fbmn",polarity_short)] = '12 not relevant'
                        continue
                    else:
                        logging.info(tab_print("Note: MZmine task for %s in %s mode has error status but has raw files. Attempting to resubmit with longer runtime."%(project_name,polarity), 2))
                        batch_submission_file = os.path.join(pathname,'%s_mzmine-sbatch.sbatch'%(parent_dir))
                        if os.path.isfile(batch_submission_file):
                            with open(batch_submission_file, 'r') as file:
                                data = file.read()
                            data = data.replace('-t 3:00:00', '-t 6:00:00')
                            with open(batch_submission_file, 'w') as file:
                                file.write(data)

                if os.path.isfile(submission_script_filename)==True:
                    with open(submission_script_filename,'r') as fid:
                        logging.info(tab_print("Submitting %s mode mzmine job for project %s"%(polarity, project_name), 2))
                        cmd = fid.read()
                        sbatch_output = subprocess.run(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE)
                        sbatch_output_str = sbatch_output.stdout.decode('utf-8').replace('\n', '')
                        logging.info(tab_print("%s"%(sbatch_output_str), 3))
                        job_id = sbatch_output_str.split()[-1]
                        job_id_filename = os.path.join(pathname,'%s_mzmine-job-id.txt'%(parent_dir))
                        with open(job_id_filename,'w') as fid:
                            fid.write("%s=%s\n"%(parent_dir,job_id))
                            logging.info(tab_print("Wrote job ID to file and setting status to ['04 running].", 3))
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
            logging.info(tab_print("MZmine submission(s) complete. Use squeue to monitor progress.", 1))
    
def set_fbmn_parameters(
    description: str,
    quant_file: str,
    spectra_file: str,
    metadata_file: str,
    raw_data: Dict
) -> None:
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

def update_mzmine_status_in_untargeted_tasks(
    background_designator: List[str],
    skip_mzmine_status: bool,
    direct_input: Optional[str] = None,
    background_ratio: int = 5,
    zero_value: float = (2/3),
    nonpolar_solvent_front: float = 0.5,
    polar_solvent_front: float = 0.8,
    nonpolar_solvent_end: float = 10,
    polar_solvent_end: float = 17.5
) -> None:
    """    
    finds initiated or running mzmine tasks
    checks if they have valid output
    updates the LIMS table with status complete if successful

    Direct_input is None (default) when all initated or running tasks haven't produced output yet. Set direct_input
    to a csv list of project names if you only want to run this function on specific untargeted_tasks
    """
    if skip_mzmine_status:
        logging.info('Skipping Step 2/7: Checking and updating status of MZmine jobs in LIMS...')
        return
    logging.info('Step 2/7: Checking and updating status of MZmine jobs in LIMS...')
    tasktype = 'mzmine'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    status_list = ['01 initiation','04 running','09 error']
    if direct_input is None:
        df = subset_df_by_status(df,tasktype,status_list)
    if df.empty:
        logging.info(tab_print("No MZmine jobs to update!", 1))
    if not df.empty:
        #logging.info(tab_print("Total of %s project(s) with MZmine status %s to attempt to update"%(df.shape[0],status_list), 1))
        index_list = []
        for i,row in df.iterrows():
            project_name = row['parent_dir']
            polarity_list = check_for_polarities(row['output_dir'],project_name)
            if polarity_list is None:
                logging.warning(tab_print("Warning! Project %s does not have a negative or a positive polarity directory. Skipping..."%(project_name), 2))
                continue
            for polarity in polarity_list:
                polarity_short = polarity[:3]
                pathname = os.path.join(row['output_dir'],'%s_%s'%(project_name,polarity))
                parent_dir = '%s_%s'%(project_name,polarity)
                job_id_filename = os.path.join(pathname,'%s_mzmine-job-id.txt'%(parent_dir))
                old_peakheight_filename = os.path.join(pathname,'%s_%s_MSMS_quant.csv'%(project_name,polarity))
                peakheight_filename = os.path.join(pathname,'%s_%s_peak-height.csv'%(project_name,polarity))
                mgf_filename = os.path.join(pathname,'%s_%s.mgf'%(project_name,polarity))
                metadata_filename = os.path.join(pathname,'%s_%s_metadata.tab'%(project_name,polarity))
                mzmine_output_filename = os.path.join(pathname,'%s_%s-mzmine.out'%(project_name,polarity))

                if direct_input is None:
                    if row['%s_%s_status'%(tasktype,polarity_short)] == '12 not relevant' or row['%s_%s_status'%(tasktype,polarity_short)] == '07 complete':
                        continue  ## This helps when one polarity is complete but the other is not, so skipped the finished one

                if os.path.isfile(mgf_filename) and os.path.isfile(metadata_filename) and \
                     os.path.isfile(mzmine_output_filename) and (os.path.isfile(peakheight_filename) or os.path.isfile(old_peakheight_filename)):
                    # MZmine is finished and status should be updated
                    logging.info(tab_print("Working on %s in %s mode"%(project_name,polarity), 2))
                    logging.info(tab_print("All MZmine output files found in %s directory, continuing..."%(polarity), 3))
                    logging.info(tab_print("Calculating feature and background counts and updating LIMS table", 3))
                    feature_count, exctrl_count, msms_count = calculate_counts_for_lims_table(peakheight_filename,mgf_filename,background_designator=background_designator)
                    df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '07 complete'
                    df.loc[i,'num_%s_features'%(polarity_short)] = int(feature_count)
                    df.loc[i,'num_%s_features_in_exctrl'%(polarity_short)] = int(exctrl_count)
                    df.loc[i,'num_%s_msms'%(polarity_short)] = int(msms_count)
                    index_list.append(i)

                    if feature_count > 0:
                        logging.info(tab_print("Filtering peak height file to remove features below background...", 4))
                        create_filtered_peakheight_file(peakheight_filename=peakheight_filename,background_designator=background_designator,output_dir=pathname, \
                                                        project_name=project_name,polarity=polarity,background_ratio=background_ratio,zero_value=zero_value, \
                                                        nonpolar_solvent_front=nonpolar_solvent_front,polar_solvent_front=polar_solvent_front,
                                                        nonpolar_solvent_end=nonpolar_solvent_end,polar_solvent_end=polar_solvent_end)
                    else:
                        logging.warning(tab_print("Warning! No features were found, not writing a filtered peak height file.", 4))

                    try:
                        recursive_chown(pathname, 'metatlas')
                    except:
                        logging.info(tab_print("Note: Could not change group ownership of %s"%(pathname), 2))

                else:
                    if row['%s_%s_status'%(tasktype,polarity_short)] == '04 running': # verify that it's actually still running, else change to error
                        if os.path.isfile(job_id_filename):
                            with open(job_id_filename,'r') as fid:
                                job_id = fid.read().strip().split('=')[-1]
                                sqs_cmd = 'squeue --job %s | wc -l'%job_id
                                sqs_output = subprocess.run(sqs_cmd, shell=True, capture_output=True, text=True)
                                if int(sqs_output.stdout.strip()) != 2: # job is not running
                                    logging.critical(tab_print("Warning! MZmine status for %s mode of %s is ['04 running'] but no job is running at NERSC. Changing to ['09 error']."%(polarity, project_name), 1))
                                    df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '09 error'    
                                    index_list.append(i)
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
        else:
            logging.info(tab_print("No MZmine jobs to update!", 1))

def create_filtered_peakheight_file(
    peakheight_filename: str,
    background_designator: List[str],
    output_dir: str,
    project_name: str,
    polarity: str,
    background_ratio: int = 3,
    zero_value: float = (2/3),
    nonpolar_solvent_front: float = 0.5,
    polar_solvent_front: float = 0.8,
    nonpolar_solvent_end: float = 10,
    polar_solvent_end: float = 17.5
) -> None:
    """
    Accepts a peakheight file and filters out features that have a max peak height in exp samples that is less than
    3 times the peak height in the control samples. Does additional filtering and reporting for EGSB
    """
    if "JGI" not in project_name.split('_')[1]:
        jgi_project = False
    else:
        jgi_project = True

    def print_empty_filter_files(data_table, filtered_filenames, jgi_project):
        empty_data_table = pd.DataFrame(columns=data_table.columns[1:])
        for filename in filtered_filenames.values():
            if jgi_project:
                if 'x_exctrl' in filename:
                    empty_data_table.to_csv(filename, index=False)
            else:
                empty_data_table.to_csv(filename, index=False)
        return

    # Read in data and create tables
    data_table = pd.read_csv(peakheight_filename, sep=',')
    filtered_filenames = {
        f'>{background_ratio}x_exctrl': os.path.join(output_dir, f'{project_name}_{polarity}_peak-height-filtered-{background_ratio}x-exctrl.csv'),
        '>1e4_exctrl': os.path.join(output_dir, f'{project_name}_{polarity}_peak-height-filtered-1e4-exctrl.csv'),
        '>1e5_exctrl': os.path.join(output_dir, f'{project_name}_{polarity}_peak-height-filtered-1e5-exctrl.csv'),
        f'>{background_ratio}x_txctrl': os.path.join(output_dir, f'{project_name}_{polarity}_peak-height-filtered-{background_ratio}x-txctrl.csv'),
        '>1e4_txctrl': os.path.join(output_dir, f'{project_name}_{polarity}_peak-height-filtered-1e4-txctrl.csv'),
        '>1e5_txctrl': os.path.join(output_dir, f'{project_name}_{polarity}_peak-height-filtered-1e5-txctrl.csv')
        }
    
    # Get background (extraction control) and sample columns
    header = data_table.columns
    exctrl_columns = [col for col in header if len(col.split('_')) > 12 and any(designator.lower() in col.split('_')[12].lower() for designator in background_designator) and 'mzml' in col.lower()]
    txctrl_columns = [col for col in header if len(col.split('_')) > 12 and 'txctrl' in col.split('_')[12].lower() and 'mzml' in col.lower()]
    sample_columns = [col for col in header if col not in exctrl_columns and col not in txctrl_columns]
    if not exctrl_columns and not txctrl_columns:
        logging.warning(tab_print(f"Warning! No columns found in peak height table with designation {background_designator} or 'TxCtrl'. Returning empty dataframes.", 5))
        print_empty_filter_files(data_table, filtered_filenames, jgi_project)
        return
    if not sample_columns:
        logging.warning(tab_print("Warning! No sample columns found in peak height table. Returning empty dataframes.", 5))
        print_empty_filter_files(data_table, filtered_filenames, jgi_project)
        return

    # Determine solvent columns
    if "_HILIC" in peakheight_filename:
        solvent_front = polar_solvent_front
        solvent_end = polar_solvent_end
    elif "_C18" in peakheight_filename:
        solvent_front = nonpolar_solvent_front
        solvent_end = nonpolar_solvent_end
    else:
        logging.warning(tab_print("Warning! Could not infer solvent front and end for filtering, check peak height files and project names. Returning empty dataframes.", 5))
        print_empty_filter_files(data_table, filtered_filenames, jgi_project)
        return
    column_limits = f'rt_in_{solvent_front}-{solvent_end}min_window'

    # Run the filtering and create summary rows
    summary_table = []
    peak_height_filtered_files = {key: [] for key in filtered_filenames.keys()}
    for i, row in data_table.iterrows():
        # Set up summary row and filtering parameters
        exctrl_max = row[exctrl_columns].max() if exctrl_columns else 1e20 # Set to a large number if no exctrl columns
        txctrl_max = row[txctrl_columns].max() if txctrl_columns else 1e20 # Set to a large number if no txctrl columns
        sample_max = row[sample_columns].max()
        rt_peak = row['row retention time']
        row_id = int(row['row ID'])
        summary_row = {
            'row ID': row_id,
            'exctrl_max': exctrl_max,
            'txctrl_max': txctrl_max,
            'sample_max': sample_max,
            'rt_peak': rt_peak
        }

        # Perform the solvent filtering
        summary_row[column_limits] = '1' if rt_peak >= solvent_front and rt_peak <= solvent_end else '0'
        
        # Perform peak height filtering
        peak_height_filters = {
            f'>{background_ratio}x_exctrl': exctrl_max * background_ratio,
            '>1e4_exctrl': exctrl_max + 1e4,
            '>1e5_exctrl': exctrl_max + 1e5,
            f'>{background_ratio}x_txctrl': txctrl_max * background_ratio,
            '>1e4_txctrl': txctrl_max + 1e4,
            '>1e5_txctrl': txctrl_max + 1e5
            }
        for filter_name, threshold in peak_height_filters.items():
            summary_row[filter_name] = '1' if sample_max >= threshold else '0'
            if summary_row[filter_name] == '1' and summary_row[column_limits] == '1':
                peak_height_filtered_files[filter_name].append(i)
        
        # Append the summary row to the summary table
        summary_table.append(summary_row)

    # Convert summary table to DataFrame
    summary_columns = ['row ID', 'exctrl_max', 'txctrl_max', 'sample_max', 'rt_peak', column_limits] + list(filtered_filenames.keys())
    summary_table = pd.DataFrame(summary_table, columns=summary_columns)
    if not exctrl_columns:
        summary_table['exctrl_max'] = np.nan # Replace large value with NA
    if not txctrl_columns:
        summary_table['txctrl_max'] = np.nan # Replace large value with NA
    if jgi_project:
        columns_to_drop = ['txctrl_max', '>1e4_exctrl', '>1e5_exctrl', f'>{background_ratio}x_txctrl', '>1e4_txctrl', '>1e5_txctrl']
        summary_table.drop(columns=[col for col in columns_to_drop if col in summary_table.columns], inplace=True)
    summary_filename = os.path.join(output_dir, f'{project_name}_{polarity}_peak-height-filtering-summary.csv')
    summary_table.to_csv(summary_filename, index=False)

    # Write the filtered peak height files
    for filter_file, indices in peak_height_filtered_files.items():
        if jgi_project:
            if filter_file in ['>1e4_exctrl', '>1e5_exctrl', f'>{background_ratio}x_txctrl', '>1e4_txctrl', '>1e5_txctrl']:
                continue
        if len(indices) == 0:  # only the header row was retained for this particular filtered file (no features retained)
            print_empty_filter_files(data_table, {filter_file: filtered_filenames[filter_file]}, jgi_project)
            continue
        else:
            data_table_filtered = data_table.loc[sorted(set([0] + indices))]

            # Replace zeros with a small value for JGI projects
            if filter_file == f'>{background_ratio}x_exctrl':
                if jgi_project:
                    values = data_table_filtered[sample_columns + exctrl_columns]
                    lowest_non_zero = values[values != 0].min().min()
                    new_value = lowest_non_zero * zero_value
                    data_table_filtered.replace(0, new_value, inplace=True)

            # Log the number of features removed and save file
            difference = data_table.shape[0] - data_table_filtered.shape[0]
            logging.info(tab_print(f"{difference} features removed by the {filter_file} background filter; {data_table_filtered.shape[0]} features retained. Saving filtered peak height table.", 5))
            data_table_filtered = data_table_filtered.drop(data_table_filtered.filter(regex="^Unnamed").columns, axis=1)
            data_table_filtered.to_csv(filtered_filenames[filter_file], index=False)

def calculate_counts_for_lims_table(
    peakheight_filename: str,
    mgf_filename: str,
    background_designator: List[str]
) -> None:
    """
    Accepts a peakheight file and mgf file and returns the number of features
    in exp and and control samples found in the peakheightfile plus number of 
    feature ids in mgf file.
    """
    data_table = pd.read_csv(peakheight_filename, sep=',')
    exctrl_count = 0
    feature_count = 0
    msms_count = 0
    exctrl_columns = [col for col in data_table.head(1) if len(col.split('_')) > 12 and any(designator.lower() in col.split('_')[12].lower() for designator in background_designator) and 'mzml' in col.lower()]
    sample_columns = [col for col in data_table.head(1) if len(col.split('_')) > 12 and any(designator.lower() not in col.split('_')[12].lower() for designator in background_designator) and 'mzml' in col.lower()]

    # Find the background designator string (default: ['ExCtrl']) in the 13th element of each project name in the peak height table
    if exctrl_columns:
        exctrl_rows = data_table[data_table[exctrl_columns].gt(0).any(axis=1)]
        exctrl_count = int(exctrl_rows.shape[0])
        logging.info(tab_print("%s features in control samples"%(exctrl_count), 4))

    if sample_columns:
        feature_rows = data_table[data_table[sample_columns].gt(0).any(axis=1)]
        feature_count = int(feature_rows.shape[0])
        logging.info(tab_print("%s features in experimental samples"%(feature_count), 4))

    if feature_count <= 1: # Exclude the header line
        logging.warning(tab_print("Warning! No features found in peak height file", 4))

    if os.path.isfile(mgf_filename):
        cmd = 'cat %s | grep FEATURE_ID| wc -l'%mgf_filename
        try:
            result = subprocess.run(cmd,stdout=subprocess.PIPE, shell=True)
            if result:
                msms_count = int(result.stdout.strip())
                logging.info(tab_print("%s msms features in mgf file"%(msms_count), 4))
            else:
                logging.warning(tab_print("Warning! Could not count msms features in mgf file", 4))
        except:
            logging.warning(tab_print("Warning! Could not count msms features in mgf file", 4))

    return feature_count, exctrl_count, msms_count

def update_fbmn_status_in_untargeted_tasks(
    skip_fbmn_status: bool,
    direct_input: Optional[str] = None
) -> None:
    """    
    finds running or waiting fbmn tasks
    checks if they have valid output
    updates the LIMS table with status complete if successful

    Direct_input is None (default) when all running or waiting tasks haven't produced output yet. Set direct_input
    to a csv list of project names if you only want to run this function on specific untargeted_tasks
    """
    if skip_fbmn_status:
        logging.info('Skipping Step 4/7: Checking and updating status of FBMN jobs in LIMS...')
        return
    logging.info('Step 4/7: Checking and updating status of FBMN jobs in LIMS...')
    tasktype='fbmn'
    df = get_table_from_lims('untargeted_tasks')
    df = filter_common_bad_project_names(df)
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
    status_list = ['04 running','13 waiting','09 error']
    if direct_input is None:
        df = subset_df_by_status(df,tasktype,status_list)
    if df.empty:
        logging.info(tab_print("No FBMN jobs to update!", 1))
    if not df.empty:
        index_list = []
        #logging.info(tab_print("Total of %s project(s) with FBMN status %s to attempt to update"%(df.shape[0],status_list), 1))
        for i,row in df.iterrows():
            polarity_list = check_for_polarities(row['output_dir'],row['parent_dir'])
            if polarity_list is None:
                logging.warning(tab_print("Warning! Project %s does not have a negative or a positive polarity directory. Skipping..."%(row['parent_dir']), 1))
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
                            logging.info(tab_print("Updating FBMN job status for %s in %s mode to complete"%(row['parent_dir'],polarity), 1))
                            df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '07 complete'
                    elif status=='FAILED':
                        if df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] != '09 error':
                            logging.info(tab_print("Updating FBMN job status for %s in %s mode to error"%(row['parent_dir'],polarity), 1))
                            df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '09 error'
                    elif status in ['RUNNING','QUEUED']:
                        if df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] != '04 running':
                            logging.info(tab_print("Updating FBMN job status for %s in %s mode to running"%(row['parent_dir'],polarity), 1))
                            df.loc[i,'%s_%s_status'%(tasktype,polarity_short)] = '04 running'
                    elif status =='N/A':
                        logging.critical(tab_print("Warning! Could not get status from GNPS2 for %s in %s mode"%(row['parent_dir'],polarity), 1))
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
        else:
            logging.info(tab_print("No FBMN jobs to update!", 1))

def write_mzmine_sbatch_and_runner(
    basepath: Union[str, bytes, os.PathLike],
    batch_filename: str,
    parent_dir: Union[str, bytes, os.PathLike],
    filelist_filename: str
) -> None:
    """
    Write the sbatch and runner files for mzmine submission via slurm
    """
    mzmine_launcher = f'{BINARY_PATH}/MZmine-3.7.2/bin/MZmine -threads auto -t /tmp'
    #mzmine_launcher = f'{BINARY_PATH}/MZmine-4.4.3/bin/mzmine -threads auto -t /tmp'

    sbatch_filename = '%s_mzmine-sbatch.sbatch'%os.path.join(basepath,parent_dir)
    runner_filename = '%s_mzmine.sh'%os.path.join(basepath,parent_dir)
    s = '%s -batch %s -input %s\necho MZMINE IS DONE\n\n\n'%(mzmine_launcher,batch_filename,filelist_filename)

    with open(sbatch_filename,'w') as fid:
        fid.write('%s\n%s\n'%(SLURM_PERLMUTTER_HEADER.replace('slurm','%s-%s'%(os.path.join(basepath,parent_dir),'mzmine')),s))
    with open(runner_filename,'w') as fid:
        fid.write('sbatch %s'%sbatch_filename)

def metadata_file_filter(data, polarity, skip_blank_filter=False, fps_files_only=False, nonstandard_filename=False):
    if nonstandard_filename and skip_blank_filter and fps_files_only:
        return data['basename'][data['basename'].apply(
                        lambda x:
                        'QC' not in x and
                        'InjBl' not in x and
                        'ISTD' not in x)].to_list()
    if nonstandard_filename and skip_blank_filter and not fps_files_only:
        return data['basename'][data['basename'].apply(
                        lambda x:
                        "_"+polarity+"_" in x and
                        'QC' not in x and
                        'InjBl' not in x and
                        'ISTD' not in x)].to_list()
    if nonstandard_filename and not skip_blank_filter and fps_files_only:
        return data['basename'][data['basename'].apply(
                        lambda x:
                        'Blank' not in x and
                        'QC' not in x and
                        'InjBl' not in x and
                        'ISTD' not in x)].to_list()
    if nonstandard_filename and not skip_blank_filter and not fps_files_only:
        return data['basename'][data['basename'].apply(
                        lambda x:
                        'Blank' not in x and
                        "_"+polarity+"_" in x and
                        'QC' not in x and
                        'InjBl' not in x and
                        'ISTD' not in x)].to_list()
    if skip_blank_filter and not nonstandard_filename and fps_files_only:
        return data['basename'][data['basename'].apply(
                        lambda x:
                        len(x.split('_')) > 9 and
                        len(x.split('_')) > 12 and
                        'QC' not in x.split('_')[12] and
                        'InjBl' not in x.split('_')[12] and
                        'ISTD' not in x.split('_')[12])].to_list()
    if skip_blank_filter and not nonstandard_filename and not fps_files_only:
        return data['basename'][data['basename'].apply(
                        lambda x:
                        len(x.split('_')) > 9 and
                        polarity in x.split('_')[9] and
                        len(x.split('_')) > 12 and
                        'QC' not in x.split('_')[12] and
                        'InjBl' not in x.split('_')[12] and
                        'ISTD' not in x.split('_')[12])].to_list()
    if not skip_blank_filter and not nonstandard_filename and fps_files_only:
        return data['basename'][data['basename'].apply(
                        lambda x:
                        'Blank' not in x and
                        len(x.split('_')) > 9 and
                        len(x.split('_')) > 12 and
                        'QC' not in x.split('_')[12] and
                        'InjBl' not in x.split('_')[12] and
                        'ISTD' not in x.split('_')[12])].to_list()
    if not skip_blank_filter and not nonstandard_filename and not fps_files_only:
        return data['basename'][data['basename'].apply(
                        lambda x:
                        'Blank' not in x and
                        len(x.split('_')) > 9 and
                        polarity in x.split('_')[9] and
                        len(x.split('_')) > 12 and
                        'QC' not in x.split('_')[12] and
                        'InjBl' not in x.split('_')[12] and
                        'ISTD' not in x.split('_')[12])].to_list()
    
def write_metadata_per_new_project(
    df: pd.DataFrame,
    background_designator: List[str],
    validate_names: bool,
    raw_data_dir: str,
    raw_data_subdir: Optional[str] = None,
    skip_blank_filter: Optional[bool] = False,
    fps_files_only: Optional[bool] = False
) -> List:
    """
    Takes a LIMS table (usually raw data from lcmsruns_plus) and creates
    mzmine metadata for writing a new untargeted_tasks directory
    """
    new_project_list = []
    # Iterate through each unique parent_dir
    for parent_dir in df['parent_dir'].unique():
        # Filter the DataFrame for the current parent_dir
        df_filtered = df[df['parent_dir'] == parent_dir]

        if raw_data_subdir is None: # This means we'll have to try to infer the locations of the mzML files from the project name
            _, validate_department, _ = vfn.field_exists(PurePath(parent_dir), field_num=1)
            try:
                if validate_department is None:
                    subdir = 'jgi' # Assume raw data location if project name is not parseable
                else:
                    subdir = validate_department.lower()
                if subdir == 'eb':
                    subdir = 'egsb'
                check_project_raw_dir = os.path.join(raw_data_dir,subdir,parent_dir)
                if os.path.exists(check_project_raw_dir):
                    full_mzml_path = check_project_raw_dir
                else:
                    possible_subdirs = ["jgi", "egsb", "akuftin", "agolini", "kblouie", "lzhan", "smkosina", "jsjordan", "mtvu", "tharwood", "bpb"]
                    for subdir in possible_subdirs:
                        check_project_raw_dir = os.path.join(raw_data_dir, subdir, parent_dir)
                        if os.path.exists(check_project_raw_dir):
                            full_mzml_path = check_project_raw_dir
                            break
                    else:
                        logging.error(f"Raw data directory for {parent_dir} could not be found after trying several possible alternatives. Not creating metadata.")
                        continue
            except Exception as e:
                logging.error(tab_print(f"Error when trying to locate mzML files on disk: {e}", 2))
                continue
        else:
            full_mzml_path = os.path.join(raw_data_dir,raw_data_subdir,parent_dir)
            if not os.path.exists(full_mzml_path):
                logging.error(tab_print(f"Raw data directory for {parent_dir} could not be found at user-supplied location of {raw_data_dir}/{raw_data_subdir}. Not creating metadata.", 2))
                continue
        
        # Initialize info dict
        project_dict = {'parent_dir': parent_dir}
        project_dict['positive'] = {}
        project_dict['negative'] = {}

        # Determine which polarities need metadata
        try:
            positive_file_subset = metadata_file_filter(df_filtered, polarity="POS", skip_blank_filter=skip_blank_filter, fps_files_only=fps_files_only, nonstandard_filename=False)
        except IndexError:
            positive_file_subset = metadata_file_filter(df_filtered, polarity="POS", skip_blank_filter=skip_blank_filter, fps_files_only=fps_files_only, nonstandard_filename=True)
        if positive_file_subset:
            positive_file_list = [os.path.join(full_mzml_path, file) for file in positive_file_subset]
            positive_file_count = len(positive_file_list)
            if validate_names is True:
                file_name_validation = [vfn.validate_file_name(PurePath(file),minimal=True,print_logger=False) for file in positive_file_list]
                try:
                    if all(file_name_valid == True for file_name_valid in file_name_validation) is False:
                        logging.warning(tab_print("Warning! Project %s has one or more mzML files that do not have valid format. Not adding project to untargeted tasks..."%(parent_dir), 2))
                        continue
                    project_name_validation = [vfn.parent_dir_num_fields(PurePath(file)) for file in positive_file_list]
                    if all(len(project_name_valid) == 0 for project_name_valid in project_name_validation) is False: # Prints empty list if project name valid
                        logging.warning(tab_print("Warning! Parent dir %s does not have valid format. Not adding project to untargeted tasks..."%(parent_dir), 2))
                        continue
                except:
                    logging.warning(tab_print("mzML file name validation failed", 2))
        else:
            positive_file_count = 0

        try:
            negative_file_subset = metadata_file_filter(df_filtered, polarity="NEG", skip_blank_filter=skip_blank_filter, fps_files_only=fps_files_only, nonstandard_filename=False)
        except IndexError:
            negative_file_subset = metadata_file_filter(df_filtered, polarity="NEG", skip_blank_filter=skip_blank_filter, fps_files_only=fps_files_only, nonstandard_filename=True)
        if negative_file_subset:
            negative_file_list = [os.path.join(full_mzml_path, file) for file in negative_file_subset]
            negative_file_count = len(negative_file_list)
            if validate_names is True:
                file_name_validation = [vfn.validate_file_name(PurePath(file),minimal=True,print_logger=False) for file in negative_file_list]
                try:
                    if all(file_name_valid == True for file_name_valid in file_name_validation) is False:
                        logging.warning(tab_print("Warning! Project %s has one or more mzML files that do not have valid format. Not adding project to untargeted tasks..."%(parent_dir), 2))
                        continue
                    project_name_validation = [vfn.parent_dir_num_fields(PurePath(file)) for file in negative_file_list]
                    if all(len(project_name_valid) == 0 for project_name_valid in project_name_validation) is False: # Prints empty list if project name valid
                        logging.warning(tab_print("Warning! Parent dir %s does not have valid format. Not adding project to untargeted tasks..."%(parent_dir), 2))
                        continue
                except:
                    logging.warning(tab_print("mzML file name validation failed", 2))
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
            ugroups = [g for g in ugroups if any(c.lower() in g for c in background_designator)]
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

def update_new_untargeted_tasks(
    background_designator: List[str],
    validate_names: bool,
    skip_sync: bool,
    output_dir: str,
    raw_data_dir: str,
    direct_input: Optional[str] = None,
    custom_mzmine_batch_params: Optional[str] = None,
    raw_data_subdir: Optional[str] = None,
    skip_blank_filter: Optional[bool] = False,
    fps_files_only: Optional[bool] = False
) -> None:
    """
    This script is called by run_mzmine.py before the untargeted pipeline kicks off

    given all directories that are in the table for potential untargeted tasks
    and all the directories in the raw data folders
    return all the directories in the raw data folders
    that are not in the untargeted tasks
    
    The strip command is because there is one folder that ends in a 
    space and labkey doesn't allow this
    """
    if skip_sync:
        logging.info('Skipping Step 1/7: Syncing LIMS and NERSC to identify new projects with raw data that are not yet in the untargeted task list...')    
        return
    logging.info('Step 1/7: Syncing LIMS and NERSC to identify new projects with raw data that are not yet in the untargeted task list...')    
    # Get lcmsrun table and subset
    df = get_table_from_lims('lcmsrun_plus')
    if direct_input is not None:
        df = df[df['parent_dir'].isin(direct_input)]
        if df.empty:
            logging.info(tab_print("Projects in direct input (%s) are not in LIMS lcmsruns table"%(direct_input), 1))
            return
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
    dirs_with_m2_files = df.groupby('parent_dir')['basename'].apply(lambda x: any(('MS2' in filename or 'MSMS' in filename) and '_MS1_' not in filename for filename in x))
    dirs_with_m2_files = dirs_with_m2_files[dirs_with_m2_files].index.tolist()
    # Print parent_dirs with files less than 3 hours old
    if time_new_folders:
        logging.info(tab_print("Note: There are data directories with files less than 3 hours old! Skipping these for now:", 2))
        for folder in time_new_folders:
            logging.info(tab_print(folder, 3))

    # Get all folders in raw data table
    if fps_files_only:
        all_folders = df.loc[df['polarity'].isin(['FPS']),'parent_dir'].unique()
        all_folders = [a.strip() for a in all_folders]
    else:
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
    new_project_info_list = write_metadata_per_new_project(df=df_new,background_designator=background_designator,raw_data_dir=raw_data_dir, fps_files_only=fps_files_only, \
                                                           skip_blank_filter=skip_blank_filter, raw_data_subdir=raw_data_subdir, validate_names=validate_names)
    new_project_info_list_subset = [d for d in new_project_info_list if d.get('polarities') is not None]

    # Check for nountargeted.txt file in new projects dirs and skip them
    new_project_info_list_subset = [
        d for d in new_project_info_list_subset 
        if not glob.glob(os.path.join(raw_data_dir, '*', d['parent_dir'], 'nountargeted.txt'))
    ]
    if not new_project_info_list_subset:
        logging.info(tab_print("No new projects to add to untargeted tasks!", 2))
        return None
    logging.info(tab_print("New projects to add to untargeted tasks:", 2))
    print_new_projects = "\n".join(["\t"*7 + d['parent_dir'] for d in new_project_info_list_subset])
    logging.info(tab_print("\n" + print_new_projects, 0))

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
            _, validate_machine_name, _ = vfn.field_exists(PurePath(project_name), field_num=6)
            logging.info(tab_print("Inferred machine name: %s"%(validate_machine_name), 2))
            if custom_mzmine_batch_params is None: # When there is not a custom input
                if validate_machine_name is None:  # Assume more lenient parameters if machine name is not going to be validated
                    logging.warning(tab_print("Warning! Could not validate machine name. Using lenient (IQX) MZmine parameters...", 2))
                    mzmine_running_parameters = mzine_batch_params_file_iqx
                    mzmine_parameter = 5
                elif any(substring in validate_machine_name.lower() for substring in ("iqx", "idx")):
                    mzmine_running_parameters = mzine_batch_params_file_iqx
                    mzmine_parameter = 5
                elif any(substring in validate_machine_name.lower() for substring in ("exp", "exploris", "qe")):
                    mzmine_running_parameters = mzine_batch_params_file
                    mzmine_parameter = 2
                else:  # Assume more lenient parameters if machine name cannot be validated
                    mzmine_running_parameters = mzine_batch_params_file_iqx
                    mzmine_parameter = 5
                logging.info(tab_print("Using MZmine parameters: %s"%(os.path.basename(mzmine_running_parameters)), 2))
                lims_untargeted_table_updater['mzmine_parameter_sheet'] = mzmine_running_parameters
                lims_untargeted_table_updater['mzmine_parameter_row'] = mzmine_parameter
            else:
                mzmine_running_parameters = ','.join(custom_mzmine_batch_params)
                mzmine_parameter = 5
                logging.info(tab_print("Using custom MZmine parameter file(s): %s"%(mzmine_running_parameters), 2))
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
                lims_untargeted_table_updater[mzmine_status_header] = '09 error' # Set default to error to indicate not passing required steps
                lims_untargeted_table_updater[fbmn_status_header] = '09 error'
                if polarity not in polarity_list: # Write blank columns to LIMS and skip writing directory/files to disk
                    lims_untargeted_table_updater[metadata_header] = ""
                    lims_untargeted_table_updater[file_count_header] = 0
                else: # Write directory/files to disk for present polarity and fill in LIMS columns
                    logging.info(tab_print("Writing MZmine submission input files...",2))
                    if os.path.exists(basepath):
                        logging.warning(tab_print("Warning! Directory %s already exists. Not writing new metadata or MZmine submission files..."%(basepath), 3))
                        continue
                    else:
                        os.mkdir(basepath)
                    try:
                        recursive_chown(basepath, 'metatlas')
                    except:
                        logging.info(tab_print("Note: Could not change group ownership of %s"%(basepath), 2))

                    logging.info(tab_print("%s metadata file (*_metadata.tab)"%(polarity), 3))
                    metadata_df = new_project_dict[polarity]['metadata_df']
                    metadata_filename = os.path.join(basepath,'%s_metadata.tab'%(parent_dir))
                    metadata_df.to_csv(metadata_filename, sep='\t', index=False)
                    
                    if custom_mzmine_batch_params is None:
                        logging.info(tab_print("%s MZmine parameter file (*_batch-params.xml)"%(polarity), 3))
                        params_filename = build_untargeted_filename(output_dir,project_name,polarity,'batch-params-mzmine')
                        with open(mzmine_running_parameters,'r') as fid:
                            orig_params = fid.read()
                        new_param_path = os.path.join(basepath,parent_dir)
                        custom_params = orig_params.replace('/Users/bpb/Downloads/mzmine_outputs',new_param_path)
                        with open(params_filename,'w') as fid:
                            fid.write('%s'%custom_params)
                    elif custom_mzmine_batch_params is not None and len(custom_mzmine_batch_params) == 1:
                        logging.info(tab_print("%s MZmine parameter file (*_batch-params.xml)"%(polarity), 3))
                        params_filename = build_untargeted_filename(output_dir,project_name,polarity,'batch-params-mzmine')
                        with open(mzmine_running_parameters,'r') as fid:
                            orig_params = fid.read()
                        new_param_path = os.path.join(basepath,parent_dir)
                        custom_params = orig_params.replace('/Users/bpb/Downloads/mzmine_outputs',new_param_path)
                        with open(params_filename,'w') as fid:
                            fid.write('%s'%custom_params)                        
                    else:
                        logging.info(tab_print("%s MZmine parameter file (see custom input above)"%(polarity), 3))
                        for custom_param in custom_mzmine_batch_params:
                            if polarity_short.upper()+"-" in custom_param:
                                mzmine_running_parameters = custom_param
                                break
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

                    if lims_untargeted_table_updater[file_count_header] is None or lims_untargeted_table_updater[file_count_header] == "NaN":
                        lims_untargeted_table_updater[file_count_header] = 0
                    if lims_untargeted_table_updater[mzmine_status_header] == "NaN":
                        lims_untargeted_table_updater[mzmine_status_header] = '12 not relevant'
                    if lims_untargeted_table_updater[fbmn_status_header] == "NaN":
                        lims_untargeted_table_updater[fbmn_status_header] = '12 not relevant'

            lims_untargeted_list.append(lims_untargeted_table_updater)

        lims_untargeted_df = pd.DataFrame(lims_untargeted_list)

        if lims_untargeted_df.shape[0] > 0:
            logging.info(tab_print("Updating LIMS table with %s new projects..."%(lims_untargeted_df.shape[0]), 1))
            update_table_in_lims(lims_untargeted_df,'untargeted_tasks',method='insert',max_size=1000)
            logging.info(tab_print("LIMS untargeted tasks table update complete.", 2))
        else:
            logging.info(tab_print("LIMS does not need updating!", 1))
    else:
        lims_untargeted_list = []
        lims_untargeted_df = pd.DataFrame(lims_untargeted_list)

    logging.info(tab_print("Exported table for %s new projects to pass to MZmine submission function."%(lims_untargeted_df.shape[0]), 1))
    return lims_untargeted_df