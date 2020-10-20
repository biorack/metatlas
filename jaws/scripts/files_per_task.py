#!/usr/bin/env python
# This script takes a task name (same name as in LAB KEY LIMS) and returns the 
# input file paths and proposed output paths. 
import sys
import os
import argparse
import re
from shutil import copyfile

# you need this!
# https://github.com/LabKey/labkey-api-python
#sys.path.insert(1,'/global/dna/shared/data/jfroula/labkey-api-python/')
#sys.path.insert(1,'/global/cscratch1/sd/jfroula/JAWS/jgi-wdl-catalog/jaws-benbowen/labkey-api-python')
import labkey as lk
import requests
import json
import pandas as pd
import hashlib
import time
import numpy as np
from datetime import datetime, time as dtime
from subprocess import check_output
import re
import argparse

# parse arguments
task_name=''
parser = argparse.ArgumentParser(description='This script prints the file paths (and output paths) for a given task, where a task is what needs to be converted (i.e. mzml -> h5).')
parser.add_argument("-t","--task",help="task names should be one of these: [mzml_to_hdf5|mzml_to_pactolus|mzml_to_spectralhits|raw_to_mzml]", type=str, required=True)
parser.add_argument("-o","--out",help="ouput file as csv. Each line is raw input file and output file with proper suffix.", type=str, required=True)
parser.add_argument("-a","--api",help="", type=str, required=True)
args = parser.parse_args()

labkey_server='metatlas-dev.nersc.gov'
project_name='/LIMS'

must_be = [
      "mzml_to_hdf5",
      "mzml_to_pactolus",
      "mzml_to_spectralhits",
      "raw_to_mzml"]

if args.task not in must_be:
    print("Please use one of the accepted task names: %s" % must_be)
    sys.exit()

def tasks_to_do(api_key):
    """
    Possible tasks are:
      mzml_to_hdf5
      mzml_to_pactolus
      mzml_to_spectralhits
      raw_to_mzml
    """
    sql = """SELECT DISTINCT task FROM file_conversion_task;"""
    con = lk.utils.create_server_context(labkey_server, project_name, use_ssl=True,api_key=api_key)
    # base execute_sql
    schema = 'lists'
    sql_result = lk.query.execute_sql(con, schema, sql,max_rows=1e6)
    if sql_result is None:
        print('execute_sql: Failed to load results from ' + schema + '.' + table)
    else:
        df = pd.DataFrame(sql_result['rows'])
        df = df[[c for c in df.columns if not c.startswith('_')]]
    return list(df['task'])

def checkInputErrors(df,task_name):
    # keep a list of all indexes of files to omit from dataframe
    myindex_to_remove=[]

    # file with omitted inputs that user needs to fix. File is named after the task.
    omitted_inputs="omitted_" + task_name
    fh = open (omitted_inputs,"w")

    # check all input files for whitespaces
    for index, row in df.iterrows():
        res = bool(re.search(r"[\t\s]", row['input_file']))
        if res == True:
            myindex_to_remove.append(index)
            fh.write("Error: file has whitespace: %s\n" % row['input_file'])
    fh.close()

    # return a new df that has bad input files removed.
    return df.drop(index=myindex_to_remove)

def get_raw_files(task_name,files_found,api_key):
    sql = "SELECT Key, file_conversion_task.input_file,file_conversion_task.output_file FROM file_conversion_task WHERE file_conversion_task.task='%s' AND file_conversion_task.status <> '09 error';" % (task_name)
    con = lk.utils.create_server_context(labkey_server, project_name, use_ssl=True,api_key=api_key)
    # base execute_sql
    schema = 'lists'
    sql_result = lk.query.execute_sql(con, schema, sql,max_rows=1e6)
    if sql_result is None:
        print('execute_sql: Failed to load results from ' + schema + '.' + table)
    else:
        df = pd.DataFrame(sql_result['rows'])
        df = df[[c for c in df.columns if not c.startswith('_')]]
    
    pd.set_option('display.max_colwidth', -1)

    # check if input files have whitespace or not.  Omit the input file from future 
    # processing but write them to a log, if they do have whitespace. The log is something like "Omitted_<task_name>".
    new_df = checkInputErrors(df,task_name)

    # returns a list of lists
    #new_df.head(100).to_csv(files_found, sep='\t', index=False, header=False)
    new_df.to_csv(files_found, sep='\t', index=False, header=False)

# grab token for LIMS
with open(args.api,'r') as fid:
    api_key = fid.read()
api_key = api_key.strip()

needed_tasks = tasks_to_do(api_key)

if args.task in needed_tasks:
    get_raw_files(args.task,args.out,api_key)
else:
    print("There were no files for task: %s" % args.task)
