#!/usr/bin/env python
# This script takes a task name (same name as in LAB KEY LIMS) and returns the 
# input file paths and proposed output paths. 
import sys
import os
import argparse
from shutil import copyfile
import sqlite_functions

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
parser.add_argument("-c","--config",help="This is a json file that will be used to create the inputs.json file for the WDL. It's format should be like a regular \"inputs.json\".", type=str, required=True)
parser.add_argument("-o","--out",help="The name of the ouput json file that the WDL will use.", type=str, required=True)
parser.add_argument("-a","--api",help="The file that contains the api key, i.e. \"apikey|23ladsf9932riadifa\".", type=str, required=True)
parser.add_argument("-w","--wdl",help="The name as shown in the first line of the WDL file. This will be used in the inputs.json file.", type=str, required=True)
args = parser.parse_args()

labkey_server='metatlas-dev.nersc.gov'
project_name='/LIMS'

must_be = [
      "mzml_to_hdf5",
      "mzml_to_pactolus",
      "mzml_to_spectralhits",
      "raw_to_mzml"]

if not os.path.exists(args.config):
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

def get_raw_files(task_name,api_key):
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

    # returns a list of lists
    ad=df.head(2).to_dict()
    tmp_dict={}

    for k in ad:
        nd=[]
        for i in ad[k]:
            nd.append(ad[k][i])
        tmp_dict[k] = nd

    # a dictionary of a dictionary
    master_dict[task_name] = tmp_dict


def create_inputs_json(master_dict):
    """ This function will create a input.json file for running the WDL in JAWS or cromwell. 
    This file is a combination of input files pulled from lims, and a config file, defining some params that the WDL expects.
    The three fields that were queried from LIMS (in the function: get_raw_files) where
        1. Key, 
        2. input_file
        3. output_file
    """
    # convert the input config json into a dictionary so we can combine it with "master_dict"
    with open(args.config) as json_file:
        config_dict = json.load(json_file)
    
    # add the name of the workflow (i.e. first line of wdl) to each dictionary key so inputs.json will have complete key names.
    inputs_only_dict = {}
    for task in master_dict:
        key = args.wdl + "." + task # add wdl name to key
        list_size = len(master_dict[task]['Key'])
        inputs_only_dict[key] = []
        for i in range(list_size):
            #inputs_only_dict[key].append([ master_dict[task]['Key'][i], master_dict[task]['input_file'][i] ])
            inputs_only_dict[key].append( {'Left': master_dict[task]['Key'][i], 'Right': master_dict[task]['input_file'][i] })
    
    
    # merge the two dictionaries
    inputs_only_dict.update(config_dict)
    
    # create the inputs.json file
    with open(args.out,'w') as fh:
        json_obj = json.dump(inputs_only_dict,fh,sort_keys=True,indent=4, separators=(',', ': '))

    return config_dict


def create_sql_db(db,master_dict):
    records=[]
    conn = sqlite_functions.sql_connect(db)
    for task in master_dict:
        array_size = len( master_dict[task]['Key'])
        for i in range(array_size):
    
            key = master_dict[task]['Key'][i]
            input_file = master_dict[task]['input_file'][i]
            output_file = master_dict[task]['output_file'][i]
            tmp_list = (key,task,input_file,output_file)
            records.append(tmp_list)
    
    sql_insert='INSERT INTO metadata(limskey,task,input,output) VALUES(?,?,?,?)'
    conn.insert_bulk(sql_insert, records)

############
### MAIN ###
############

# grab token for LIMS
with open(args.api,'r') as fid:
    api_key = fid.read()
api_key = api_key.strip()

needed_tasks = tasks_to_do(api_key)

# for each task, grab files that need to be processed from LIMS KEY
master_dict = {}
for task in must_be:
    # only try to get files if there were any files to get.
    if task in needed_tasks:
        get_raw_files(task,api_key)
    else:
        print("There were no files for task: %s" % task)

###                             ###
### Create the inputs.json file ###
###                             ###
config_dict = create_inputs_json(master_dict)

###               ###
### create sql db ###
###               ###
k=args.wdl + '.db' # the name of the key in the config file
create_sql_db(config_dict[k], master_dict)
