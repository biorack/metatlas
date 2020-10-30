# import shutil
import os
import subprocess

import argparse
import labkey as lk
import pandas as pd
from labkey.api_wrapper import APIWrapper



key_file = '/global/cfs/cdirs/metatlas/labkey_user.txt'
with open(key_file,'r') as fid:
    api_key = fid.read().strip()
labkey_server='metatlas-dev.nersc.gov'
project_name='LIMS/'
api = APIWrapper(labkey_server, project_name, use_ssl=True,api_key=api_key)



def get_tasks_from_lims():
    sql = """SELECT input_file,output_file FROM file_conversion_task WHERE task='raw_to_mzml';"""
#     context_path = "labkey"
    
    sql_result = api.query.execute_sql('lists', sql,max_rows=1e6)
    if sql_result is None:
        print(('execute_sql: Failed to load results from ' + schema + '.' + table))
        return None
    else:
        df = pd.DataFrame(sql_result['rows'])
        df = df[[c for c in df.columns if not c.startswith('_')]]
        return df
    
def convert_raw(input_file):
    """
    shifter --volume=$DATADIR:/mywineprefix --image=biocontainers/pwiz:phenomenal-v3.0.18205_cv1.2.54 mywine msconvert --32 --mzML $DATAFILE -o /global/cfs/cdirs/metatlas/projects/fake_rawdata/fakeuser/fake_experiment_for_tests

    a good temp_basedir is /global/cscratch1/sd/bpb/temp_raw
    you probably need to mkdir on it to be sure
    
    helpful tips are here:
    https://docs.nersc.gov/development/shifter/issues/#invalid-volume-map
    
    user nobody needs permission all the way up
    setfacl -R -m u:nobody:x /global/cfs/cdirs/metatlas/raw_data
    
    """
    command = 'shifter --volume=DATADIR:/mywineprefix --image=biocontainers/pwiz:phenomenal-v3.0.18205_cv1.2.54 mywine msconvert --32 --mzML DATAFILE -o OUTDIR'
    if os.path.isfile(input_file):
        print('given file',input_file)
        input_file = os.path.realpath(input_file)
        if "/project/projectdirs" in input_file:
            input_file = input_file.replace("/project/projectdirs","/cfs/cdirs")
        wd = os.path.dirname(input_file)
        command = command.replace('DATADIR',wd)
        command = command.replace('OUTDIR',wd)
        command = command.replace('DATAFILE',input_file)
        print('\n%s\n'%command)
        subprocess.run(command, shell=True,check=True)
    
       
def main():
    parser = argparse.ArgumentParser(description='a command line tool for converting raw to mzml at nersc.')
    # query the lims for files to convert
    parser.add_argument('-input_file','--input_file', help='Full path to input raw file', default=False)
#     parser.add_argument('-output_file','--output_file', help='Full path to output mzML file', default=False)
#     parser.add_argument('-temp_dir','--temp_dir', help='Full path to directory where files are actually converted', default='/global/cscratch1/sd/bpb/temp_raw')

    args = vars(parser.parse_args())
    df = get_tasks_from_lims()

    for i,row in df.iterrows():
#         args['input_file']
        try:
            convert_raw(row['input_file'])
        except:
            pass
    
    # update tasks in the lims that are completed
    
    # catch any errors and update task as "error" in the lims
    
    
if __name__ == "__main__":
    main()
