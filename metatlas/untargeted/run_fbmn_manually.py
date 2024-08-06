
import os
import sys
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
from pathlib2 import PurePath
from metatlas.tools.validate_filenames import parent_dir_num_fields

##### Kick off the script
mzm.kickoff()

##### Parse the 'projects' parameter
if len(sys.argv) < 2:
    print("Usage: python run_fbmn_at_gnps2.py <projects>\n")
    print("<projects> should be a comma-separated list of valid project names")
    sys.exit(1)

##### Get project list from CLI standard input and validate
projects = sys.argv[1]
project_list = projects.split(',')
for project in project_list:
    validate = parent_dir_num_fields(PurePath(project))
    if not validate:
        print(f'{project} is not a valid project name, exiting.')
        sys.exit(1)

##### Check status of fbmn jobs
print('Checking status of fbmn jobs at GNPS2 and updating LIMS...')
mzm.update_fbmn_status_in_untargeted_tasks(polarity='positive',polarity_short='pos')
mzm.update_fbmn_status_in_untargeted_tasks(polarity='negative',polarity_short='neg')

##### Pick up secret file and extract password
with open('/global/homes/m/msdata/gnps2/gnps2_bpbowen.txt','r') as fid:
    t = fid.read()
t = t.split('\n')[0]
var, pword = t.split('=')
pword = pword.replace('"', '').replace("'", '').strip()
username = "bpbowen"
password = pword

##### Submit jobs
polarities = ['positive','negative']
task_list = []
print("Attempting to submit %s new FBMN jobs at GNPS2:"%(len(project_list)*2))
print("\t%s",project_list)
for exp in project_list:
    for pol in polarities:    
        department = exp.split('_')[1].lower()
        if department =='eb':
            department = 'egsb'
        if not department in ['jgi','egsb']:
            print(f'{exp} does not have a valid department name in the second field')
            continue
        description = f'{exp}_{pol}'
        spectra_file = f'USERUPLOAD/bpbowen/untargeted_tasks/{exp}_{pol}/{exp}_{pol}.mgf'
        quant_file = f'USERUPLOAD/bpbowen/untargeted_tasks/{exp}_{pol}/{exp}_{pol}_quant.csv'
        metadata_file = f'USERUPLOAD/bpbowen/untargeted_tasks/{exp}_{pol}/{exp}_{pol}_metadata.tab'
        raw_data = f'USERUPLOAD/bpbowen/raw_data/{department}/{exp}'
        if not os.path.exists(raw_data.replace('USERUPLOAD/bpbowen/','/global/cfs/cdirs/metatlas/')):
            print(f'{raw_data} does not exist')
            continue
        username = username
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
        job_id = mzm.submit_quickstart_fbmn(params, username, password)
        task_list.append({'experiment':exp,'polarity':pol,'response':job_id})

##### Write GNPS2 task IDs to file for reference
mzm.write_tasks_to_file(task_list)


