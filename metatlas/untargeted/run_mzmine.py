import sys
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
from metatlas.tools.validate_filenames import parent_dir_num_fields
from pathlib2 import PurePath

##### Kick off the script
mzm.start_script(script="run_mzmine.py")

##### Set up the download and output directories
download_folder = '/global/cfs/cdirs/metatlas/projects/untargeted_outputs'
output_dir = '/global/cfs/cdirs/metatlas/projects/untargeted_tasks'
raw_data_dir = '/global/cfs/cdirs/metatlas/raw_data'

##### Get project list from CLI standard input and validate
if len(sys.argv) > 1:
    projects = sys.argv[1]
    project_list = projects.split(',')
    for project in project_list:
        validate = parent_dir_num_fields(PurePath(project))
        if not validate:
            print(f'{project} is not a valid project name, exiting.')
            sys.exit(1)
else:
    project_list = []

##### Submit new mzmine jobs to perlmutter
if project_list:
    print("\nNotice! Not syncing LIMS and NERSC before executing script!\n")
    print("Notice! Only executing script with user-selected projects!\n")
    print('Step 1/2: Checking and updating status of MZmine jobs in LIMS...')
    mzm.update_mzmine_status_in_untargeted_tasks(direct_input=project_list)
    print('\nStep 2/2: Submitting new MZmine jobs that are in the command line input...')
    mzm.submit_mzmine_jobs(direct_input=project_list)
else:
    print('\nStep 1/3: Syncing LIMS and NERSC to identify new projects with raw data that are not yet in the untargeted task list...')
    new_projects = mzm.update_new_untargeted_tasks(update_lims=True,output_dir=output_dir,raw_data_dir=raw_data_dir)
    print('\nStep 2/3: Checking and updating status of MZmine jobs in LIMS...')
    mzm.update_mzmine_status_in_untargeted_tasks()
    print('\nStep 3/3: Submitting new MZmine jobs that are "initialized"...')
    mzm.submit_mzmine_jobs(new_projects=new_projects)

##### Wrap up the script
mzm.end_script(script="run_mzmine.py")