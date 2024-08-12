import sys
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
from metatlas.tools.validate_filenames import parent_dir_num_fields
from pathlib2 import PurePath

##### Kick off the script
mzm.start_script(script="download_fbmn_results.py")

##### Set up the download and output directories
download_folder = '/global/cfs/cdirs/metatlas/projects/untargeted_outputs'
output_dir = '/global/cfs/cdirs/metatlas/projects/untargeted_tasks'

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

##### Export untargeted results if mzmine and fbmn are complete
if project_list:
    print("\nNotice! Only executing script with user-selected projects...\n")
    print('Step 1/3: Checking and updating status of MZmine jobs in LIMS...')
    mzm.update_mzmine_status_in_untargeted_tasks(direct_input=project_list)
    print('\nStep 2/3: Checking and updating status of FBMN jobs in LIMS...')
    mzm.update_fbmn_status_in_untargeted_tasks(direct_input=project_list)
    print('\nStep 3/3: Checking for completed FBMN jobs and downloading results...')
    mzm.download_fbmn_results(output_dir=output_dir,overwrite=True,direct_input=project_list)
else:
    print('\nStep 1/3: Checking and updating status of MZmine jobs in LIMS...')
    mzm.update_mzmine_status_in_untargeted_tasks()
    print('\nStep 2/3: Checking and updating status of FBMN jobs in LIMS...')
    mzm.update_fbmn_status_in_untargeted_tasks()
    print('\nStep 3/3: Checking for completed FBMN jobs and downloading results...')
    mzm.download_fbmn_results(output_dir=output_dir,overwrite=False)

##### Wrap up the script
mzm.end_script(script="download_fbmn_results.py")