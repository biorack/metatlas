import sys
from pathlib2 import PurePath
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
from metatlas.tools.validate_filenames import parent_dir_num_fields

##### Kick off the script
mzm.start_script(script="check_untargeted_status.py")

##### Check for CLI input
if len(sys.argv) < 1:
    print('Please provide a project name or list of project names separated by commas.')
    sys.exit(1)

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

##### Run updaters before retrieving status from LIMS and printing to stdout
if project_list:
    print('\nStep 1/3: Checking and updating status of MZmine jobs in LIMS...')
    mzm.update_mzmine_status_in_untargeted_tasks(direct_input=project_list)
    print('\nStep 2/3: Checking and updating status of FBMN jobs in LIMS...')
    mzm.update_fbmn_status_in_untargeted_tasks(direct_input=project_list)
    print('\nStep 3/3: Printing job status for project queries...')
    mzm.get_untargeted_status(direct_input=project_list)
else:
    print('\nStep 1/3: Checking and updating status of MZmine jobs in LIMS...')
    mzm.update_mzmine_status_in_untargeted_tasks()
    print('\nStep 2/3: Checking and updating status of FBMN jobs in LIMS...')
    mzm.update_fbmn_status_in_untargeted_tasks()
    print('\nStep 3/3: Printing job status for project queries...')
    mzm.get_untargeted_status()

##### Wrap up the script
mzm.end_script(script="check_untargeted_status.py")