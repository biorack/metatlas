import sys
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
from metatlas.tools.validate_filenames import parent_dir_num_fields
from pathlib2 import PurePath

##### Kick off the script
mzm.start_script(script="run_fbmn.py")

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

##### Submit new fbmn jobs to perlmutter
if project_list:
    print("\nNotice: Only executing script with user-selected projects...\n")
    print('Checking and updating status of MZmine jobs in LIMS...')
    mzm.update_mzmine_status_in_untargeted_tasks(direct_input=project_list)
    print('\nChecking and updating status of FBMN jobs in LIMS...')
    mzm.update_fbmn_status_in_untargeted_tasks(direct_input=project_list)
    print('\nSubmitting new FBMN jobs to GNPS2...')
    mzm.submit_fbmn_jobs(direct_input=project_list)
else:
    print('\nChecking and updating status of MZmine jobs in LIMS...')
    mzm.update_mzmine_status_in_untargeted_tasks()
    print('\nChecking and updating status of FBMN jobs in LIMS...')
    mzm.update_fbmn_status_in_untargeted_tasks()
    print('\nSubmitting new FBMN jobs to GNPS2...')
    mzm.submit_fbmn_jobs()

##### Wrap up the script
mzm.end_script(script="run_fbmn.py")