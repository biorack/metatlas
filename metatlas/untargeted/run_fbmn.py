import sys
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
from metatlas.tools.validate_filenames import parent_dir_num_fields
from pathlib2 import PurePath

##### Kick off the script
mzm.kickoff()

##### Check status of mzmine jobs and update LIMS untargeted tasks table (critical because fbmn
##### job submission depends on mzmine jobs being at completed status)
print('Updating status of completed mzmine jobs in LIMS...')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='positive',polarity_short='pos')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='negative',polarity_short='neg')

##### Check and updating status of fbmn jobs so they're not run multiple times
print('Checking status of fbmn jobs at GNPS2 and updating LIMS...')
mzm.update_fbmn_status_in_untargeted_tasks(polarity='positive',polarity_short='pos')
mzm.update_fbmn_status_in_untargeted_tasks(polarity='negative',polarity_short='neg')

##### Get project list from CLI standard input and validate
project_list = []
if len(sys.argv) > 1:
    projects = sys.argv[1]
    project_list = projects.split(',')
    for project in project_list:
        validate = parent_dir_num_fields(PurePath(project))
        if not validate:
            print(f'{project} is not a valid project name, exiting.')
            sys.exit(1)

##### Submit jobs to GNPS2
if project_list:
    print("\tNotice: Submitting fbmn jobs for only user-defined projects...")
    pos_task_list = mzm.submit_fbmn_jobs(polarity='positive',direct_input=project_list)
    neg_task_list = mzm.submit_fbmn_jobs(polarity='negative',direct_input=project_list)
    task_list = pos_task_list + neg_task_list
else:
    pos_task_list = mzm.submit_fbmn_jobs(polarity='positive')
    neg_task_list = mzm.submit_fbmn_jobs(polarity='negative')
    task_list = pos_task_list + neg_task_list

##### Write GNPS2 task IDs to file for reference
mzm.write_fbmn_tasks_to_file(task_list)

