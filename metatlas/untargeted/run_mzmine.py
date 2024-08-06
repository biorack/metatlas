import sys
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
import subprocess
from metatlas.tools.validate_filenames import parent_dir_num_fields
from pathlib2 import PurePath

##### Kick off the script
mzm.kickoff()

##### Find new projects in raw data that aren't in LIMS
print('Finding projects in raw data that are not already in the LIMS untargeted tasks table')
new_projects = mzm.update_new_untargeted_tasks()

##### Check status of mzmine jobs and update LIMS untargeted tasks table
print('Checking and updating status of initialized or running mzmine jobs...')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='positive',polarity_short='pos')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='negative',polarity_short='neg')

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

##### Submit new mzmine jobs
print('Submitting new mzmine jobs that are intitialized...')
if project_list:
    print("\tNotice: Submitting fbmn jobs for only user-defined projects...")
    mzm.submit_mzmine_jobs(polarity='positive',polarity_short='pos',direct_input=project_list)
    mzm.submit_mzmine_jobs(polarity='negative',polarity_short='neg',direct_input=project_list)
else:
    mzm.submit_mzmine_jobs(polarity='positive',polarity_short='pos')
    mzm.submit_mzmine_jobs(polarity='negative',polarity_short='neg')

##### End
print("Mzmine submission script complete. Monitor jobs with sqs:")
subprocess.run(['sqs'], check=True)