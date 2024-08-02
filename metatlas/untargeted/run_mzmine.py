import sys
import os
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
import subprocess

##### Kick off the script
mzm.kickoff()

##### Find new projects in raw data that aren't in LIMS
print('Finding projects in raw data that are not already in the LIMS untargeted tasks table')
new_projects = mzm.update_new_untargeted_tasks()

##### Add some guardrails for number of new projects that are identified
if len(new_projects['parent_dir']) == 0:
    print('No new projects found in raw data that are not already in untargeted tasks. Exiting script.')
    sys.exit(0)
if len(new_projects['parent_dir']) > 20:
    print('There are more than 20 new projects in raw data that are not in untargeted tasks, please check if this is accurate. Exiting script.')
    sys.exit(0)

##### Check status of mzmine jobs and update LIMS untargeted tasks table
print('Checking and updating status of initialized or running mzmine jobs...')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='positive',polarity_short='pos')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='negative',polarity_short='neg')

##### Submit new mzmine jobs
print('Submitting new mzmine jobs that are intitialized...')
mzm.submit_mzmine_jobs(polarity='positive',polarity_short='pos')
mzm.submit_mzmine_jobs(polarity='negative',polarity_short='neg')

##### End
print("Mzmine submission script complete. Monitor jobs with sqs:")
subprocess.run(['sqs'], check=True)