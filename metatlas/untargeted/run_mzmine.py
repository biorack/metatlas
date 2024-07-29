import sys
import os
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
import subprocess

##### Kick off the script
mzm.kickoff()

##### Find new projects in raw data that aren't in LIMS
print('\nFinding projects in raw data that are not already in the LIMS untargeted tasks table')
new_folders = mzm.update_new_untargeted_tasks()

if len(new_folders['parent_dir']) == 0:
    print('\nNo new projects found in raw data. Exiting script.')
    sys.exit(0)

##### Remove any DNA SIP jobs from untargeted tasks in LIMS
print('\nRemoving any DNA SIP jobs from untargeted tasks in LIMS...')
df = mzm.get_table_from_lims('untargeted_tasks',columns=['Key','parent_dir','mzmine_pos_status','mzmine_neg_status','fbmn_pos_status','fbmn_neg_status'])
df = df[(df['parent_dir'].str.contains('DNA-SIP')) & (df['parent_dir'].str.contains('_RM'))]
cols = [c for c in df.columns if c.endswith('status')]
df[cols] = '12 not relevant'
mzm.update_table_in_lims(df,'untargeted_tasks',method='update')

##### Check status of mzmine jobs
print('\nChecking and updating status of initialized or running mzmine jobs...')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='positive',polarity_short='pos')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='negative',polarity_short='neg')

##### Submit new mzmine jobs
print('\nSubmitting new mzmine jobs that are pending intitialization...\n')
mzm.submit_mzmine_jobs(polarity='positive',polarity_short='pos')
mzm.submit_mzmine_jobs(polarity='negative',polarity_short='neg')

##### End
print("\nMzmine submission script complete. Monitor jobs with sqs:\n")
subprocess.run(['sqs'], check=True)