import sys
from pathlib2 import PurePath
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
from metatlas.tools.validate_filenames import parent_dir_num_fields

# Parse the 'projects' parameter
if len(sys.argv) < 2:
    print("Usage: python check_untargeted_status.py <projects>\n")
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

##### Update status of mzmine jobs in LIMS untargeted tasks table
print('\nUpdating status of completed mzmine jobs in LIMS...')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='positive',polarity_short='pos',direct_input=project_list)
mzm.update_mzmine_status_in_untargeted_tasks(polarity='negative',polarity_short='neg',direct_input=project_list)

##### Update status of fbmn jobs in LIMS untargeted tasks table
print('\nChecking status of fbmn jobs at GNPS2 and updating LIMS...')
mzm.update_fbmn_status_in_untargeted_tasks(polarity='positive',polarity_short='pos',direct_input=project_list)
mzm.update_fbmn_status_in_untargeted_tasks(polarity='negative',polarity_short='neg',direct_input=project_list)

# Export status of mzmine jobs
print('\nPrinting job status for project queries...\n')
status_df = mzm.get_untargeted_statuses(project_list)

print(status_df.to_string(),"\n")