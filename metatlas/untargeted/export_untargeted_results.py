import sys
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
from metatlas.tools.validate_filenames import parent_dir_num_fields
from pathlib2 import PurePath

# ##### Kick off the script
mzm.kickoff()

# ##### Get ready for exporting results
output_dir = '/global/cfs/cdirs/metatlas/projects/untargeted_tasks'
download_folder = '/global/cfs/cdirs/metatlas/projects/untargeted_outputs'

##### Check status of mzmine jobs
print('Checking and updating status of initialized or running mzmine jobs in LIMS...')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='positive',polarity_short='pos')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='negative',polarity_short='neg')

##### Check status of fbmn jobs
print('Checking and updating status of initialized or running fbmn jobs in LIMS...')
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

##### Download completed fbmn results to untargeted_tasks dir
print('Checking for completed fbmn jobs and downloading results...')
if project_list:
    print("\tNotice: Downloading and overwriting fbmn results only for user-defined projects...")
    mzm.download_fbmn_results(polarity='positive',polarity_short='pos',output_dir=output_dir,overwrite=True,direct_input=project_list)
    mzm.download_fbmn_results(polarity='negative',polarity_short='neg',output_dir=output_dir,overwrite=True,direct_input=project_list)
else:
    mzm.download_fbmn_results(polarity='positive',polarity_short='pos',output_dir=output_dir,overwrite=False)
    mzm.download_fbmn_results(polarity='negative',polarity_short='neg',output_dir=output_dir,overwrite=False)

##### Zip up the output folders and upload to google drive
print("Zipping up and (optionally) uploading output folders to gdrive for recently completed projects...")
if project_list:
    print("\tNotice: Zipping and exporting results only for user-defined projects...")
    mzm.zip_and_upload_untargeted_results(download_folder=download_folder,output_dir=output_dir, \
                                  upload=True,overwrite=True,direct_input=project_list)
else:
    mzm.zip_and_upload_untargeted_results(download_folder=download_folder,output_dir=output_dir, \
                                  upload=True,overwrite=False)



