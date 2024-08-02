import sys
import os
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
import subprocess

##### Kick off the script
mzm.kickoff()

##### Get ready for exporting results
output_dir = '/global/cfs/cdirs/metatlas/projects/untargeted_tasks'
download_folder = '/global/cfs/cdirs/metatlas/projects/untargeted_outputs'

# #### Remove contaminants from recent mgf files (202 ions for fluoranthene)
print("Removing fluoranthene from recent mgf files in mzmine results...")
recent_files = mzm.get_recent_mgf_files(time_back=30)
for f in recent_files:
    if f:
        mzm.remove_contaminant_from_mgf(f)

##### Check status of mzmine jobs
print('Checking and updating status of initialized or running mzmine jobs...')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='positive',polarity_short='pos')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='negative',polarity_short='neg')

##### Zip up the output folders
print("\nZipping up mzmine output folders for any recently completed projects...")
df = mzm.get_table_from_lims('untargeted_tasks')
df = mzm.filter_common_bad_project_names(df)
df = df[(df['mzmine_pos_status'] == '07 complete') & (df['mzmine_neg_status'] == '07 complete')]
override=False
new_results = []
for i,row in df.iterrows():
    output_zip_archive = os.path.join(download_folder,'%s.zip'%row['parent_dir'])
    input_folder = os.path.join(output_dir,'%s*'%row['parent_dir'])
    if override==True or not os.path.isfile(output_zip_archive):
        new_results.append(output_zip_archive)
        cmd = 'zip -rjuq %s %s'%(output_zip_archive,input_folder)
        os.system(cmd)
        mzm.recursive_chown(output_zip_archive, 'metatlas')
        print("\tNew pos and neg mzmine results written to %s\n"%output_zip_archive)

##### Copy files to Google Drive
if len(new_results) > 0:
    print("Copying zipped mzmine project directory to untargeted folder on Google Drive...")
    for project in new_results:
        if os.path.exists(project):
            cmd = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone copy --size-only {project} ben_lbl_gdrive:/untargeted_outputs'
            subprocess.check_output(cmd, shell=True)
            project_folder = os.path.basename(project)
            check_upload = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone ls ben_lbl_gdrive:/untargeted_outputs | grep "{project_folder}"'
            check_upload_out = subprocess.check_output(check_upload, shell=True)
            if check_upload_out.decode('utf-8').strip():
                print("\tCopied %s to Google Drive.\n"%project)
            else:
                print("\tWarning: Copying %s to Google Drive may have failed.\n"%project)
        else:
            print("\tWarning: Local ath %s does not exist.\n"%project)
else:
    print("No new mzmine results to copy to Google Drive.\n")
