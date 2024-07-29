import sys
import os
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm
import subprocess
import re

##### Kick off the script
mzm.kickoff()

##### Get ready for exporting results
output_dir = '/global/cfs/cdirs/metatlas/projects/untargeted_tasks'
download_folder = '/global/cfs/cdirs/metatlas/projects/untargeted_outputs'

# #### Remove contaminants from recent mgf files (202 ions for fluoranthene)
print("\nRemoving fluoranthene from recent mgf files in mzmine results...")
recent_files = mzm.get_recent_mgf_files()
for f in recent_files:
    if f:
        mzm.remove_contaminant_from_mgf(f)

##### Zip up the output folders
print("\nZipping up the output folders...")
df = mzm.get_table_from_lims('untargeted_tasks')
df = df[~(df['parent_dir'].str.contains('DNA-SIP')) & ~(df['parent_dir'].str.contains('_RM'))]
df = df[~(df['parent_dir'].str.contains(' '))]
df['parent_dir'] = df['parent_dir'].str.replace(r'&', r'\&', regex=True) # escape ampersands in some project names
override=False
new_results = []
for i,row in df.iterrows():
    output_zip_archive = os.path.join(download_folder,'%s.zip'%row['parent_dir'])
    input_folder = os.path.join(output_dir,'%s*'%row['parent_dir'])
    if override==True or not os.path.isfile(output_zip_archive.replace(r'\&', '&')):
        new_results.append(output_zip_archive.replace(r'\&', '&'))
        cmd = 'zip -rjuq %s %s'%(output_zip_archive,input_folder)
        os.system(cmd)
        print("\tNew pos and neg mzmine results written to %s\n"%output_zip_archive)

##### Copy files to Google Drive
if len(new_results) > 0:
    print("Copying files to untargeted folder on Google Drive...")
    for project in new_results:
        if os.path.exists(project):
            cmd = f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone copy --size-only {project} bpk_lbl_gdrive:/untargeted_outputs'
            subprocess.check_output(cmd, shell=True)
            print("\tCopied %s to Google Drive.\n"%project)
        else:
            print("\tWarning: Path %s does not exist.\n"%project)
else:
    print("No new mzmine results to copy to Google Drive.\n")
