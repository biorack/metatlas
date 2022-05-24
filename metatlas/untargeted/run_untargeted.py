import sys
import os
sys.path.insert(0,'/global/homes/b/bpb/repos/metatlas')
from metatlas.untargeted import tools as mzm
import subprocess


new_folders = mzm.update_new_untargeted_tasks()
print('\nChecking status of mzmine jobs')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='positive',polarity_short='pos')
mzm.update_mzmine_status_in_untargeted_tasks(polarity='negative',polarity_short='neg')
print('\nChecking status of fbmn jobs')
mzm.update_fbmn_status_in_untargeted_tasks(polarity='positive',polarity_short='pos')
mzm.update_fbmn_status_in_untargeted_tasks(polarity='negative',polarity_short='neg')
print('\nSubmitting new mzmine jobs')
mzm.submit_mzmine_jobs(polarity='positive',polarity_short='pos')
mzm.submit_mzmine_jobs(polarity='negative',polarity_short='neg')
print('\nSubmitting new fbmn jobs')
mzm.submit_fbmn_jobs(polarity='positive',polarity_short='pos',N=20)
mzm.submit_fbmn_jobs(polarity='negative',polarity_short='neg',N=20)


mzm.update_num_features()
mzm.update_num_msms()

# import datetime
# datetime.datetime.now()

# update gnps zip file downloads


output_dir = '/global/cfs/cdirs/metatlas/projects/untargeted_tasks'
download_folder = '/global/cfs/cdirs/metatlas/projects/untargeted_outputs'

df = mzm.get_table_from_lims('untargeted_tasks')
print(df.shape)
for i,row in df.iterrows():
    for polarity in ['positive','negative']:
        mzm.get_gnps_zipfile(row['parent_dir'],output_dir,polarity,row['fbmn_%s_status'%polarity[:3]],override=False) 



override=True
for i,row in df.iterrows():
    output_zip_archive = os.path.join(download_folder,'%s.zip'%row['parent_dir'])
    input_folder = os.path.join(output_dir,'%s*'%row['parent_dir'])
    if (override==True) | (not os.path.isfile(output_zip_archive)):
        cmd = 'zip -rju %s %s'%(output_zip_archive,input_folder)
        os.system(cmd)

os.system("chgrp -R metatlas %s"%download_folder)
os.system("chmod o+rx %s"%download_folder)
os.system("chmod o+rx %s"%download_folder)

# cmd = """/global/cfs/cdirs/m342/USA/shared-repos/rclone/bin/rclone copy --update /global/cfs/cdirs/metatlas/projects/untargeted_outputs ben_lbl_gdrive:/untargeted_outputs"""
cmd = """/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone copy --size-only /global/cfs/cdirs/metatlas/projects/untargeted_outputs ben_lbl_gdrive:/untargeted_outputs"""
subprocess.check_output(cmd, shell=True)
