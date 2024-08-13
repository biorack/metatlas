import sys
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm

##### Kick off the script
mzm.start_script(script="run_untargeted_pipeline.py")

##### Set up specific paths and filenames
download_folder = '/global/cfs/cdirs/metatlas/projects/untargeted_outputs'
output_dir = '/global/cfs/cdirs/metatlas/projects/untargeted_tasks'
raw_data_dir = '/global/cfs/cdirs/metatlas/raw_data'
gnps2_doc_name="Untargeted_metabolomics_GNPS2_Guide.docx"

##### Run the pipeline steps
print('\nStep 1/7: Syncing LIMS and NERSC to identify new projects with raw data that are not yet in the untargeted task list...')
new_projects = mzm.update_new_untargeted_tasks(update_lims=True, validate_names=True, output_dir=output_dir, raw_data_dir=raw_data_dir)

print('\nStep 2/7: Checking and updating status of MZmine jobs in LIMS...')
mzm.update_mzmine_status_in_untargeted_tasks()

print('\nStep 3/7: Submitting new MZmine jobs that are "initialized"...')
mzm.submit_mzmine_jobs(new_projects=new_projects, submit_only_new=True)

print('\nStep 4/7: Checking and updating status of FBMN jobs in LIMS...')
mzm.update_fbmn_status_in_untargeted_tasks()

print('\nStep 5/7: Submitting new FBMN jobs to GNPS2...')
mzm.submit_fbmn_jobs(output_dir=output_dir, overwrite=False, validate_names=True)

print('\nStep 6/7: Checking for completed FBMN jobs and downloading results...')
mzm.download_fbmn_results(output_dir=output_dir, overwrite=False)

print("\nStep 7/7: Zipping up and (optionally) uploading output folders to gdrive...")
mzm.zip_and_upload_untargeted_results(download_folder=download_folder, output_dir=output_dir, upload=True, overwrite_zip=True, \
                                        overwrite_drive=True, min_features_admissible=0, add_documentation=True, doc_name=gnps2_doc_name)

##### End the script
mzm.end_script(script="run_untargeted_pipeline.py")