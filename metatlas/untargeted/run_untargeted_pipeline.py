import os
import sys
import argparse
import logging
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm

def main():
    ##### Parse command line arguments
    parser = argparse.ArgumentParser(description='Run untargeted pipeline.')
    parser.add_argument('--download_dir', type=str, default='/global/cfs/cdirs/metatlas/projects/untargeted_outputs', help='Path to the download folder')
    parser.add_argument('--output_dir', type=str, default='/global/cfs/cdirs/metatlas/projects/untargeted_tasks', help='Path to the output directory')
    parser.add_argument('--raw_data_dir', type=str, default='/global/cfs/cdirs/metatlas/raw_data', help='Path to the raw data directory')
    parser.add_argument('--update_lims', type=bool, default=True, help='Update LIMS with new untargeted tasks')
    parser.add_argument('--validate_names', type=bool, default=True, help='Validate filenames and project names')
    parser.add_argument('--overwrite_zip',type=bool, default=False, help='Overwrite existing zip files in download folder')
    parser.add_argument('--overwrite_drive', type=bool, default=False, help='Overwrite existing zip files on google drive')
    parser.add_argument('--overwrite_fbmn', type=bool, default=False, help='Overwrite existing fbmn results files that are already in the output directory')
    parser.add_argument('--gdrive_upload', type=bool, default=True, help='Upload to google drive')
    parser.add_argument('--min_features', type=int, default=0, help='Set minimum number of MZmine features for a project polarity to be zipped')
    parser.add_argument('--add_gnps2_documentation', type=bool, default=True, help='File name of the GNPS2 documentation to add to project zips')
    parser.add_argument('--gnps2_doc_name', type=str, default='Untargeted_metabolomics_GNPS2_Guide.docx', help='File name of the GNPS2 documentation to add to project zips')
    parser.add_argument('--log_file', type=str, default='/global/homes/b/bkieft/metatlas/metatlas/untargeted/untargeted_pipeline_log.txt', help='Log file name with full path')
    parser.add_argument('--log_level', type=str, default='INFO', help='Logger level. One of [DEBUG, INFO, WARNING, ERROR, or CRITICAL]')
    parser.add_argument('--log_format', type=str, default='%(asctime)s - %(levelname)s - %(message)s', help='Logger format')
    parser.add_argument('--direct_input', type=str, default=None, help='Input project names from command line as a comma-separated list')
    parser.add_argument('--background_designator', type=str, default="TxCtrl,ExCtrl", help='Input control/background sample names from command line as a comma-separated list')

    ##### Configure logging
    mzm.call_logger(log_filename=args.log_file, log_level=args.log_level, log_format=args.log_format)

    ##### Set arguments to pass to the pipeline steps and run some checks
    args = parser.parse_args()
    if args.direct_input:
        args.direct_input = args.direct_input.split(',')
    if args.background_designator:
        args.background_designator = args.background_designator.split(',')
    if args.overwrite_drive is True and args.gdrive_upload is False:
        logging.error('Incompatible flags. Cannot overwrite google drive if not uploading to google drive.')
        sys.exit(1)
    if args.gnps2_doc_name is None and args.add_gnps2_documentation is True:
        logging.error('Incompatible flags. Must provide GNPS2 document name if you want to add it to the zip.')
        sys.exit(1)
    if not os.path.exists(args.download_dir):
        logging.error('Download directory does not exist. Please check flag.')
        sys.exit(1)
    if not os.path.exists(args.output_dir):
        logging.error('Output directory does not exist. Please check flag.')
        sys.exit(1)
    if not os.path.exists(args.raw_data_dir):
        logging.error('Raw data directory does not exist. Please check flag.')
        sys.exit(1)

    ##### Kick off the script
    mzm.start_script(script="run_untargeted_pipeline.py")

    ##### Run the pipeline steps
    logging.info('Step 1/7: Syncing LIMS and NERSC to identify new projects with raw data that are not yet in the untargeted task list...')
    new_projects = mzm.update_new_untargeted_tasks(update_lims=args.update_lims, validate_names=args.validate_names, \
                                                   output_dir=args.output_dir, raw_data_dir=args.raw_data_dir,
                                                   background_designator=args.background_designator)

    logging.info('Step 2/7: Checking and updating status of MZmine jobs in LIMS...')
    mzm.update_mzmine_status_in_untargeted_tasks(direct_input=args.direct_input,background_designator=args.background_designator)

    logging.info('Step 3/7: Submitting new MZmine jobs that are "initialized"...')
    mzm.submit_mzmine_jobs(new_projects=new_projects, direct_input=args.direct_input)

    logging.info('Step 4/7: Checking and updating status of FBMN jobs in LIMS...')
    mzm.update_fbmn_status_in_untargeted_tasks(direct_input=args.direct_input)

    logging.info('Step 5/7: Submitting new FBMN jobs to GNPS2...')
    mzm.submit_fbmn_jobs(output_dir=args.output_dir, overwrite_fbmn=args.overwrite_fbmn, validate_names=args.validate_names,direct_input=args.direct_input)

    logging.info('Step 6/7: Checking for completed FBMN jobs and downloading results...')
    mzm.download_fbmn_results(output_dir=args.output_dir, overwrite_fbmn=args.overwrite_fbmn,direct_input=args.direct_input)

    logging.info("Step 7/7: Zipping up and (optionally) uploading output folders to gdrive...")
    mzm.zip_and_upload_untargeted_results(download_folder=args.download_dir,output_dir=args.output_dir,upload=args.gdrive_upload,overwrite_zip=args.overwrite_zip, \
                                          overwrite_drive=args.overwrite_drive, min_features_admissible=args.min_features, \
                                          add_documentation=args.add_gnps2_documentation,doc_name=args.gnps2_doc_name,direct_input=args.direct_input)

    ##### End the script
    mzm.end_script(script="run_untargeted_pipeline.py")

if __name__ == "__main__":
    main()