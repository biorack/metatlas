import os
import sys
import argparse
import logging
sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/')
from metatlas.untargeted import tools as mzm

##### Parse command line arguments
def main():

    ##### Set arguments to pass to the pipeline steps and run some checks
    parser = argparse.ArgumentParser(description='Run untargeted pipeline.')
    add_arguments(parser)
    args = parser.parse_args()
    check_args(args)
    step_bools = check_skipped_steps(args)

    ##### Configure logging
    mzm.call_logger(log_filename=args.log_file, log_level=args.log_level, log_format=args.log_format)

    ##### Kick off the script
    start_message = mzm.start_script(script="run_untargeted_pipeline.py")
    logging.info('\n')
    logging.info(start_message)

    ##### Write args to log for reference
    logging.info(f'Arguments used: {args}')

    ##### Step 1/7: Syncing LIMS and NERSC to identify new projects with raw data that are not yet in the untargeted task list
    new_projects = mzm.update_new_untargeted_tasks(update_lims=args.update_lims, validate_names=args.validate_names, direct_input=args.direct_input,\
                                                   output_dir=args.output_dir, raw_data_dir=args.raw_data_dir, raw_data_subdir=args.raw_data_subdir, \
                                                   background_designator=args.background_designator,skip_sync=step_bools[0])
    
    ##### Step 2/7: Checking and updating status of MZmine jobs in LIMS
    mzm.update_mzmine_status_in_untargeted_tasks(direct_input=args.direct_input,background_designator=args.background_designator, \
                                                 skip_mzmine_status=step_bools[1],background_ratio=args.background_ratio, \
                                                 zero_value=args.zero_value,nonpolar_solvent_front=args.nonpolar_solvent_front, \
                                                 polar_solvent_front=args.polar_solvent_front)

    ##### Step 3/7: Submitting new MZmine jobs that are "initialized"
    mzm.submit_mzmine_jobs(new_projects=new_projects,direct_input=args.direct_input,skip_mzmine_submit=step_bools[2], \
                           overwrite_mzmine=args.overwrite_mzmine)

    ##### Step 4/7: Checking and updating status of FBMN jobs in LIMS
    mzm.update_fbmn_status_in_untargeted_tasks(direct_input=args.direct_input,skip_fbmn_status=step_bools[3])

    ##### Step 5/7: Submitting new FBMN jobs to GNPS2
    mzm.submit_fbmn_jobs(output_dir=args.output_dir, overwrite_fbmn=args.overwrite_fbmn, direct_input=args.direct_input, skip_fbmn_submit=step_bools[4], \
                        mirror_raw_data=args.mirror_raw_data, mirror_mzmine_results=args.mirror_mzmine_results, raw_data_dir=args.raw_data_dir,raw_data_subdir=args.raw_data_subdir)

    ##### Step 6/7: Checking for completed FBMN jobs and downloading results
    mzm.download_fbmn_results(output_dir=args.output_dir, overwrite_fbmn=args.overwrite_fbmn,direct_input=args.direct_input, \
                              skip_fbmn_download=step_bools[5])

    ##### Step 7/7: Zipping up and (optionally) uploading output folders to gdrive
    mzm.zip_and_upload_untargeted_results(download_folder=args.download_dir,output_dir=args.output_dir, \
                                          raw_data_subdir=args.raw_data_subdir, upload=args.gdrive_upload,overwrite_zip=args.overwrite_zip, \
                                          overwrite_drive=args.overwrite_drive, min_features_admissible=args.min_features, skip_zip_upload=step_bools[6], \
                                          add_documentation=args.add_gnps2_documentation,doc_name=args.gnps2_doc_name,direct_input=args.direct_input, \
                                          abridged_filenames=args.abridged_filenames)

    ##### End the script
    end_message = mzm.end_script(script="run_untargeted_pipeline.py")
    logging.info(end_message)

def add_arguments(parser):
    ## Multiple functions
    parser.add_argument('--output_dir', type=str, default='/global/cfs/cdirs/metatlas/projects/untargeted_tasks', help='Path to the output directory')
    parser.add_argument('--overwrite_fbmn', type=bool, default=False, help='Overwrite existing fbmn results files that are already in the output directory')
    parser.add_argument('--direct_input', type=str, default=None, help='Input project names from command line as a CSV list')
    parser.add_argument('--background_designator', type=str, default="ExCtrl", help='Input control/background sample names from command line as a CSV list')
    parser.add_argument('--skip_steps', type=str, default=None, help='Skip pipeline step(s) with --direct_input. CSV list with elements [Sync,MZmine_Status,MZmine_Submit,FBMN_Status,FBMN_Submit,FBMN_Download,Zip_and_Upload]')
    parser.add_argument('--raw_data_dir', type=str, default='/global/cfs/cdirs/metatlas/raw_data', help='Path to the raw data directory')
    parser.add_argument('--raw_data_subdir', type=str, default=None, help='Name of the raw_data subdirectory (e.g., jgi, egsb) to use for tasks that require raw data files. If not given, will try to infer from the project name.')
    ## Step 1 only
    parser.add_argument('--update_lims', type=bool, default=True, help='Update LIMS with new untargeted tasks')
    parser.add_argument('--validate_names', type=bool, default=False, help='Validate filenames and project names')
    ## Step 2 only
    parser.add_argument('--background_ratio', type=float, default=5, help='Ratio of background to sample intensity for filtering features')
    parser.add_argument('--zero_value', type=float, default=(2/3), help='Proportion of the lowest intensity value from the experiment to use as replacement zero value')
    parser.add_argument('--polar_solvent_front', type=float, default=0.8, help='Retention time to use as HILIC solvent front (mins) for filtering features')
    parser.add_argument('--nonpolar_solvent_front', type=float, default=0.5, help='Retention time to use as C18/LIPID solvent front (mins) for filtering features')
    ## Step 3 only
    parser.add_argument('--overwrite_mzmine', type=bool, default=False, help='Overwrite existing mzmine results files that are already in the output directory')
    ## Step 5 only
    parser.add_argument('--mirror_raw_data', type=bool, default=True, help='Mirror raw data files to GNPS2')
    parser.add_argument('--mirror_mzmine_results', type=bool, default=True, help='Mirror mzmine results files to GNPS2')
    ## Step 7 only
    parser.add_argument('--download_dir', type=str, default='/global/cfs/cdirs/metatlas/projects/untargeted_outputs', help='Path to the download folder')
    parser.add_argument('--overwrite_zip',type=bool, default=False, help='Overwrite existing zip files in download folder')
    parser.add_argument('--overwrite_drive', type=bool, default=False, help='Overwrite existing zip files on google drive')
    parser.add_argument('--gdrive_upload', type=bool, default=True, help='Upload to google drive')
    parser.add_argument('--min_features', type=int, default=0, help='Set minimum number of MZmine features for a project polarity to be zipped')
    parser.add_argument('--add_gnps2_documentation', type=bool, default=True, help='File name of the GNPS2 documentation to add to project zips')
    parser.add_argument('--gnps2_doc_name', type=str, default='Untargeted_metabolomics_GNPS2_Guide.docx', help='File name of the GNPS2 documentation to add to project zips')
    parser.add_argument('--abridged_filenames', type=bool, default=True, help='Use abridged filenames in the zipped folders. Can be set to False for non-conforming project names')
    ## Logger/helper
    parser.add_argument('--log_file', type=str, default='/global/cfs/cdirs/m2650/untargeted_logs/untargeted_pipeline.log', help='Log file name with full path')
    parser.add_argument('--log_level', type=str, default='INFO', help='Logger level. One of [DEBUG, INFO, WARNING, ERROR, or CRITICAL]')
    parser.add_argument('--log_format', type=str, default='%(asctime)s - %(levelname)s - %(message)s', help='Logger format')

def check_args(args):
    ##### Check if the input arguments are valid
    if args.direct_input:
        args.direct_input = args.direct_input.split(',')
    if args.background_designator:
        args.background_designator = args.background_designator.split(',')
    if args.skip_steps:
        args.skip_steps = args.skip_steps.split(',')
    if args.skip_steps is not None and args.direct_input is None:
        logging.error('Incompatible flags. Must provide direct input if you want to skip steps.')
        sys.exit(1)
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

def check_skipped_steps(args):
    ##### Check for steps to skip during direct input
    if args.skip_steps is not None and args.direct_input is not None:
        step_skip_bool_list = [False] * 7
        if 'Sync' in args.skip_steps:
            step_skip_bool_list[0] = True
        if 'MZmine_Status' in args.skip_steps:
            step_skip_bool_list[1] = True
        if 'MZmine_Submit' in args.skip_steps:
            step_skip_bool_list[2] = True
        if 'FBMN_Status' in args.skip_steps:
            step_skip_bool_list[3] = True
        if 'FBMN_Submit' in args.skip_steps:
            step_skip_bool_list[4] = True
        if 'FBMN_Download' in args.skip_steps:
            step_skip_bool_list[5] = True
        if 'Zip_and_Upload' in args.skip_steps:
            step_skip_bool_list[6] = True
    else:
        step_skip_bool_list = [False] * 7
    return step_skip_bool_list
    
if __name__ == "__main__":
    main()