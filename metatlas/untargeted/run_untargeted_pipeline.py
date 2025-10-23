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
    mzm.call_logger(log_filename=args.log_file, log_level=args.log_level, log_format=args.log_format, log_to_stdout=args.log_to_stdout)

    ##### Kick off the script
    start_message = mzm.start_script(script="run_untargeted_pipeline.py")
    logging.info('\n')
    logging.info(start_message)

    ##### Write args to log for reference
    logging.info(f'Arguments used: {args}')

    ##### Step 1/7: Syncing LIMS and NERSC to identify new projects with raw data that are not yet in the untargeted task list
    new_projects = mzm.update_new_untargeted_tasks(
        background_designator=args.background_designator,
        validate_names=args.validate_names,
        skip_sync=step_bools[0],
        output_dir=args.output_dir,
        raw_data_dir=args.raw_data_dir,
        direct_input=args.direct_input,
        custom_mzmine_batch_params=args.custom_mzmine_batch_params,
        raw_data_subdir=args.raw_data_subdir,
        skip_blank_filter=args.skip_blank_filter,
        fps_files_only=args.fps_files_only,
        project_tag=args.project_tag
    )

    ##### Step 2/7: Checking and updating status of MZmine jobs in LIMS
    mzm.update_mzmine_status_in_untargeted_tasks(
        background_designator=args.background_designator,
        skip_mzmine_status=step_bools[1],
        direct_input=args.direct_input,
        background_ratio=args.background_ratio,
        zero_value=args.zero_value,
        nonpolar_solvent_front=args.nonpolar_solvent_front,
        polar_solvent_front=args.polar_solvent_front,
        nonpolar_solvent_end=args.nonpolar_solvent_end,
        polar_solvent_end=args.polar_solvent_end,
        project_tag=args.project_tag
    )

    ##### Step 3/7: Submitting new MZmine jobs that are "initialized"
    mzm.submit_mzmine_jobs(
        new_projects=new_projects,
        direct_input=args.direct_input,
        skip_mzmine_submit=step_bools[2],
        overwrite_mzmine=args.overwrite_mzmine,
        project_tag=args.project_tag
    )

    ##### Step 4/7: Checking and updating status of FBMN jobs in LIMS
    mzm.update_fbmn_status_in_untargeted_tasks(
        direct_input=args.direct_input,
        skip_fbmn_status=step_bools[3],
        project_tag=args.project_tag
    )

    ##### Step 5/7: Submitting new FBMN jobs to GNPS2
    mzm.submit_fbmn_jobs(
        output_dir=args.output_dir,
        overwrite_fbmn=args.overwrite_fbmn,
        direct_input=args.direct_input,
        skip_fbmn_submit=step_bools[4],
        skip_mirror_raw_data=args.skip_mirror_raw_data,
        skip_mirror_mzmine_results=args.skip_mirror_mzmine_results,
        raw_data_dir=args.raw_data_dir,
        raw_data_subdir=args.raw_data_subdir,
        project_tag=args.project_tag
    )

    ##### Step 6/7: Checking for completed FBMN jobs and downloading results
    mzm.download_fbmn_results(
        output_dir=args.output_dir,
        overwrite_fbmn=args.overwrite_fbmn,
        direct_input=args.direct_input,
        skip_fbmn_download=step_bools[5],
        project_tag=args.project_tag
    )

    ##### Step 7/7: Zipping up and (optionally) uploading output folders to gdrive
    mzm.zip_and_upload_untargeted_results(
        download_folder=args.download_folder,
        output_dir=args.output_dir,
        doc_name=args.doc_name,
        add_documentation=args.add_documentation,
        skip_zip_upload=step_bools[6],
        abridged_filenames=args.abridged_filenames,
        upload=args.upload,
        overwrite_zip=args.overwrite_zip,
        overwrite_drive=args.overwrite_drive,
        direct_input=args.direct_input,
        min_features_admissible=args.min_features_admissible,
        project_tag=args.project_tag
    )

    ##### End the script
    end_message = mzm.end_script(script="run_untargeted_pipeline.py")
    logging.info(end_message)

def add_arguments(parser):
    ## Multiple functions
    parser.add_argument('--output_dir', type=str, default='/global/cfs/cdirs/metatlas/projects/untargeted_tasks', help='Path to the output directory')
    parser.add_argument('--direct_input', type=str, default=None, help='Input project names from command line as a CSV list')
    parser.add_argument('--project_tag', type=str, default=None, help='Tag to append to project names for running multiple parameter sets (e.g., with flag --custom_mzmine_batch_params). If given, project names will be modified as PROJECTNAME_TAG')
    parser.add_argument('--background_designator', type=str, default="ExCtrl", help='Input control/background sample names from command line as a CSV list')
    parser.add_argument('--skip_steps', type=str, default=None, help='Skip pipeline step(s) with --direct_input. CSV list with elements [Sync,MZmine_Status,MZmine_Submit,FBMN_Status,FBMN_Submit,FBMN_Download,Zip_and_Upload]')
    parser.add_argument('--raw_data_dir', type=str, default='/global/cfs/cdirs/metatlas/raw_data', help='Path to the raw data directory')
    parser.add_argument('--raw_data_subdir', type=str, default=None, help='Name of the raw_data subdirectory (e.g., jgi, egsb) to use for tasks that require raw data files. If not given, will try to infer from the second field of the project name.')
    parser.add_argument('--overwrite_fbmn', action='store_true', help='Overwrite existing fbmn results files that are already in the output directory')
    ## Step 1 only
    parser.add_argument('--validate_names', action='store_true', help='Validate filenames and project names')
    parser.add_argument('--custom_mzmine_batch_params', type=str, default=None, help='Full path to custom mzmine batch parameters xml. If using FPS only mode, supply a csv list of pos and neg parameter files')
    parser.add_argument('--skip_blank_filter', action='store_true', help='Do not filter out files with "Blank" in the name from the untargeted task list')
    parser.add_argument('--fps_files_only', action='store_true', help='Only FPS files will be input, so do not check for polarity in file name and use custom mzmine batch parameters')
    ## Step 2 only
    parser.add_argument('--background_ratio', type=float, default=3, help='Ratio of background to sample intensity for filtering features')
    parser.add_argument('--zero_value', type=float, default=(2/3), help='Proportion of the lowest intensity value from the experiment to use as replacement zero value')
    parser.add_argument('--polar_solvent_front', type=float, default=0.8, help='Retention time to use as HILIC solvent front (mins) for filtering features')
    parser.add_argument('--nonpolar_solvent_front', type=float, default=0.5, help='Retention time to use as C18/LIPID solvent front (mins) for filtering features')
    parser.add_argument('--polar_solvent_end', type=float, default=17.5, help='Retention time to use as HILIC solvent end (mins) for filtering features')
    parser.add_argument('--nonpolar_solvent_end', type=float, default=10, help='Retention time to use as C18/LIPID solvent end (mins) for filtering features')
    ## Step 3 only
    parser.add_argument('--overwrite_mzmine', action='store_true', help='Overwrite existing mzmine results files that are already in the output directory')
    ## Step 5 only
    parser.add_argument('--skip_mirror_raw_data', action='store_true', help='Skip the default mirror of raw data files to GNPS2')
    parser.add_argument('--skip_mirror_mzmine_results', action='store_true', help='Skip the default mirror of mzmine results files to GNPS2')
    ## Step 7 only
    parser.add_argument('--download_dir', type=str, default='/global/cfs/cdirs/metatlas/projects/untargeted_outputs', help='Path to the download folder')
    parser.add_argument('--min_features', type=int, default=0, help='Set minimum number of MZmine features for a project polarity to be zipped')
    parser.add_argument('--gnps2_doc_name', type=str, default='Untargeted_metabolomics_GNPS2_Guide.docx', help='File name of the GNPS2 documentation to add to project zips')
    parser.add_argument('--overwrite_zip', action='store_true', help='Overwrite existing zip files in download folder')
    parser.add_argument('--overwrite_drive', action='store_true', help='Overwrite existing zip files on google drive')
    parser.add_argument('--gdrive_upload', action='store_true', help='Upload to google drive')
    parser.add_argument('--no_gdrive_upload', action='store_false', dest='gdrive_upload', help='Do not upload to google drive')
    parser.add_argument('--add_gnps2_documentation', action='store_true', help='Add GNPS2 documentation to project zips')
    parser.add_argument('--no_add_gnps2_documentation', action='store_false', dest='add_gnps2_documentation', help='Do not add GNPS2 documentation to project zips')
    parser.add_argument('--abridged_filenames', action='store_true', help='Use abridged filenames in the zipped folders. Can be set to False for non-conforming project names')
    parser.add_argument('--no_abridged_filenames', action='store_false', dest='abridged_filenames', help='Do not use abridged filenames in the zipped folders')
    parser.set_defaults(gdrive_upload=True)
    parser.set_defaults(add_gnps2_documentation=True)
    parser.set_defaults(abridged_filenames=True)
    ## Logger/helper
    parser.add_argument('--log_file', type=str, default='/global/cfs/cdirs/m2650/untargeted_logs/untargeted_pipeline.log', help='Log file name with full path')
    parser.add_argument('--log_level', type=str, default='INFO', help='Logger level. One of [DEBUG, INFO, WARNING, ERROR, or CRITICAL]')
    parser.add_argument('--log_format', type=str, default='%(asctime)s - %(levelname)s - %(message)s', help='Logger format')
    parser.add_argument('--log_to_stdout', action='store_true', help='Log to stdout instead of a file')

def check_args(args):
    ##### Check if the input arguments are valid
    if args.fps_files_only and args.custom_mzmine_batch_params is None:
        logging.error('FPS files only flag requires custom mzmine batch parameters list (e.g., POS-, NEG- prefixes). Please check flags.')
        sys.exit(1)
    if args.direct_input:
        args.direct_input = args.direct_input.split(',')
    if args.background_designator:
        args.background_designator = args.background_designator.split(',')
    if args.skip_steps:
        args.skip_steps = args.skip_steps.split(',')
    if args.custom_mzmine_batch_params:
        args.custom_mzmine_batch_params = args.custom_mzmine_batch_params.split(',')
    if args.overwrite_drive is True and args.gdrive_upload is False:
        logging.error('Incompatible flags. Cannot overwrite google drive if not uploading to google drive.')
        sys.exit(1)
    if not os.path.exists(args.download_dir):
        logging.error('Download directory does not exist. Please check flag and that you are at NERSC.')
        sys.exit(1)
    if not os.path.exists(args.output_dir):
        logging.error('Output directory does not exist. Please check flag and that you are at NERSC.')
        sys.exit(1)
    if not os.path.exists(args.raw_data_dir):
        logging.error('Raw data directory does not exist. Please check flag and that you are at NERSC.')
        sys.exit(1)

def check_skipped_steps(args):
    ##### Check for steps to skip during direct input
    if args.skip_steps is not None:
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