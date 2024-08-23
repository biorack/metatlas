# Untargeted Analysis Workflow

## Overview

The metatlas untargeted analysis pipeline currently operates as a single cronjob that is run
on `perlmutter` by the NERSC pseudo user `msdata`. The `msdata` cronjob is named `run_untargeted_pipeline`
and is (currently) initiated every 6 hours (4x daily) on the 45th minute. This pipeline takes
raw metabolite data (`*.mzML`) and runs several third-party software workflows to create a suite of
untargeted metabolomics results for JGI users.

## Glossary and Acronyms

- NERSC (National Energy Research Scientific Computing Center): used for some computational analysis
  in the untargeted pipeline; we use the perlmutter nodes and the `msdata` pseudo user
- GNPS2 (Global Natural Product Social Molecular Networking): a community for natural product researchers 
  working with mass spectrometry data; we currently use version 2
- FBMN (feature-based molecular networking): a module of GNPS2 for discovering the identity of
  poorly-characterized molecules
- MZmine: open-source software for mass spectrometry data processing; we currently use version 3

## Pipeline

The untargeted pipeline is split into several discrete steps which all get run during the
`run_untargeted_pipeline` job. While the pipeline gets run in its entirety each cycle, a given project
will take 3 full cycles to complete. Since the pipeline is (currently) run every 6 hours, a project
should have complete untargeted results within about 24 hours, with built-in [flags](#pipeline-flags) 
that allow for each cycle to be re-run in case of an error.

In summary, during the first cycle a results directory will be made on disk at NERSC and 
an MZmine job will be submitted via slurm; during the second cycle an FBMN job will be submitted at GNPS2; during
the third cycle the FBMN files will be downloaded and the full results directory will be zipped and uploaded to Google
Drive for the user. The functions that define each pipeline step are in the `tools.py` script at NERSC in
`/global/common/software/m2650/metatlas-repo/metatlas/untargeted/`. The basic summary and some important 
checkpoints/features of each step (function) are described below.

### Step 1 (Cycle 1)

`update_new_untargeted_tasks()`

1. Find all projects that have raw data in the LIMS that are over 3 hours old (indicating that
   the transfer from NERSC to LIMS is complete) but do not have a row in the LIMS `untargeted_tasks`
   table yet. Use raw data to determine which polarity(ies) are present. Upstream of Step 1, another `msdata`
   cronjob called `untargeted_lims_update` runs every 30 minutes and ensures that raw data files
   are synced between NERSC disk and LIMS before initiating a new untargeted task.
2. Create a row in the LIMS `untargeted_tasks` table and create a matching `untargeted_tasks` 
   directory at NERSC for each polarity. Into the latter directory, write files necessary for MZmine submission
   at NERSC (i.e., `*mzmine.sh`, `*.sbatch`, `*metadata.tab`, `*batch-params.xml`, `*filelist.txt`).
3. Set the status of the new projects/polarities in the LIMS `untargeted_tasks` table to `initiation` for
   MZmine and `waiting` for FBMN.

### Step 2 (Cycle 1)

`update_mzmine_status_in_untargeted_tasks()`

1. Find all projects in the LIMS `untargeted_tasks` table that have an MZmine status of `initation`, `running`, or `error`.
2. For each polarity in each of these projects, look inside its `untargeted_tasks` directory at NERSC
   for MZmine results.
3. If it has them all, update status to `complete` in the LIMS `untargeted_tasks` table . Do some other checks
   (e.g., if the MZmine results aren't there yet and its status is `running`, use an `squeue` command to
   ensure it's actually still running at NERSC).
4. If the project/polarity updates to `complete` during the function call, calculate number of features and 
   number of msms hits and update the LIMS`untargeted_tasks`  table. Then, create a filtered peak height table 
   if samples of the background (control, default to `ExCtrl`) condition exist.

### Step 3 (Cycle 1)

`submit_mzmine_jobs()`

1. By default, submit *new* projects which were initated during Step 1. You can change this by
   passing a list of project names to the `--direct_input` flag when running the `run_untargeted_pipeline.py` script
   via command line.
2. When the newly initated projects (default) are passed to the function, use the `*mzmine.sh` and `*.sbatch` files
   created in Step 1 to submit an MZmine slurm job at NERSC.
3. If submitted successfully, update the LIMS `untargeted_tasks` table to `running` for that project/polarity.
4. Analysts can use `sqs` or `squeue -u msdata` to monitor the progress of submitted MZmine jobs.

### Step 4 (Cycle 2)

`update_fbmn_status_in_untargeted_tasks()`

1. Find all projects in the LIMS `untargeted_tasks` table that have an FBMN status of `waiting`, `running`, or `error`.
2. Query the GNPS2 task that is running FBMN analysis (using a `*_gnps2-fbmn-task.txt` file in the `untargeted_tasks`
   directory at NERSC for a project/polarity) and parse the HTML output to find its status at GNPS2.
3. Update the the LIMS `untargeted_tasks` table according to the GNPS2 status.

### Step 5 (Cycle 2)

`submit_fbmn_jobs()`

1. Find all projects in the LIMS `untargeted_tasks` table that have an FBMN status of `waiting`, or `error` and
   an MZmine status of `complete`.
2. Before submitting FBMN jobs to GNPS2, check that the MZmine results (which are needed to run FBMN) have been
   successfully mirrored to GNPS2. This mirror is accomplished by an `msdata` cronjob at NERSC called `sync_untargeted_mzmine_results_to_gnps2`,
   which runs every 30 minutes. While this cronjob is set to run 15 mintues before the untargeted pipeline, in case
   Step 5 is initated before the data are mirrored the `submit_fbmn_jobs()` function will run the mirror manually.
3. Submit an FBMN task at GNPS2 and capture the task ID as a file (`*_gnps2-fbmn-task.txt`) in the project/polarity `untargeted_tasks`
   directory at NERSC. The FBMN task requires the mirrored `*quant.csv` (peak area), `*mgf` (msms data), and `*metadata.tab` files.
4. Update the the LIMS `untargeted_tasks` table to `running` if the FBMN task was submitted successfully.

### Step 6 (Cycle 3)

`download_fbmn_results()`

1. Check for projects in the LIMS `untargeted_tasks` table which have an MZmine and FBMN status of `complete` for all
   relevant polarties.
2. Download the FBMN file (`*_gnps2-fbmn-network.graphml`) and the results table file (`*_gnps2-fbmn-library-results.tsv`)
   from GNPS2 using the task ID from `*_gnps2-fbmn-task.txt`. After successfully downloading these files, create a
   `*_gnps2-page-link.txt` file for the user to get to the GNPS2 page with their full results.

### Step 7 (Cycle 3)

`zip_and_upload_untargeted_results()`

1. Check for projects in the LIMS `untargeted_tasks` table which have an MZmine and FBMN status of `complete` for all
   relevant polarties, and also checks that the MZmine and FBMN files have been generated and contain data. For example,
   the `*.mgf` file must not be empty and the `*peak-height.csv` file must have at least one feature (minimum number of
   features allowed to pass can be changed with a flag).
2. Download the latest GNPS2 documentation from Google Drive and places it into the project/polarity
   `untargeted_tasks` directory at NERSC.
3. Zip all polarity directories together along with the GNPS2 documentation (can be turned off with a flag) and place
   the archive (named after the project) into the `untargeted_outputs` directory at NERSC.
4. Automatically upload the archive to Google Drive if the zip is successful (can be turned off with a flag).

## Pipeline Flags

There are over a dozen flags that can be passed to the `run_untargeted_pipeline.py` script, but they all
have default values which get passed during a typical run with the automated cronjob. If the untargeted pipeline - 
or any individual steps therein - needs to be run outside the cronjob (e.g., to fasttrack a project or to fix
an error that occurred during the automated job), the flags come in handy. Here are some examples:

### Example 1
```
/global/common/software/m2650/python3-matchms/bin/python /global/homes/b/bkieft/metatlas/metatlas/untargeted/run_untargeted_pipeline.py
--direct_input <project_1>,<project_2>,<project_3> --overwrite_mzmine True --overwrite_fbmn True --log_file </my/home/directory/custom_untargeted_run.log>
```
This command will run all steps of the pipeline on just the three provided projects using the `--direct_input` flag and providing
a comma-separated list of valid project names. Since you know that MZmine and FBMN have already been run on one or more of these projects,
the `--overwrite_mzmine` and `--overwrite_fbmn` flags are used to replace any existing MZmine and FBMN data on disk at NERSC with new data. 
**An important note: since the untargeted pipeline is run in three cycles, you may have to run this command multiple times to get all steps**
**to properly run**. For this reason, the `--log_file` flag is passed to monitor the progress of the pipeline. Make sure that the group
permissions for the log file are set to `metatlas`.

### Example 2

```
/global/common/software/m2650/python3-matchms/bin/python /global/homes/b/bkieft/metatlas/metatlas/untargeted/run_untargeted_pipeline.py
--direct_input <project_1> --overwrite_zip True --overwrite_drive True --skip_steps Sync,MZmine_Status,MZmine_Submit,FBMN_Status,FBMN_Submit,FBMN_Download
```
This command will run only the final step (Step 7: Zip_and_Upload) on `project_1` data. This command might be run if you want new files
or results to show up for the user in the zip folder on Google Drive (hence the `--overwrite_zip` and `--overwrite_drive` flags).

### Example 3

```
/global/common/software/m2650/python3-matchms/bin/python /global/homes/b/bkieft/metatlas/metatlas/untargeted/run_untargeted_pipeline.py
--direct_input <project_1> --min_features 1000 --log_file </my/home/directory/custom_untargeted_run.log>
```
This command will run all steps of the pipeline on just `project_1` data, and only successfully zip and upload the results folder if
at least 1,000 features were identified in the MZmine peak height table. **As in Example 1, you would need to monitor the log file**
**and run this command until all the results/uploads you want are complete**

## Warning/Error Notifications for Automated Pipeline

Because the untargeted pipeline runs as an `msdata` cronjob, errors and warnings that may need attention are emailed directly to an analyst
(currently `bkieft@lbl.gov`). To proactively catch warnings and errors early that might break the pipeline, the pipeline log file
(currently at `/global/cfs/cdirs/m2650/untargeted_logs/untargeted_pipeline.log`) is queried every *N* number of pipeline cycles.
The warning/error query happens via an `msdata` cronjob called `untargeted_warnings-errors_email_report` that currently runs once per day
and is set to 4 cycles because the pipeline runs every 6 hours - so, a check every 24 hours. ERROR and WARNING logging lines are
reported in the email for the analyst to double-check.

## Progress Notifications for Automated Pipeline

Once per week, the number of successful and failed projects running through the automated untargeted pipeline for the previous *N*
days will be counted up and emailed to a subset of the metabolomics team. This report is generated via the `msdata` cronjob
`untargeted_project_email_report`. The job uses the script `check_untargeted_status.py`.