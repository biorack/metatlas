# Untargeted Analysis Workflow

## Overview

The metatlas untargeted analysis pipeline currently operates as a single cronjob that is run
on `perlmutter` by the NERSC pseudo user `msdata`. The `msdata` cronjob is named `run_untargeted_pipeline`
and is (currently) initiated every 6 hours (4x daily) on the 45th minute. This pipeline takes
raw metabolite data (`*.mzML`) and runs third-party software workflows to create a suite of
untargeted metabolomics results for JGI users.

## Glossary and Acronyms

- NERSC (National Energy Research Scientific Computing Center): used for some computational analysis
  in the untargeted pipeline; we use the perlmutter nodes and the `msdata` pseudo user
- GNPS2 (Global Natural Product Social Molecular Networking): a community for natural product researchers 
  working with mass spectrometry data; we currently use version 2
- MZmine: open-source software for mass spectrometry data processing; we currently use version 3 intalled
  on perlmutter. This is the first major tool used in the pipeline.
- FBMN (feature-based molecular networking): a module of GNPS2 which runs molecular networking, a method to
  visualize and annotate the chemical space in non-targeted mass spectrometry data. This is the second
  major tool used in the pipeline.

## Pipeline

The untargeted pipeline is split into several discrete steps which all get run during every
`run_untargeted_pipeline` cronjob. While the pipeline by default gets run in its entirety each time it's called,
a given project will take 3 full runs to complete (see details of each [step](#step-1-cycle-1) below). This
is because steps that submit compute jobs (Steps 3 and 5, MZmine and FBMN submissions) don't continually
check whether the job is finished, but rather wait until the next scheduled run to check and update
the project's status accordingly. Since the pipeline is run every 6 hours (current cronjob settings), a project
will have complete untargeted results within about 24 hours. Built-in [flags](#pipeline-flags) allow for any step
to be re-run in case of an error.

In summary, for each project, during the first cycle (Steps 1-3) a results directory will be made on disk at NERSC and 
an MZmine job will be submitted via slurm; during the second cycle (Steps 4-5) an FBMN job will be submitted at GNPS2
using the results of a successful MZmine run; and during the third cycle (Steps 6-7) all results files will be compiled
and the results directory will be zipped and uploaded to Google Drive for the user. Each cycle, project status will
be also be updated in the LIMS `untargeted_tasks` table to keep track of successes/failures and notify analysts of any
behavior that might require intervention (see [notifications](#email-notifications) section below). The full pipeline
code is in the `run_untargeted_pipeline.py` script and functions that define each pipeline step are in the `tools.py`
script, both on NERSC in `/global/common/software/m2650/metatlas-repo/metatlas/untargeted/`. The overview and some 
important checkpoints and features of each step are described below.

### Step 1 (Cycle 1)

`update_new_untargeted_tasks()`

0. Upstream of Step 1, an `msdata` cronjob called `untargeted_lims_update` runs every 30 minutes and ensures
   that raw data files (`*.mzML`) are synced between NERSC disk and LIMS before initiating a new untargeted task.
1. Find all projects that have raw data in the LIMS that are over 3 hours old (indicating that
   the transfer from NERSC to LIMS is complete) but do not have a row in the LIMS `untargeted_tasks`
   table yet. Use raw data to determine which polarity(ies) should be initiated.
2. Create a row in the LIMS `untargeted_tasks` table and create a matching `untargeted_tasks` 
   directory at NERSC for each polarity. Into the latter directory, write files necessary for MZmine submission
   at NERSC (i.e., `*mzmine.sh`, `*.sbatch`, `*metadata.tab`, `*batch-params.xml`, `*filelist.txt`).
3. Set the status of the new projects/polarities in the LIMS `untargeted_tasks` table to `initiation` for
   MZmine and `waiting` for FBMN.

### Step 2 (Cycle 1)

`update_mzmine_status_in_untargeted_tasks()`

1. Find all projects in the LIMS `untargeted_tasks` table that have an MZmine status of `initation`, `running`, or `error`
   for each relevant polarity
2. Look inside its `untargeted_tasks` directory at NERSC for MZmine results.
3. If it has all expected results files, update status to `complete` in the LIMS `untargeted_tasks` table.
4. If the project gets moved from `initation`, `running`, or `error` to `complete` (or if the `--direct_input` flag is used),
   calculate the number of features and number of msms hits and update the relevant columns of the LIMS`untargeted_tasks` table.
   Then, create a filtered peak height table with various rules for retaining features (see [flags](#pipeline-flags) for options).

### Step 3 (Cycle 1)

`submit_mzmine_jobs()`

1. By default, submit *only new* projects which were initated during Step 1. If you want to submit an MZmine job manually for
   a given project, you can the `run_untargeted_pipeline.py` script manually via the NERSC command line and pass a comma-separated
   list of project names with the `--direct_input` flag (see [examples](#example-1) below).
2. Check for existing MZmine results in the NERSC `untargeted_tasks` directory for each relevant polarity and do not re-submit
   MZmine if they're found. Note that this condition can be superceeded with a `--overwrite_mzmine` flag (see [flags](#pipeline-flags)).
3. Use the `*mzmine.sh` and `*.sbatch` files created in Step 1 to submit an MZmine slurm job at NERSC.
4. If submitted successfully, update the LIMS `untargeted_tasks` table to `running` for that project/polarity.
5. Analysts can use `sqs` or `squeue -u msdata` to manually monitor the progress of submitted MZmine jobs.

### Step 4 (Cycle 2)

`update_fbmn_status_in_untargeted_tasks()`

1. Find all projects in the LIMS `untargeted_tasks` table that have an FBMN status of `waiting`, `running`, or `error`
   for each relevant polarity.
2. Query the GNPS2 task that is running FBMN analysis (using a `*_gnps2-fbmn-task.txt` file in the `untargeted_tasks`
   directory at NERSC for a project/polarity) and parse the HTML output to find its status at GNPS2.
3. Update the the LIMS `untargeted_tasks` table according to the GNPS2 status.

### Step 5 (Cycle 2)

`submit_fbmn_jobs()`

1. Find all projects in the LIMS `untargeted_tasks` table that have an FBMN status of `waiting` or `error`, plus
   an MZmine status of `complete`, for each relevant polarity.
2. Before submitting FBMN jobs to GNPS2, check that the MZmine results (which are required to run FBMN) have been
   successfully mirrored to GNPS2. This mirror is accomplished by an `msdata` cronjob at NERSC called `sync_untargeted_mzmine_results_to_gnps2`
   that runs every 30 minutes. While this cronjob is set to run 15 mintues before the untargeted pipeline, in case
   Step 5 is initated before the data are mirrored the `submit_fbmn_jobs()` function will run the mirror manually.
   FBMN requires the mirrored `*quant.csv` (peak area), `*mgf` (msms data), and `*metadata.tab` (experimental design) files.
3. Submit an FBMN task at GNPS2 and write the task ID to a file (`*_gnps2-fbmn-task.txt`) in the project/polarity `untargeted_tasks`
   directory at NERSC.
4. Update the the LIMS `untargeted_tasks` table to `running` if the FBMN task was submitted successfully.

### Step 6 (Cycle 3)

`download_fbmn_results()`

1. Check for projects in the LIMS `untargeted_tasks` table which have an MZmine and FBMN status of `complete` for all
   relevant polarities, but that do not yet have FBMN results in the NERSC `untargeted_tasks` directory. Note that
   the latter condition can be superceeded with a `--overwrite_fbmn` flag (see [flags](#pipeline-flags)).
2. Download the FBMN file (`*_gnps2-fbmn-network.graphml`) and the results table file (`*_gnps2-fbmn-library-results.tsv`)
   from GNPS2 using an API call with the task ID (`*_gnps2-fbmn-task.txt`). After successfully downloading these files, create a
   `*_gnps2-page-link.txt` file so the user can navitage to the GNPS2 page with their full results.

### Step 7 (Cycle 3)

`zip_and_upload_untargeted_results()`

1. Check for projects in the LIMS `untargeted_tasks` table which have an MZmine and FBMN status of `complete` for all
   relevant polarties, and also check that the MZmine and FBMN files have been generated and contain data. For example,
   the `*.mgf` file must not be empty and the `*peak-height.csv` file must have at least one feature (minimum number of
   features allowed to pass can be changed with a flag; this can be modified with the `--min_features` flag).
2. Optionally (but by default), download the latest GNPS2 documentation from Google Drive and place it into the
   `untargeted_tasks` directory at NERSC.
3. Zip all polarity directories together along with the GNPS2 documentation and place the archive (named after the project)
   into the `untargeted_outputs` directory at NERSC. If the archive for a project already exists it will not be overwritten
   unless the `--overwrite_zip` flag is used.
4. Optionally (but by default), upload the archive to Google Drive if the zip is successful (automatic uploading can be
   disabled with the `--gdrive_upload` flag).

## Pipeline Flags

Use the `--help` flag to print out all possible arguments that can be used to modulate/customize the pipeline:

```/global/common/software/m2650/python3-matchms/bin/python /global/homes/b/bkieft/metatlas/metatlas/untargeted/run_untargeted_pipeline.py --help```

There are over a dozen flags that can be passed to the `run_untargeted_pipeline.py` script, but they all
have default values which get passed during a typical run with the automated cronjob. If the untargeted pipeline - 
or any individual steps therein - needs to be run outside the cronjob (e.g., to fasttrack a project or to fix
an error that occurred during the automated job), the flags come in handy. Here are some examples:

### Example 1
```
/global/common/software/m2650/python3-matchms/bin/python /global/homes/b/bkieft/metatlas/metatlas/untargeted/run_untargeted_pipeline.py
--direct_input <project_1>,<project_2>,<project_3> --overwrite_mzmine True --overwrite_fbmn True --log_file </my/home/directory/custom_untargeted_run.log>
```
This command will run all steps of the pipeline on just the three provided projects using the `--direct_input` flag with
a comma-separated list of valid project names. Since you know that MZmine and FBMN have already been run on one or more of these projects,
the `--overwrite_mzmine` and `--overwrite_fbmn` flags set to True are used to replace any existing MZmine and FBMN data on disk at NERSC with new data. 
**An important note: since the untargeted pipeline is run in three cycles, you may have to run this command multiple times to get all steps**
**to properly run since you're resubmitting MZmine and FBMN jobs**. For this reason, the `--log_file` flag is passed to monitor the progress of the pipeline
in a personal directory. Make sure that the group permissions for the log file are set to `metatlas` so the pipeline can write to the log.

### Example 2

```
/global/common/software/m2650/python3-matchms/bin/python /global/homes/b/bkieft/metatlas/metatlas/untargeted/run_untargeted_pipeline.py
--direct_input <project_1> --min_features 1000 --log_file </my/home/directory/custom_untargeted_run.log> --overwrite_zip True --overwrite_drive True
--skip_steps Sync,MZmine_Status,MZmine_Submit,FBMN_Status,FBMN_Submit,FBMN_Download
```
This command will run only the final step (Step 7) on project_1 data, and only upload the results to Google Drive if the number of
features in the MZmine `*peak-height.csv` is greater than 1,000. This command might be run if you want a new updated
results archive to show up on Google Drive for the user (hence the `--overwrite_zip` and `--overwrite_drive` flags), but only if it meets the
minimum threshold.

### Example 3

```
/global/common/software/m2650/python3-matchms/bin/python /global/homes/b/bkieft/metatlas/metatlas/untargeted/run_untargeted_pipeline.py
--direct_input <project_1,project_2,project_3> --background_ratio 10 --polar_solvent_front 1.0 --overwrite_zip True --overwrite_drive True
--min_features 100 --skip_steps Sync,MZmine_Submit,FBMN_Status,FBMN_Submit,FBMN_Download
```
This command will run Steps 2 and 7 on project_1, project_2, and project_3. During Step 2, a new `*peak-height-filtered.csv`
file will be produced using the new (non-default) parameters provided in the `--background_ratio` and `--polar_solvent_front` flags. Then,
a new archive (with the new filtered peak height table) will be zipped and uploaded to Google Drive. This more complex
command is useful if you want to regenerate new results based on some observations you made about the existing data. In this example, the
background ratio is set very high and the solvent front is scaled up to be conservative about which features remain in the filtered table,
then the `--min_features` flag is used to check that at least 100 features passed these new thresholds.

## Email Notifications

### Warning/Error Notifications for Automated Pipeline

Because the untargeted pipeline runs as an `msdata` cronjob, errors and warnings that may need attention are emailed directly to an analyst
(default: `bkieft@lbl.gov`). To proactively catch warnings and errors early that might break the pipeline, the pipeline log file
(currently at `/global/cfs/cdirs/m2650/untargeted_logs/untargeted_pipeline.log`) is queried every *N* number of pipeline cycles.
The warning/error query happens via an `msdata` cronjob called `untargeted_warnings-errors_email_report` that currently runs once per day
and is set to 4 cycles because the pipeline runs every 6 hours - so, a check every 24 hours. ERROR and WARNING logging lines are
reported in the email for the analyst to double-check.

### Progress Notifications for Automated Pipeline

Once per week, the number of successful and failed projects running through the automated untargeted pipeline for the previous *N*
days will be counted up and emailed to a subset of the metabolomics team. This report is generated via the `msdata` cronjob
`untargeted_project_email_report`. The job uses the script `check_untargeted_status.py`.