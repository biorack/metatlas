# Targeted HILIC/Polar Analysis Workflow

## One-time setup

### RClone configuration

#### For MacOS/Linux

Open a terminal on Cori and run:
`
mkdir -p ~/.config/rclone
`

Open a terminal on your local machine and run:
```
curl --silent --show-error https://rclone.org/install.sh | sudo bash > /dev/null
# You will be prompted to enter your password, this allows the installation of rclone
rclone config create metabolomics drive root_folder_id 0B-ZDcHbPi-aqZzE5V3hOZFc0dms
# You will be prompted in your web browser to grant rclone access to Google Drive
scp $(rclone config file | tail -1) dtn01.nersc.gov:~/.config/rclone/rclone.conf
```

Now you can verify that rclone is working correctly by running this on Cori:
```
/global/cfs/cdirs/m342/USA/shared-repos/rclone/bin/rclone lsd metabolomics:Analysis_uploads
```

which should yield a listing of metabolomics experiment names similar to:

```
          -1 2021-08-30 10:01:06        -1 20210323_JGI-AK_SS_504264_GEBA_Pantoea-final_QE-HF_HILICZ_USHXG01602
          -1 2021-08-30 12:32:39        -1 20210518_JGI-AK_IG-SS_503256_BETO_Pceleri_QE-HF_HILICZ_USHXG01602
          -1 2021-09-13 16:39:15        -1 20210721_JGI-AK_JB_504782_PseudoOphi_final_QE-139_HILICZ_USHXG01490
          -1 2021-09-13 17:40:55        -1 20210723_JGI-AK_DB-TM_506963_LemCreek_final_QE-HF_HILICZ_USHXG01494
          -1 2021-09-13 16:39:15        -1 20210728_JGI-AK_MD_507130_Bioscales_pilot2_QE-139_HILICZ_USHXG01490
          -1 2021-09-10 16:05:18        -1 20210804_JGI-AK_PA-CT_507784_Frtlzr_Set1_QE-139_HILICZ_USHXG01490
          -1 2021-09-13 16:34:45        -1 20210819_JGI-AK_MK_506588_SoilWaterRep_final_QE-139_HILICZ_USHXG01490
```

#### For Windows

1. Download and unzip [rclone-current-windows-amd64.zip](https://downloads.rclone.org/rclone-current-windows-amd64.zip).
2. Open a powershell window
3. Run `rclone.exe config create metabolomics drive root_folder_id 0B-ZDcHbPi-aqZzE5V3hOZFc0dms`
4. Find your RClone configuration file location by running `rclone config file`
5. Transfer the RClone configuration file to `~/.config/rclone/rclone.conf` on cori

### Checkout Metatlas code from git repository

On Cori run:
```
cd
git clone https://github.com/biorack/metatlas.git
cd ~/metatlas
git checkout oo_mads3
```

### Make directory to store work in progress

On Cori run
```
mkdir -p ~/metabolomics_data
```

## Per-project workflow

### Get updates to the metatlas code

On Cori run:
```
cd ~/metatlas
git pull
```

### Perform RT correction

#### Parameters
The `experiment_name` parameter can be found in the [Sample Tracking and QC Checkpoints - Northen Lab](https://docs.google.com/spreadsheets/d/126t1OeXQnCCgP6e-6Pac_Ku_A1R7MQLm_tl_Dkqsv_w/edit#gid=1548851545) Google Sheet. The experiment names can be found on the 'New Extraction' sheet in either column 'M' or 'N' depending on the type of chromotography that was performed. This value will be something like `20210723_JGI-AK_DB-TM_506963_LemCreek_final_QE-HF_HILICZ_USHXG01494`.

The `analysis_number` parameter is an integer that you'll need to increment if you redo an analysis. It should be set to 0 initially.

The `project_directory` is where you want to store the analysis while working on it. You should use `~/metabolomics_data`.

#### Execution

On Cori run (where you substitute the 3 parameters described above):
```
cd ~/metatlas/papermill
./launch_rt_prediction.sh experiment_name analysis_number project_directory
```

This will submit a slurm job. You will receive an email when the job starts executing and when it has completed.

#### Outputs

Once the job has completed, you should check the files generated to make sure the RT correction models look acceptable. You can find the output files at `~/metabolomics_data/EXPERIMENT_NAME/${USER}ANALYSIS_NUMBER/`. One easy way to view these files is to open them from the [jupyter](https://jupyter.nersc.gov/) file browser.


### Perform ISTDsEtc Analysis

1. Launch [jupyter.nersc.gov](https://jupyter.nersc.gov/) in your web browser and start a 'Shared CPU Node'
2. Open `~/metabolomics_data/EXPERIMENT_NAME/${USER}ANALYSIS_NUMBER/PROJECT_ID_ISTDsEtc_POS.ipynb` within JupyterLab (you no longer need to use the Classic Notebook interface). If you are prompted to select a kernel, select 'Metatlas Targeted'.
3. The first code cell of the notebook contains descriptions of the parameters and their default values. The second code cell of the notebook contain parameter values that were auto-populated from the RT correction slurm job. These values in the second code block will override the default values from the first code block. The third code block validates your parameter values and also validates that your environment is correctly configured. Execute the first 3 code cells and see if there are any errors. If you get an error message (usually error messages will be in red), you will need to correct the issue so that the cell executes without giving an error before moving on. The error messages commonly see at this point in the workflow generally include some description of what action is needed to correct the problem.
4. Execute the code blocks 4 and 5 to read in data and bring up the Annotation GUI.
5. For each of the compound-adduct pairs in your atlas, set the RT min and RT max boundaries to just contain the EIC peak that corresponds to the compound you are currently evaluating. For each compound-adduct pair, you must either select one of the MSMS-quality descriptors (upper set of radio buttons) or use the bottom set of radio buttons to mark the compound-adduct pair for removal. Failure to set either MSMS-quality descriptors or the remove state for each compound-adduct pair will result in the subsequent step throwing an error.
6. Execute the 6th code block to generate output files and upload them to Google Drive.
7. Review your output files, which will be under [this Google Drive folder](https://drive.google.com/drive/folders/19Ofs5AHB3O8-NYApJUwj4YvH8TbKCGJW?usp=sharing).
8. Repeat steps 1-7 for the the corresponding NEG mode notebook.
9. Move your output folder on Google Drive into the location indicated in column 'M' of the 'New Extraction' sheet in [Sample Tracking and QC Checkpoints - Northen Lab](https://docs.google.com/spreadsheets/d/126t1OeXQnCCgP6e-6Pac_Ku_A1R7MQLm_tl_Dkqsv_w/edit#gid=1548851545) Google Sheet.
10. Email Katherine a link to the output folder so that she can review your outputs.


### Perform FinalEMA-HILIC Analysis

1. Follow the same steps as the ISTDsEtc analysis except use the notebook name `PROJECT_ID_FinalEMA-HILIC_POS.ipynb`.
2. Open the `POS_PROJECT_ID_Final_Identifications.xlsx` file in the output directory on Google Drive.
3. Make sure everything looks okay in the spreadsheet.
4. If there are any compound-adduct pairs that need to be removed at this point (because they are duplicated or you can now determine a similar compound was a better match for a given peak), you can place 'REMOVE' in columns B, M, and N. In columns B and N you should also include some description such as 'REMOVE - duplicate' or 'REMOVE - other isomer prefered (tryptophan matches MSMS reference)' or 'REMOVE - other isomer prefered (tryptophan matches reference RT)'.
5. If you are able to resolve some overlapping identifications at this point, then update the value in column B for the prefered match to no longer include the name of the molecule that is no longer considered a possible match.
6. Repeat steps 1-5 for the corresponding NEG mode notebook.
7. Move your output folder on Google Drive into the location indicated in column 'M' of the 'New Extraction' sheet in [Sample Tracking and QC Checkpoints - Northen Lab](https://docs.google.com/spreadsheets/d/126t1OeXQnCCgP6e-6Pac_Ku_A1R7MQLm_tl_Dkqsv_w/edit#gid=1548851545) Google Sheet.
8. Email Katherine a link to the output folder so that she can review your outputs.

## Using the Annotation GUI

### Key Bindings

`l` or right arrow - next compound-adduct pair

`h` or left arrow - previous compound-adduct pair

`k` or up arrow - next MSMS reference for this compound-adduct pair

`j` or down arrow - previous MSMS reference for this compound-adduct pair

`x` - sets the remove radio button

`m` - steps through the similar compound-adduct pairs and matches the RT bounds to those of the similar compound-adduct pair

`z` - steps though zoom levels of 1x, 5x, 25x on the MSMS mirror plot by 5

`s` - toggles on/off the red and blue shading of EIC plot that show RT ranges for similar compounds
