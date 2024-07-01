# Create metatlas atlas from CSV and deposit into database

## Purpose

This documentation will describe how to convert a CSV containing metabolite
information (for ISTD, QC, EMA analyses, etc.) into a metatlas MySQL database
object that can be referenced by name or unique ID.

## Gather atlas CSV

The metatlas compound atlas you want to create and deposit should be a table in 
comma-separated value (CSV) format. It may contain many different fields, but
must contain the following required columns:

`rt_max` - float

`rt_min` - float

`rt_peak` - float

`mz` - float

Additional suggested columns are:

`label` - free ASCII text. For internal stardards atlases, the label column should include a 
parenthetical indicating 'unlabeled' or the isotopic labeling pattern (e.g. 'U - 13C, 15N')

`inchi_key` - XXXXXXXXXXXXXX-XXXXXXXXXX-X

`polarity` - negative or positive

`adduct` - the adduct of the ionized compound

Once the compound table is in the correct format, move it to a directory on perlmutter
at NERSC where you have read/write permissions.

## Run the Add_Atlas_to_DB.ipynb notebook

The [Add_Atlas_to_DB.ipynb](
https://github.com/biorack/metatlas/blob/main/notebooks/reference/Add_Atlas_to_DB.ipynb) notebook will deposit your CSV into the metatlas MySQL database.
A useful way to run the notebook is to use the [JupyterLab ](https://jupyter.nersc.gov/) server in your web browser
by starting a 'Login Node' session on Perlmutter. You can read more about how to use the 
service at NERSC's [Jupyter Overview page](https://docs.nersc.gov/services/jupyter/).

To do this, ensure you have cloned or pulled the most up-to-date version of the [metatlas github repo](https://github.com/biorack/metatlas/tree/main)
main branch on perlmutter and navigate there within the JupyterLab ecosystem. Open up the notebook at
``` <REPO_DIR>/notebooks/reference/Add_Atlas_to_DB.ipynb ``` and choose the "Metatlas Targeted" kernel if prompted.

You can also run the notebook on perlmutter from an IDE such as VSCode, ensuring you are using an
appropriate kernel, [e.g. Metatlas Targeted](https://github.com/biorack/metatlas/blob/main/docker/shifter.kernel.json).

### Set parameters in notebook cells

#### Parameters

The first cell contains several parameters that should be set by the user.

`csv_atlas_file_name` - direct path (from root) to the atlas CSV file on perlmutter

`atlas_name` - suggested: **C18/HILIC/LIPID**_**date-as-Ymd**_**TPL**_**Workflow**_**Lab/Unlab**_**NEG/POS**
(e.g., `C18_20240624_TPL_EMA_Unlab_NEG`)

`polarity` - negative or positive

`mz_tolerance` - integer (default is 5), other common values are 10 and 20

`sort_atlas` - True or False (default is True, which sorts atlas ascending by RT and MZ 
[and optionally compound labeling, see below])

`istd_atlas` - True or False (default is False, meaning this is not an ISTD atlas and
no isotopically labeled compounds are present in CSV. If True, CSV should have the strings
"(unlabeled)" or "(_labeling pattern_)" e.g. "(U - 1C, 15N)" in the label column, and
the notebook will sort the atlas so labeled compunds appear first)

`run_retrieval_check` - True or False (default is False because it takes a while, this tests
if the deposited atlas has the same number of compounds as the input CSV)

Remaining parameters can be set based on specific needs and usage is described in the notebook.

#### Initialization

The second cell ininitalizes all the variables set above, performs system checks (e.g., on 
the kernel environment), and import various functions from the metatlas universe.

#### Atlas Sort

The third cell performs a convenience sort if `sort_atlas` is True by reading, sorting, and 
writing the atlas to an updated CSV.

Sorting helps the analyst who will be using the template notebook for 
metabolomics workflows and should generally be kept as True. The default behavior 
of the notebook is to (1) sort your atlas CSV first by retention time (RT)
and second by mass/charge (MZ), with the lowest values of each appearing first in the 
final atlas, and (2) write the sorted CSV to a file at the same location as 
the original CSV (with the suffix "sorted" appended before the ".csv").  If you are
inputting a file which is in CSV format (required) but whose name does not end in ".csv",
you will need to ammend this line or update your filename. 

If you are depositing an internal standards atlas with both unlabeled and isotopically labeled
compounds, also ensure that `istd_atlas` is set to True so the sorting function
knows to put the labeled compounds first in the final atlas.

Check that the logger info printed in this cell matches your expectations about which CSV (sorted
or unsorted) will be deposited to the database.

#### Atlas Generation

This is the main call of the notebook and will:

1. Grab the designated atlas CSV
2. Run checks on fields/formatting
3. Transform it to a MySQL table
4. Deposit it to the local database at NERSC
5. Assign it the provided name and an auto-generated unique ID
6. Transform the deposited atlas to a pandas dataframe for optional inspection

This cell is a timed block, and should take on the order of minutes to complete. Once
completed, it will print to standard output the atlas name and unique IDs for use in 
metabolomics analysis configs with metatlas workflows.

#### Deposited Atlas Check

The final cell runs an optional check on whether the number of compounds in the input CSV matches the
number of compounds in the deposited atlas. This is not run by default because retrieval of
the atlas from MySQL by unique ID is time consuming (>15 minutes). Set the `run_retrieval_check`
parameter to `True` in the notebook Parameters cell to run this check.

## Document the atlas name and unique ID

Important note: To use the newly deposited atlas for metatlas workflows you will need to 
update the config file(s) with the atlas name and the unique ID (printed after the Atlas
Generation cell).
