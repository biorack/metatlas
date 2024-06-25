# Creating atlas from CSV table and depositing it to MySQL database

### Gather atlas CSV

The metatlas compound atlas you want to create and deposit should be a table in 
comma-separated value (CSV) format. It may contain any number of fields, but
must contain the following required columns:

`rt_max` - float

`rt_min` - float

`rt_peak` - float

`mz` - float

Additional suggested columns are:

`inchi_key` - XXXXXXXXXXXXXX-XXXXXXXXXX-X

`name` - free text with only ASCII characters (check for special characters)

`polarity` - negative or positive

Once the compound table is in the correct format, move it to a directory on perlmutter
at NERSC where you have read/write permissions.

### Run the Add_Atlas_to_DB.ipynb notebook

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

### Set parameters of Add_Atlas_to_DB.ipynb notebook

#### Cell 1

The first cell contains several parameters that should be set by the user.

`csv_atlas_file_name` - direct path (from root) to the atlas CSV file

`atlas_name` - suggested: `<C18/HILIC/LIPID>_<date>_<TPL>_<Workflow>_<Lab/Unlab>_<NEG/POS>`

`polarity` - negative or positive

`mz_tolerance` - integer (default is 5)

`sort_atlas` - True or False (default is True)

`run_retrieval_check` - True or False (default is False because it taks a while)

Remaining parameters can be set based on specific needs and useage is described.

#### Cell 2

No user-defined settings are required for this cell. Ensure that the metatlas repo commit 
you expect (likely, most recent) is printed in the info log below the cell.

#### Cells 3 and 4

The default behavior of the notebook is to sort your atlas CSV first by retention time (RT)
and second by mass/charge (MZ), then write the sorted CSV to a file at the same location as 
the original csv (with the suffix "sorted" appended before the ".csv"). Sorting helps the
analyst who will be using the template notebook for metabolomics workflows.

Cell 3 defines a function to perform the sorting; Cell 4 reads, sorts, and writes the
atlas CSV, then defines whether the atlas to be deposited will be sorted or unsorted (depending
on the `sort_atlas` parameter in Cell 1).

Ensure the logger info printed under Cell 4 matches your expectations about which atlas will be
deposited to the database.

#### Cell 5

This is the main call of the notebook and will:

1. Grab the designated atlas CSV
2. Run checks on fields/formatting
3. Transform it to a MySQL table
4. Deposit it to the local database at NERSC
5. Assign it the provided name and an auto-generated unique ID
6. Transform the deposited atlas to a pandas dataframe for optional inspection

Cell 5 is a timed block, and should take on the order of single-digit minutes to complete. Once
completed, it will print the atlas name and unique IDs for use in metabolomics analysis with
metatlas workflows.

#### Cell 6

Cell 6 runs a check on whether the number of compounds in the input atlas CSV matches the
number of compounds in the deposited atlas. This is not run by default because retrieval of
the atlas from MySQL by unique ID is time consuming (>15 minutes). Set the `run_retrieval_check`
parameter to `True` in Cell 1 to run this check.

### Document the atlas name and unique ID

To use the newly deposited atlas for metatlas workflows, you will need to update the config file(s)
with the atlas name (user-defined in Cell 1) and the unique ID (auto-generated in Cell 5). Both of
these strings are displayed in the output of the notebook. 
