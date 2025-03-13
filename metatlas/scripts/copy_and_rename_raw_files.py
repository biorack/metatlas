import os
import pandas as pd
import shutil
from pathlib import Path
import argparse

from metatlas.tools.validate_filenames import validate_file_name
from metatlas.tools.logging import activate_module_logging

class ValidationError(Exception):
    pass

def check_file_name(file_name):
    valid_name = validate_file_name(file_name, minimal=True, print_logger=False)
    if valid_name is not True:
        raise ValidationError(f"New filename '{file_name}' failed required validation check.")

def check_directory_name(directory_name):
    base_name = directory_name.name
    num_fields = len(base_name.split("_"))
    if num_fields < 9:
        raise ValidationError(f"New project name '{directory_name}' failed required validation check.")

def check_file_name_uniqueness(table):
    if len(table.iloc[:, 1].unique()) != len(table):
        raise ValidationError(f"New filenames provided in input naming table are not unique.")

def default_logger():
    return activate_module_logging(
            __name__,
            console_level="INFO",
            console_format="{color}{levelname:8}{reset} {message}",
            file_level="INFO",
            filename="/global/cfs/cdirs/m2650/copy_rename_logs/copy_and_rename_raw_files.log")

def copy_and_rename(table_file, current_path, new_path, checks, logger=None):

    logger = logger if logger is not None else default_logger()

    if not os.path.exists(table_file):
        raise ValidationError(f"Table file not found: {table_file}")
    else:
        rename_table = pd.read_csv(table_file, header=None, index_col=False)

    if not os.path.exists(new_path):
        if checks is True:
            check_directory_name(Path(new_path))
        try:
            os.makedirs(new_path)
            shutil.chown(new_path, group='metatlas')
        except Exception as e:
            logger.error(f"Failed to create new directory: {new_path}")
            logger.error(e)

    if checks is True:
        check_file_name_uniqueness(rename_table)

    for index, row in rename_table.iterrows():
        current_name, new_name = row.iloc[0], row.iloc[1]
        if checks is True:
            check_file_name(Path(new_name))
        current_file = os.path.join(current_path, current_name)
        new_file = os.path.join(new_path, new_name)
        if os.path.exists(current_file):
            shutil.copy2(current_file, new_file)
            shutil.chown(new_file, group='metatlas')
            logger.info(f"Copied: {current_file} -> {new_file}")
        else:
            logger.error(f"File not found: {current_file}")

# def main():
#     parser = argparse.ArgumentParser(description="Rename raw files based on a CSV table and copy to new location.")
#     parser.add_argument("table_file", type=str, help="CSV file with two columns: current_name, new_name")
#     parser.add_argument("current_path", type=str, help="Path to the current files")
#     parser.add_argument("new_path", type=str, help="Path to the new files")
#     parser.add_argument("checks", type=bool, help="True for checks on new names, else False")

#     args = parser.parse_args()

#     logger = activate_module_logging(
#             __name__,
#             console_level="INFO",
#             console_format="{color}{levelname:8}{reset} {message}",
#             file_level="INFO",
#             filename="/global/cfs/cdirs/m2650/copy_rename_logs/copy_and_rename_raw_files.log")

#     copy_and_rename(args.table_file, args.current_path, args.new_path, checks=args.checks, logger=logger)

# if __name__ == "__main__": 
#     main()