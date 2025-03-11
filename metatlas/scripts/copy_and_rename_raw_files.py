import sys
import os
import csv
import shutil
from pathlib import Path
import argparse

from metatlas.tools.validate_filenames import validate_file_name
from metatlas.tools.logging import activate_module_logging

def activate_logging():
    logger = activate_module_logging(
        __name__,
        console_level="INFO",
        console_format="{color}{levelname:8}{reset} {message}",
        file_level="INFO",
        filename=Path("/global/cfs/cdirs/m2650/copy_rename_logs"),
    )
    return logger

def read_in_rename_table(table_file):
    with open(table_file, mode='r') as file:
        reader = csv.reader(file)
        return list(reader)

def check_file_name(file_name, logger):
    valid_name = validate_file_name(file_name, minimal=True)
    if valid_name is not True:
        logger.critical(f"New filename {file_name} failed required validation check. Check that file name conforms.")
        sys.exit(1)

def check_directory_name(directory_name, logger):
    base_name = Path(directory_name).name
    num_fields = len(base_name.split("_"))
    if num_fields < 9:
        logger.critical(f"New directory name {directory_name} failed required validation check. It has {num_fields} fields, but at least 9 are required.")
        sys.exit(1)

def copy_and_rename_raw_files(table, current_path, new_path, logger):

    if not os.path.exists(new_path): # Create new project directory if it doesn't exist
        check_directory_name(new_path, logger)
        os.makedirs(new_path)
        shutil.chown(new_path, group='metatlas')
        logger.info(f"Notice! Had to make new project directory: {new_path}")

    for current_name, new_name in table:
        check_file_name(new_name, logger)
        current_file = os.path.join(current_path, current_name)
        new_file = os.path.join(new_path, new_name)
        if os.path.exists(current_file):
            shutil.copy2(current_file, new_file)
            logger.info(f"Copied: {current_file} -> {new_file}")
        else:
            logger.warning(f"File not found: {current_file}")

def main():
    parser = argparse.ArgumentParser(description="Rename raw files based on a CSV table and copy to new location.")
    parser.add_argument("table_file", help="CSV file with two columns: current_name, new_name")
    parser.add_argument("current_path", help="Path to the current files")
    parser.add_argument("new_path", help="Path to the new files")

    args = parser.parse_args()

    logger = activate_logging()

    rename_table = read_in_rename_table(args.table_file)
    copy_and_rename_raw_files(rename_table, args.current_path, args.new_path, logger)

if __name__ == "__main__":
    main()