import sys
import os
import csv
import shutil

if len(sys.argv) != 4:
    print("Usage: python script.py <table.csv> <current_path> <new_path>")
    sys.exit(1)

def read_in_rename_table(table_file):
    with open(table_file, mode='r') as file:
        reader = csv.reader(file)
        return list(reader)

def copy_and_rename_raw_files(table, current_path, new_path):
    if not os.path.exists(new_path): # Create new project directory if it doesn't exist
        os.makedirs(new_path)
        shutil.chown(new_path, group='metatlas')
        print(f"Notice! Had to make new project directory: {new_path}")
    for current_name, new_name in table:
        current_file = os.path.join(current_path, current_name)
        new_file = os.path.join(new_path, new_name)
        if os.path.exists(current_file):
            shutil.copy2(current_file, new_file)
            shutil.chown(new_file, group='metatlas')
            print(f"Copied: {current_file} -> {new_file}")
        else:
            print(f"File not found: {current_file}")

rename_table_file = sys.argv[1]
current_files_basepath = sys.argv[2]
new_files_basepath = sys.argv[3]

rename_table = read_in_rename_table(rename_table_file)
copy_and_rename_raw_files(rename_table, current_files_basepath, new_files_basepath)
