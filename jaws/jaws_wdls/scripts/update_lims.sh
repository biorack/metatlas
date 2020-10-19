#!/bin/bash
fullpath=$(realpath $0)
path_to_script=$(dirname $fullpath)
$path_to_script/update_lcmsfiles_in_lims.py --update_file_tables True --update_lcmsrun_names True --update_lcmsrun_items True
$path_to_script/update_lcmsfiles_in_lims.py --update_fileconversion_tasks True
