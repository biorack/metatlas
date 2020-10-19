#!/bin/bash
id=$1

if [[ ! $id ]]; then 
	echo "Usage: $0 <cromwell id only>"
	exit 1
fi

rawfiles=cromwell-executions/metatlas_wf/$id/call-files_for_mzml_alias/execution/raw_files.tsv
rawfiles=$(realpath $rawfiles)

dir=cromwell-executions/metatlas_wf/$id/call-raw_to_mzml
dir=$(realpath $dir)

find $dir -name meta.json -exec grep 'warning' {} \; | sort -u | awk '{print $4}' > mzml_too_small.files

for i in `cat mzml_too_small.files`; do 
	grep $i $rawfiles | awk '{print $1,$3}'
done  > mzml_too_small_fullpath.files 
