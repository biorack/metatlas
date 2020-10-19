#!/bin/bash
id=$1

set -e

if [[ ! $id ]]; then 
	echo "Usage: $0 <cromwell id>"
	exit 1
fi

rawfiles=cromwell-executions/metatlas_wf/$id/call-files_for_spectralhits_alia/sexecution/raw_files.tsv
dir=cromwell-executions/metatlas_wf/$id/call-mzml_to_spectralhits

function find_err_raws {
	err=$1
	find $dir -name stderr -exec grep -l "$err" {} \; > spec_hits.err

	for i in `cat spec_hits.err`; do
  		dname=$(dirname $i)
		filename=$(grep outpath $dname/meta.json | awk '{print $2}')
 		filename=$(echo ${filename%_spectral-hits.tab.gz}.mzML)
		if [[ -e $filename ]]; then
			echo $filename
		else
			echo Error, file does not exist: $file
			exit 1
		fi
	done > spec_hits.files
	rm spec_hits.err
}

find_err_raws "not convertable to dataframe"

