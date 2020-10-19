#!/bin/bash
set -e 
id=$1
if [[ ! $id ]]; then
	echo "Usage: $0 <cromwell id>"
	exit
fi

echo mzML
find cromwell-executions/metatlas_wf/${id}/call-raw_to_mzml -name "*.mzML" | grep -v glob | wc -l

echo hdf5
find cromwell-executions/metatlas_wf/${id}/call-mzml_to_hdf5/ -name "*.h5" | grep -v glob | wc -l

echo pactolus
find cromwell-executions/metatlas_wf/${id}/call-mzml_to_pactolus -name "*.pactolus.gz" | grep -v glob | wc -l

echo spectralhits
find cromwell-executions/metatlas_wf/${id}/call-mzml_to_spectralhits -name "*_spectral-hits.tab.gz" | grep -v glob | wc -l
