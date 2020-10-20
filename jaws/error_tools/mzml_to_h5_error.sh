#!/bin/bash
# rawfiles='../call-files_for_hdf5_alias/execution/raw_files.tsv' 
id=$1

if [[ ! $id ]]; then 
	echo "Usage: $0 <cromwell id>"
	exit 1
fi

rawfiles=cromwell-executions/metatlas_wf/$id/call-files_for_hdf5_alias/execution/raw_files.tsv
rawfiles=$(realpath $rawfiles)

dir=cromwell-executions/metatlas_wf/$id/call-mzml_to_hdf5
dir=$(realpath $dir)

function find_err_raws {
	err=$1
	find $dir -name stderr -exec grep -l "$err" {} \; > $err.err

	for i in `cat $err.err`; do
  		dname=$(dirname $i)
  		filename=$(ls $dname/*.h5)
  		bname=$(basename $filename)
		fullname=$(grep $bname $rawfiles)
		echo "$fullname $err"
	done > h5_${err}.files
	rm $err.err
}

expected=$(wc -l $rawfiles | cut -f1 -d' ')
created=$(find $dir -name "*.h5" |  grep -v 'glob-' | wc -l)
>2& echo there were $created h5 files found and expected $expected

find_err_raws ms2_neg
find_err_raws ParseError


find $dir -name "*.h5" -exec ls -lh {} \; | grep -v glob- > m
for i in `awk '$5 !~ "M" {print $5":"$9}' m`; do 
	size=$(echo $i | cut -f1 -d:)
	bname=$(echo $i | cut -f2 -d:)
	bname=$(basename $bname)
	fullname=$(grep $bname $rawfiles)
	echo "$fullname too_small_[${size}]"
done > h5_too_small_errors
rm m

awk '{print $1}' h5_too_small_errors | sort -n > t
awk '{print $1}' h5_ParseError.files | sort -n > p
comm -23 p t > unique_to_parse_error
comm -13 p t > unique_to_too_small
comm -12 p t > in_common
>2& echo unique_to_parse_error: $(wc -l unique_to_parse_error)
>2& echo unique_to_too_small: $(wc -l unique_to_too_small)
>2& echo in common: $(wc -l in_common)
for i in `cat in_common`; do 
	grep $i h5_too_small_errors h5_*.files | cut -f2 -d: | awk '{print $1,$3,$4}'
done

