#!/bin/bash
#set -e	# to exit program if one command fails
#set -x

INPUT_DIR=
OUTPUT_DIR=
USER="bpb"

usage()
{
cat<<EOF
  Usage: $0 
    <-i <inputdir> (required)> 
    <-o <outdir> (optional)>

    default outdir is same path as original path of the input file as recorded in LIMS
EOF
exit 1
}


# check for all required args
ARGS=$*
if [[ $ARGS =~ "-i" ]]; then 
    echo; 
else
    echo -e "Missing some required arguments\n"
    usage
fi


while getopts 'i:o:' OPTION
do 
  case $OPTION in 
  i)    INPUT_DIR="$OPTARG"
        ;;
  o)    OUTPUT_DIR="$OPTARG"
        ;;
  ?)    echo usage
        exit 1
        ;;
  esac
done

if [ ! -e $INPUT_DIR ]; then
  echo
  echo "Error: No input file was found."
  echo
  exit
fi

# verify user defined output dir exists
if [[ $OUTPUT_DIR ]]; then
	if [[ ! -d $OUTPUT_DIR ]]; then
		mkdir -p $OUTPUT_DIR || echo Failed to create $OUTPUT_DIR
	fi
fi

function copy_files {
	NAME=$1

    #
    # find all mzML files and copy them from the jaws repo to final destination
    #
    SOURCE_FILES=($(find $INPUT_DIR -name "$NAME" | grep -v glob-))
    num=${#SOURCE_FILES[@]}
	tmplog=copy.$$.log
    for file in ${SOURCE_FILES[@]}; do
    	if [[ $OUTPUT_DIR ]]; then
    		cp $file $OUTPUT_DIR
    		echo cp $file $OUTPUT_DIR >> $tmplog
    	else
			# grab the destination file path. This path has been stored in meta.json.
    		dname=$(dirname $file)
    		outpath=$(grep outpath $dname/meta.json | awk '{print $2}')
    		cp $file $outpath
    		echo cp $file $outpath >> $tmplog 

			# give USER full permissions on copied file
    		setfacl -m "u:$USER:rwx" $outpath
    
    	fi
    done
    echo "there were $num files"
}

echo -n "## copying mzml files ..."
copy_files '*.mzML'

echo -n "## copying h5 files ..."
copy_files '*.h5'

echo -n "## copying pactolus files ..."
copy_files '*.pactolus.gz'

echo -n "## copying spectral-hits files ..."
copy_files '*_spectral-hits.tab.gz'

