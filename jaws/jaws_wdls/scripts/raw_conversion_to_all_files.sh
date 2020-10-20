#!/bin/bash

usage()
{
cat<<EOF
  Usage: $0 <options> <desired files to produce, as list separated by space [mzml h5 spec pactolus]>
    <-r raw file (required)>
    <-l tree lookup [tree_lookup.npy] (optional)>
    <-c config [spectral_hits.config] (optional)>
    <-t threads [1] (optional)>

    defaults are shown in square brackets
EOF
exit 1
}

# check for all required args
ARGS=$*
if [[ $ARGS =~ "-r" ]]; then
    echo;
else
    echo -e "Missing some required arguments\n"
    usage
fi


while getopts 'r:l:c:t:' OPTION
do
  case $OPTION in
  c)    CONFIG="$OPTARG"
        ;;
  l)    TREE_LOOKUP="$OPTARG"
        ;;
  r)    RAW_FILE="$OPTARG"
        ;;
  t)    THREADS="$OPTARG"
        ;;
  ?)    echo usage
        exit 1
        ;;
  esac
done

shift $((OPTIND - 1))
files_to_create="$*"
read -r -a files_to_create <<< "$*"

if [[ ${#files_to_create[@]} == 0 ]]; then
	echo "You need to choose one or more files to create. Use one or more of the arguments [mzml h5 spec pactolus] and list them after any command line options."
	exit 1
fi

# validate that user included only acceptable values
acceptable=("mzml" "h5" "spec" "pactolus")

for value in ${files_to_create[@]}; do
	if [[ ! " ${acceptable[@]} " =~ " ${value} " ]]; then
		echo "${value} is not an exceptable file type. Please use only (${acceptable[@]})"
		exit 1;
	fi
done

# defaults
CONFIG=${CONFIG:-"spectral_hits.config"}
TREE_LOOKUP=${TREE_LOOKUP:-"tree_lookup.npy"}
THREADS=${THREADS:-1}

bname=$(basename $RAW_FILE)
bname=${bname%.raw}
if [[ ! -s $RAW_FILE ]]; then
	>&2 echo "Error: missing or empty file: $RAW_FILE"
	exit 1
fi

#
# raw -> mzml
#
if [[ " ${files_to_create[@]} " =~ " mzml " ]]; then
	echo "### creating mzml"
    shifter --volume=$(pwd):/mywineprefix --image=biocontainers/pwiz:phenomenal-v3.0.18205_cv1.2.54 mywine msconvert --32 --mzML $RAW_FILE 
    if [[ ! -s "$bname.mzML" ]]; then
    	>&2 echo "Error: no mzml file was created."
    	exit 1
    fi
fi

#
# mzml -> h5 conversion
#
if [[ " ${files_to_create[@]} " =~ " h5 " ]]; then
	echo -e "\n\n### creating h5"
    shifter --image=jfroula/jaws-pactolus-spectral:1.1.0 mzml_loader_jaws.py \
    -i $bname.mzML \
    -o $bname.h5
    if [[ ! -s "$bname.h5" ]]; then
        >&2 echo "Error: no h5 file was created."
        exit 1
    fi
fi
    
#
# create commands for creating spectral_hits.gz files
#
if [[ " ${files_to_create[@]} " =~ " spec " ]]; then
	echo -e "\n\n### creating spectral hits"
    shifter --image=jfroula/jaws-pactolus-spectral:1.1.0 \
    /usr/local/bin/spectral_hits_jaws.py -f -n -c ${CONFIG} -w spectral_hits_cmds $(pwd)
    if [[ ! -s "spectral_hits_cmds" ]]; then
        >&2 echo "Error: no spectral_hits_cmds was created."
        exit 1
    fi
    
    # creates spectral-hits.tab.gz from spectral_hits_cmd
    spectral_single_cmd=$(cat spectral_hits_cmds)
    shifter --image=jfroula/jaws-pactolus-spectral:1.1.0 $spectral_single_cmd -c ${CONFIG} || echo "warning: no spectral-hits.tab.gz produced"
fi
    
#
# create pactolus.gz files
#
if [[ " ${files_to_create[@]} " =~ " pactolus " ]]; then
	echo -e "\n\n### creating pactolus"
    shifter --image=jfroula/jaws-pactolus-spectral:1.1.0 python /root/pactolus/pactolus/score_mzmlfile.py \
    --infile $bname.mzML \
    --ms2_tolerance 0.0100 \
    --ms1_tolerance 0.0100 \
    --ms1_pos_neutralizations 1.007276 18.033823 22.989218 \
    --ms2_pos_neutralizations -1.00727646677 -2.01510150677 0.00054857990946 \
    --ms1_neg_neutralizations -1.007276 59.013851 \
    --ms2_neg_neutralizations 1.00727646677 2.01510150677 -0.00054857990946 \
    --tree_file ${TREE_LOOKUP} \
    --num_cores $THREADS
    
    if [[ ! "$bname.pactolus.gz" ]]; then
        >&2 echo "Error: no pactolus file created."
        exit 1
    fi
fi    
