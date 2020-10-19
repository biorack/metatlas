workflow metatlas_wf {
    Int threads
    File config
    File tree_lookup
	Array[Pair[Int,File]] raw_to_mzml
	Array[Pair[Int,File]] mzml_to_hdf5
	Array[Pair[Int,File]] mzml_to_pactolus
	Array[Pair[Int,File]] mzml_to_spectralhits

    #
    # convert raw_to_mzml
    #
	scatter (path in raw_to_mzml){
      call raw_mzml {
           input: path=path
      }
    }

    #
    # convert mzml_to_hdf5
    #
    scatter (path in mzml_to_hdf5) {
      call mzml_hdf5 {
          input: path=path
      }
    }

#    #
#    # convert mzml_to_pactolus
#    #
#    scatter (paths in files_for_pactolus_alias.file_pairs) {
#      call mzml_pactolus {
#          input: path=path,
#                 config = config,
#                 tree_lookup = tree_lookup,
#                 threads = threads
#          }
#        }

    #
    # convert mzml_to_spectralhits
    #
#    scatter (paths in mzml_to_spectralhits) {
#      call mzml_spectralhits {
#          input: path=path,
#				 config=config
#          }
#	}
}


### ------------------------------------------ ###

task raw_mzml {

	Pair[Int,File] path
    Int limskey = path.left
    File input_path = path.right

    command <<<
        SECONDS=0
        required_mzml_size=1000000 # 1Mb
	    dollar='$'


        echo "### creating mzml"
        shifter --volume=$(pwd):/mywineprefix --image=biocontainers/pwiz:phenomenal-v3.0.18205_cv1.2.54 mywine msconvert --32 --mzML ${input_path} 
        mzml_file=$(ls *.mzML)

        if [[ ! -s $mzml_file ]]; then
            >&2 echo "Warning: no mzml file was created."
            exit 0
        fi  

        # run a little file size test to make sure the newly created mzml file is at least 1Mb.
        my_mzml=$(ls *.mzML)
        size=$(ls -l $my_mzml | awk '{print $5}')
        warning=
        if [[ $size -lt $required_mzml_size ]]; then
          warning="warning: file $my_mzml is less than $required_mzml_size."
        fi

        # write output path to file so we can copy the resulting 
        # file from this task to the proper destination after jaws is completed.

		# Write some metadata about the conversion process to a file for parsing downstream.
		# Info from this table will be used to update LIMS.

		cat <<EOF > meta.json
		{
		  limskey: ${limskey}
		  filesize: $size
		  warning: $warning
		}
		EOF

        # get time this task took
        hrs=$(( SECONDS/3600 )); mins=$(( (SECONDS-hrs*3600)/60)); secs=$(( SECONDS-hrs*3600-mins*60 ))
        printf 'Time spent: %02d:%02d:%02d\n' $hrs $mins $secs
    >>> 

    runtime {
        time: "01:00:00"
        mem: "5G"
        poolname: "metatlas"
        node: 1
        nwpn: 8
    }

    output {
      Array[File]? mzml_file = glob("*.mzML")
    }
}

task mzml_hdf5 {
	Pair[Int,File] path
    Int limskey = path.left
    File input_path = path.right
    String dollar='$'

    command <<<
        SECONDS=0
        bname=$(basename ${input_path})
        filename="${dollar}{bname%.*}"
        required_h5_size=1000000  # 1Mb

        shifter --image=jfroula/jaws-pactolus-spectral:1.2.0 mzml_loader_jaws.py \
        -i $input_path \
        -o $filename.h5

        if [[ -s "$filename.h5" ]]; then
            # validate that h5 file is greater than 5M
            filesize=$(stat --printf="%s" $filename.h5)
            if [[ $filesize -lt $required_h5_size ]]; then
              echo "File size for $filename.h5 is less than the minimum size: $required_h5_size."
            fi
        else
            echo "Warning: h5 file missing or empty."
        fi  


        # run a little file size test to make sure the newly created mzml file is at least 1Mb.
        my_h5=$(ls *.h5)
        size=$(ls -l $my_h5 | awk '{print $5}')
        warning=
        if [[ $size -lt $required_h5_size ]]; then
           warning="warning: file $my_h5 is less than $required_h5_size."
        fi

        # write output path to file so we can copy the resulting 
        # file from this task to the proper destination after jaws is completed.
        #echo $key $output_path $size $warning> outpath

		# Write some metadata about the conversion process to a file for parsing downstream.
        # Info from this table will be used to update LIMS.
        cat <<EOF > meta.json
        {
          key: ${limskey}
          filesize: $size
          warning: $warning
        }
		EOF

        # get time this task took
        hrs=$(( SECONDS/3600 )); mins=$(( (SECONDS-hrs*3600)/60)); secs=$(( SECONDS-hrs*3600-mins*60 ))
        printf 'Time spent: %02d:%02d:%02d\n' $hrs $mins $secs
    >>> 

    runtime {
        time: "01:00:00"
        mem: "5G"
        poolname: "metatlas"
        node: 1
        nwpn: 8
    }

    output {
      Array[File]? h5 = glob("*.h5")
    }
}

task mzml_pactolus {
	Pair[Int,File] path
    Int limskey = path.left
    File input_path = path.right
    File config
    File tree_lookup
    Int threads

    command <<<
        SECONDS=0

        echo -e "\n\n### creating pactolus"
        shifter --image=jfroula/jaws-pactolus-spectral:1.2.0 python /root/pactolus/pactolus/score_mzmlfile.py \
        --infile ${input_path} \
        --ms2_tolerance 0.0100 \
        --ms1_tolerance 0.0100 \
        --ms1_pos_neutralizations 1.007276 18.033823 22.989218 \
        --ms2_pos_neutralizations -1.00727646677 -2.01510150677 0.00054857990946 \
        --ms1_neg_neutralizations -1.007276 59.013851 \
        --ms2_neg_neutralizations 1.00727646677 2.01510150677 -0.00054857990946 \
        --tree_file ${tree_lookup} \
        --num_cores ${threads}

        if [[ ! "$lcms_filename.pactolus.gz" ]]; then
            >&2 echo "Warning: no pactolus file created."
            exit 0
        fi

        # write output path to file so we can copy the resulting 
        # file from this task to the proper destination after jaws is completed.
        #echo $key $output_path > outpath

		# Write some metadata about the conversion process to a file for parsing downstream.
        # Info from this table will be used to update LIMS.
        cat <<EOF > meta.json
        {
          key: ${limskey}
        }
		EOF

        # get time this task took
        hrs=$(( SECONDS/3600 )); mins=$(( (SECONDS-hrs*3600)/60)); secs=$(( SECONDS-hrs*3600-mins*60 ))
        printf 'Time spent: %02d:%02d:%02d\n' $hrs $mins $secs
    >>>

    runtime {
        time: "01:00:00"
        mem: "5G"
        poolname: "metatlas"
        node: 1
        nwpn: 8
    }

    output {
      Array[File]? pactolus = glob("*.pactolus.gz")
    }
}

task mzml_spectralhits {
	Pair[Int,File] path
    Int limskey = path.left
    File input_path = path.right
    File config

    command <<<
        SECONDS=0

        # creates spectral-hits.tab.gz from spectral_hits_cmd
        shifter --image=jfroula/jaws-pactolus-spectral:1.2.0 \
        /usr/local/bin/spectral_hits_jaws.py -f -m ${input_path} -l "myfile_spectral-hits.tab.gz" -c ${config} 2> myerr.log
        if [[ -s "myerr.log" ]]; then
          >&2 cat myerr.log
        fi

        if [[ ! "myfile_spectral-hits.tab.gz" ]]; then
            >&2 echo "Warning: no spectralhits file created."
            exit 0
        fi

        # write output path to file so we can copy the resulting 
        # file from this task to the proper destination after jaws is completed.
        #echo $key $output_path > outpath

		# Write some metadata about the conversion process to a file for parsing downstream.
        # Info from this table will be used to update LIMS.
        cat <<EOF > meta.json
        {
          key: ${limskey}
          filesize: $size
          warning: $warning
        }
		EOF

        # get time this task took
        hrs=$(( SECONDS/3600 )); mins=$(( (SECONDS-hrs*3600)/60)); secs=$(( SECONDS-hrs*3600-mins*60 ))
        printf 'Time spent: %02d:%02d:%02d\n' $hrs $mins $secs
    >>> 

    runtime {
        time: "01:00:00"
        mem: "5G"
        poolname: "metatlas"
        node: 1
        nwpn: 8
    }

    output {
      Array[File]? spectralhits =  glob("*_spectral-hits.tab.gz")
    }
}

