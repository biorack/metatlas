workflow metatlas_wf {
    File mzml_dir
    File config
    File tree_lookup
    Int threads

    call find_mzml_files {
        input: mzml_dir=mzml_dir
    }

    scatter (file in find_mzml_files.mzml_files) {
        call mzml_file_conversion {
            input: mzml_file=file
        }
    }

    call spectral_hits_create_cmds {
        input: mzml_dir=mzml_dir,
            config=config
    }

    Array[String] cmds = read_lines(spectral_hits_create_cmds.spectral_cmds)

    scatter (cmd in cmds) {
        call spectral_hits_run_cmds {
            input: spectral_single_cmd = cmd,
                config = config
        }
    }

    scatter (file in find_mzml_files.mzml_files) {
        call pactolus {
            input: mzml_file = file,
                tree_lookup = tree_lookup,
                threads = threads
        }
    }
}

### ------------------------------------------ ###
task find_mzml_files {
    File mzml_dir
    command {
        find ${mzml_dir} -name "*.mzML" > f.tmp
    }
    output {
        Array[File] mzml_files = read_lines("f.tmp")
    }
}

task mzml_file_conversion {
    File mzml_file
    String output_file = basename(mzml_file) 
    String h5_name = sub(output_file,".mzML",".h5")

    command {
        shifter --image=jfroula/jaws-pactolus-spectral:1.0.2 mzml_loader_jaws.py \
                -i ${mzml_file} \
                -o ${h5_name}
    }
    output {
        File h5_file = h5_name
    }
}

task spectral_hits_create_cmds {
    File mzml_dir
    File config

    command {
        shifter --image=jfroula/jaws-pactolus-spectral:1.0.2 spectral_hits_jaws.py -f -n -c ${config} -w spectral_hits_cmds ${mzml_dir}
    }
    output {
        File spectral_cmds = "spectral_hits_cmds"
    }
}

task spectral_hits_run_cmds {
   String spectral_single_cmd
   File config

   command {
       shifter --image=jfroula/jaws-pactolus-spectral:1.0.2 ${spectral_single_cmd} -c ${config}
   }
   output {
        String out = read_string(stdout())
   }
}

task pactolus {
    File mzml_file
    File tree_lookup
    Int threads = 2

    command {
        shifter --image=jfroula/jaws-pactolus-spectral:1.0.2 score_mzmlfile.py \
        --infile ${mzml_file} \
        --ms2_tolerance 0.0100 \
        --ms1_tolerance 0.0100 \
        --ms1_pos_neutralizations 1.007276 18.033823 22.989218 \
        --ms2_pos_neutralizations -1.00727646677 -2.01510150677 0.00054857990946 \
        --ms1_neg_neutralizations -1.007276 59.013851 \
        --ms2_neg_neutralizations 1.00727646677 2.01510150677 -0.00054857990946 \
        --tree_file ${tree_lookup} \
        --num_cores ${threads}
    }
    output {
    }
}

