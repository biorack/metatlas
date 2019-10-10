workflow raw_to_mzml_wf {
    File raw_dir
    String out_dir
    Int threads

    call find_raw_files {
        input: raw_dir=raw_dir
    }

    scatter (file in find_raw_files.raw_files) {
        call raw_to_mzml_conversion {
            input: raw_file=file
        }
    }

    call copy_mzml_files {
        input: mzml_files=raw_to_mzml_conversion.mzml_file,
			   out_dir = out_dir
    }

}

### ------------------------------------------ ###
task find_raw_files {
    File raw_dir

    command {
        find ${raw_dir} -name "*.raw" 
    }

    output {
        Array[File] raw_files = read_lines(stdout())
    }
}

task raw_to_mzml_conversion {
    File raw_file
	String bname = basename(raw_file)
	String mzml_file_name = sub(bname,".raw$",".mzML")

    command {
	    shifter --volume=$(pwd):/mywineprefix \
		  --image=biocontainers/pwiz:phenomenal-v3.0.18205_cv1.2.54 \
		  mywine msconvert --32 --mzML ${raw_file}
    }

    output {
        File mzml_file = mzml_file_name
    }
}

task copy_mzml_files {
    Array[File] mzml_files
	String out_dir

    command {
        cp ${sep=' ' mzml_files} ${out_dir}
    }
}
