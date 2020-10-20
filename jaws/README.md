# Metabolomics Pipeline (Ben Bowen)

## How to Run Using Test Data
```
1) # on cori
   ssh cori.nersc.gov

2) git clone https://gitlab.com/jfroula/jaws-benbowen.git

3) cd jaws-benbowen && module load python # you need python3

4) cusomize your inputs.json file to point to the input & output dirs. (see below)

5) # I'm running /usr/bin/java v11.0.6 
   java -jar /global/cfs/projectdirs/jaws/cromwell/cromwell.jar run raw_to_mzml.wdl -i inputs.json
```

## Understanding the inputs.json
You need to customize the inputs.json with paths to your .raw (or .mzML) files (i.e. replace <full-path-to>). Include paths to your spectral_hits.config and tree_lookup.npy file.
Also, indicate with a true or false whether the input is mzML on this line "raw_to_mzml_wf.mzml_inputs_bool" (default is 'false' if this line is missing).

```
{
  "raw_to_mzml_wf.mzml_dir": "<full-path-to>/raw_data_test",
  "raw_to_mzml_wf.dest_dir": "<full-path-to>/raw_data_test_out",
  "raw_to_mzml_wf.config": "./spectral_hits.config",
  "raw_to_mzml_wf.tree_lookup": "./tree_lookup.npy",
  "raw_to_mzml_wf.threads": 2,
  "raw_to_mzml_wf.mzml_inputs_bool": false
}
```

## Docker environment
There are only two docker images required.  They are hardcoded in the WDL and they may need updating from time to tome. Furthermore, 
they are maintained on docker hub. For the Dockerfile, see the docker hub link (except for pwiz).

image for msconvert (no Dockerfile on hub.docker)
[biocontainers/pwiz:phenomenal-v3.0.18205_cv1.2.54](https://hub.docker.com/r/biocontainers/pwiz)

image for everything else
[jfroula/jaws-pactolus-spectral](https://hub.docker.com/repository/docker/jfroula/jaws-pactolus-spectral)


## Notes
Python scripts and libraries from these repositories were used in the new docker images
```
git clone https://github.com/biorack/metatlas.git
git clone https://github.com/biorack/pactolus.git
```

Ben's original document defining environment for this pipeline
Master environment doc  
`https://github.com/biorack/nersc_python_install/blob/master/install_python.sh`

