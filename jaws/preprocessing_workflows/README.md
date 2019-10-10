# Metabolomics Pipeline (Ben Bowen)

## Docker environment
These are links to the docker environments (in hub.docker) that are required for this workflow.

* [mzml -> h5 conversion](https://cloud.docker.com/repository/docker/jfroula/mzml-convert)
* [spectral hits](https://cloud.docker.com/repository/docker/jfroula/spectral-hits)
* [pactolus](https://cloud.docker.com/repository/docker/jfroula/pactolus)
* [mzmine](https://cloud.docker.com/repository/docker/jfroula/mzmine)

## Notes
Python scripts and libraries from these repositories were used in the new docker images
```
git clone https://github.com/biorack/metatlas.git
git clone https://github.com/biorack/pactolus.git
```

The mzml_to_h5 conversion step uses this script  
`<this_repository>/metatlas/mzml_loader_jaws.py`

Ben's original document defining environment for this pipeline
Master environment doc  
`https://github.com/biorack/nersc_python_install/blob/master/install_python.sh`

his python with correct env  
`/global/common/software/m2650/python-cori/bin/python`


