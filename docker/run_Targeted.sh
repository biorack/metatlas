#!/bin/bash

if [[ $# == 1 ]]; then
	IMAGE=$1
elif [[ $# == 0 ]]; then
	IMAGE=metatlas
else
	echo "Usage: $0 [image_name]"
	exit 1
fi

docker run -it --rm -v ~/metatlas-dev:/src -v $(pwd)/out:/out $IMAGE papermill -p experiment 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583 -p metatlas_repo_path /src -p project_directory /out -p max_cpus 2 /src/notebooks/reference/Targeted.ipynb /out/Targeted.ipynb
