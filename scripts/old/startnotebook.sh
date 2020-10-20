#!/bin/sh
module load python/2.7-anaconda
# Get the HOME directory after startup
IP=$(ip addr show dev ib0|grep 128.55|awk '{print $2}'|sed 's|/.*||')
echo "*****************************************************"
echo "From your laptop do: "
echo "ssh -L 8888:$IP:8888 cori.nersc.gov"
echo "*****************************************************"
echo ""

HOME=$PWD 
jupyter-notebook --ip='*' --no-browser $@
