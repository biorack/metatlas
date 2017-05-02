These notes are a rambling work in progress.  Both how to run container, run at NERSC, and learning Docker.

from a denovo login node
ssh denovo
module load shifter

to create the jupyter session:
shifter --image=scanon/jup:test /chos/global/common/shared/das/start-jupyter-shifter.sh

from terminal:
ssh -L 8888:<IP>:8888 genepool.nersc.gov
leave this terminal session open while jupyter session is open

The script will tell you something like this:
Copy/paste this URL into your browser when you connect for the first time,
to login with a token:
http://localhost:8889/?token=65b21ad135....

the note that tells you how to open the tunnel can be wrong port.  Change tunnel port number to match url port number.

TODO: test tunnel with windows user. If I can't get it to work, send them to office hours Thursday.

also the tunnel command written only works if your local username is the same as nersc it will be different if you have a config file or none of the above
safest to use is username@genepool.nersc.gov
just keep that in mind when showing people how to use it

Docker stacks
https://github.com/jupyter/docker-stacks

TODO: configure system wide jupyter config with CSS:
http://jupyter.readthedocs.io/en/latest/projects/jupyter-directories.html

On my laptop, build a container:
cd to where the Dockerfile is
docker build -t metatlas Dockerfile
docker images 
docker tag metatlas benbowen/metatlas:latest
docker login
docker push benbowen/metatlas:latest


ssh denovo
module load shifter
to create the jupyter session:
shifterimg pull benbowen/metatlas:latest
started at 1:40 PM.
shifterimg images | grep metatlas

```
Ctrl-C out 
```

#shifter --image=benbowen/metatlas:latest --volume=/global/common/shared/das/:/das /das/start-jupyter-shifter.sh
#shifter --image=benbowen/metatlas:latest --volume=/global/project/projectdirs/metatlas/:/metatlas /das/start-jupyter-shifter.sh
shifter --image=benbowen/metatlas:latest ./start-jupyter-shifter.sh

### start-jupyter-shifter.sh

```
#!/bin/sh

# Get the HOME directory after startup
IP=$(ip addr show dev ib0|grep 128.55|awk '{print $2}'|sed 's|/.*||')
echo "*****************************************************"
echo "From your laptop do: "
echo "ssh -L 8888:$IP:8888 genepool.nersc.gov"
echo "*****************************************************"
echo ""

HOME=$PWD 
jupyter-notebook --ip='*' --no-browser $@
```

# for interactive bash experience with shifter
shifter --image=benbowen/metatlas:latest /bin/bash

### notebook running scripts are here:
* /usr/local/bin/start.sh
* /usr/local/bin/start-notebook.sh
* /usr/local/bin/start-singleuser.sh
* /home/$NB_USER/.jupyter/jupyter_notebook_config.py
* chown -R $NB_USER:users /home/$NB_USER/.jupyter

shifter doesn't allow you to volume-mount /global/common. You'll have to copy the start-jupyter-shifter.sh script somewhere it will allow mounting from (e.g. /global/projectb) and mount it from there.

# Shifter at nersc

You are still your nersc username inside the container.  Therefore, you are using the config files in the $USERNAME are used by the notebook by default:

```
$USERNAME/.ipython
README  extensions  kernels  nbextensions  profile_default
$USERNAME/.jupyter
custom  migrated  nbconfig
```

We want everyone to have the same experience in the container and for these not to be used or to be optionally added after invoking a default experience.




from terminal:
ssh -L 8888:<IP>:8888 genepool.nersc.gov
leave this terminal session open while jupyter session is open

TODO: Add my R packages to the R conda packages

ssh -L 8888:<IP>:8888 genepool.nersc.gov
128.55.161.43

http://localhost:8888/?token=2d23167b7594....

docker run -p <port>:<port> <imagename>

Docker 'run' command to start an interactive BaSH session
docker run -it <image> /bin/bash


This file needs to be deleted sometimes:
  623  ls -lah ~/Library/Containers/com.docker.docker/Data/com.docker.driver.amd64-linux/Docker.qcow2 
  624  rm /Users/bpb/Library/Containers/com.docker.docker/Data/com.docker.driver.amd64-linux/Docker.qcow2