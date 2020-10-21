#!/bin/bash

#as username pasteur
#give genome group access to metatlas files as user pasteur 
cd /project/projectdirs/metatlas
setfacl -Rdm g:genome:r .
#from metatlas directory

#NOTE THAT THIS IS ALSO RUN AS BPB
