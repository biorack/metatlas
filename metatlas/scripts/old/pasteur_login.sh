#!/bin/bash

module load globus
myproxy-logon -s nerscca.nersc.gov
gsissh localhost -p 2222 -l pasteur
