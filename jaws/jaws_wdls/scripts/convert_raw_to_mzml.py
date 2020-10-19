#!/usr/bin/env python
import os,sys
import subprocess

raw_dir_tmp = sys.argv[1]
bname=os.path.basename(raw_dir_tmp)

# go into each subdirectory and check for .raw files.  Process these files
# on the spot (i.e. write the mzML files to this subdir).
rootDir = os.path.abspath('.')
for dirName, subdirList, fileList in os.walk(raw_dir_tmp):
    for fname in fileList:
          if fname.endswith(".raw"):
              os.chdir(dirName)
              print("processing raw file: %s" % fname)
              command = "shifter --volume=%s:/mywineprefix --image=biocontainers/pwiz:phenomenal-v3.0.18205_cv1.2.54 mywine msconvert --32 --mzML %s" % (os.path.abspath("."),fname)
              print(command)

              process=subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)
              output = process.stdout
              stderror = process.stderr
              thereturncode = process.returncode
              print(output)

              if thereturncode:
                  exit("There was an error running msconvert: \n%s" % (stderror))

          os.chdir(rootDir)

