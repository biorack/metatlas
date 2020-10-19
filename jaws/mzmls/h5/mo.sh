for i in `cat ../list`; do 
	f=${i%.mzML}; 
	shifter --image=jfroula/jaws-pactolus-spectral:1.2.0 mzml_loader_jaws.py -i $i -o $f.h5; 
done
