cd /root/
ls
./spectral_hits_jaws.py -h
./mzml_loader_jaws.py -h
mzml_loader_jaws.py -i 
ls
ls ../
ls /dat
ls
cd /root/
ls
mzml_loader_jaws.py -i /dat/*.mzML -o frog
ls 
spectral_hits_jaws.py -f -n -c spectral_hits.config -w cmds MZML
mkdir MZML
spectral_hits_jaws.py -f -n -c spectral_hits.config -w cmds MZML
cat cmds
cat spectral_hits.config
ls
vi metatlas.wdl
less metatlas.wdl
score_mzmlfile.py         --infile ${mzml_file}         --ms2_tolerance 0.0100         --ms1_tolerance 0.0100         --ms1_pos_neutralizations 1.007276 18.033823 22.989218         --ms2_pos_neutralizations -1.00727646677 -2.01510150677 0.00054857990946         --ms1_neg_neutralizations -1.007276 59.013851         --ms2_neg_neutralizations 1.00727646677 2.01510150677 -0.00054857990946         --tree_file ${tree_lookup}         --num_cores ${threads}
