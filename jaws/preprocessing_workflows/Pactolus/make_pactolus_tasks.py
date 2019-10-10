from __future__ import print_function
import sys
import os
import fnmatch
import glob as glob
import argparse

def make_sbatch_filename(file):
    return ''.join(os.path.basename(file).split('.')[:-1])+'.sbatch'

def create_job_script(infile=None,
               ms1_tolerance=0.01,
               ms2_tolerance=0.01,
               ms1_pos_neutralizations=[1.007276, 18.033823, 22.989218],
               ms1_neg_neutralizations=[-1.007276, 59.013851],
              ms2_pos_neutralizations=[-1.00727646677,-2.0151015067699998,0.00054857990946],
              ms2_neg_neutralizations=[1.00727646677,2.0151015067699998,-0.00054857990946],
               tree_file='/project/projectdirs/metatlas/projects/clean_pactolus_trees/tree_lookup.npy',
               num_cores=16,
              project_id='m1541',
              time='03:59:00',
              job_name='pactolus',
              outlog = '/project/projectdirs/metatlas/projects/logs/pactolus_logs/out.out',
               errlog='/project/projectdirs/metatlas/projects/logs/pactolus_logs/err.out',
               python='python',
               pactolus_script='score_mzmlfile.py',
              partition='realtime',
              cores=8,
              use_haswell=True,
              license_project=True,
              license_scratch=False):
    """
    creates a slurm srun string to search pactolus and make an output.
    outfile is a csv.
    
    infile is an mzml file
    
    Your job string should look like this:
    srun -A m2650 -p debug -N 1 -t 00:20:00 ...
    -o '/project/projectdirs/metatlas/projects/logs/pactolus_logs/out.out' ...
    -e '/project/projectdirs/metatlas/projects/logs/pactolus_logs/err.out' ...
    python score_pactolus.py --infile 'None' --outfile 'None' ...
    --polarity positive --ms2_tolerance 0.0100 --ms1_tolerance 0.0100 ...
    --ms1_neutralizations 1.007276 18.033823 22.989218 ...
    --ms2_neutralizations -1.00727646677 -2.01510150677 0.00054857990946 ...
    --tree_file '/project/projectdirs/metatlas/projects/clean_pactolus_trees/tree_lookup.npy' --num_cores 64 &
    
    For negative mode, ms2_neutralizations can be the positive ones multiplied by -1
    For negative mode, ms1_neutralizations would likely by: [-1.007276 59.013851 34.969402]
    
    Consider replacing python with '/project/projectdirs/metatlas/anaconda/bin/python'
    
    pactolus_script with path to where script is:
    /global/homes/b/bpb/repos/pactolus/pactolus/sc
    
    """
    hdf_flag = 'export HDF5_USE_FILE_LOCKING=FALSE'

    #add trailing spaces to each row.
    header_str = ['#!/bin/bash -l',
                  '#SBATCH --account=%s'%project_id,
                  '#SBATCH --job-name="%s"'%job_name,
                  '#SBATCH --time=%s'%time,
                  '#SBATCH --nodes=1',
                  '#SBATCH --output="%s"'%outlog,
                  '#SBATCH --error="%s"'%errlog,
                  '#SBATCH --partition=%s'%partition]
    if use_haswell:
        header_str.append('#SBATCH -C haswell')
    if license_project:
        header_str.append('#SBATCH -L project')
    if license_scratch:
        header_str.append('#SBATCH -L scratch')

    # cmd_str = ("%s && %s %s "%(hdf_flag,python,pactolus_script),
    cmd_str = ("%s %s "%(python,pactolus_script),
              '--infile "%s" '%(infile),
               "--ms2_tolerance %.4f --ms1_tolerance %.4f "%(ms2_tolerance,ms1_tolerance),
            "--ms1_pos_neutralizations %s "%(' '.join(map(str,ms1_pos_neutralizations))),
            "--ms2_pos_neutralizations %s "%(' '.join(map(str,ms2_pos_neutralizations))),
            "--ms1_neg_neutralizations %s "%(' '.join(map(str,ms1_neg_neutralizations))),
            "--ms2_neg_neutralizations %s "%(' '.join(map(str,ms2_neg_neutralizations))),
            "--tree_file '%s' "%(tree_file),
            "--num_cores %d\n"%(cores) )
    # header_str = '\n'.join(header_str)
    cmd_str = ''.join(cmd_str)
    # job_str = '%s\n%s'%(header_str,cmd_str)
    
    return cmd_str


def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


def main():
    parser = argparse.ArgumentParser(description='a command line tool for making pactolus tasks.  outputs one line for each file to run pactolus on')
    parser.add_argument('-i','--indir', help='mzml folder to search recursively for files', required=True)
    parser.add_argument('-rr','--rerun',help='rerun already processed files',default=False)
    parser.add_argument('-o','--outfile', help='shell script to output', required=True)
    parser.add_argument('-m2t','--ms2_tolerance', help='tolerance in Daltons for ms2', type=float,default=0.01)
    parser.add_argument('-m1t','--ms1_tolerance', help='tolerance in Daltons for ms1', type=float,default=0.01)
    parser.add_argument('-m1pn','--ms1_pos_neutralizations', help='adducts to neutralize for in ms1: 1.007276,18.033823,22.989218', type=float,nargs='+',default=[1.007276,18.033823,22.989218])
    parser.add_argument('-m2pn','--ms2_pos_neutralizations', help='ionization states to neutralize for in ms2: -1.00727646677,-2.0151015067699998,0.00054857990946', type=float,nargs='+',default=[-1.00727646677,-2.0151015067699998,0.00054857990946])
    parser.add_argument('-m1nn','--ms1_neg_neutralizations', help='adducts to neutralize for in ms1: -1.007276, 59.013851', type=float,nargs='+',default=[-1.007276, 59.013851])
    parser.add_argument('-m2nn','--ms2_neg_neutralizations', help='ionization states to neutralize for in ms2: 1.00727646677,2.0151015067699998,-0.00054857990946', type=float,nargs='+',default=[1.00727646677,2.0151015067699998,-0.00054857990946])
    parser.add_argument('-t','--tree_file', help='tree file: /project/projectdirs/metatlas/projects/clean_pactolus_trees/tree_lookup.npy', default='/project/projectdirs/metatlas/projects/clean_pactolus_trees/tree_lookup.npy')
    parser.add_argument('-n','--num_cores', help='number of cores to use for multiprocessing', type=int,default=8)
    parser.add_argument('-pid','--project_id', help='project id', type=str,default='m1541')
    parser.add_argument('-partition','--partition', help='nersc partition (realtime, debug, norma)', type=str,default='realtime')
    parser.add_argument('-time','--wall_time', help='job wall time', type=str,default='03:59:00')
    parser.add_argument('-job_name','--job_name', help='job name', type=str,default='pactolus')
    parser.add_argument('-log_path','--log_path', help='path to log directory', type=str,default='/project/projectdirs/metatlas/projects/logs/pactolus_logs/')
    parser.add_argument('-overwrite','--overwrite', help='Overwrite pre-existing file(s): True/False', type=bool,default=False)
    parser.add_argument('-python_binary','--python_binary', help='path to the python binary you want to use', type=str,default='/global/common/software/m2650/python-cori/bin/python')
    parser.add_argument('-scoring','--scoring', help='full path to scoring function (/global/homes/b/bpb/repos/pactolus/pactolus/score_mzmlfile.py)', type=str,default='/global/homes/b/bpb/repos/pactolus/pactolus/score_mzmlfile.py')
    args = vars(parser.parse_args())
    files = find_files(args['indir'],'*.mzML')
    files = list(files)
    print(len(files),'files with mzML')
    done_files = find_files(args['indir'],'*.pactolus.gz') #find all the completed files
    done_files = list(done_files)
    no_extension_files = list(set([f.split('.')[0] for f in files]) - set([f.split('.')[0] for f in done_files]))
    prescreened_files = ['%s.mzML'%f for f in no_extension_files]

    paths = [os.path.split('%s.mzML'%f)[0] for f in no_extension_files]
    paths = list(set(paths))
    bad_paths = []
    for p in paths:
        if not os.access(p, os.W_OK):
            print('not writable',p)
            bad_paths.append(p)

    files = []
    for f in prescreened_files:
        if not os.path.split(f)[0] in bad_paths:
            files.append(f)

    # files.extend(glob.glob(os.path.join(args['indir'],'*.mzml')))
    # print(len(files),'files with mzml')

    print(len(files),'files to process')
    with open(args['outfile'],'w') as script_fid:
        for file in files:

            job_script = create_job_script(infile=file,
                                            errlog=os.path.join(args['log_path'],os.path.basename(file)+'.err'),
                                            outlog=os.path.join(args['log_path'],os.path.basename(file)+'.out'),
                                            ms1_tolerance=args['ms1_tolerance'],
                                            ms2_tolerance=args['ms2_tolerance'],
                                            ms1_pos_neutralizations=args['ms1_pos_neutralizations'],
                                            ms1_neg_neutralizations=args['ms1_neg_neutralizations'],
                                            ms2_pos_neutralizations=args['ms2_pos_neutralizations'],
                                            ms2_neg_neutralizations=args['ms2_neg_neutralizations'],
                                            tree_file=args['tree_file'],
                                            num_cores=args['num_cores'],
                                            project_id=args['project_id'],
                                            time=args['wall_time'],
                                            job_name=args['job_name'],
                                            python=args['python_binary'],
                                            pactolus_script=args['scoring'],
                                            partition=args['partition'],
                                            cores=args['num_cores'],
                                            use_haswell=True,
                                            license_project=True,
                                            license_scratch=False)
            script_fid.write('%s'%job_script)

if __name__ == "__main__":
    main()


