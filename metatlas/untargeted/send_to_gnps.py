from __future__ import absolute_import
from __future__ import print_function
import sys
import os
import argparse
# repos/metatlas/metatlas/untargeted/gnps_quickstart_util.py
sys.path.insert(0,'/global/homes/b/bpb/repos/metatlas/metatlas/untargeted')
exec(open('/global/homes/b/bpb/repos/metatlas/metatlas/untargeted/gnps_quickstart_util.py').read())

import uuid

with open('/global/homes/b/bpb/gnps_password.txt','r') as fid:
    gnps_password = fid.read().strip()

# execfile('/global/homes/b/bpb/repos/GNPS_quickstart/util.py')

def copy_and_submit_to_gnps(basedir,basename,override=False):
    new_basename = os.path.join(basedir,basename)
    taskid_filename = '%s_%s.txt'%(new_basename,'gnps-uuid')
    if (os.path.isfile(taskid_filename)) & (override==False):
        print('Already submitted')
        with open(taskid_filename,'r') as fid:
            print(fid.read())
        return False
    # files_filename = '%s_filelist.txt'%new_basename
    # if os.path.isfile(files_filename):
    #     with open(files_filename,'r') as fid:
    #         files = fid.read()
    #     files = files.split('\n')
    #     for f in files:
    #         upload_to_gnps(f,basename,'rawdata','bpbowen',gnps_password)
    # else:
    #     import tempfile
    #     with tempfile.NamedTemporaryFile() as tmp:
    #         print(tmp.name)
    #         tmp.write('no raw data')
    #         upload_to_gnps(tmp.name,basename,'rawdata','bpbowen',gnps_password)
    
    metadata_filename = '%s_%s.tab'%(new_basename,'metadata')
    if os.path.isfile(metadata_filename):
        upload_to_gnps(metadata_filename,basename,'samplemetadata','bpbowen',gnps_password)
        
    else:
        print('METADATA NOT FOUND %s'%metadata_filename)
        return False
        
    mgf_filename = '%s_%s.mgf'%(new_basename,'MSMS')
    if os.path.isfile(mgf_filename):
        upload_to_gnps(mgf_filename,basename,'featurems2','bpbowen',gnps_password)
    else:
        print('SPECTRA NOT FOUND %s'%mgf_filename)
        return False
        
    quant_filename = '%s_%s.csv'%(new_basename,'peak-area')
    if os.path.isfile(quant_filename):
        upload_to_gnps(quant_filename,basename,'featurequantification','bpbowen',gnps_password)
    else:
        print('QUANT NOT FOUND %s'%quant_filename)
        return False
    
    task_id = launch_GNPS_featurenetworking_workflow( basename,basename, 
                                           'bpbowen', gnps_password, 'ben.bowen@gmail.com',
                                           'MZMINE2', [], 'HIGHRES',
                                           os.path.basename(metadata_filename),
                                           os.path.basename(mgf_filename),
                                           os.path.basename(quant_filename))
#                                           uuid.uuid1()) #I don't think this uuid does anything
    if task_id is not None:
        with open(taskid_filename,'w') as fid:
            fid.write('https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=%s'%task_id)
        # SUCCESS!!!
        return task_id
    else:
        return False

    
      
def main():
    # print command line arguments
    parser = argparse.ArgumentParser(description='a tool for submitting jobs to GNPS at nersc.')
    # parser.add_argument('-m2t','--ms2_tolerance', help='tolerance in Daltons for ms2', type=float,default=0.01)
    # parser.add_argument('-m1pn','--ms1_pos_neutralizations', help='adducts to neutralize for in ms1: 1.007276,18.033823,22.989218', type=float,nargs='+',default=[1.007276,18.033823,22.989218])
    parser.add_argument('-basedir','--basedir', help='Directory with mzmine outputs', default='')
    parser.add_argument('-basename','--basename', help='Prefix Filename of all outputs', default='')
    parser.add_argument('-override','--override', help='Override GNPS job', default=False)
   
    args = vars(parser.parse_args())
    # trees = np.load(args['tree_file'])
    print(args)
    taskid = copy_and_submit_to_gnps(args['basedir'],args['basename'],override=args['override'])
    print(taskid)

if __name__ == "__main__":
    main()

    
    
