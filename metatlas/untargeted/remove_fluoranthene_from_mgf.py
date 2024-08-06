import sys
#sys.path.insert(0,'/global/common/software/m2650/metatlas-repo/metatlas')
sys.path.insert(0,'/global/homes/b/bkieft/metatlas/')
from metatlas.untargeted import tools as mzm

#### Remove contaminants from recent mgf files (202 ions for fluoranthene)
print("Removing fluoranthene from recent mgf files in mzmine results...")
recent_files = mzm.get_recent_mgf_files(time_back=2)
for f in recent_files:
    if f:
        mzm.remove_contaminant_from_mgf(f)