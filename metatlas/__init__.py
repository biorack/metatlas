__version__ = '0.4'
import os
import sys
# set up the NERSC enviroment
if os.path.exists('/global/project'):
    # NOTE: this line needs to be added to $HOME/.ipython/default_profile/startup/ipython_startup.py:
    #  sys.path.insert(0, '/global/project/projectdirs/metatlas/python_pkgs/')
    os.environ['R_LIBS_USER'] = '/global/project/projectdirs/metatlas/r_pkgs/'
    """
    try:
        sys.path.remove('/anaconda/lib/python2.7/site-packages')
    except ValueError:
        pass
    sys.path.append('/global/project/projectdirs/openmsi/jupyterhub_libs/anaconda/lib/python2.7/site-packages')
    """

from .mzml_loader import mzml_to_hdf
from .h5_query import plot_heatmap, plot_spectrogram, plot_XIC
from .h5_query import get_data, get_XIC, get_heatmap, get_spectrogram
from .metatlas_objects import (
    Method, Sample, LcmsRun, ReferenceDatabase, FunctionalSet,
    Compound, Reference, IdentificationGrade, CompoundId, Atlas,
    Group, MzIntensityPair, FragmentationReference, RtReference,
    MzReference, queryDatabase
)



# plan of attack:
# make our own anaconda environment at /metatlas/anaconda
# let Ben manage it using conda install commands
# then, all we have to do is put that guy first
