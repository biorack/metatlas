__version__ = '0.4'
import os
# set up the NERSC enviroment
if os.path.exists('/global/project'):
    os.environ['R_LIBS_USER'] = '/global/project/projectdirs/metatlas/r_pkgs/'

from .mzml_loader import mzml_to_hdf
from .h5_query import plot_heatmap, plot_spectrogram, plot_XIC
from .h5_query import get_data, get_XIC, get_heatmap, get_spectrogram, get_info
from .metatlas_objects import (
    Method, Sample, LcmsRun, ReferenceDatabase, FunctionalSet,
    Compound, Reference, IdentificationGrade, CompoundId, Atlas,
    Group, MzIntensityPair, FragmentationReference, RtReference,
    MzReference, queryDatabase
)
