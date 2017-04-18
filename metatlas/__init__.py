__version__ = '0.7'
import os
# set up the NERSC enviroment
if os.path.exists('/global/project'):
    os.environ['R_LIBS_USER'] = '/global/project/projectdirs/metatlas/r_pkgs/'

from .mzml_loader import mzml_to_hdf
from .plotting import plot_heatmap, plot_spectrogram, plot_chromatogram
from .h5_query import (
    get_data, get_chromatogram, get_heatmap, get_spectrogram, get_info)
from .metatlas_objects import (
    Method, Sample, LcmsRun, FunctionalSet,
    Compound, Reference, IdentificationGrade, CompoundIdentification, Atlas,
    Group, MzIntensityPair, FragmentationReference, RtReference,
    MzReference, retrieve, store, remove, remove_objects, workspace,
    to_dataframe
)
#ReferenceDatabase, 
from .gui import show_experiments, show_lcms_run, edit_objects
