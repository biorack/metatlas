from __future__ import print_function
import os

import numpy as np
from numpy.testing.decorators import skipif
import tables

from metatlas.mzml_loader import mzml_to_hdf, get_test_data
from metatlas.h5_query import (
    get_chromatogram, get_data, get_spectrogram, get_heatmap,
    get_info)
from metatlas.plotting import (
    plot_heatmap, plot_spectrogram, plot_chromatogram
)


from metatlas.helpers import dill2plots as dp

from metatlas.helpers import metatlas_get_data_helper_fun as ma_data
from metatlas import metatlas_objects as metob
from metatlas.object_helpers import ON_NERSC
from ipywidgets import interact, interactive, fixed



@skipif(not ON_NERSC)
def test_interact_get_metatlas_files():
    experiment = '%violacein%'
    name = '%_%'
    most_recent = False
    files = dp.get_metatlas_files(experiment=experiment, name=name, most_recent=most_recent)
    assert len(files) == 1350

    experiment = '%violacein%'
    name = '%_%'
    most_recent = True
    files = dp.get_metatlas_files(experiment=experiment, name=name, most_recent=most_recent)
    assert len(files) == 1350

    experiment = '%violacein%'
    name = '%_%'
    most_recent = False
    files = dp.get_metatlas_files(experiment=experiment, name=name, most_recent=most_recent)
    assert len(files) == 1350

    experiment = '%violacein%'
    name = '%_%'
    most_recent = False
    files = dp.get_metatlas_files(experiment=experiment, name=name, most_recent=most_recent)
    assert len(files) == 1350

    experiment = '%violacein%'
    name = '%_%'
    most_recent = False
    files = dp.get_metatlas_files(experiment=experiment, name=name, most_recent=most_recent)
    assert len(files) == 1350

    experiment = '%violacein%'
    name = '%_%'
    most_recent = False
    files = dp.get_metatlas_files(experiment=experiment, name=name, most_recent=most_recent)
    assert len(files) == 1350



@skipif(not ON_NERSC)
def test_make_empty_fileinfo_sheet():
    experiment = '%violacein%'
    name = '%_%'
    most_recent = False
    f = dp.get_metatlas_files(experiment=experiment, name=name, most_recent=most_recent)
    temp_file_name = '/global/homes/b/bpb/Downloads/empty_violacein_384_finfo.tab'
    dp.make_empty_fileinfo_sheet(temp_file_name, f)

    # see if the file exits
    assert os.path.exists(temp_file_name)



