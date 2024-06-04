# pylint: disable=missing-function-docstring, missing-module-docstring, line-too-long

"""unit tests of atlas prefilter functions"""

from metatlas.tools import atlas_prefilter as pf
import numpy as np
import pandas as pd


def test_order_ms2_spectrum():

    unordered_array = np.array([[1e+01,1e+03,1e+05,1e+04,1e+02],[1e+05,1e+03,1e+01,1e+02,1e+04]])

    ordered_array = pf.order_ms2_spectrum(unordered_array)

    assert np.all(ordered_array == np.array([[1e+01,1e+02,1e+03,1e+04,1e+05],[1e+05,1e+04,1e+03,1e+02,1e+01]]))

def test_get_sample_file_paths(analysis_ids, lcmsrun):

    sample_file_paths = pf.get_sample_file_paths(analysis_ids)
    
    assert sample_file_paths == [lcmsrun.hdf5_file]