# pylint: disable=missing-function-docstring, missing-module-docstring, line-too-long

"""unit tests of atlas prefilter functions"""

from metatlas.tools import atlas_prefilter as pf
from metatlas.io.metatlas_get_data_helper_fun import make_atlas_df
import numpy as np
import pandas as pd


def test_order_ms2_spectrum():
    unordered_array = np.array([[1e+01, 1e+03, 1e+05, 1e+04, 1e+02], [1e+05, 1e+03, 1e+01, 1e+02, 1e+04]])
    ordered_array = pf.order_ms2_spectrum(unordered_array)

    assert np.all(ordered_array == np.array([[1e+01, 1e+02, 1e+03, 1e+04, 1e+05], [1e+05, 1e+04, 1e+03, 1e+02, 1e+01]]))


def test_get_sample_file_paths(analysis_ids, lcmsrun):
    sample_file_paths = pf.get_sample_file_paths(analysis_ids)

    assert sample_file_paths == [lcmsrun.hdf5_file]


def test_score_ms2_data_01(mocker, msms_refs, atlas_with_2_cids):
    mocker.patch('pandas.read_csv', return_value=msms_refs)
    atlas_df = make_atlas_df(atlas_with_2_cids)
    ms2_data = pd.DataFrame(columns=['label', 'spectrum', 'rt', 'precursor_mz', 'precursor_peak_height'])

    ms2_data_scored = pf.score_ms2_data(ms2_data, atlas_df, 'positive', '/dummy/path/to/msms_refs', 0.02)

    assert ms2_data_scored.empty


def test_score_ms2_data_02(mocker, msms_refs, atlas_with_2_cids, pf_ms2_pos):
    mocker.patch('pandas.read_csv', return_value=pd.DataFrame(columns=msms_refs.columns))
    atlas_df = make_atlas_df(atlas_with_2_cids)

    ms2_data_scored = pf.score_ms2_data(pf_ms2_pos, atlas_df, 'positive', '/dummy/path/to/msms_refs', 0.02)

    assert ms2_data_scored.empty
