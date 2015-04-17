from __future__ import print_function
import os

import numpy as np
import tables

from metatlas.mzml_loader import mzml_to_hdf, get_test_data
from metatlas.h5_query import (
    get_XIC, get_data, get_spectrogram, get_HeatMapRTMZ)

fid = None


def setup():
    global fid
    path = get_test_data()
    out_file = 'test.h5'

    mzml_to_hdf(path, out_file_name=out_file)
    fid = tables.open_file('test.h5')


def rmse(targets, predictions):
    targets = targets / targets.max()
    predictions = predictions / predictions.max()

    return np.sqrt(((predictions - targets) ** 2).mean())


def test_xicof():
    x, y = get_XIC(fid, 1, 1000, 1, 0)
    dname = os.path.dirname(__file__)
    xicof_scidb = np.load(os.path.join(dname, 'xic_scidb.npy'))

    assert rmse(y, xicof_scidb[:, 1]) < 0.06


def test_spectrogram():
    x, y = get_spectrogram(fid, 1, 5, 1, 0)

    assert np.allclose(x.mean(), 855.718857765)
    assert y.min() >= 0
    assert y.max() <= 100


def test_heatmap():
    data = get_HeatMapRTMZ(fid, 1000, 1000, 1, 0)

    assert np.allclose(data['arr'].mean(), 8743.73010776)
    assert np.allclose(data['mz_bins'][0], 30.838549386)
    assert np.allclose(data['rt_bins'][-1], 19.2570870609)

    data = get_HeatMapRTMZ(fid, 1000, 1000, 1, 0, min_mz=50)
    assert np.allclose(data['mz_bins'][0], 50.8247537537)


def test_get_data():
    data = get_data(fid, 1, 0, min_rt=5, min_mz=100, precursor_MZ=0)
    assert np.allclose(data['i'].mean(), 7825.55387233)
    assert np.allclose(data['mz'][0], 100.979026794)
    assert np.allclose(data['rt'][-1], 5.00666666031)
