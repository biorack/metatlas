from __future__ import print_function
from __future__ import absolute_import
import os

import numpy as np
import tables

from metatlas.mzml_loader import mzml_to_hdf, get_test_data
from metatlas.h5_query import (
    get_chromatogram, get_data, get_spectrogram, get_heatmap,
    get_info)
from metatlas.plotting import (
    plot_heatmap, plot_spectrogram, plot_chromatogram
)

fid = None


def setup():
    global fid
    paths = get_test_data()
    out_file = 'test_query.h5'

    mzml_to_hdf(paths['basic'], out_file_name=out_file)
    fid = tables.open_file('test_query.h5')


def teardown():
    fid.close()


def rmse(targets, predictions):
    targets = targets / targets.max()
    predictions = predictions / predictions.max()

    return np.sqrt(((predictions - targets) ** 2).mean())


def test_XIC():
    x, y = get_chromatogram(fid, 1, 1000)
    dname = os.path.dirname(__file__)
    xicof_scidb = np.load(os.path.join(dname, 'xic_scidb.npy'))

    assert rmse(y, xicof_scidb[:, 1]) < 0.06
    plot_chromatogram(x, y)


def test_BPC():
    x, y = get_chromatogram(fid, 1, 1000, np.max)

    assert y.max() > 2.5e+06
    assert y.max() < 2.6e+06
    plot_chromatogram(x, y, title='BPC for Sample')


def test_spectrogram():
    x, y = get_spectrogram(fid, 1, 5)

    assert np.allclose(x.mean(), 855.718857765)
    assert y.min() >= 0

    plot_spectrogram(x, y)


def test_heatmap():
    data = get_heatmap(fid, 1000)

    assert np.allclose(data['arr'].mean(), 3790.08673939)
    assert np.allclose(data['mz_bins'][0], 30.004)
    assert np.allclose(data['rt_bins'][-1], 19.2667)

    data = get_heatmap(fid, 1000, min_mz=50)
    assert np.allclose(data['mz_bins'][0], 50.0002)
    plot_heatmap(data['arr'], data['rt_bins'], data['mz_bins'])


def test_get_data():
    data = get_data(fid, min_rt=5, min_mz=100)
    assert np.allclose(data['i'].mean(), 7825.55387233)
    assert np.allclose(data['mz'][0], 100.979026794)
    assert np.allclose(data['rt'][0], 5.00666666031)


def test_get_info():
    data = get_info(fid)
    assert data['ms2_pos']['nrows'] == 0
    assert data['ms1_pos']['nrows'] == 0
    assert data['ms1_neg']['nrows'] == 933367
    assert data['ms2_neg']['nrows'] == 0
