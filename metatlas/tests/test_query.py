from __future__ import print_function

from metatlas.mzml_loader import mzml_to_hdf, get_test_data
from metatlas.h5_query import get_XICof


def rmse(target, predictions):
    target = target / target.max()
    predictions = predictions / predictions.max()

    return np.sqrt(((predictions - targets) ** 2).mean())


def test_xicof():
    return
    fid = tables.open_file('140808_1_RCH2_neg.h5')
    x, y = get_XICof(fid, 1, 1000, 1, 0)
    
    xicof_scidb = np.load('xicof_scidb.npy')

    assert rmse(y, xicof_scidb[:, 1]) < 0.01

