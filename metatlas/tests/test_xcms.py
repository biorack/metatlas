from __future__ import print_function

from metatlas.xcms import get_xmcs_set, group, peak_table
from metatlas.mzml_loader import get_test_data


def test_xcms():
    path = get_test_data()['basic']
    xset = get_xmcs_set([path], 'neg')
    xset = group(xset)
    df = peak_table(xset)
    assert 'mz' in df.columns
    assert 'rtmin' in df.columns
    assert df.shape == (279, 13)
