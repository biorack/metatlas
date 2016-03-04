from __future__ import print_function

import sys

import tables
from numpy.testing import assert_raises

from metatlas.mzml_loader import mzml_to_hdf, get_test_data, main
from metatlas.h5_query import get_info


def test_loader():
    paths = get_test_data()
    out_file = 'test_loader.h5'

    mzml_to_hdf(paths['basic'], out_file_name=out_file)

    check_file(out_file)


def check_file(out_file):
    fid = tables.open_file(out_file)
    table = fid.root.ms1_neg

    assert table.nrows == 933367
    assert table[0][0] == 59.01387023925781
    assert table[-1][0] == 1666.9520263671875
    scan_time = table.read()['rt']
    assert len(scan_time) == 933367


def test_loader_main():
    paths = get_test_data()
    out_file = 'test_loader_main.h5'

    sys.argv = ['dummy', '--input', paths['basic'],
                '--output', out_file,
                '--debug']
    main()
    check_file(out_file)


def test_invalid_file():
    paths = get_test_data()
    out_file = 'test_invalid.h5'

    assert_raises(TypeError, mzml_to_hdf, paths['wrong_fmt'],
                  out_file_name=out_file)


def test_32_64():
    paths = get_test_data()
    fname = 'mix32_64.h5'

    mzml_to_hdf(paths['mix32_64'], out_file_name=fname)

    fid = tables.open_file(fname)
    print(fid)
    table = fid.root.ms1_pos

    assert table.nrows == 1041183, table.nrows

    table2 = fid.root.ms1_neg

    assert table2.nrows == 762699, table2.nrows


def test_ms_ms():
    paths = get_test_data()
    file2 = 'ms_ms.h5'

    mzml_to_hdf(paths['ms_ms'], out_file_name=file2)

    fid = tables.open_file(file2)

    table = fid.root.ms2_pos

    assert table.nrows == 604775, table.nrows


