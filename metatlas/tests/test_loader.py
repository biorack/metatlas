from __future__ import print_function

import sys

import tables

from metatlas.mzml_loader import mzml_to_hdf, get_test_data, main


def test_loader():
    path = get_test_data()
    out_file = 'test_loader.h5'

    mzml_to_hdf(path, out_file_name=out_file)

    fid = tables.open_file(out_file)
    table = fid.root.spectra

    assert table.nrows == 933367
    assert table[0][0] == 59.01387023925781
    assert table[-1][0] == 1666.9520263671875
    scan_time = [y['rt'] for y in table.where('(ms_level==1)')]
    assert len(scan_time) == 933367

    table.close()


def test_loader_main():
    path = get_test_data()
    out_file = 'test_loader.h5'

    sys.argv = ['dummy', '--input', path,
                '--output', out_file,
                '--debug']
    main()

    fid = tables.open_file(out_file)
    table = fid.root.spectra

    assert table.nrows == 933367
    assert table[0][0] == 59.01387023925781
    assert table[-1][0] == 1666.9520263671875
    scan_time = [y['rt'] for y in table.where('(ms_level==1)')]
    assert len(scan_time) == 933367

    table.close()
