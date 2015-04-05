from __future__ import print_function
import tables

from metatlas.mzml_loader import mzml_to_hdf, get_test_data


def test_loader():
    path = get_test_data()

    mzml_to_hdf(path)

    out_file = path.replace('.mzML', '.h5')

    fid = tables.open_file(out_file)
    table = fid.root.spectra

    assert table.nrows == 933367
    assert table[0][0] == 59.01387023925781
    assert table[-1][0] == 1666.9520263671875
    scan_time = [y['scan_time'] for y in table.where('(ms_level==1)')]
    assert len(scan_time) == 5082
