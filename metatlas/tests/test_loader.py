from __future__ import print_function
import os
import sys

import requests
import tables

from metatlas.mzml_loader import mzml_to_hdf


def test_loader():
    dname = os.path.dirname(__file__)
    path = os.path.join(dname, 'test.mzML')

    url = ("https://drive.google.com/uc?"
           "export=download&id=0B2pT935MmTv2TlZxcTBkdGczWHM")

    if not os.path.exists(path):
        # NOTE the stream=True parameter
        print('Downloading: %s\n' % url, file=sys.stderr)
        r = requests.get(url, stream=True)
        with open(path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)
                    f.flush()
        print('Download complete\n', file=sys.stderr)

    mzml_to_hdf(path)

    out_file = path.replace('.mzML', '.h5')

    fid = tables.open_file(out_file)
    table = fid.root.spectra

    assert table.nrows == 933367
    assert table[0][0] == 59.01387023925781
    assert table[-1][0] == 1666.9520263671875
    scan_time = [y['scan_time'] for y in table.where('(ms_level==1)')]
    assert len(scan_time) == 5082

    os.chmod(out_file, 777)
    os.remove(out_file)
