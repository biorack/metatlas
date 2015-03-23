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
    assert list(table[0]) == (59.0, 0.03933333232998848, 9, 0.0, 1,
                              0.0, 0.0, 0.0)
    assert list(table[-1]) == (1666.0, 19.266700744628906, 225, 0.0, 1,
                               0.0, 0.0, 0.0)
    ms_1 = [y['scan_time'] for y in table.where('(ms_level==1)')]
    assert len(ms_1) == 5082

    os.chmod(out_file, 777)
    os.remove(out_file)
