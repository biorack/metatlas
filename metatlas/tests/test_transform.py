from __future__ import print_function
import os
import sys
import simplejson

import tables

from metatlas.trns_transform_mzML_LCMS_to_MetaboliteAtlas2_MAFileInfo import transform
from metatlas.mzml_loader import get_test_data


def test_transform():
    path = get_test_data()

    setup = dict(polarity='negative',
                      group='test',
                      retention_correction='none',
                      normalization_factor='11')

    transform(input_directory=os.path.dirname(path),
                      working_directory=os.path.dirname(path),
                      **setup)

    assert os.path.exists(path.replace('.mzML', '_finfo.json'))

    with open(path.replace('.mzML', '_finfo.json')) as fid:
        data = simplejson.load(fid)

    print(data)
    assert data['name'] == os.path.basename(path).replace('.mzML',  '')
    assert data['atlases'] == []
    assert data['run_file_id'] == 'dummy_shock_id'

    for (key, value) in setup.items():
        assert data[key] == value
