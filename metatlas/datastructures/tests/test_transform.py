from __future__ import print_function
from __future__ import absolute_import
import os
import sys

import simplejson

from metatlas.kbase.trns_transform_mzML_LCMS_to_MetaboliteAtlas2_MAFileInfo \
    import transform, main
from metatlas.mzml_loader import get_test_data


def test_transform():
    path = prep()

    setup = dict(polarity='negative',
                 group='test',
                 retention_correction='none',
                 normalization_factor='11')

    transform(input_directory=os.path.dirname(path),
              working_directory=os.path.dirname(path),
              **setup)
    check_dir(path, **setup)


def prep():
    path = get_test_data()['basic']
    if os.path.exists(path.replace('.mzML', '_finfo.json')):
        os.remove(path.replace('.mzML', '_finfo.json'))
    return path


def check_dir(path, **setup):
    assert os.path.exists(path.replace('.mzML', '_finfo.json'))

    with open(path.replace('.mzML', '_finfo.json')) as fid:
        data = simplejson.load(fid)

    print(data)
    assert data['mzml_file_name'] == os.path.basename(path).replace('.mzML',
                                                                    '')
    assert data['atlases'] == []

    for (key, value) in setup.items():
        assert data[key] == value


def test_transform_main():
    return
    path = prep()
    sys.argv = ['transform', '--shock_service_url', '',
                '--input_directory', os.path.dirname(path),
                '--working_directory', os.path.dirname(path),
                '--debug']

    main()
    check_dir(path)
