import os
import shutil
import glob
import simplejson as json
import tempfile

import tables

from .trns_transform_mzML_LCMS_to_MetaboliteAtlas2_MAFileInfo import \
    transform

from .kbase_utils import (
    get_object_uid, SHOCK_URL, save_ws_object, getHandles, get_object,
    get_object_from_ref, download_file_from_shock)


def create_compound(name, formula, adducts, mz, mz_threshold, rt_min,
                    rt_max, rt_peak, neutral_mass, pubchem_id):
    """Create a MetAtlas compound.

    Parameters
    ----------
    name : str
        Common name of the compound.
    formula : str
        Chemical forumla.
    adducts : str
        Adduct ions.
    mz : float
        Mass-to-charge ratio.
    mz_threshold : float
        Threshold in ppm.
    rt_min : float
        Min retention time (minutes).
    rt_max : float
        Max retention time (minutes).
    rt_peak : float
        Peak retention time (minutes).
    pubchem_id : str
        Pubchem ID for compound.

    Returns
    -------
    name, id : tuple
        Name of compound and kbase id.
    """
    dictData = dict(name=name, formula=formula, adducts=adducts, mz=mz,
                    mz_threshold=mz_threshold, rt_min=rt_min,
                    rt_max=rt_max, rt_peak=rt_peak,
                    neutral_mass=neutral_mass, pubchem_id=pubchem_id)
    dict_save_params = dict(type='MetaboliteAtlas2.MACompound',
                            data=dictData, name=name, hidden=0)
    save_ws_object(dict_save_params)
    return name, get_object_uid(name)


def create_atlas(name, compounds, sample, method):
    """Create a MetAtlas Dictionary.

    Parameters
    ----------
    name : str
        Common name of the atlas.
    compounds : list of strings
        Names of compounds (names of kbase workspace objects).
    sample : str
        Description of sample.
    method : str
        Description of method.

    Returns
    -------
    name, id : tuple
        Name of atlas and kbase id.
    """
    compounds = [get_object_uid(c) for c in compounds]
    dictData = dict(name=name, compounds=compounds, sample=sample,
                    method=method)
    dict_save_params = dict(type='MetaboliteAtlas2.MADictionary',
                            data=dictData, name=name, hidden=0)
    save_ws_object(dict_save_params)
    return name, get_object_uid(name)


def create_experiment(name, description, reconstitution_method,
                      quenching_method, extraction_method,
                      chromatography_method, atlases):
    """Create a MetAtlas Experiment.

    Parameters
    ----------
    name : str
        Name of the experiment.
    description : str
        Description of experiment.
    reconstitution_method : str
        Reconstitution method used in experiment.
    quenching_method : str
        Quenching method used in experiment.
    extraction_method : str
        Extraction method used in experiment.
    chromatography_method : str
        Chromatography method used in experiment.
    atlases : list of strings
        Names of atlases (names of kbase workspace objects)

    Returns
    -------
    name, id : tuple
        Name of experiment and kbase id.
    """
    atlases = [get_object_uid(a) for a in atlases]
    dictData = dict(name=name, description=description,
                    reconstitution_method=reconstitution_method,
                    quenching_method=quenching_method,
                    extraction_method=extraction_method,
                    chromatography_method=chromatography_method,
                    atlas_ids=atlases)
    dict_save_params = dict(type='MetaboliteAtlas2.MAExperiment',
                            data=dictData, name=name, hidden=0)
    save_ws_object(dict_save_params)
    return name, get_object_uid(name)


def create_fileinfo(input_file, mzml_file_name=None, polarity=None,
                    atlases=None, group=None, inclusion_order=None,
                    normalization_factor=None, retention_correction=None):
    """Create a MetAtlas FileInfo object.

    Parameters
    ----------
    input_file : str
        Path to input_file.
    mzml_file_name : str, optional
        Name of the file object (defaults to root of file name).
    atlases : list of strings, optional
        Names of atlas dictionaries (kbase workspace names).
    polarity : int, optional
        Polarity of ions in file.
    group : str, optional
        Group.
    inclusion_order : int, optional
        Inclusion order.
    normalization_factor : float, optional
        Normalization factor.
    retention_correction : float, optional
        Retention correction factor

    Returns
    -------
    name, id : tuple
        Name of finfo object and kbase id.
    """
    tempdir = tempfile.mkdtemp()
    # create the file in a directory
    shutil.copy2(input_file, tempdir)
    # pass the directory to transform
    transform(input_directory=tempdir,
              shock_service_url=SHOCK_URL,
              working_directory=tempdir,
              mzml_file_name=mzml_file_name, polarity=polarity,
              atlases=atlases,
              group=group, inclusion_order=inclusion_order,
              normalization_factor=normalization_factor,
              retention_correction=retention_correction)
    # get the resulting json
    output_file = glob.glob('%s/*.json' % tempdir)[0]
    with open(output_file) as fid:
        data = json.load(fid)

    if not SHOCK_URL:
        # move the HDF file to the input directory
        input_dir = os.path.dirname(mzml_file_name)
        shutil.copy2(data['run_file_id'], input_dir)
        base = os.path.basename(data['run_file_id'])
        data['run_file_id'] = os.path.join(input_dir, base)

    shutil.rmtree(tempdir)

    data['atlases'] = [get_object_uid(a) for a in data['atlases']]

    dict_save_params = dict(type='MetaboliteAtlas2.MAFileInfo',
                            data=data, name=data['mzml_file_name'], hidden=0)
    save_ws_object(dict_save_params)

    return data['mzml_file_name'], get_object_uid(data['mzml_file_name'])


def get_compound(name):
    """Get a workspace compound by name as a dictionary"""
    return get_object(name)


def get_atlas(name):
    """Get a workspace atlas by name as a nested dictionary"""
    obj = get_object(name)
    print(obj['compounds'])
    obj['compounds'] = [get_object_from_ref(c) for c in obj['compounds']]
    return obj


def get_experiment(name):
    """Get a workspace experiment by name as a nested dictionary"""
    obj = get_object(name)
    obj['atlas_ids'] = [get_object_from_ref(a) for a in obj['atlas_ids']]
    for a in obj['atlas_ids']:
        a['compounds'] = [get_object_from_ref(c) for c in a['compounds']]
    return obj


def get_finfo(name):
    """Get a workspace fileinfo by name as a nested dictionary.

    Contains a "fid", which is an open pytables object.
    """
    obj = get_object(name)
    for a in obj['atlases']:
        a['compounds'] = [get_object_from_ref(c) for c in a['compounds']]
    if SHOCK_URL:
        shock_id = getHandles(handle_ids=[obj['run_file_id']])[0]['id']
        filename = obj['name'] + '.h5'
        download_file_from_shock(shock_id=shock_id, filename=filename)
    else:
        filename = obj['run_file_id']
    obj['fid'] = tables.open_file(filename)
    return obj
