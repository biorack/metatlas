import os
import sys
import shutil
import glob
import simplejson as json
import tempfile

from .kbase_utils import get_object_uid, SHOCK_URL, save_ws_object
from .trns_transform_mzML_LCMS_to_MetaboliteAtlas2_MAFileInfo import \
    transform


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
    dict_save_params = dict(type='MetaboliteAtlas2.MACompound-0.1',
                            data=dictData, name='test_dict2', hidden=0)
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
    dict_save_params = dict(type='MetaboliteAtlas2.MADictionary-2.0',
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
                    atlases=atlases)
    dict_save_params = dict(type='MetaboliteAtlas2.MAExperiment-2.0',
                            data=dictData, name=name, hidden=0)
    save_ws_object(dict_save_params)
    return name, get_object_uid(name)


def create_fileinfo(input_file, name='', polarity='',
                    atlases=None, group='', inclusion_order='',
                    normalization_factor='', retention_correction=''):
    """Create a MetAtlas FileInfo object.

    Parameters
    ----------
    input_file : str
        Path to input_file.
    name : str, optional
        Name of the file object (defaults to root of file name).
    atlases : list of strings, optional
        Names of atlas dictionaries (kbase workspace names).
    polarity : str, optional
        Polarity of ions in file.
    group : str, optional
        Group.
    inclusion_order : str, optional
        Inclusion order.
    normalization_factor : str, optional
        Normalization factor.
    retention_correction : str, optional
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
              name=name, polarity=polarity, atlases=atlases,
              group=group, inclusion_order=inclusion_order,
              normalization_factor=normalization_factor,
              retention_correction=retention_correction)
    # get the resulting json
    output_file = glob.glob('%s/*.json' % tempdir)[0]
    with open(output_file) as fid:
        data = json.load(fid)

    shutil.rmtree(tempdir)

    data['atlases'] = [get_object_uid(a) for a in data['atlases']]
    dict_save_params = dict(type='MetaboliteAtlas2.MAFileInfo-0.2',
                            data=data, name=data['name'], hidden=0)
    save_ws_object(dict_save_params)

    return data['name'], get_object_uid(data['name'])


def get_from_nersc(user, relative_path):
    """Load a data file from NERSC to an H5 file

    Parameters
    ----------
    user : str
        NERSC user account
    relative_path : str
        Path to file from "/project/projectdirs/metatlas/original_data/<user>/"
    """
    import pexpect
    from IPython.display import clear_output

    cmd = 'scp -o StrictHostKeyChecking=no '
    path = "/project/projectdirs/metatlas/original_data/%s/%s"
    path = path % (user, relative_path)
    cmd += '%s@edisongrid.nersc.gov:%s . && echo "Download Complete"'
    cmd = cmd % (user, path)
    print(cmd)
    return
    proc = pexpect.spawn(cmd)
    proc.expect("assword:*")
    if sys.version.startswith('3'):
        passwd = input()
    else:
        passwd = raw_input()
    clear_output()
    proc.send(passwd)
    proc.send('\r')
    proc.expect('Download Complete')
    proc.close()
    return os.path.abspath(os.path.basename(relative_path))
