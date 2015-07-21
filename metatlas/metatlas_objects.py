import os
import dataset
import pickle
from .mzml_loader import mzml_to_hdf
from IPython.utils.traitlets import (
    HasTraits, CUnicode, CFloat, List, CInt, Bool, Instance)

NERSC_WORKSPACE = '/project/projectdirs/metatlas/workspace'


class Compound(HasTraits):

    """
    name: str
      Name of the compound
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
    """
    name = CUnicode()
    formula = CUnicode()
    adducts = CUnicode()
    mz = CFloat()
    mz_threshold = CFloat()
    rt_min = CFloat()
    rt_max = CFloat()
    rt_peak = CFloat()
    neutral_mass = CFloat()
    pubchem_id = CFloat()


class Atlas(HasTraits):

    """
    name : str
        Common name of the atlas.
    compounds : list of strings
        Names of compounds (names of kbase workspace objects).
    sample : str
        Description of sample.
    method : str
        Description of method.
    """
    name = CUnicode()
    compounds = List(Instance(Compound))
    sample = CUnicode()
    method = CUnicode()


class FileSpec(HasTraits):

    """
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
    """
    polarity = CInt()
    group = CUnicode()
    inclusion_order = CInt()
    normalization_factor = CFloat()
    retention_correction = CFloat()


class FInfo(FileSpec):

    """
    mzml_file : str
        Path to input_file.
    hdf_file : str
        Name of the hdf_file.
    """ + FileSpec.__doc__

    mzml_file = CUnicode()
    hdf_file = CUnicode()

    def parse(self):
        """Parse a file info spec"""
        self.hdf_file = mzml_to_hdf(self.mzml_file, self.hdf_file or None)


class Experiment(HasTraits):

    """
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
    atlases : list of atlases
        List of atlas objects.
    finfos: list of finfos
        List of finfo objects.
    """
    name = CUnicode()
    description = CUnicode()
    reconstitution_method = CUnicode()
    quenching_method = CUnicode()
    extraction_method = CUnicode()
    chromatography_method = CUnicode()
    atlases = List(Instance(Atlas))
    finfos = List(Instance(FInfo))
    _shared = Bool()

    def save(self):
        """Save the experiment to the workspace"""
        _WORKSPACE.save(self)

    def share(self):
        """Share the experiment"""
        pass

    def load_files(self, files, filespec):
        """Load a list of files into the experiment

        Parameters
        ----------
        files : list of str
            List of file paths to load.
        filespec: FileSpec
            The file spec for the files.
        """
        state = filespec.__getstate__()
        for fname in files:
            finfo = FInfo(mzml_file=fname, **state)
            finfo.parse()
            self.finfos.append(finfo)


def get_experiment(name, username=None):
    """Get an experiment by name.

    Parameters
    ----------
    name : str
        Name of the experiment.
    username: str, optional
        Username of original experimenter.

    Returns
    -------
    out : Experiment
        Experiment object.
    """
    if username:
        raise NotImplemented
    return _WORKSPACE.load(name)


class _Workspace(object):

    def __init__(self):
        path = os.environ['USER'] + '.db'
        # allow for fallback when not on NERC
        if os.path.exists(NERSC_WORKSPACE):
            path = os.path.join(NERSC_WORKSPACE, path)
        self.db = dataset.connect('sqlite:///%s' % path)
        os.chmod(path, 0o775)

    def save(self, experiment):
        top = pickle.dumps(experiment)
        finfos = [pickle.dumps(f) for f in experiment.finfos]
        finfos = pickle.dumps(finfos)
        atlases = []
        for atlas in experiment.atlases:
            compounds = [pickle.dumps(c) for c in atlas.compounds]
            atlases.append([pickle.dumps(atlas)] + compounds)
        atlases = pickle.dumps(atlases)
        self.db[experiment.name].insert(dict(top=top, finfos=finfos,
                                             atlases=atlases))

    def load(self, name):
        objects = self.db[name].all()
        if objects:
            obj = list(objects)[-1]
            experiment = pickle.loads(obj['top'])
            finfos = pickle.loads(obj['finfos'])
            for f in finfos:
                finfo = pickle.loads(f)
                experiment.finfos.append(finfo)
            atlas_blobs = pickle.loads(obj['atlases'])
            for blob in atlas_blobs:
                atlas = pickle.loads(blob[0])
                experiment.atlases.append(atlas)
                for b in blob[1:]:
                    compound = pickle.loads(b)
                    atlas.compounds.append(compound)
            return experiment


# Singleton Workspace object
_WORKSPACE = _Workspace()
