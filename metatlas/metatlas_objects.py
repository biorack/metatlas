import os
import dataset
import pickle
import random
import pandas as pd
from .mzml_loader import mzml_to_hdf
from IPython.utils.traitlets import (
    HasTraits, CUnicode, CFloat, List, CInt, Bool, Instance)

NERSC_WORKSPACE = '/project/projectdirs/metatlas/workspace'


class _MetatlasObject(HasTraits):
    name = CUnicode()
    unique_id = CInt()

    def _unique_id_default(self):
        return random.randint(0, 1e6)


class Compound(_MetatlasObject):

    """
    name: str
      Name of the compound
    formula : str
        Chemical formula of the intact, neutral form of the molecule.
    adducts : str
        Adduct accounting form the ionization and potential in source degradation.
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
    neutral_mass : float
        The intact, neutral mass of the molecule
    pubchem_id : str
        Pubchem ID for compound.
    """
    formula = CUnicode()
    adducts = CUnicode()
    mz = CFloat()
    mz_threshold = CFloat()
    rt_min = CFloat()
    rt_max = CFloat()
    rt_peak = CFloat()
    neutral_mass = CFloat()
    pubchem_id = CUnicode()


class Atlas(_MetatlasObject):

    """
    name : str
        Common name of the atlas.
    compounds : list of strings
        Names of compounds (names of kbase workspace objects).
    sample : str
        Descri 84t4ption of sample.
    method : str
        Description of method.
    """
    name = CUnicode()
    compounds = List(Instance(Compound))
    sample = CUnicode()
    method = CUnicode()

    def edit(self):
        import qgrid
        qgrid.nbinstall(overwrite=False)
        from IPython.html.widgets import Button, HBox
        from IPython.display import display

        # create a visualization for the dataframe
        if not self.compounds:
            self.compounds.append(Compound())
        data = [c._trait_values for c in self.compounds]
        dataframe = pd.DataFrame(data)
        cols = ['name', 'formula', 'adducts', 'mz', 'mz_threshold',
                'rt_min', 'rt_max', 'rt_peak', 'neutral_mass',
                'pubchem_id']
        grid = qgrid.QGridWidget(df=dataframe[cols], remote_js=True)

        add_row = Button(description="Add Row")
        add_row.on_click(grid.add_row)

        rem_row = Button(description="Remove Row")
        rem_row.on_click(grid.remove_row)

        display(HBox((add_row, rem_row)), grid)

        grid.on_msg(self._incoming_msg)

    def _incoming_msg(self, widget, msg):
        if msg['type'] == 'add_row':
            self.compounds.append(Compound())
        elif msg['type'] == 'remove_row':
            self.compounds.pop(msg['row'])
        elif msg['type'] == 'cell_change':
            c = self.compounds[msg['row']]
            setattr(c, msg['column'], msg['value'])


class FileSpec(_MetatlasObject):

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


class FileInfo(FileSpec):

    """
    mzml_file : str
        Path to input_file.
    hdf_file : str
        Name of the hdf_file.
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

    mzml_file = CUnicode()
    hdf_file = CUnicode()

    def parse(self):
        """Parse a file info spec"""
        self.hdf_file = mzml_to_hdf(self.mzml_file, self.hdf_file or None)


class Experiment(_MetatlasObject):

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
    finfos = List(Instance(FileInfo))
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
        keys = filespec.trait_names()
        state = {}
        for key in keys:
            state[key] = getattr(filespec, key)
        for fname in files:
            finfo = FileInfo(mzml_file=fname, **state)
            finfo.parse()
            self.finfos.append(finfo)

    def edit(self):
        import qgrid
        qgrid.nbinstall(overwrite=False)
        from IPython.html.widgets import Button
        from IPython.display import display

        # create a visualization for the dataframe
        if not self.finfos:
            raise TypeError('Please load files first')
        data = [f._trait_values for f in self.finfos]
        dataframe = pd.DataFrame(data)
        cols = ['name', 'polarity', 'group', 'inclusion_order',
                'normalization_factor', 'retention_correction']

        grid = qgrid.QGridWidget(df=dataframe[cols], remote_js=True)

        rem_row = Button(description="Remove Finfo")
        rem_row.on_click(grid.remove_row)

        display(rem_row, grid)

        grid.on_msg(self._incoming_msg)

    def _incoming_msg(self, widget, msg):
        if msg['type'] == 'remove_row':
            self.finfos.pop(msg['row'])
        elif msg['type'] == 'cell_change':
            f = self.finfos[msg['row']]
            setattr(f, msg['column'], msg['value'])


def list_experiments(username=None):
    """Get a list of available experiments.

    Parameters
    ----------
    username: str, optional
        Username of original experimenter.

    Returns
    -------
    out : Experiment
        Experiment object.
    """
    if username:
        raise NotImplemented
    experiments = _WORKSPACE.load_all()
    return [e.name for e in experiments if isinstance(e, Experiment)]


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

    def save(self, obj):
        self.db[obj.name].insert(dict(data=pickle.dumps(obj)))

    def load(self, name):
        try:
            objects = self.db[name].all()
        except Exception:
            errmsg = 'Object called "%s" not found in database"' % name
            raise ValueError(errmsg)
        if objects:
            obj = list(objects)[-1]
            return pickle.loads(obj['data'])

    def load_all(self):
        objs = []
        for t in self.db.tables:
            try:
                objs.append(self.load(t))
            except Exception:
                pass
        return objs

# Singleton Workspace object
_WORKSPACE = _Workspace()
