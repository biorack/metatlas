import getpass
import uuid
import time
import os
import json
import pprint
import tables
import numpy as np
from pwd import getpwuid
from tabulate import tabulate
import pandas as pd

from .utils import (
    set_docstring, Workspace, format_timestamp, MetList,
    MetUnicode, MetFloat, MetInstance, MetInt, MetEnum,
    edit_traits, HasTraits, CBool, POLARITY, Stub, List
)


# Whether to fetch stubs automatically, disabled when we want to display
# a large number of objects.
FETCH_STUBS = True


def retrieve(object_type, **kwargs):
    """Get objects from the Metatlas object database.

    This will automatically select only objects created by the current
    user unless `username` is provided. Use `username='*'` to search
    against all users.

    Parameters
    ----------
    object_type: string
      The type of object to search for (i.e. "Groups").
    **kwargs
      Specific search queries (i.e. name="Sargasso").
      Use '%' for wildcard patterns (i.e. description='Hello%').
      If you want to match a '%' character, use '%%'.

    Returns
    -------
    objects: list
      List of Metatlas Objects meeting the criteria.  Will return the
      latest version of each object.
    """
    return WORKSPACE.retrieve(object_type, **kwargs)


def remove(object_type, **kwargs):
    """Remove objects from the Metatlas object database.

    Parameters
    ----------
    object_type: string
      The type of object to remove (i.e. "Groups").
    **kwargs
      Specific search queries (i.e. name="Sargasso").
      Use '%' for wildcard patterns (i.e. description='Hello%').
      If you want to match a '%' character, use '%%'.
    """
    return WORKSPACE.remove(object_type, **kwargs)


def remove_objects(objects, all_versions=True, **kwargs):
    """Remove objects from the database.

    Parameters
    ----------
    all_versions: boolean, optional
        If True, remove all versions of the object sharing the current
        head_id.
    """
    return WORKSPACE.remove_objects(objects, all_versions, **kwargs)


def store(objects, **kwargs):
    """Store Metatlas objects in the database.

    Parameters
    ----------
    objects: Metatlas object or list of Metatlas Objects
        Object(s) to store in the database.
    """
    WORKSPACE.save_objects(objects, **kwargs)


@set_docstring
class MetatlasObject(HasTraits):

    name = MetUnicode('Untitled', help='Name of the object')
    description = MetUnicode('No description',
                             help='Description of the object')
    unique_id = MetUnicode(help='Unique identifier for the object',
                           readonly=True)
    creation_time = MetInt(help='Unix timestamp at object creation',
                           readonly=True)
    username = MetUnicode(help='Username who created the object',
                          readonly=True)
    last_modified = MetInt(help='Unix timestamp at last object update',
                           readonly=True)
    prev_uid = MetUnicode(help='Unique id of previous version', readonly=True)
    head_id = MetUnicode(help='Unique id of most recent version of this object', readonly=True)
    _loopback_guard = CBool(False, readonly=True)
    _changed = CBool(False, readonly=True)

    def __init__(self, **kwargs):
        """Set the default attributes."""
        kwargs.setdefault('unique_id', uuid.uuid4().hex)
        kwargs.setdefault('head_id', kwargs['unique_id'])
        kwargs.setdefault('username', getpass.getuser())
        kwargs.setdefault('creation_time', int(time.time()))
        kwargs.setdefault('last_modified', int(time.time()))
        super(MetatlasObject, self).__init__(**kwargs)
        self._changed = True
        self.on_trait_change(self._on_update)

    def _update(self, override_user=False):
        """Store the object in the workspace, including child objects.

        Child objects are stored in their own tables, and Lists of
        Child objects are captured in link tables.
        """
        self._loopback_guard = True
        # see if we need to create a new object
        if not override_user and self.username != getpass.getuser():
            self._changed = True
            self.prev_uid = self.unique_id
            self.unique_id = uuid.uuid4().hex
            self.head_id = self.unique_id
            self.username = getpass.getuser()
            self.last_modified = time.time()
            return True, None
        else:
            changed, prev_uid = self._changed, self.prev_uid
            if changed:
                self.last_modified = time.time()
                self.prev_uid = uuid.uuid4().hex
        self._changed = False
        self._loopback_guard = False
        return changed, prev_uid

    def clone(self, recursive=False):
        """Create a new version of this object.

        Parameters
        ----------
        recursive: boolean, optional
            If true, clone all of the descendant objects as well.

        Returns
        -------
        obj: MetatlasObject
            Cloned object.
        """
        obj = self.__class__()
        for (tname, trait) in self.traits().items():
            if tname.startswith('_') or trait.get_metadata('readonly'):
                continue
            val = getattr(self, tname)
            if recursive and isinstance(trait, List):
                val = [v.clone(True) for v in val]
            elif recursive and isinstance(trait, MetInstance) and val:
                val = val.clone(True)
            setattr(obj, tname, val)
        obj.prev_uid = self.unique_id
        obj.head_id = self.unique_id
        obj.unique_id = uuid.uuid4().hex
        return obj

    def edit(self):
        """Create an editor for the object

        Creates labels and editors for all traits.

        If the trait is marked readyonly or is an instance of another Trait,
        it will be disabled.
        If the trait is an Enum, there will be a dropdown.
        Otherwise, there will be a text editor.
        """
        edit_traits(self)

    def show_diff(self, unique_id=None):
        """Show a diff of what has changed between this and previous version.

        Parameters
        ----------
        unique_id: optional, string
            Unique id to compare against (defaults to current entry in db).
        """
        if unique_id is None:
            unique_id = self.unique_id
        obj = retrieve(self.__class__.__name__, unique_id=unique_id)
        if len(obj) != 1:
            print('No change!')
            return
        obj = obj[0]
        msg = []
        for (tname, trait) in self.traits().items():
            if tname.startswith('_') or trait.get_metadata('readonly'):
                continue
            val = getattr(self, tname)
            other = getattr(obj, tname)
            if isinstance(trait, MetInstance):
                if val.unique_id != other.unique_id:
                    msg.append((tname, other.unique_id, val.unique_id))
            elif isinstance(trait, MetList):
                if not len(val) == len(other):
                    msg.append((tname, '%s items' % len(other),
                                '%s items' % len(val)))
                else:
                    for (v, o) in zip(val, other):
                        if not v.unique_id == o.unique_id:
                            msg.append((tname, 'objects changed', ''))
                            break
            elif val != other:
                msg.append((tname, str(other), str(val)))
        print(tabulate(msg))

    def _on_update(self, name):
        """When the model changes, set the update fields.
        """
        if self._loopback_guard or name.startswith('_'):
            return
        self._changed = True

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        names = sorted(self.trait_names())
        names.remove('name')
        names = ['name'] + [n for n in names if not n.startswith('_')]
        state = dict([(n, getattr(self, n)) for n in names])
        state['creation_time'] = format_timestamp(self.creation_time)
        state['last_modified'] = format_timestamp(self.last_modified)
        return pprint.pformat(state)

    def __getattribute__(self, name):
        """Automatically resolve stubs on demand.
        """
        value = super(MetatlasObject, self).__getattribute__(name)
        if isinstance(value, Stub) and FETCH_STUBS:
            value = value.retrieve()
            setattr(self, name, value)
        elif isinstance(value, list) and value and FETCH_STUBS:
            new = []
            changed = False
            for subvalue in value:
                if isinstance(subvalue, Stub):
                    new.append(subvalue.retrieve())
                    changed = True
                else:
                    new.append(subvalue)
            if changed:
                setattr(self, name, new)
                value = new
        return value


@set_docstring
class Method(MetatlasObject):
    """
    For each LCMS run, a Method is a consistent description of
    how the sample was prepared and LCMS data was collected.
    """
    protocol_ref = MetUnicode(help='Reference to a published protocol: ' +
        'identical to the protocol used.')
    quenching_method = MetUnicode(help='Description of the method used to ' +
        'stop metabolism.')
    extraction_solvent = MetUnicode(help='Solvent or solvent mixture used to ' +
        'extract metabolites.')
    reconstitution_method = MetUnicode(help='Solvent or solvent mixture ' +
        'the extract is reconstituted in prior to injection for an LCMS Run.')
    mobile_phase_a = MetUnicode(help='Solvent or solvent mixture.')
    mobile_phase_b = MetUnicode(help='Solvent or solvent mixture.')
    temporal_parameters = MetUnicode('List of temporal changes to the' +
        'mixing of mobile_phase_a and mobile_phase_b.')
    column_model = MetUnicode(help='Brand and catalog number of the column')
    column_type = MetUnicode(help='Class of column used.')
    scan_mz_range = MetUnicode(help='Minimum and ' +
        'maximum mz recorded for a run.')
    instrument = MetUnicode(help='Brand and catalog number for the ' +
        'mass spectrometer.')
    ion_source = MetUnicode(help='Method for ionization.')
    mass_analyzer = MetUnicode(help='Method for detection.')
    polarity = MetEnum(POLARITY, 'positive', help='polarity for the run')


@set_docstring
class Sample(MetatlasObject):
    """A Sample is the material that is processed with a Method."""
    pass


def load_lcms_files(mzml_files):
    """Parse mzML files and load them into LcmsRun objects.

    Note: This should be done automatically for runs in
    /project/projectdirs/metatlas/raw_data/<username>

    Parameters
    ----------
    mzml_files: list of str
       List of paths to mzml_files.

    Returns
    -------
    runs: list
       List of LcmsRun objects.
    """
    runs = []
    for fname in mzml_files:
        hdf5_file = fname.replace('.mzML', '.h5')
        if os.path.exists(hdf5_file):
            print('File already exists: %s' % hdf5_file)
            continue
        try:
            from metatlas import LcmsRun, mzml_to_hdf
            hdf_file = mzml_to_hdf(fname)
            user = getpwuid(os.stat(fname).st_uid).pw_name
            filename = os.path.splitext(os.path.basename(fname))[0]
            dirname = os.path.dirname(fname)
            experiment = os.path.basename(dirname)
            description = experiment + ' ' + filename
            ctime = os.stat(fname).st_ctime
            run = LcmsRun(name=filename, description=description,
                          created_by=user,
                          modified_by=user,
                          created=ctime, last_modified=ctime,
                          mzml_file=fname, hdf5_file=hdf_file)
            runs.append(run)
        except Exception as e:
            print(e)
    store(runs)
    return runs


@set_docstring
class LcmsRun(MetatlasObject):
    """An LCMS run is the reference to a file prepared with liquid
    chromatography and mass spectrometry.

    The msconvert program is used to convert raw data (centroided is prefered)
    to mzML.

    Note: These objects are not intented to be created directly, but by putting
    the files in /project/projectdirs/metatlas/raw_data/<username> or by
    running `load_lcms_files()`.
    """
    method = MetInstance(Method)
    hdf5_file = MetUnicode(help='Path to the HDF5 file at NERSC')
    mzml_file = MetUnicode(help='Path to the MZML file at NERSC')
    sample = MetInstance(Sample)

    def interact(self, min_mz=None, max_mz=None, polarity=None, ms_level=1):
        """Interact with LCMS data - XIC linked to a Spectrogram plot.

        Parameters
        ----------
        min_mz: float
            Minimum m/z (defaults to file min)
        max_mz: float
            Maximum m/z (defaults to file max)
        polarity: {0, 1}
            Polarity (defaults to neg. if present in file, else pos.)
        ms_level: {0, 1}
            The ms level.
        """
        import matplotlib.pyplot as plt
        from metatlas import get_chromatogram, get_spectrogram, get_info
        fid = tables.open_file(self.hdf5_file)

        info = get_info(fid)
        if polarity is None:
            if info['ms%s_neg' % ms_level]['nrows']:
                polarity = 0
            else:
                polarity = 1
        if polarity == 0:
            table_name = 'ms%s_neg' % ms_level
        else:
            table_name = 'ms%s_pos' % ms_level
        if min_mz is None:
            min_mz = info[table_name]['min_mz']
        if max_mz is None:
            max_mz = info[table_name]['max_mz']

        rt, irt = get_chromatogram(fid, min_mz, max_mz, 1, polarity)
        mz, imz = get_spectrogram(fid, rt[0], rt[1], 1, polarity)

        fig, (ax1, ax2) = plt.subplots(ncols=2)
        ax1.plot(rt, irt)
        ax1.set_title('XIC: %0.1f - %0.1f m/z' % (min_mz, max_mz))
        ax1.set_xlabel('Time (min)')
        ax1.set_ylabel('Intensity')
        ax1._vline = ax1.axvline(rt[0])

        ax2.vlines(mz, 0, imz)
        ax2.set_xlabel('Mass (m/z)')
        ax2.set_title('MS%s Spectrogram at %0.1f min' % (ms_level, rt.min()))

        def callback(event):
            if event.inaxes == ax1:
                rt_event = event.xdata
                # get the closest actual RT
                idx = (np.abs(rt - rt_event)).argmin()
                mz, imz = get_spectrogram(fid, rt[idx], rt[idx], 1, polarity)

                ax1._vline.remove()
                ax1._vline = ax1.axvline(rt_event, color='k')

                ax2.clear()
                ax2.vlines(mz, 0, imz)
                ax2.set_xlabel('Mass (m/z)')
                ax2.set_title('Spectrogram at %0.1f min' % rt_event)
                fig.canvas.draw()
        fig.canvas.mpl_connect('button_press_event', callback)
        fig.canvas.draw()


@set_docstring
class ReferenceDatabase(MetatlasObject):
    """External databases (PubChem, KEGG, MetaCyc, KBase, etc)."""
    enabled = CBool(True)


@set_docstring
class FunctionalSet(MetatlasObject):
    """Functional sets of compounds.  For example set called "hexose" would
    include "glucose, galactose, etc".  Functional sets can be sets-of-sets.
    "Sugars" would be a set that contains "Hexoses".
    """
    enabled = CBool(True)
    members = MetList(MetInstance(MetatlasObject))


@set_docstring
class Compound(MetatlasObject):
    """A Compound as a name and description to capture incomplete IDs:
    for example "hexose", "32:0 diacylglycerol".
    For IDs that have high enough confidence for structural assignments an
    InChi string is the ID.
    """
    inchi = MetUnicode(help='IUPAC International Chemical Identifier, optional')
    formula = MetUnicode()
    mono_isotopic_molecular_weight = MetFloat()
    synonyms = MetUnicode()
    url = MetUnicode(help='Reference database table url')
    permanent_charge = MetInt()
    inchi_key = MetUnicode()
    number_components = MetInt(help='Must be one or greater')
    neutralized_inchi = MetUnicode()
    neutralized_inchi_key = MetUnicode()
    neutralized_2d_inchi = MetUnicode()
    neutralized_2d_inchi_key = MetUnicode()
    reference_xrefs = MetList(MetInstance(ReferenceDatabase),
                           help='Tag a compound with compound ids from ' +
                                'external databases')
    functional_sets = MetList(MetInstance(FunctionalSet))


@set_docstring
class Reference(MetatlasObject):
    """Place holder for future reference sources.
    We expect many in silico methods will soon be robust enough to suggest
    retention times, m/z, and fragmentation.
    MIDAS is a great example of this.
    """
    lcms_run = MetInstance(LcmsRun)
    enabled = CBool(True)
    ref_type = MetUnicode(help='The type of reference')


@set_docstring
class IdentificationGrade(MetatlasObject):
    """
    Each CompoundIdentification will have an identification_grade
    Identification Grades:
    1) High intensity and verifiable by MSMS and RT
    authentic standard
    2) Verifiable by MSMS from database or
    publication
    3) Has fragment ion or neutral loss characteristic
    of a class of compounds
    4) definitive chemical formula and adduct
    5) Significant changing metabolite with MSMS
    suggestion from MIDAS
    6) Significant changing metabolite
    7) Not Significant changing metabolite with
    MSMS suggestion from MIDAS
    8) Not Significant changing metabolite
    """
    pass


ID_GRADES = dict()


class _IdGradeTrait(MetInstance):

    klass = IdentificationGrade

    def validate(self, obj, value):
        global ID_GRADES
        if not value:
            return
        if isinstance(value, self.klass):
            return value
        elif isinstance(value, str):
            if value.upper() in ID_GRADES:
                return ID_GRADES[value.upper()]
            objects = WORKSPACE.retrieve('identificationgrade', name=value.upper())
            if objects:
                ID_GRADES[value.upper()] = objects[-1]
                return objects[-1]
            else:
                self.error(obj, value)
        else:
            self.error(obj, value)


@set_docstring
class Group(MetatlasObject):
    """A Group can be a:
    file,
    group of files, or
    group of groups
    """
    items = MetList(MetInstance(MetatlasObject),
                 help='Can contain other groups or LCMS Runs')


@set_docstring
class MzIntensityPair(MetatlasObject):
    mz = MetFloat()
    intensity = MetFloat()


@set_docstring
class FragmentationReference(Reference):
    polarity = MetEnum(POLARITY, 'positive')
    precursor_mz = MetFloat()
    mz_intensities = MetList(MetInstance(MzIntensityPair),
                          help='list of [mz, intesity] tuples that describe ' +
                               ' a fragmentation spectra')


@set_docstring
class RtReference(Reference):
    rt_peak = MetFloat()
    rt_min = MetFloat()
    rt_max = MetFloat()
    rt_units = MetEnum(('sec', 'min'), 'sec')


@set_docstring
class MzReference(Reference):
    """Source of the assertion that a compound has a given m/z and
    other properties directly tied to m/z.
    """
    mz = MetFloat()
    mz_tolerance = MetFloat()
    mz_tolerance_units = MetEnum(('ppm', 'Da'), 'ppm')
    detected_polarity = MetEnum(POLARITY, 'positive')
    adduct = MetUnicode(help='Optional adduct')
    modification = MetUnicode(help='Optional modification')
    observed_formula = MetUnicode(help='Optional observed formula')


@set_docstring
class CompoundIdentification(MetatlasObject):
    """A CompoundIdentification links multiple sources of evidence about a
    compound's identity to an Atlas."""
    compound = MetList(MetInstance(Compound))
    identification_grade = _IdGradeTrait(
        help='Identification grade of the id (can be specified by a letter A-H'
    )
    mz_references = MetList(MetInstance(MzReference))
    rt_references = MetList(MetInstance(RtReference))
    frag_references = MetList(MetInstance(FragmentationReference))


@set_docstring
class Atlas(MetatlasObject):
    """An atlas contains many compound_ids."""
    compound_identifications = MetList(
        MetInstance(CompoundIdentification),
        help='List of Compound Identification objects')


def find_invalid_runs(**kwargs):
    """Find invalid runs.
    """
    override = kwargs.pop('_override', False)
    if not override:
        kwargs.setdefault('username', getpass.getuser())
    all_runs = retrieve('lcmsruns', **kwargs)
    invalid = []
    for run in all_runs:
        if (not os.path.exists(run.hdf5_file) or
            not os.path.exists(run.mzml_file) or
            not os.stat(run.hdf5_file).st_size):
                invalid.append(run)
    return invalid


# Singleton Workspace object
# Must be instantiated after all of the Metatlas Objects
# are defined so we can get all of the subclasses.
WORKSPACE = Workspace()


def _create_qgrid(objects):
    """Create a qgrid from a list of metatlas objects.
    """
    global FETCH_STUBS
    import qgrid
    qgrid.nbinstall(overwrite=False)

    # we want to handle dates, enums, and use ids for objects
    FETCH_STUBS = False
    objs = [o._trait_values.copy() for o in objects if o.__class__ == objects[0].__class__]
    FETCH_STUBS = True
    enums = []
    cols = []
    # remove lists and use strings for objects
    for (tname, trait) in objects[0].traits().items():
        if tname.startswith('_') or tname in ['head_id', 'prev_uid']:
            continue
        cols.append(tname)
        if isinstance(trait, MetList):
            [o.__setitem__(tname, str([i.unique_id for i in o[tname]]))
             for o in objs]
        if isinstance(trait, MetInstance):
            [o.__setitem__(tname, getattr(o[tname], 'unique_id', 'None')) for o in objs]
        if isinstance(trait, MetEnum):
            enums.append(tname)

    dataframe = pd.DataFrame(objs)[sorted(cols)]
    for col in enums:
        dataframe[col] = dataframe[col].astype('category')
    for col in ['last_modified', 'creation_time']:
        dataframe[col] = pd.to_datetime(dataframe[col], unit='s')

    options = qgrid.grid.defaults.grid_options
    grid = qgrid.grid.QGridWidget(precision=6,
                       grid_options=json.dumps(options),
                       remote_js=True)
    grid.df = dataframe

    def handle_msg(widget, content, buffers=None):
        if content['type'] == 'cell_change':
            obj = objects[content['row']]
            try:
                setattr(obj, content['column'], content['value'])
            except Exception:
                pass

    grid.on_msg(handle_msg)
    return grid


def edit_objects(objects):
    """Edit a set of metatlas object in a QGrid.
    """
    from IPython.display import display
    grid = _create_qgrid(objects)
    display(grid)


if __name__ == '__main__':
    m1 = Group(name='spam')
    store(m1)
    m1.description = 'baz'
    store(m1)
    print(retrieve('group', name='spam'))
