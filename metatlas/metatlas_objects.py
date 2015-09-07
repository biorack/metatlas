import getpass
import uuid
import time
import os
import pickle
import functools

try:
    from traitlets import (
        HasTraits, CUnicode, List, CInt, Instance, Enum,
        CFloat, TraitError, CBool)
except ImportError:
    from IPython.utils.traitlets import (
        HasTraits, CUnicode, List, CInt, Instance, Enum,
        CFloat, TraitError, CBool)
from IPython.display import display
try:
    from ipywidgets import Text, Dropdown, HBox, VBox
except ImportError:
    from IPython.html.widgets import Text, Dropdown, HBox, VBox
import dataset

from metatlas import mzml_to_hdf

POLARITY = ('positive', 'negative', 'alternating')
NERSC_WORKSPACE = '/project/projectdirs/metatlas/workspace'


class _Workspace(object):

    def __init__(self):
        # allow for fallback when not on NERSC
        #if os.path.exists(NERSC_WORKSPACE):
        #    path = 'mysql+mysqldb://%s/
        #else:
        self.path = getpass.getuser() + '_workspace.db'
        self._db = None
        # handle circular references
        self.seen = dict()

    @property
    def db(self):
        if self._db:
            return self._db
        self._db = dataset.connect('sqlite:///%s' % self.path)
        os.chmod(self.path, 0o775)

    def insert(self, name, state):
        if name not in self.db:
            self.db.create_table(name, primary_id='unique_id',
                                 primary_type='String(32)')
        self.db[name].insert(state)

    def find_one(self, name, **kwargs):
        return self.db[name].find_one(**kwargs)

    def find(self, name, **kwargs):
        return self.db[name].find(**kwargs)

    def save(self, obj, name=None):
        if obj.unique_id in self.seen:
            return
        name = name or obj.__class__.__name__
        if self.find_one(name, unique_id=obj.unique_id):
            return
        self.seen[obj.unique_id] = ''
        state = obj.traits().copy()
        for (tname, trait) in obj.traits().items():
            if tname.startswith('_'):
                del state[tname]
                continue
            value = getattr(obj, tname)
            if isinstance(trait, List) and value:
                # if it is not composed of other objects, just pickle it
                if not isinstance(value[0], MetatlasObject):
                    state[tname] = pickle.dumps(value)
                    continue
                # handle a list of objects by using a Link table
                # create the link table if necessary
                table_name = '_'.join([name, tname])
                # create an entry in the table for each item
                # store the item in its own table
                for subvalue in value:
                    subvalue.store()
                    link = dict(unique_id=obj.unique_id,
                                target_id=subvalue.unique_id,
                                target_table=subvalue.__class__.__name__)
                    if not self.find_one(table_name, unique_id=obj.unique_id,
                                         target_id=subvalue.unique_id):
                        self.insert(table_name, link)
                # do not store this entry in our own table
                del state[tname]
            elif isinstance(trait, Instance):
                # handle a sub-object
                # if it is not assigned, use and empty unique_id
                if value is None:
                    state[tname] = ''
                # otherwise, store the uid and allow the object to store
                # itself
                else:
                    state[tname] = value.unique_id
                    value.store()
            else:
                # store the raw value in this table
                state[tname] = value
        self.insert(name, state)
        del self.seen[obj.unique_id]

    def load(self, obj, name=None):
        if obj.unique_id in self.seen:
            return self.seen[obj.unique_id]
        name = name or obj.__class__.__name__
        # get our table entry
        if name not in self.db:
            return
        entry = self.find_one(name, unique_id=obj.unique_id)
        if not entry:
            return
        state = obj.traits()
        self.seen[obj.unique_id] = obj
        for (tname, trait) in state.items():
            if tname.startswith('_'):
                continue
            if isinstance(trait, List):
                # check for a pickled object
                if tname in entry:
                    setattr(obj, tname, entry[tname])
                    continue
                # get the items using the link table
                # make the items retrieve themselves by unique id
                items = []
                table_name = '_'.join([name, tname])
                rows = self.find(table_name, unique_id=obj.unique_id)
                for row in rows:
                    klass = eval(row['target_table'])
                    target = klass(unique_id=row['target_id'])
                    target = target.retrieve()
                    if target:
                        items.append(target)
                setattr(obj, tname, items)
            elif isinstance(trait, Instance):
                # handle a sub-object
                # if it was not assigned, use None
                if not entry[tname]:
                    setattr(obj, tname, None)
                # allow the object to retrieve itself by unique id
                else:
                    subobject = trait.klass(unique_id=entry[tname])
                    setattr(obj, tname, subobject.retrieve())
            else:
                # set the attribute directly
                setattr(obj, tname, entry[tname])
        del self.seen[obj.unique_id]
        # return object for parent usage
        return obj

# Singleton Workspace object
workspace = _Workspace()


def set_docstring(cls):
    """Set the docstring for a MetatlasObject object"""
    doc = cls.__doc__
    if not doc:
        doc = cls.__name__ + ' object.'
    doc += '\n\nParameters\n----------\n'
    for (tname, trait) in sorted(cls.class_traits().items()):
        if tname.startswith('_'):
            continue
        descr = trait.__class__.__name__.lower()
        if descr.startswith('c'):
            descr = descr[1:]
        elif descr == 'enum':
            descr = '{' + ', '.join(trait.values) + '}'
        doc += '%s: %s\n' % (tname, descr)
        help_text = trait.get_metadata('help')
        if not help_text:
            help_text = '%s value.' % tname
        help_text = help_text.strip()
        if help_text.endswith('.'):
            help_text = help_text[:-1]
        if trait.get_metadata('readonly'):
            help_text += ' (read only)'
        help_text += '.'
        doc += '    %s\n' % help_text
    cls.__doc__ = doc
    return cls


@set_docstring
class MetatlasObject(HasTraits):

    name = CUnicode('Untitled', help='Name of the object')
    description = CUnicode('No description', help='Description of the object')
    unique_id = CUnicode(help='Unique identifier for the object',
                         readonly=True)
    created = CInt(help='Unix timestamp at object creation',
                   readonly=True)
    created_by = CUnicode(help='User who created the object',
                          readonly=True)
    last_modified = CInt(help='Unix timestamp at last object update',
                         readonly=True)
    modified_by = CUnicode(help='User who last updated the object',
                           readonly=True)
    prev_unique_id = CUnicode(help='Unique id of previous version',
                              readonly=True)
    _loopback_guard = CBool(False, readonly=True)
    _changed = CBool(False, readonly=True)

    def __init__(self, **kwargs):
        """Set the default attributes."""
        kwargs.setdefault('unique_id', uuid.uuid4().hex)
        kwargs.setdefault('created_by', getpass.getuser())
        kwargs.setdefault('modified_by', getpass.getuser())
        kwargs.setdefault('created', int(time.time()))
        kwargs.setdefault('last_modified', int(time.time()))
        super(MetatlasObject, self).__init__(**kwargs)
        self.on_trait_change(self._on_update)

    def store(self, name=None):
        """Store the object in the workspace, including child objects.

        Child objects are stored in their own tables, and Lists of
        Child objects are captured in link tables.
        """
        self._loopback_guard = True
        if self._changed:
            self.prev_unique_id = self.unique_id
            self.unique_id = uuid.uuid4().hex
            self.last_modified = time.time()
            self.modified_by = getpass.getuser()
        workspace.save(self, name)
        self._changed = False
        self._loopback_guard = False

    def retrieve(self, name=None):
        """Retrieve the object from the workspace, including child objects.
        """
        self._loopback_guard = True
        value = workspace.load(self, name)
        self._changed = False
        self._loopback_guard = False
        return value

    def edit(self):
        """Create an editor for the object

        Creates labels and editors for all traits.

        If the trait is marked readyonly or is an instance of another Trait,
        it will be disabled.
        If the trait is an Enum, there will be a dropdown.
        Otherwise, there will be a text editor.
        """
        edit_traits(self)

    def _on_update(self, name):
        """When the model changes, set the update fields.
        """
        if self._loopback_guard or name.startswith('_'):
            return
        print('changed', name)
        self._changed = True

    def __str__(self):
        return self.name + ' (%s)' % self.unique_id


@set_docstring
class Method(MetatlasObject):
    """
    For each LCMS run, a Method is a consistent description of
    how the sample was prepared and LCMS data was collected.
    """
    protocol_ref = CUnicode(help='Reference to a published protocol: ' +
        'identical to the protocol used.')
    quenching_method = CUnicode(help='Description of the method used to ' +
        'stop metabolism.')
    extraction_solvent = CUnicode(help='Solvent or solvent mixture used to ' +
        'extract metabolites.')
    reconstitution_method = CUnicode(help='Solvent or solvent mixture ' +
        'the extract is reconstituted in prior to injection for an LCMS Run.')
    mobile_phase_a = CUnicode(help='Solvent or solvent mixture.')
    mobile_phase_b = CUnicode(help='Solvent or solvent mixture.')
    temporal_parameters = CUnicode('List of temporal changes to the' +
        'mixing of mobile_phase_a and mobile_phase_b.')
    column_model = CUnicode(help='Brand and catalog number of the column')
    column_type = CUnicode(help='Class of column used.')
    scan_mz_range = CUnicode(help='Minimum and ' +
        'maximum mz recorded for a run.')
    instrument = CUnicode(help='Brand and catalog number for the ' +
        'mass spectrometer.')
    ion_source = CUnicode(help='Method for ionization.')
    mass_analyzer = CUnicode(help='Method for detection.')
    polarity = Enum(POLARITY, 'positive', help='polarity for the run')


@set_docstring
class Sample(MetatlasObject):
    """A Sample is the material that is processed with a Method."""
    pass


@set_docstring
class LcmsRun(MetatlasObject):
    """An LCMS run is the reference to a file prepared with liquid 
    chromatography and mass spectrometry.

    The msconvert program is used to convert raw data (centroided is prefered)
    to mzML.  The mzML file can be converted to HDF5 by the metatlas method below (parse).
    """
    method = Instance(Method, allow_none=True)
    hdf5_file = CUnicode(help='Path to the HDF5 file at NERSC')
    mzml_file = CUnicode(help='Path to the MZML file at NERSC')
    sample = Instance(Sample, allow_none=True)

    def _name_default(self):
        return self.mzml_file.replace('.mzml', '')

    def parse(self):
        """Parse a file info spec"""
        self.hdf_file = mzml_to_hdf(self.mzml_file, self.hdf_file or None)


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
    members = List(Instance(MetatlasObject))


@set_docstring
class Compound(MetatlasObject):
    """A Compound as a name and description to capture incomplete IDs:
    for example "hexose", "32:0 diacylglycerol".
    For IDs that have high enough confidence for structural assignments an
    InChi string is the ID.
    """
    InChl = CUnicode(help='IUPAC International Chemical Identifier, optional')
    reference_xrefs = List(Instance(ReferenceDatabase),
                           help='Tag a compound with compound ids from ' +
                                'external databases')
    functional_sets = List(Instance(FunctionalSet))


@set_docstring
class Reference(MetatlasObject):
    """Place holder for future reference sources.
    We expect many in silico methods will soon be robust enough to suggest
    retention times, m/z, and fragmentation.
    MIDAS is a great example of this.
    """
    lcms_run = Instance(LcmsRun, allow_none=True)
    enabled = CBool(True)


@set_docstring
class IdentificationGrade(MetatlasObject):
    """
    Each CompoundId will have an identification_grade
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


@set_docstring
class CompoundId(MetatlasObject):
    """A CompoundId links multiple sources of evidence about a compound's 
    identity to an Atlas."""
    compound = Instance(Compound, allow_none=True)
    identification_grade = Instance(IdentificationGrade, allow_none=True)
    references = List(Instance(Reference))


@set_docstring
class Atlas(MetatlasObject):
    """An atlas contains many compound_ids."""
    compounds_ids = List(Instance(CompoundId),
                         help='List of Compound Identification objects')


@set_docstring
class Group(MetatlasObject):
    """A Group can be a:
    file,
    group of files, or
    group of groups
    """
    items = List(Instance(MetatlasObject),
                 help='Can contain other groups or LCMS Runs')


@set_docstring
class MzIntensityPair(MetatlasObject):
    mz = CFloat()
    intensity = CFloat()


@set_docstring
class FragmentationReference(Reference):

    polarity = Enum(POLARITY, 'positive')
    precursor_mz = CFloat()
    mz_intensities = List(Instance(MzIntensityPair),
                          help='list of [mz, intesity] tuples that describe ' +
                               ' a fragmentation spectra')


@set_docstring
class RtReference(Reference):

    RTpeak = CFloat()
    RTmin = CFloat()
    RTmax = CFloat()
    RTUnits = Enum(('sec', 'min'), 'sec')


@set_docstring
class MzReference(Reference):
    """Source of the assertion that a compound has a given m/z and
    other properties directly tied to m/z.
    """ 
    mz = CFloat()
    mz_tolerance = CFloat()
    mz_tolerance_units = Enum(('ppm', 'Da'), 'ppm')
    detected_polarity = Enum(POLARITY, 'positive')
    adduct = CUnicode(help='Optional adduct')
    modification = CUnicode(help='Optional modification')
    observed_formula = CUnicode(help='Optional observed formula')


def edit_traits(obj):
    """Create an IPython widget editor for a Traits object"""
    names = sorted(obj.trait_names())
    names.remove('name')
    names = ['name'] + names
    items = [Text('', disabled=True)]
    for name in names:
        if name.startswith('_'):
            continue
        try:
            value = getattr(obj, name)
        except TraitError:
            value = None
        trait = obj.traits()[name]
        if (trait.get_metadata('readonly') or
                isinstance(trait, Instance) or value is None):
            items.append(Text(str(value), disabled=True))

        elif isinstance(trait, Enum):
            # create a closure around "name" for the on_trait_change
            # callback
            def create_dropdown(name):
                dd = Dropdown(value=value, options=trait.values)

                def callback(dummy, value):
                    setattr(obj, name, value)
                dd.on_trait_change(callback, 'value')
                items.append(dd)

            create_dropdown(name)
        else:
            def callback(name, textbox):
                try:
                    setattr(obj, name, textbox.value)
                except Exception:
                    textbox.color = 'red'

            text = Text(str(value))
            text.on_submit(functools.partial(callback, name))
            items.append(text)
    labels = [Text(name, disabled=True) for name in names]
    labels = [Text(obj.__class__.__name__, disabled=True)] + labels
    display(HBox(children=[VBox(children=labels), VBox(children=items)]))


# TODO:
# make identification grade a custom trait
# go through the workflow and determine how best to utilize the objects

if __name__ == '__main__':
    m = Group()
    m.__doc__
