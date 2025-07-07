from __future__ import absolute_import
from __future__ import print_function
import getpass
import logging
import os
import pprint
import time
import uuid
import re

from rdkit import Chem

from copy import deepcopy
from typing import Dict, Optional

import pandas as pd
from tabulate import tabulate

from .object_helpers import (
    set_docstring, Workspace, format_timestamp, MetList,
    MetUnicode, MetFloat, MetInstance, MetInt, MetEnum, MetBool, HasTraits,
    Stub
)
from six.moves import zip

logger = logging.getLogger(__name__)

#Making a new table means adding a new class to metatlas_objects.py.
#Floats are set as single precision by default, unfortunately, so here is the best way to create a table containing floats:
#Create a new table
#1) Create a new class in metatlas_objects.py.
#2) Create a new object of that class and store it.
#3) Log into database
#4) Run alter table TABLE_NAME modify COLUMN_NAME double; for each float column
#Add a floating point object to a new table
#1) Update the class in metatlas_objects.py.
#2) Create an object with the updated class and store it.
#3) Log into database
#4) Run alter table TABLE_NAME modify COLUMN_NAME double; for each new float column


# Whether to fetch stubs automatically, disabled when we want to display
# a large number of objects.
FETCH_STUBS = True
#ADDUCTS = ('','[M]+','[M+H]+','[M+H]2+','[M+2H]2+','[M+H-H2O]2+','[M+K]2+','[M+NH4]+','[M+Na]+','[M+H-H2O]+','[M-H]-','[M-2H]-','[M-H+Cl]-','[M-2H]2-','[M+Cl]-','[2M+H]+','[2M-H]-','[M-H+Na]+','[M+K]+','[M+2Na]2+','[M-e]+','[M+acetate]-','[M+formate]-','[M-H+Cl]2-','[M-H+2Na]+','[M+3H]3+','[M-3H]3-')
ADDUCTS = ('', '[2M+ACN+H]+', '[2M+ACN+Na]+', '[2M+FA-H]-', '[2M+H]+', '[2M+Hac-H]-', '[2M+K]+', '[2M+NH4]+', '[2M+Na]+', '[2M-2H+Na]-', '[2M-H]-', '[3M-H]-', '[M+2ACN+2H]2+', '[M+2ACN+H]+', '[M+2H+Na]3+', '[M+2H]2+', '[M+2K-H]+', '[M+2Na-H]+', '[M+2Na]2+', '[M+3ACN+2H]2+', '[M+3H]3+', '[M+3Na]3+', '[M+ACN+2H]2+', '[M+ACN+H]+', '[M+ACN+Na]+', '[M+Br]-', '[M+CH3COO]-', '[M+CH3OH+H]+', '[M+Cl]-', '[M+DMSO+H]+', '[M+FA-H]-', '[M+H+2Na]3+', '[M+H+K]2+', '[M+H+NH4]2+', '[M+H+Na]2+', '[M+H-CH3NH2]+', '[M+H-H2O]+', '[M+H-H2O]2+', '[M+H-NH3]+', '[M+H]+', '[M+H]2+', '[M+Hac-H]-', '[M+IsoProp+H]+', '[M+IsoProp+Na+H]+', '[M+K-2H]-', '[M+K]+', '[M+K]2+', '[M+NH4]+', '[M+Na-2H]-', '[M+Na]+', '[M+OAc]-', '[M+TFA-H]-', '[M+acetate]-', '[M+formate]-', '[M-2H2O+H]+', '[M-2H]-', '[M-2H]2-', '[M-3H]3-', '[M-CH3]-', '[M-H+2Na]+', '[M-H+Cl]-', '[M-H+Cl]2-', '[M-H+Na]+', '[M-H-CO2-2HF]-', '[M-H2O+H]+', '[M-H2O-H]-', '[M-H]-', '[M-e]+', '[M]+', '[M]-')
POLARITY = ('positive', 'negative', 'alternating')
FRAGMENTATION_TECHNIQUE = ('hcd','cid','etd','ecd','irmpd')


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
    workspace = Workspace.get_instance()
    return workspace.retrieve(object_type, **kwargs)


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
    if not isinstance(object_type, str):
        print('remove() expects a string argument, use remove_objects() to'
              'delete actual objects.')
    workspace = Workspace.get_instance()
    workspace.remove(object_type, **kwargs)


def remove_objects(objects, all_versions=True, **kwargs):
    """Remove objects from the database.

    Parameters
    ----------
    all_versions: boolean, optional
        If True, remove all versions of the object sharing the current
        head_id.
    """
    if isinstance(objects, str):
        print('remove_objects() expects actual objects, use remove() to'
              'remove objects by type.')
    workspace = Workspace.get_instance()
    workspace.remove_objects(objects, all_versions, **kwargs)


def store(objects, **kwargs):
    """Store Metatlas objects in the database.

    Parameters
    ----------
    objects: Metatlas object or list of Metatlas Objects
        Object(s) to store in the database.
    """
    workspace = Workspace.get_instance()
    workspace.save_objects(objects, **kwargs)


@set_docstring
class MetatlasObject(HasTraits):

    name = MetUnicode('Untitled', help='Name of the object')
    description = MetUnicode('No description', help='Description of the object')
    unique_id = MetUnicode(help='Unique identifier for the object').tag(readonly=True)
    creation_time = MetInt(help='Unix timestamp at object creation').tag(readonly=True)
    username = MetUnicode(help='Username who created the object').tag(readonly=True)
    last_modified = MetInt(help='Unix timestamp at last object update').tag(readonly=True)
    prev_uid = MetUnicode(help='Unique id of previous version').tag(readonly=True)
    head_id = MetUnicode(help='Unique id of most recent version of this object').tag(readonly=True)
    _loopback_guard = MetBool(False).tag(readonly=True)
    _changed = MetBool(False).tag(readonly=True)

    def __init__(self, **kwargs):
        """Set the default attributes."""
        #logger.debug('Creating new instance of %s with parameters %s', self.__class__.__name__, kwargs)
        kwargs.setdefault('unique_id', uuid.uuid4().hex)
        kwargs.setdefault('head_id', kwargs['unique_id'])
        kwargs.setdefault('username', getpass.getuser())
        kwargs.setdefault('creation_time', int(time.time()))
        kwargs.setdefault('last_modified', int(time.time()))
        super(MetatlasObject, self).__init__(**kwargs)
        self._changed = True
        self.observe(self._on_update, type='change')

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
        logger.debug('Started cloning instance of %s with recursive=%s', self.__class__.__name__, recursive)
        obj = self.__class__()
        for (tname, trait) in self.traits().items():
            if tname.startswith('_') or trait.metadata.get('readonly', False):
                continue
            val = getattr(self, tname)
            if recursive and isinstance(trait, MetList):
                val = [v.clone(True) for v in val]
            elif recursive and isinstance(trait, MetInstance) and val:
                val = val.clone(True)
            setattr(obj, tname, val)
        obj.prev_uid = self.unique_id
        obj.head_id = self.unique_id
        obj.unique_id = uuid.uuid4().hex
        logger.debug('Finished cloning instance of %s object (%s to %s)', self.__class__.__name__, self.unique_id, obj.unique_id)
        return obj

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
            if tname.startswith('_') or trait.metadata['readonly']:
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
        print((tabulate(msg)))

    def _on_update(self, change):
        """When the model changes, set the update fields.
        """
        if self._loopback_guard or change['name'].startswith('_'):
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
        return pprint.pformat(state) #str(state)#

    def __getattribute__(self, name):
        """Automatically resolve stubs on demand.
        """
#         value = super(MetatlasObject, self).__getattribute__(name)
        value = super().__getattribute__(name)

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
    #Time = MetList()
    #Flow = MetList()
    #A_percent = MetList()
    #B_percent = MetList()
    #Chromatography Stack
    #Model
    #Serial Number
    #Modules
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


# def load_lcms_files(mzml_files):
#     """Parse mzML files and load them into LcmsRun objects.

#     Note: This should be done automatically for runs in
#     /project/projectdirs/metatlas/raw_data/<username>

#     Parameters
#     ----------
#     mzml_files: list of str
#        List of paths to mzml_files.

#     Returns
#     -------
#     runs: list
#        List of LcmsRun objects.
#     """
#     runs = []
#     for fname in mzml_files:
#         hdf5_file = fname.replace('.mzML', '.h5')
#         if os.path.exists(hdf5_file):
#             print('File already exists: %s' % hdf5_file)
#             continue
#         try:
#             from metatlas import LcmsRun, mzml_to_hdf
#             hdf_file = mzml_to_hdf(fname)
#             user = getpwuid(os.stat(fname).st_uid).pw_name
#             filename = os.path.splitext(os.path.basename(fname))[0]
#             dirname = os.path.dirname(fname)
#             experiment = os.path.basename(dirname)
#             description = experiment + ' ' + filename
#             ctime = os.stat(fname).st_ctime
#             run = LcmsRun(name=filename, description=description,
#                           created_by=user,
#                           modified_by=user,
#                           created=ctime, last_modified=ctime,
#                           mzml_file=fname, hdf5_file=hdf_file)
#             runs.append(run)
#         except Exception as e:
#             print(e)
#     store(runs)
#     return runs


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
    experiment = MetUnicode(help='The name of the experiment')
    hdf5_file = MetUnicode(help='Path to the HDF5 file at NERSC')
    mzml_file = MetUnicode(help='Path to the MZML file at NERSC')
    acquisition_time = MetInt(help='Unix timestamp when data was acquired creation')
    injection_volume = MetFloat()
    injection_volume_units = MetEnum(('uL', 'nL'), 'uL')
    pass_qc = MetBool(help= 'True/False for if the LCMS Run has passed a quality control assessment')
    sample = MetInstance(Sample)


@set_docstring
class FunctionalSet(MetatlasObject):
    """Functional sets of compounds.  For example set called "hexose" would
    include "glucose, galactose, etc".  Functional sets can be sets-of-sets.
    "Sugars" would be a set that contains "Hexoses".
    """
    enabled = MetBool(True)
    members = MetList(MetInstance(MetatlasObject))


@set_docstring
class Compound(MetatlasObject):
    """A Compound is a structurally distinct entry.  The majority of MetAtlas compounds are from a merge
    of WikiData, miBIG, HMDB, ChEBI, LipidMaps, MetaCyc, GNPS, ENZO-Library, MSMLS-Library.  Compounds that
    had an unparseable structural identifier by RDKIT June, 2016, were ignored.  Distinct molecules are found
    by inchi-key of neutralized and de-salted molecules.
    """
    #name is inherited by all metatlas objects and is the most commonly used name for each compound
    #Description is a short text description of the compound
    iupac_name = MetUnicode(help='IUPAC International Chemical Identifier, optional')
    synonyms = MetUnicode()
    source=MetUnicode()
    chebi_id=MetUnicode()
    hmdb_id=MetUnicode()
    img_abc_id=MetUnicode()
    kegg_id=MetUnicode()
    lipidmaps_id=MetUnicode()
    metacyc_id=MetUnicode()
    pubchem_compound_id=MetUnicode()
    pubchem_url = MetUnicode(help='Reference database table url')
    wikipedia_url = MetUnicode(help='Reference database table url')
    kegg_url = MetUnicode(help='Reference database table url')
    hmdb_url = MetUnicode(help='Reference database table url')
    chebi_url = MetUnicode(help='Reference database table url')
    lipidmaps_url = MetUnicode(help='Reference database table url')

    #RDKIT Calculates these with some helper functions
    formula = MetUnicode()
    mono_isotopic_molecular_weight = MetFloat()
    permanent_charge = MetInt()
    number_components = MetInt(help='Must be one or greater')
    num_free_radicals = MetInt()
    inchi = MetUnicode()
    inchi_key = MetUnicode()
    neutralized_inchi = MetUnicode()
    neutralized_inchi_key = MetUnicode()
    neutralized_2d_inchi = MetUnicode()
    neutralized_2d_inchi_key = MetUnicode()

    #reference_xrefs = MetList(MetInstance(ReferenceDatabase),
    #                       help='Tag a compound with compound ids from ' +
    #                            'external databases')
    #functional_sets = MetList(MetInstance(FunctionalSet))


@set_docstring
class Reference(MetatlasObject):
    """Place holder for future reference sources.
    We expect many in silico methods will soon be robust enough to suggest
    retention times, m/z, and fragmentation.
    Pactolus is a great example of this.
    """
    lcms_run = MetInstance(LcmsRun)
    enabled = MetBool(True)
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


ID_GRADES: Dict[str, IdentificationGrade] = dict()


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
            objects = Workspace.get_instance().retrieve('identificationgrade', name=value.upper())
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
    short_name = MetUnicode()




@set_docstring
class MzIntensityPair(MetatlasObject):
    mz = MetFloat()
    intensity = MetFloat()


@set_docstring
class FragmentationReference(Reference):
    #This is specific for storing MS2 fragmentation spectra
    #A Fragmentation Tree will be added as a datatype when MS^n is deposited
    polarity = MetEnum(POLARITY, 'positive')
    precursor_mz = MetFloat()
    isolation_window = MetFloat(-1.0,help='width of the isolation window in Daltons')
    collision_energy = MetUnicode()#MetFloat()
    # adduct = MetEnum(ADDUCTS,'',help='Adduct')
    technique = MetEnum(FRAGMENTATION_TECHNIQUE,'cid')
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
    adduct = MetEnum(ADDUCTS,'',help='Adduct')
    #add when needed: charge = MetFloat(help='the charge on the m/z feature')
    modification = MetUnicode(help='Optional modification')
    observed_formula = MetUnicode(help='Optional observed formula')

@set_docstring
class IntensityReference(Reference):
    """Source of the assertion that a compound has a given m/z and
    other properties directly tied to m/z.
    """
    peak_height = MetFloat()
    peak_area = MetFloat()
    amount = MetFloat()
    amount_units = MetEnum(('nmol', 'mol'), 'nmol')

@set_docstring
class CompoundIdentification(MetatlasObject):
    """A CompoundIdentification links multiple sources of evidence about a
    compound's identity to an Atlas."""
    compound = MetList(MetInstance(Compound))
    identification_grade = _IdGradeTrait(
        help='Identification grade of the id (can be specified by a letter A-H'
    )
    identification_notes = MetUnicode('',
                             help='notes about this identifiation')
    ms1_notes = MetUnicode('',
                             help='notes about ms1 peaks')
    ms2_notes = MetUnicode('',
                             help='notes about ms2 matches')
    mz_references = MetList(MetInstance(MzReference))
    rt_references = MetList(MetInstance(RtReference))
    frag_references = MetList(MetInstance(FragmentationReference))
    intensity_references = MetList(MetInstance(IntensityReference))
    internal_standard_id = MetUnicode(help='Freetext identifier for an internal standard')
    do_normalization = MetBool(False)
    internal_standard_to_use = MetUnicode(help='identifier of which internal standard to normalize by')


@set_docstring
class Atlas(MetatlasObject):
    """An atlas contains many compound_ids."""
    compound_identifications = MetList(
        MetInstance(CompoundIdentification),
        help='List of Compound Identification objects')


# @set_docstring
# class SampleSet(MetatlasObject):
#     lcmsruns = MetList(MetInstance(LcmsRun))



@set_docstring
class MZMineTask(MetatlasObject):
    """
    For a collection of lcms runs, perform untargeted analysis with a scriptable binary
    Store the run parameters in the database for reuse later
    """
    lcmsruns = MetList(MetInstance(LcmsRun))
    output_csv = MetUnicode(help='Path to the output csv file at NERSC')
    output_project = MetUnicode(help='Path to the output project file at NERSC')
    input_xml = MetUnicode(help='Path to the input xml file at NERSC')
    input_xml_text = MetUnicode(help='Text of the batch xml file')

    mz_tolerance = MetFloat(8.0)
    mz_tolerance_units = MetEnum(('ppm', 'Da'), 'ppm')
    polarity = MetEnum(POLARITY, 'positive')

    min_peak_duration = MetFloat(0.015)
    max_peak_duration = MetFloat(30)
    rt_tol_perfile = MetFloat(0.015)
    rt_tol_multifile = MetFloat(0.15)
    rt_units = MetEnum(('sec', 'min'), 'min')

    noise_floor = MetFloat(40000.0,help='Signals below this value are not considered')
    min_peak_height = MetFloat(100000.0,help='max of eic must be at least this big')

    mzmine_launcher = MetUnicode(help='Path to the shell script that launches mzmine file at NERSC')

# @set_docstring
# class PactolusTask(MetatlasObject):
#     """
#     For an LCMS Run, search its msms spectra against pactolus trees
#     """
#     lcmsrun = MetInstance(LcmsRun)
#     output_file = MetUnicode(help='Path to the output hdf5 file at NERSC')
#     input_file = MetUnicode(help='Path to the input hdf5 container file at NERSC')
#     mz_tol = MetFloat(help='mz tolerance in Daltons')
#     polarity = MetEnum(POLARITY, 'positive')
#     min_intensity = MetFloat(help='minimum precursor ion intensity')
#     pactolus_tree_directory = MetUnicode(help='Path to the directory containing pactolus trees at NERSC')
#     pactolus_launcher = MetUnicode(help='Path to the shell script that launches pactolus search at NERSC')

def find_invalid_runs(**kwargs):
    """Find invalid runs.
    """
    override = kwargs.pop('_override', False)
    if not override:
        kwargs.setdefault('username', getpass.getuser())
    else:
        kwargs.setdefault('username', '*')
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
# workspace = Workspace()


def to_dataframe(objects):
    """
    Convert a set of Metatlas objects into a dataframe.
    """
    global FETCH_STUBS
    # we want to handle dates, enums, and use ids for objects
    FETCH_STUBS = False
    objs = [o._trait_values.copy() for o in objects if o.__class__ == objects[0].__class__]
    if not objs:
        FETCH_STUBS = True
        return pd.DataFrame()
    FETCH_STUBS = True
    enums = []
    cols = []
    # remove lists and use strings for objects
    for (tname, trait) in objects[0].traits().items():
        if tname.startswith('_') or tname in ['head_id', 'prev_uid']:
            continue
        cols.append(tname)
        if isinstance(trait, MetList):
            [o.__setitem__(tname, [i.unique_id for i in o[tname]])
             for o in objs]
        elif isinstance(trait, MetInstance):
            for obj_id,o in enumerate(objs):
                if tname not in o:
                    o[tname] = 'None'
        elif isinstance(trait, MetEnum):
            enums.append(tname)
        else:
            for obj_id,o in enumerate(objs):
                if tname not in o:
                    o[tname] = 'None'

    dataframe = pd.DataFrame(objs)[sorted(cols)]
#     for col in enums:
#         dataframe[col] = dataframe[col].astype('category')
    for col in ['last_modified', 'creation_time']:
        dataframe[col] = pd.to_datetime(dataframe[col], unit='s')
    return dataframe


def adjust_atlas_rt_range(
    in_atlas: Atlas, rt_min_delta: Optional[float], rt_max_delta: Optional[float]
) -> Atlas:
    """Reset the rt_min and rt_max values by adding rt_min_delta or rt_max_delta to rt_peak"""
    if rt_min_delta is None and rt_max_delta is None:
        return in_atlas
    out_atlas = deepcopy(in_atlas)
    for cid in out_atlas.compound_identifications:
        rts = cid.rt_references[0]
        rts.rt_min = rts.rt_min if rt_min_delta is None else rts.rt_peak + rt_min_delta
        rts.rt_max = rts.rt_max if rt_max_delta is None else rts.rt_peak + rt_max_delta
    return out_atlas

def is_inchi_key(s: str) -> bool:
    return bool(re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$").match(s))

def smiles_to_inchi_key(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    return Chem.inchi.MolToInchiKey(mol)

def ensure_inchi_keys(compounds):
    result = []
    for c in compounds:
        if is_inchi_key(c):
            result.append(c)
        else:
            try:
                inchi_key = smiles_to_inchi_key(c)
                result.append(inchi_key)
            except Exception as e:
                raise ValueError(f"Could not convert '{c}' to InChIKey: {e}")
    return result

def subset_atlas_by_compounds(
    in_atlas: Atlas, compounds: Optional[MetList] = None
) -> Atlas:
    """Return a new atlas with only the specified compounds or inchi_keys."""
    if compounds is None or not compounds:
        logger.info(f"No compounds specified for subsetting {in_atlas.name} ({in_atlas.unique_id}), returning the original atlas.")
        return in_atlas
    out_atlas = deepcopy(in_atlas)
    # If compounds is a list of strings, ensure all are inchi_keys (convert SMILES if needed)
    if isinstance(compounds, list) and all(isinstance(c, str) for c in compounds):
        try:
            inchi_keys = set(ensure_inchi_keys(compounds))
        except Exception as e:
            logger.error(f"Error converting compounds to InChIKeys: {e}")
            raise
        out_atlas.compound_identifications = [
            cid for cid in out_atlas.compound_identifications
            if any(getattr(comp, "inchi_key", None) in inchi_keys for comp in cid.compound)
        ]
    else:
        out_atlas.compound_identifications = [
            cid for cid in out_atlas.compound_identifications
            if any(c in compounds for c in cid.compound)
        ]
    logger.info(f"Compounds were specified for subsetting atlas {in_atlas.name} ({in_atlas.unique_id}): {set(compounds)}")
    logger.info(f"Filtered atlas from {len(in_atlas.compound_identifications)} to {len(out_atlas.compound_identifications)} compounds.")
    return out_atlas
