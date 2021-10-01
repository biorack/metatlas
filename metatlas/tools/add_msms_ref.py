""" For minipulating msms_refs files """
import json
import logging
import math
import os
import uuid

from typing import Any, cast, Dict, Optional, List, Sequence, Tuple, TypedDict

import ipysheet
import ipywidgets as widgets
import matchms
import numpy as np
import pandas as pd
import traitlets

from pandas.api.types import CategoricalDtype
from rdkit import Chem
from traitlets import Float, HasTraits, Instance, Int, TraitError, TraitType, Unicode, validate

from metatlas.datastructures import metatlas_objects as metob
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.plots import dill2plots as dp
from metatlas.tools import cheminfo

logger = logging.getLogger(__name__)

NUM_STARTING_ROWS = 10

COLUMN_WIDTH = 3

POLARITIES = ["negative", "positive"]

INSTRUMENT_TYPES = [
    "ESI-ITFT",
    "ESI-ITTOF",
    "ESI-QTOF",
    "LC-ESI-IT",
    "LC-ESI-ITFT",
    "LC-ESI-ITTOF",
    "LC-ESI-QFT",
    "LC-ESI-QIT",
    "LC-ESI-QQ",
    "LC-ESI-QTOF",
    "Orbitrap",
]

FRAG_METHODS = [
    "CID",
    "EBEQ",
    "FT-ICR/FTMS",
    "FT-ICR/Fourier transform ion cyclotron resonance",
    "HCD",
    "IT-FT/ion trap with FTMS",
    "IT/ion trap",
    "LOW-ENERGY CID",
    "Q-TOF",
    "QqIT",
    "QqLIT",
    "QqQ",
]

polarity_type = CategoricalDtype(categories=POLARITIES, ordered=True)
frag_method_type = CategoricalDtype(categories=FRAG_METHODS, ordered=False)
instrument_type_type = CategoricalDtype(categories=INSTRUMENT_TYPES, ordered=False)

REFS_TYPES = {  # these values are pandas dtypes
    "database": "string",
    "id": "string",
    "name": "string",
    "spectrum": "string",
    "decimal": np.ushort,
    "precursor_mz": np.float64,
    "polarity": polarity_type,
    "adduct": "string",
    "fragmentation_method": frag_method_type,
    "collision_energy": "string",
    "instrument": "string",
    "instrument_type": instrument_type_type,
    "formula": "string",
    "exact_mass": np.float64,
    "inchi_key": "string",
    "inchi": "string",
    "smiles": "string",
}

REFS_DEFAULTS = {
    "database": "",
    "id": "",
    "name": "",
    "spectrum": "[[],[]]",
    "decimal": 4,
    "precursor_mz": 0,
    "polarity": "positive",
    "adduct": "",
    "fragmentation_method": FRAG_METHODS[0],
    "collision_energy": "0eV",
    "instrument": "",
    "instrument_type": INSTRUMENT_TYPES[0],
    "formula": "",
    "exact_mass": 0.0,
    "inchi_key": "InChi=",
    "inchi": "",
    "smiles": "",
}


class Input():
    def __init__(self, identifier, label, basic_type, validator):
        self.identifier = identifier
        self.label = label
        self.basic_type = basic_type
        self.validator = validator


INPUTS = [
    Input("name", "Name", "text", lambda x: not is_bad_name(x["name"], x["inchi"])),
    Input("molecule_id", "Inchi or Smiles", "text", lambda x: inchi_or_smiles_to_molecule(x["molecule_id"]) is not None),
    Input("adduct", "Adduct", "text", lambda x: valid_adduct(x["adduct"])),
    Input("instrument", "Instrument", "text", lambda x: len(x["instrument"]) > 0,),
    Input("instrument_type", "Instrument Type", INSTRUMENT_TYPES, lambda x: x["instrument_type"] in INSTRUMENT_TYPES),
    Input("fragmentation_method", "Fragmentation Method", FRAG_METHODS, lambda x: x["fragmentation_method"] in FRAG_METHODS),
    Input("mz_tolerance", "m/z Tolerance", "numeric", lambda x: is_pos_number(x["mz_tolerance"])),
    Input("rt_min", "Min RT", "numeric", lambda x: is_valid_rt_min(x["rt_min"], x["rt_max"])),
    Input("rt_max", "Max RT", "numeric", lambda x: is_valid_rt_max(x["rt_min"], x["rt_max"])),
    Input("h5_file_name", "File Name (.h5)", "text", lambda x: is_readable_file(x["h5_file_name"])),
]

REFS_V3_FILE_NAME = "/global/project/projectdirs/metatlas/projects/spectral_libraries/msms_refs_v3.tab"


def is_number(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def is_pos_number(value):
    return is_number(value) and float(value) >= 0


def is_inchi(test_inchi: str) -> bool:
    """True if input can be parsed as an inchi string"""
    return Chem.inchi.MolFromInchi(test_inchi) is not None


def is_valid_inchi_pair(test_inchi: str, test_inchi_key: str) -> bool:
    """True if if test_inchi has the inchi key test_inchi_key"""
    if not is_inchi(test_inchi):
        return False
    return test_inchi_key == Chem.inchi.InchiToInchiKey(test_inchi)


def is_valid_inchi_smiles_pair(test_inchi: str, test_smiles: str) -> bool:
    """
    True if test_inchi and test_smiles have the same structure.
    Also True if test_smiles is None and test_inchi is a valid inchi
    """
    if pd.isna(test_smiles):
        return is_inchi(test_inchi)
    mol_from_inchi = Chem.inchi.MolFromInchi(test_inchi)
    mol_from_smiles = Chem.MolFromSmiles(test_smiles)
    if mol_from_inchi is None or mol_from_smiles is None:
        return False
    return are_equal(mol_from_inchi, mol_from_smiles)


def inchi_or_smiles_to_molecule(molecule_id: str) -> Optional[Chem.rdchem.Mol]:
    return Chem.inchi.MolFromInchi(molecule_id) or Chem.MolFromSmiles(molecule_id)


def get_compound(inchi_key: str) -> Optional[metob.Compound]:
    """
    Returns first compound from database matching inchi_key and with username pasteur
    or None if not found
    """
    try:
        return metob.retrieve("Compounds", inchi_key=inchi_key, username="*")[0]
    except IndexError:
        return None


def are_equal(molecule1: Chem.rdchem.Mol, molecule2: Chem.rdchem.Mol) -> bool:
    """True if both molecules are substructures of each other"""
    return molecule1.HasSubstructMatch(molecule2) and molecule2.HasSubstructMatch(molecule1)


def is_synonym(name: str, synonym_string: str) -> bool:
    """
    Inputs:
        name: string to check for within synonym_string
        synonym_string: string with /// between names
    Returns True if case insensitive match of name to full name in synonym_string
    """
    return name.lower() in [x.lower() for x in synonym_string.split("///")]


def is_bad_name(name: str, inchi: str) -> bool:
    """Returns true if the molecule is in the database but name is not a synonym"""
    if len(name) == 0:
        return True
    inchi_key = Chem.inchi.InchiToInchiKey(inchi)
    compound_result = get_compound(inchi_key)
    if compound_result is None:
        return False
    return not is_synonym(name, compound_result.name)


class Proposal(TypedDict):
    """for use with traitlets.validate"""

    owner: HasTraits
    value: object
    trait: TraitType


class Spectrum(HasTraits):
    # pylint: disable=too-few-public-methods
    """List of intensities with list of corresponding MZ values"""
    intensities: List[float] = traitlets.List(trait=Float())
    mzs: List[float] = traitlets.List(trait=Float())

    def __init__(self, mzs: Sequence[float], intensities: Sequence[float], **kwargs) -> None:
        """required fields are inputs"""
        with self.hold_trait_notifications():
            super().__init__(**kwargs)
            self.intensities = intensities
            self.mzs = mzs

    def __repr__(self) -> str:
        """Return representation of data"""
        nested_list_form = [[f"{m:.5f}" for m in self.mzs], [f"{x:.3f}" for x in self.intensities]]
        return str(nested_list_form).replace('\'', '')

    def __str__(self) -> str:
        """Return string representation of data"""
        return self.__repr__()

    @validate("intensities")
    def _valid_intensities(self, proposal: Proposal) -> List[float]:
        """validate positive values, not empty, and same length as mzs list"""
        value = cast(List[float], proposal["value"])
        if len(value) != len(self.mzs):
            raise TraitError("length of intensities and mzs must be equal")
        if any(x <= 0 for x in value):
            raise TraitError("intensities must be positive")
        return value

    @validate("mzs")
    def _valid_mzs(self, proposal: Proposal) -> List[float]:
        """validate positive values, not empty, and same length as intensities list"""
        value = cast(List[float], proposal["value"])
        if len(value) != len(self.intensities):
            raise TraitError("length of intensities and mzs must be equal")
        if value != sorted(value):
            raise TraitError("mzs values must be monotonically increasing")
        if any(x <= 0 for x in value):
            raise TraitError("mzs values must be positive")
        return value


def str_to_spectrum(spectrum_str: str) -> Spectrum:
    """Converts a spectrum string into a Spectrum class instance"""
    if spectrum_str is None or spectrum_str == '':
        return Spectrum(mzs=[], intensities=[])
    try:
        decoded = json.loads(spectrum_str)
    except (TypeError, json.JSONDecodeError):
        logger.error("Cannot convert '%s' to a Spectrum object, setting to empty spectrum", spectrum_str)
        return Spectrum(mzs=[], intensities=[])
    if len(decoded) != 2:
        logger.error("Invalid specturm '%s'. Truncating elements after first two lists.", spectrum_str)
    return Spectrum(mzs=decoded[0], intensities=decoded[1])


def _valid_enum(proposal, name, values_list):
    """generic validation for enumerated type"""
    if proposal["value"] not in values_list:
        raise TraitError(f"{name} must be one of {', '.join(values_list)}")
    return proposal["value"]


def _valid_not_len_zero(proposal, name):
    """generic validation for length greater than 0"""
    if len(proposal["value"]) == 0:
        raise TraitError(f"{name} cannot have a length of zero")
    return proposal["value"]


def _valid_positive(proposal, name):
    """generic validation for positive value"""
    if proposal["value"] < 0:
        raise TraitError(f"{name} must be positive")
    return proposal["value"]


def valid_adduct(value):
    adducts = matchms.importing.load_adducts_dict()
    return matchms.utils.looks_like_adduct(value) and value in adducts


def is_readable_file(value):
    return os.path.isfile(value) and os.access(value, os.R_OK)


def is_valid_rt_min(rt_min: Any, rt_max: Any) -> bool:
    if not is_pos_number(rt_min):
        return False
    return is_pos_number(rt_max) and rt_min < rt_max


def is_valid_rt_max(rt_min: Any, rt_max: Any) -> bool:
    if not is_pos_number(rt_max):
        return False
    return is_pos_number(rt_min) and rt_min < rt_max


class MsmsRef(HasTraits):
    # pylint: disable=too-few-public-methods,too-many-instance-attributes
    """one line from msms_refs file"""
    database: str = Unicode(allow_none=True)
    id: str = Unicode(default=uuid.uuid4(), allow_none=True)
    name: str = Unicode(allow_none=True)
    spectrum: Spectrum = Instance(klass=Spectrum, allow_none=True)
    decimal: np.ushort = Int(default_value=4, allow_none=True)
    precursor_mz: np.float64 = Float(allow_none=True)
    polarity: str = Unicode(allow_none=True)
    adduct: str = Unicode(allow_none=True)
    fragmentation_method: str = Unicode(allow_none=True)
    collision_energy: str = Unicode(allow_none=True)
    instrument: str = Unicode(allow_none=True)
    instrument_type: str = Unicode(allow_none=True)
    formula: str = Unicode(allow_none=True)
    exact_mass: np.float64 = Float(allow_none=True)
    inchi_key: str = Unicode(allow_none=True)
    inchi: str = Unicode(allow_none=True)
    smiles: str = Unicode(allow_none=True)

    # pylint: disable=no-self-use,too-many-arguments
    def __repr__(self) -> str:
        not_na_values = []
        for k in REFS_TYPES:
            value = getattr(self, k) if self.trait_has_value(k) else ''
            not_na_values.append('' if value is None else str(value))
        return ';'.join(not_na_values)

    def __str__(self) -> str:
        return self.__repr__()

    @validate("database")
    def _valid_database(self, proposal):
        """valid if database string has positive length"""
        return _valid_not_len_zero(proposal, "database")

    @validate("id")
    def _valid_id(self, proposal):
        """valid if id string has positive length"""
        return _valid_not_len_zero(proposal, "id")

    @validate("name")
    def _valid_name(self, proposal):
        """valid if name string has positive length"""
        return _valid_not_len_zero(proposal, "name")

    @validate("decimal")
    def _valid_decimal(self, proposal):
        """valid if decimal is positive"""
        return _valid_positive(proposal, "decimal")

    @validate("precursor_mz")
    def _valid_precursor_mz(self, proposal):
        """valid if precursor_mz is positive"""
        return _valid_positive(proposal, "precursor_mz")

    @validate("polarity")
    def _valid_polarity(self, proposal):
        """valid if polarity is in POLARITIES"""
        return _valid_enum(proposal, "polarity", POLARITIES)

    @validate("adduct")
    def _valid_adduct(self, proposal):
        """valid if adduct string has positive length"""
        return _valid_not_len_zero(proposal, "adduct")

    @validate("fragmentation_method")
    def _valid_fragmentation_method(self, proposal):
        """valid if fragmentation_method in FRAG_METHODS"""
        return _valid_enum(proposal, "fragmentation_method", FRAG_METHODS)

    @validate("collision_energy")
    def _valid_collision_energy(self, proposal):
        """valid if collision_energy has positive length"""
        return _valid_not_len_zero(proposal, "collision_energy")

    @validate("instrument")
    def _valid_instrument(self, proposal):
        """valid if instrument has positive length"""
        return _valid_not_len_zero(proposal, "instrument")

    @validate("instrument_type")
    def _valid_instrument_type(self, proposal):
        """valid if instrument_type is in INSTRUMENT_TYPES"""
        return _valid_enum(proposal, "instrument_type", INSTRUMENT_TYPES)

    @validate("formula")
    def _valid_formula(self, proposal):
        """valid if formula has positive length"""
        return _valid_not_len_zero(proposal, "formula")

    @validate("exact_mass")
    def _valid_exact_mass(self, proposal):
        """valid if exact_mass is positive"""
        return _valid_positive(proposal, "exact_mass")

    @validate("inchi_key")
    def _valid_inchi_key(self, proposal):
        """valid if inchi_key has positive length"""
        return _valid_not_len_zero(proposal, "inchi_key")

    @validate("inchi")
    def _valid_inchi(self, proposal):
        """valid if inchi matches with inchi_key"""
        if not is_inchi(proposal["value"]):
            raise TraitError("not valid inchi")
        if not is_valid_inchi_pair(proposal["value"], self.inchi_key):
            raise TraitError("inchi and inchi_key do not represent the same molecule")
        return proposal["value"]

    @validate("smiles")
    def _valid_smiles(self, proposal):
        """valid if smiles matches with inchi"""
        if not is_valid_inchi_smiles_pair(self.inchi, proposal["value"]):
            raise TraitError("inchi and smiles do not represent the same molecule")
        return proposal["value"]

    @traitlets.default("smiles")
    def _get_default_smiles(self):
        """generate smiles from inchi"""
        if self.inchi is not None and self.inchi != '':
            return Chem.MolToSmiles(Chem.inchi.MolFromInchi(self.inchi))
        return None

    def has_missing_fields(self) -> bool:
        """Returns True if there are fields with None values, logs an error message for each field missing"""
        out = False
        for name in REFS_TYPES:
            value = getattr(self, name, None)
            if value is None or value == '':
                out = True
                logger.error("No '%s' field in %s", name, str(self))
        return out

    def is_bad(self) -> bool:
        """
        If returns True, then the inputs are bad, but if returns False do not assume the inputs are good
        returning False only means that there is no evidence the inputs are bad. Conclusively saying
        the inputs are good for unusual chemicals that are not in databases is hard.
        """
        bad = self.has_missing_fields()
        if self.fragmentation_method not in FRAG_METHODS:
            logger.error('Invalid fragmentation method "%s" for %s.', self.fragmentation_method, self.name)
            bad = True
        if not is_valid_inchi_pair(self.inchi, self.inchi_key):
            logger.error("Invalid inchi/inchi_key pair for %s.", self.name)
            bad = True
        if not is_valid_inchi_smiles_pair(self.inchi, self.smiles):
            logger.error("Invalid inchi/smiles pair for %s.", self.name)
            bad = True
        results = metob.retrieve("compounds", username="*", inchi_key=self.inchi_key)
        if len(results) == 0:
            logger.warning("Could not find inchi_key=%s in database (name=%s), so skipping some tests.", self.inchi_key, self.name)
            return bad
        ref_compound = results[0]
        if self.formula != ref_compound.formula:
            logger.error(
                'Formula "%s" for %s does not match value "%s" in database.',
                self.formula,
                self.name,
                ref_compound.formula,
            )
            bad = True
        if not math.isclose(self.exact_mass, ref_compound.mono_isotopic_molecular_weight, rel_tol=1e-9):
            logger.error(
                "Exact mass %s for %s does not match value %s in database.",
                self.exact_mass,
                self.name,
                ref_compound.mono_isotopic_molecular_weight,
            )
            bad = True
        if not is_synonym(self.name, ref_compound.synonyms):
            logger.error("The entry with inchi_key=%s does not contain name '%s' in database.", self.inchi_key, self.name)
            bad = True
        return bad


def read_msms_refs(file_name: str, sep="\t", **kwargs) -> pd.DataFrame:
    """Read in msms refs from file with correct types"""
    file_df = pd.read_csv(file_name, sep=sep, dtype=REFS_TYPES, **kwargs)
    logger.info("Read in %d existing references from %s", len(file_df), file_name)
    return file_df


def get_empty_refs() -> pd.DataFrame:
    """Returns an empty MSMS refs DataFrame with the correct columns and types"""
    return pd.DataFrame(data={k: [] for k, v in REFS_TYPES.items()}).astype(REFS_TYPES)


def df_row_to_ref(data: dict) -> MsmsRef:
    """ converts a row from df.to_dict(orient='records') to a MsmsRef instance"""
    data_minus_na = {k: v for k, v in data.items() if pd.notna(v)}
    if 'spectrum' in data_minus_na:
        data_minus_na['spectrum'] = str_to_spectrum(data_minus_na['spectrum'])
    return MsmsRef(**data_minus_na)


def get_num_bad_refs(refs_df: pd.DataFrame) -> int:
    """Return number of rows that fail validation in refs_df. Info on failures to logger"" """
    return sum([0] + [0 if df_row_to_ref(row).is_bad() else 1 for row in refs_df.to_dict(orient='records')])


def in_rt_mz_ranges(rt, rt_min, rt_max, mz, mz_target, mz_tol):
    return dp.within_tolerance(mz, mz_target, mz_tol) and (rt_min <= rt <= rt_max)


def extract_most_intense(in_df, rt_min, rt_max, mz_target, mz_tol):
    group_cols = ['rt', 'polarity', 'precursor_MZ', 'precursor_intensity', 'collision_energy']
    in_tol_df = in_df.groupby(group_cols).filter(lambda x: in_rt_mz_ranges(x.iloc[0]['rt'], rt_min, rt_max, x.iloc[0]['precursor_MZ'], mz_target, mz_tol))
    precursor_intensity_max = in_tol_df['precursor_intensity'].max()
    most_intense_df = in_tol_df.groupby(group_cols).filter(lambda x: precursor_intensity_max == x.iloc[0]['precursor_intensity'])
    spectrum = Spectrum(tuple(most_intense_df['mz']), tuple(most_intense_df['i']))
    most_intense = most_intense_df.iloc[0]
    return (spectrum, most_intense['rt'], most_intense['precursor_MZ'], most_intense['collision_energy'])


def extract_spectrum(h5_file_name, molecule_id, adduct, rt_min, rt_max, mz_tolerance) -> Spectrum:
    if matchms.utils.is_valid_inchi(molecule_id):
        mol = Chem.MolFromInchi(molecule_id)
    elif matchms.utils.is_valid_smiles(molecule_id):
        mol = Chem.MolFromSmiles(molecule_id)
    else:
        raise ValueError(f"molecule_id '{molecule_id}' is not a valid inchi or smiles string")
    h5_df = ma_data.df_container_from_metatlas_file(h5_file_name)
    parent_mass = Chem.Descriptors.ExactMolWt(mol)
    precursor_mz = cheminfo.get_precursor_mz(parent_mass, adduct)
    return extract_most_intense(h5_df, rt_min, rt_max, precursor_mz, mz_tolerance)


def sheet_row_to_spectrum(input_sheet, input_defs, row_num) -> Spectrum:
    row_dict = row_list_to_dict(input_sheet.cells[0].value[row_num], input_defs)
    return extract_spectrum(row_dict["h5_file_name"], row_dict["inchi"], row_dict["adduct"], float(row_dict["rt_min"]), float(row_dict["rt_max"]), float(row_dict["mz_tolerance"]))


def row_col_to_cell_num(in_sheet: ipysheet.sheet, row_num: int, col_num: int) -> int:
    return in_sheet.columns * row_num + col_num


def row_list_to_dict(values: List[Any], input_defs: List[Input]) -> Dict[str, Any]:
    return dict(zip([x.identifier for x in input_defs], values))


def get_invalid_cells(input_sheet: ipysheet.sheet, input_defs: List[Input]) -> List[Tuple[int, str]]:
    bad_cells = []
    for row_num, values in enumerate(input_sheet.cells[0].value):
        row_dict = row_list_to_dict(input_defs, values)
        for column_num, current_def in enumerate(input_defs):
            try:
                is_good = current_def.validators(row_dict)
            except Exception:
                is_good = False
            if not is_good:
                bad_cells.append((row_num, current_def.label))
    return bad_cells


def spectrums_from_sheet(input_sheet):
    for row_num in range(input_sheet.rows):
        sheet_row_to_spectrum(input_sheet, INPUTS, row_num)
    pass


def display_inputs_ui(num_rows: int) -> widgets.Box:
    """Display spreadsheet for entering input values"""
    col_headers = [x.label for x in INPUTS]
    input_sheet = ipysheet.sheet(rows=num_rows, columns=len(INPUTS), column_headers=col_headers, column_width=COLUMN_WIDTH, column_resizing=False)
    ipysheet.easy.cell_range([['']*len(INPUTS)]*num_rows)
    extract = widgets.Button(description="Extract Spectrums")
    log_box = widgets.Output()

    def on_extract_clicked(_):
        log_box.clear_output()
        invalid = get_invalid_cells(input_sheet, [x.validator for x in INPUTS])
        with log_box:
            for row_num, col_name in invalid:
                logger.error("In row %d, invalid value for '%s'.", row_num+1, col_name)
            if len(invalid) > 0:
                logger.error("All inputs must pass validation before spectrum extraction")
                return
        spectrums_from_sheet(input_sheet)

    extract.on_click(on_extract_clicked)
    return widgets.VBox([input_sheet, extract, log_box])
