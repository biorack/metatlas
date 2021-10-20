""" For minipulating msms_refs files """
# pylint: disable=too-few-public-methods,missing-function-docstring,too-many-arguments

import functools
import logging
import math
import os

from enum import Enum
from pathlib import Path
from typing import Any, Callable, cast, Dict, Optional, List, Tuple, TypedDict, Union

import ipysheet
import ipywidgets as widgets
import matchms
import numpy as np
import pandas as pd
import traitlets

from pandas.api.types import CategoricalDtype
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import RDLogger
from traitlets import Float, HasTraits, Instance, Int, TraitError, TraitType, Unicode, validate
from tqdm.notebook import tqdm

from metatlas.datastructures import metatlas_objects as metob
from metatlas.datastructures.spectrum import Spectrum, str_to_spectrum
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.plots import dill2plots as dp
from metatlas.tools import cheminfo

logger = logging.getLogger(__name__)

COLUMN_WIDTH = 5

H5_FILE_NAME_SYSTEM_POS = 6
H5_FILE_NAME_OPTIONAL_POS = 14

NEW_REFS_DB_NAME = "metatlas"

DEFAULT_COLLISION_ENERGY = "CE102040"

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

OUTPUT_COLUMNS = {
    "Name": "name",
    "Chemical source": "id",
    "Adduct": "adduct",
    "Polarity": "polarity",
    "Exact mass": "exact_mass",
    "Precursor m/z": "precursor_mz",
    "Spectra": "spectrum",
    "m/z by intensity": "most_intense_mzs",
    "Collision energy": "collision_energy",
    "Frag. method": "fragmentation_method",
    "Instrument": "instrument",
    "Instrument type": "instrument_type",
    "Formula": "formula",
    "Inchi key": "inchi_key",
}

HELP_TEXT = (
    "Compound Name must be in the synonym list for the corresponding database entry\n"
    "Inchi or Smiles cannot be an Inchi Key\n"
    "The supported adducts can be found at "
    "https://github.com/matchms/matchms/blob/master/matchms/data/known_adducts_table.csv\n"
    f"Allowed values for Instrument Type are {', '.join(INSTRUMENT_TYPES)}.\n"
    f"Allowed values for Frag. Method are {', '.join(FRAG_METHODS)}.\n"
    "m/z Tol. is a relative tolerance value in expressed units of parts per million\n"
    "File Name should contain an absolute path.\n"
    "\n"
    "All rows must be filled in. In the parameter block, set num_rows_to_add to change the number of rows."
)

GROUP_SPECTRUM_COLS = ["rt", "polarity", "precursor_MZ", "precursor_intensity", "collision_energy"]


class Input:
    """Properties of an input to the spectrum extraction"""

    def __init__(self, identifier, label, basic_type, validator):
        self.identifier: str = identifier
        self.label: str = label
        self.basic_type: str = basic_type
        self.validator: Callable = validator


INPUTS = [
    Input("name", "Compound name", "text", lambda x: not is_bad_name(x["name"], x["molecule_id"])),
    Input(
        "molecule_id",
        "Inchi or Smiles",
        "text",
        lambda x: cheminfo.inchi_or_smiles_to_molecule(x["molecule_id"]) is not None,
    ),
    Input("chemical_source", "Chemical source", "text", lambda x: len(x["chemical_source"]) > 0),
    Input("adduct", "Adduct", "text", lambda x: cheminfo.valid_adduct(x["adduct"])),
    Input(
        "instrument_type",
        "Instrument type",
        INSTRUMENT_TYPES,
        lambda x: x["instrument_type"] in INSTRUMENT_TYPES,
    ),
    Input(
        "fragmentation_method",
        "Frag. method",
        FRAG_METHODS,
        lambda x: x["fragmentation_method"] in FRAG_METHODS,
    ),
    Input("mz_tolerance", "m/z tol. [ppm]", "numeric", lambda x: is_pos_number(x["mz_tolerance"])),
    Input("rt_min", "Min RT [min.]", "numeric", lambda x: is_valid_rt_min(x["rt_min"], x["rt_max"])),
    Input("rt_max", "Max RT [min.]", "numeric", lambda x: is_valid_rt_max(x["rt_min"], x["rt_max"])),
    Input("h5_file_name", "File name (.h5)", "text", lambda x: is_readable_file(x["h5_file_name"])),
]


class InputRecord:
    """Type for holding one row from input sheet"""

    # pylint: disable=too-many-instance-attributes

    def __init__(
        self,
        name: str,
        chemical_source: str,
        molecule_id: str,
        adduct: str,
        instrument_type: str,
        fragmentation_method: str,
        h5_file_name: str,
        mz_tolerance: float,
        rt_min: float,
        rt_max: float,
    ) -> None:
        self.name = name
        self.chemical_source = chemical_source
        self.molecule_id = molecule_id
        self.adduct = adduct
        self.instrument_type = instrument_type
        self.fragmentation_method = fragmentation_method
        self.h5_file_name = h5_file_name
        self.mz_tolerance = mz_tolerance
        self.rt_min = rt_min
        self.rt_max = rt_max

    def __repr__(self) -> str:
        return (
            f"{self.name};{self.chemical_source};{self.molecule_id};{self.adduct};"
            f"{self.instrument_type};{self.fragmentation_method};{self.h5_file_name};"
            f"{self.mz_tolerance};{self.rt_min};{self.rt_max}"
        )

    def __hash__(self) -> int:
        return hash(repr(self))


class LayoutPosition(Enum):
    """Define vertical ordering of element in GUI"""

    NAME_INPUT = 0
    SEARCH_OUTPUT = 1
    SHEET_INPUT = 2
    EXTRACT_HELP = 3
    SHEET_OUTPUT = 4
    SAVE = 5
    LOG = 6


def is_number(value: Any) -> bool:
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


def to_float(value: str) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def is_pos_number(value: Any) -> bool:
    return is_number(value) and float(value) >= 0


def get_compound(inchi_key: str) -> Optional[metob.Compound]:
    """
    Returns first compound from database matching inchi_key
    or None if not found
    """
    try:
        return metob.retrieve("Compounds", inchi_key=inchi_key, username="*")[0]
    except IndexError:
        return None


def is_bad_name(name: str, inchi: str) -> bool:
    """Returns true if the molecule is in the database but name is not a synonym"""
    if len(name) == 0:
        return True
    inchi_key = Chem.inchi.InchiToInchiKey(inchi)
    compound_result = get_compound(inchi_key)
    if compound_result is None:
        return False
    return not cheminfo.is_synonym(name, compound_result.name)


class Proposal(TypedDict):
    """for use with traitlets.validate"""

    owner: HasTraits
    value: str
    trait: TraitType


def _valid_enum(proposal: Proposal, name: str, values_list: List[str]) -> Any:
    """generic validation for enumerated type"""
    if proposal["value"] not in values_list:
        raise TraitError(f"{name} must be one of {', '.join(values_list)}")
    return proposal["value"]


def _valid_not_len_zero(proposal: Proposal, name: str) -> Any:
    """generic validation for length greater than 0"""
    if len(proposal["value"]) == 0:
        raise TraitError(f"{name} cannot have a length of zero")
    return proposal["value"]


def _valid_positive(proposal: Proposal, name: str) -> Any:
    """generic validation for positive value"""
    if float(proposal["value"]) < 0:
        raise TraitError(f"{name} must be positive")
    return proposal["value"]


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
    id: str = Unicode(allow_none=True)
    name: str = Unicode(allow_none=True)
    spectrum: Spectrum = Instance(klass=Spectrum, allow_none=True)
    decimal: np.ushort = Int(default_value=5, allow_none=True)
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
            value = getattr(self, k) if self.trait_has_value(k) else ""
            not_na_values.append("" if value is None else str(value))
        return ";".join(not_na_values)

    def __str__(self) -> str:
        return self.__repr__()

    @validate("database")
    def _valid_database(self, proposal: Proposal) -> str:
        """valid if database string has positive length"""
        return _valid_not_len_zero(proposal, "database")

    @validate("id")
    def _valid_id(self, proposal: Proposal) -> str:
        """valid if id string has positive length"""
        return _valid_not_len_zero(proposal, "id")

    @validate("name")
    def _valid_name(self, proposal: Proposal) -> str:
        """valid if name string has positive length"""
        return _valid_not_len_zero(proposal, "name")

    @validate("decimal")
    def _valid_decimal(self, proposal: Proposal) -> int:
        """valid if decimal is positive"""
        return _valid_positive(proposal, "decimal")

    @validate("precursor_mz")
    def _valid_precursor_mz(self, proposal: Proposal) -> float:
        """valid if precursor_mz is positive"""
        return _valid_positive(proposal, "precursor_mz")

    @validate("polarity")
    def _valid_polarity(self, proposal: Proposal) -> str:
        """valid if polarity is in POLARITIES"""
        return _valid_enum(proposal, "polarity", POLARITIES)

    @validate("adduct")
    def _valid_adduct(self, proposal: Proposal) -> str:
        """valid if adduct string has positive length"""
        return _valid_not_len_zero(proposal, "adduct")

    @validate("fragmentation_method")
    def _valid_fragmentation_method(self, proposal: Proposal) -> str:
        """valid if fragmentation_method in FRAG_METHODS"""
        return _valid_enum(proposal, "fragmentation_method", FRAG_METHODS)

    @validate("collision_energy")
    def _valid_collision_energy(self, proposal: Proposal) -> str:
        """valid if collision_energy has positive length"""
        return _valid_not_len_zero(proposal, "collision_energy")

    @validate("instrument")
    def _valid_instrument(self, proposal: Proposal) -> str:
        """valid if instrument has positive length"""
        return _valid_not_len_zero(proposal, "instrument")

    @validate("instrument_type")
    def _valid_instrument_type(self, proposal: Proposal) -> str:
        """valid if instrument_type is in INSTRUMENT_TYPES"""
        return _valid_enum(proposal, "instrument_type", INSTRUMENT_TYPES)

    @validate("formula")
    def _valid_formula(self, proposal: Proposal) -> str:
        """valid if formula has positive length"""
        return _valid_not_len_zero(proposal, "formula")

    @validate("exact_mass")
    def _valid_exact_mass(self, proposal: Proposal) -> float:
        """valid if exact_mass is positive"""
        return _valid_positive(proposal, "exact_mass")

    @validate("inchi_key")
    def _valid_inchi_key(self, proposal: Proposal) -> str:
        """valid if inchi_key has positive length"""
        return _valid_not_len_zero(proposal, "inchi_key")

    @validate("inchi")
    def _valid_inchi(self, proposal: Proposal) -> str:
        """valid if inchi matches with inchi_key"""
        new_value = str(proposal["value"])
        if not matchms.utils.is_valid_inchi(new_value):
            raise TraitError("not valid inchi")
        if not cheminfo.is_valid_inchi_pair(new_value, self.inchi_key):
            raise TraitError("inchi and inchi_key do not represent the same molecule")
        if not cheminfo.is_valid_inchi_smiles_pair(new_value, self.smiles):
            raise TraitError("inchi and smiles do not represent the same molecule")
        return new_value

    @validate("smiles")
    def _valid_smiles(self, proposal: Proposal) -> str:
        """valid if smiles matches with inchi"""
        new_value = str(proposal["value"])
        if not matchms.utils.is_valid_smiles(new_value):
            raise TraitError("Invalid smiles values")
        if not cheminfo.is_valid_inchi_smiles_pair(self.inchi, new_value):
            raise TraitError("inchi and smiles do not represent the same molecule")
        return new_value

    @traitlets.default("formula")
    def _get_default_formula(self) -> Optional[str]:
        """generate formula from inchi"""
        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None:
            return mol
        return Chem.rdMolDescriptors.CalcMolFormula(mol)

    @traitlets.default("smiles")
    def _get_default_smiles(self) -> Optional[str]:
        """generate smiles from inchi"""
        if self.inchi is not None and self.inchi != "":
            return Chem.MolToSmiles(Chem.inchi.MolFromInchi(self.inchi))
        return None

    @traitlets.default("inchi")
    def _get_default_inchi(self) -> Optional[str]:
        """generate inchi from smiles"""
        if self.smiles is not None and self.smiles != "":
            return Chem.inchi.MolToInchi(Chem.MolFromSmiles(self.smiles))
        return None

    def has_missing_fields(self) -> bool:
        """Returns True if there are fields with None values, logs an error message for each field missing"""
        out = False
        for name in REFS_TYPES:
            value = getattr(self, name, None)
            if value is None or value == "":
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
        if not cheminfo.is_valid_inchi_pair(self.inchi, self.inchi_key):
            logger.error("Invalid inchi/inchi_key pair for %s.", self.name)
            bad = True
        if not cheminfo.is_valid_inchi_smiles_pair(self.inchi, self.smiles):
            logger.error("Invalid inchi/smiles pair for %s.", self.name)
            bad = True
        ref_compound = get_compound(self.inchi_key)
        if ref_compound is None:
            logger.warning(
                "Could not find inchi_key=%s in database (name=%s), so skipping some tests.",
                self.inchi_key,
                self.name,
            )
            return bad
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
        if not cheminfo.is_synonym(self.name, ref_compound.synonyms):
            logger.error(
                "The entry with inchi_key=%s does not contain name '%s' in database.",
                self.inchi_key,
                self.name,
            )
            bad = True
        return bad


def read_msms_refs(file_name: str, sep: str = "\t", **kwargs) -> pd.DataFrame:
    """Read in msms refs from file with correct types"""
    file_df = pd.read_csv(file_name, sep=sep, dtype=REFS_TYPES, **kwargs)
    logger.info("Read in %d existing references from %s", len(file_df), file_name)
    return file_df


def get_empty_refs() -> pd.DataFrame:
    """Returns an empty MSMS refs DataFrame with the correct columns and types"""
    return pd.DataFrame(data={k: [] for k, v in REFS_TYPES.items()}).astype(REFS_TYPES)


def df_row_to_ref(data: dict) -> MsmsRef:
    """converts a row from df.to_dict(orient='records') to a MsmsRef instance"""
    data_minus_na = {k: v for k, v in data.items() if pd.notna(v)}
    if "spectrum" in data_minus_na:
        data_minus_na["spectrum"] = str_to_spectrum(data_minus_na["spectrum"])
    return MsmsRef(**data_minus_na)


def get_num_bad_refs(refs_df: pd.DataFrame) -> int:
    """Return number of rows that fail validation in refs_df. Info on failures to logger"" """
    return sum([0] + [1 if df_row_to_ref(row).is_bad() else 0 for row in refs_df.to_dict(orient="records")])


def get_invalid_cells(sheet: ipysheet.sheet, input_defs: List[Input]) -> List[Tuple[int, str]]:
    bad_cells = []
    for row_num, values in enumerate(sheet.cells[0].value):
        row_dict = row_list_to_dict(values, input_defs)
        for current_def in input_defs:
            is_good = current_def.validator(row_dict)
            if not is_good:
                bad_cells.append((row_num, current_def.label))
    return bad_cells


def row_list_to_dict(values: List[str], input_defs: List[Input]) -> Dict[str, Any]:
    return {x.identifier: v for x, v in zip(input_defs, values)}


def row_list_to_rec(values: List[str], input_defs: List[Input]) -> InputRecord:
    parameters = {x.identifier: v for x, v in zip(input_defs, values)}
    return InputRecord(
        name=parameters["name"],
        chemical_source=parameters["chemical_source"],
        molecule_id=parameters["molecule_id"],
        adduct=parameters["adduct"],
        instrument_type=parameters["instrument_type"],
        fragmentation_method=parameters["fragmentation_method"],
        h5_file_name=parameters["h5_file_name"],
        mz_tolerance=float(parameters["mz_tolerance"]),
        rt_min=float(parameters["rt_min"]),
        rt_max=float(parameters["rt_max"]),
    )


def in_rt_mz_ranges(
    rtime: float, rt_min: float, rt_max: float, m_z: float, mz_target: float, mz_tol: float
) -> bool:
    """
    Inputs:
        rt: measure retention time in minutes
        rt_min: lower bound of passing RT values
        rt_max: upper bound of passing RT values
        m_z: measured mass-charge ratio
        mz_target: passing within mz_tol of this value
        mz_tol: mz_tolerance in units of ppm
    """
    return dp.within_tolerance(m_z, mz_target, mz_tol * 1e-6) and (rt_min <= rtime <= rt_max)


def refs_list_to_df(refs: List[MsmsRef]) -> pd.DataFrame:
    data: Dict[str, List[Any]] = {k: [] for k in REFS_TYPES}
    for ref in refs:
        for key in REFS_TYPES:
            data[key].append(getattr(ref, key))
    return pd.DataFrame(data=data)


def extract_most_intense(
    h5_file_name: str, molecule_id: str, adduct: str, rt_min: float, rt_max: float, mz_tol: float
) -> Tuple[Spectrum, float, float, float, str]:
    """
    Inputs:
        molecule_id: either inchi or smiles string
        mz_tol: mz_tolerance in units of ppm
    Returns Spectrum, RT, parent_mass, precursor_mz, collision_energy
    """
    mol = cheminfo.inchi_or_smiles_to_molecule(molecule_id)
    parent_mass = ExactMolWt(mol)
    precursor_mz = cheminfo.get_precursor_mz(parent_mass, adduct)
    h5_data = ma_data.df_container_from_metatlas_file(h5_file_name)
    msms_df = h5_data["ms2_pos"] if cheminfo.is_positive_mode(adduct) else h5_data["ms2_neg"]
    in_tol_df = msms_df.groupby(GROUP_SPECTRUM_COLS).filter(
        lambda x: in_rt_mz_ranges(
            x.iloc[0]["rt"], rt_min, rt_max, x.iloc[0]["precursor_MZ"], precursor_mz, mz_tol
        )
    )
    precursor_intensity_max = in_tol_df["precursor_intensity"].max()
    most_intense_df = in_tol_df.groupby(GROUP_SPECTRUM_COLS).filter(
        lambda x: precursor_intensity_max == x.iloc[0]["precursor_intensity"]
    )
    most_intense = most_intense_df.iloc[0]
    return (
        Spectrum(tuple(most_intense_df["mz"]), tuple(most_intense_df["i"])),
        most_intense["rt"],
        parent_mass,
        float(most_intense["precursor_MZ"]),
        get_collision_energy(h5_file_name),
    )


@functools.lru_cache(maxsize=None)
def build_msms_ref(in_rec: InputRecord) -> Optional[MsmsRef]:
    """MsmsRef factory"""
    ref_keys = MsmsRef().class_trait_names()
    ref_dict = {key: getattr(in_rec, key) for key in ref_keys if hasattr(in_rec, key)}
    try:
        (
            ref_dict["spectrum"],
            _,
            ref_dict["exact_mass"],
            ref_dict["precursor_mz"],
            ref_dict["collision_energy"],
        ) = extract_most_intense(
            in_rec.h5_file_name,
            in_rec.molecule_id,
            in_rec.adduct,
            in_rec.rt_min,
            in_rec.rt_max,
            in_rec.mz_tolerance,
        )
    except IndexError:
        logger.error("Matching spectrum not found for %s.", in_rec.name)
        return None
    ref_dict["polarity"] = "positive" if cheminfo.is_positive_mode(str(ref_dict["adduct"])) else "negative"
    mol = cheminfo.normalize_molecule(cheminfo.inchi_or_smiles_to_molecule(in_rec.molecule_id))
    ref_dict["inchi"] = Chem.inchi.MolToInchi(mol)
    ref_dict["smiles"] = Chem.MolToSmiles(mol)
    ref_dict["inchi_key"] = Chem.inchi.InchiToInchiKey(ref_dict["inchi"])
    ref_dict["database"] = NEW_REFS_DB_NAME
    ref_dict["id"] = in_rec.chemical_source
    ref_dict["instrument"] = get_instrument(in_rec.h5_file_name)
    return MsmsRef(**ref_dict)


def save_msms_refs(existing_refs_df: pd.DataFrame, output_file_name: str, layout: widgets.Box) -> None:
    """Create CSV file containing old and new MSMS refs"""
    with get_new_log_box(layout):
        if not is_valid_input_sheet():
            return
        new_df = generate_msms_refs_df(ipysheet.sheet("input"))
        if new_df.empty:
            logger.error("Incomplete new MSMS reference definitions. Not writing an output file.")
            return
        out_df = pd.concat([existing_refs_df, new_df])
        out_df.to_csv(output_file_name, sep="\t", index=False)
        logger.info("New MSMS references file with %d records written to %s.", len(out_df), output_file_name)


def generate_msms_refs_df(sheet: ipysheet.sheet) -> pd.DataFrame:
    """Create DataFrame containing the new MSMS refs"""
    new_refs = [build_msms_ref(row_list_to_rec(row, INPUTS)) for row in tqdm(sheet.cells[0].value)]
    if None in new_refs:
        new_refs = []
    return refs_list_to_df(cast(List[MsmsRef], new_refs))


def load_msms_refs(file_name: Optional[str], validate_existing: bool = False) -> List[MsmsRef]:
    """Read in existing MSMS refs file and valid if desired"""
    if file_name is not None and not is_readable_file(file_name):
        try:
            raise FileNotFoundError(f"{file_name} does not exist or is not readable.")
        except ValueError as err:
            logger.exception(err)
            raise err
    if file_name is None:
        refs_df = get_empty_refs()
        num_bad = 0
    else:
        refs_df = read_msms_refs(file_name)
        num_bad = get_num_bad_refs(refs_df)
        logger.info("Number of existing msms reference records not passing validation is %d", num_bad)
    if validate_existing and num_bad > 0:
        try:
            raise ValueError("All existing MSMS references must pass validation before spectrum extraction")
        except ValueError as err:
            logger.exception(err)
            raise err
    return refs_df


def create_input_sheet(
    inputs: List[Input], num_rows: int, data: Optional[List[Any]] = None
) -> ipysheet.sheet:
    input_sheet = ipysheet.sheet(
        key="input",
        rows=num_rows,
        columns=len(inputs),
        column_headers=[x.label for x in inputs],
        column_resizing=False,
        column_width=COLUMN_WIDTH,
    )
    ipysheet.easy.cell_range(data or [[""] * len(inputs)] * num_rows)
    return input_sheet


def get_log_box(layout: widgets.Box) -> widgets.Text:
    return layout.children[LayoutPosition.LOG.value]


def get_new_log_box(layout: widgets.Box) -> widgets.Text:
    log_box = get_log_box(layout)
    log_box.clear_output()
    return log_box


def display_to_log_box(layout: widgets.Box, message: str) -> None:
    with get_new_log_box(layout):
        print(message)


def is_valid_input_sheet() -> bool:
    """Validate the input sheet, logs problems"""
    input_sheet = ipysheet.sheet("input")
    invalid = get_invalid_cells(input_sheet, INPUTS)
    for row_num, col_name in invalid:
        logger.error("In row %d, invalid value for '%s'.", row_num + 1, col_name)
    if len(invalid) > 0:
        logger.error("All inputs must pass validation before spectrum extraction")
        return False
    logger.info("All inputs fields for new references have passed validation.")
    return True


def extract_all(layout: widgets.Box) -> None:
    """launch msms refs extraction and export"""
    with get_new_log_box(layout):
        if is_valid_input_sheet():
            logger.info("Extracting MSMS reference spectrums....")
            input_sheet = ipysheet.sheet("input")
            new_refs_df = generate_msms_refs_df(input_sheet)
            if new_refs_df.empty:
                return
            sheet = create_refs_sheet(new_refs_df)
            layout.children = swap_layout(layout.children, LayoutPosition.SHEET_OUTPUT.value, sheet)


def add_msms_refs(
    existing_refs_file_name: Optional[str], output_file_name: str, validate_existing: bool, num_rows: int
) -> widgets.Box:
    existing_refs_df = load_msms_refs(existing_refs_file_name, validate_existing)
    return display_ui(existing_refs_df, output_file_name, num_rows)


def display_ui(existing_refs_df: pd.DataFrame, output_file_name: str, num_rows: int) -> widgets.VBox:
    """Display spreadsheet for entering input values"""
    # Layout:
    #    Row 0: molecule name text sub-string search input box and submission button
    #    Row 1: molecule name search results: sheet or status message
    #    Row 2: main input sheet for collecting most fields
    #    Row 3: Execute and Help buttons
    #    Row 4: Output sheet
    #    Row 5: Save to File button
    #    Row 6: logging output
    elements = {}
    name_input = widgets.Text(description="Name")
    min_mw = widgets.FloatText(value=0, description="Min. MW")
    max_mw = widgets.FloatText(value=99999, description="Max. MW")
    search_button = widgets.Button(description="Molecule Search")
    elements[LayoutPosition.NAME_INPUT.value] = widgets.HBox([name_input, min_mw, max_mw, search_button])
    elements[LayoutPosition.SEARCH_OUTPUT.value] = widgets.HTML(value="")
    elements[LayoutPosition.SHEET_INPUT.value] = create_input_sheet(INPUTS, num_rows)
    extract_button = widgets.Button(description="Extract spectrums")
    help_button = widgets.Button(description="Help")
    elements[LayoutPosition.EXTRACT_HELP.value] = widgets.HBox([extract_button, help_button])
    elements[LayoutPosition.SHEET_OUTPUT.value] = widgets.HTML(value="No new MSMS refs have been added")
    elements[LayoutPosition.SAVE.value] = widgets.Button(description="Save to file")
    elements[LayoutPosition.LOG.value] = widgets.Output()
    element_list = [elements[k.value] for k in LayoutPosition]
    layout = widgets.VBox(element_list)
    search_button.on_click(lambda _: search(name_input.value, min_mw.value, max_mw.value, layout))
    extract_button.on_click(lambda _: extract_all(layout))
    help_button.on_click(lambda _: display_to_log_box(layout, HELP_TEXT))
    elements[LayoutPosition.SAVE.value].on_click(
        lambda _: save_msms_refs(existing_refs_df, output_file_name, layout)
    )
    return layout


def swap_layout(existing: List[widgets.Box], index: int, update: widgets.Box) -> List[widgets.Box]:
    """Replace existing[index] with update"""
    return [x if i != index else update for i, x in enumerate(existing)]


def is_valid_num_results(num: int, input_value: str, layout: widgets.Box, max_valid: int = 100) -> bool:
    if 0 < num <= max_valid:
        return True
    if num == 0:
        message = f"<b>No molecule names containing '{input_value}' were found in the database.</b>"
    else:
        message = f"""<b>Too many matches (>{max_valid}).
                      {num} matches of '{input_value}' were found in the database.</b>"""
    message_widget = widgets.HTML(value=message)
    layout.children = swap_layout(layout.children, LayoutPosition.SEARCH_OUTPUT.value, message_widget)
    return False


def get_synonym_matches(query: str) -> List[metob.Compound]:
    """
    Search DB for all molecules where query is a substring match within the synonym or name
    fields and then filter out duplicates by inchi_key
    """
    workspace = metob.Workspace.get_instance()
    workspace.get_connection()
    # query based on from http://mysql.rjweb.org/doc.php/groupwise_max
    sql = f"""\
            SELECT
                inchi, name -- The desired columns
            FROM
              ( SELECT  @prev := '' ) init
            JOIN
              ( SELECT  inchi_key != @prev AS first,  -- the 'GROUP BY'
                        @prev := inchi_key,           -- the 'GROUP BY'
                        inchi, name -- Also the desired columns
                    FROM  compounds -- The table
                    WHERE name LIKE '%{query}%' or synonyms LIKE '%{query}%'
                    ORDER BY inchi_key --  need to order for similar to be together
                    LIMIT 999999  -- kludge to keep the ORDER BY from being ignored
              ) x
            WHERE first;"""
    if workspace.path.startswith("sqlite:"):
        sql = f"SELECT inchi, name FROM compounds WHERE name LIKE '%{query}%' or synonyms LIKE '%{query}%'"
    out = list(workspace.db.query(sql))
    workspace.close_connection()
    return out


def filter_to_norm_inchi_in_db(dicts: List[metob.Compound]) -> List[Dict[str, str]]:
    inchi_list = [x["norm_inchi"] for x in dicts]
    results = metob.retrieve("Compound", inchi=inchi_list, username="*")
    in_db = {x.inchi for x in results}
    return [x for x in dicts if x["norm_inchi"] in in_db]


def filter_by_mw(dicts: List[Dict[str, str]], min_mw, max_mw) -> List[Dict[str, str]]:
    return [x for x in dicts if min_mw <= x["MW"] <= max_mw]


def clear_search_output(layout: widgets.Box) -> None:
    blank = widgets.HTML(value="")
    layout.children = swap_layout(layout.children, LayoutPosition.SEARCH_OUTPUT.value, blank)


def search(query: str, min_mw: float, max_mw: float, layout: widgets.Box) -> None:
    with get_new_log_box(layout):
        clear_search_output(layout)
        results = get_synonym_matches(query)
        for cur in results:
            RDLogger.DisableLog("rdApp.*")  # hide rdkit warnings
            cur["mol"] = cheminfo.normalize_molecule(Chem.inchi.MolFromInchi(cur["inchi"]))
            RDLogger.EnableLog("rdApp.*")
            cur["norm_inchi"] = Chem.inchi.MolToInchi(cur["mol"])
            cur["MW"] = ExactMolWt(cur["mol"])
        filtered = filter_by_mw(filter_to_norm_inchi_in_db(results), min_mw, max_mw)
        logger.debug("Found %d matches to %s.", len(filtered), query)
        if not is_valid_num_results(len(filtered), query, layout):
            return
        final = sorted(filtered, key=lambda x: x["MW"])
        logger.debug("Num mols: %d", len(final))
        column_names = ["", "Name", "MW", "Structure"]
        sheet = ipysheet.sheet(
            rows=len(final),
            columns=len(column_names),
            column_headers=column_names,
            column_resizing=False,
            column_width=[1, 4, 2, 10],
        )
        buttons = [widgets.Button(description="use", layout=widgets.Layout(width="100%")) for x in final]
        for button in buttons:
            button.on_click(lambda current: on_use_button_clicked(current, final, layout))
        ipysheet.column(0, buttons)
        ipysheet.column(1, [x["name"] for x in final])
        ipysheet.column(2, [ExactMolWt(x["mol"]) for x in final])
        ipysheet.column(3, [cheminfo.mol_to_image(x["mol"]) for x in final])
        layout.children = swap_layout(layout.children, LayoutPosition.SEARCH_OUTPUT.value, sheet)


def on_use_button_clicked(
    current: widgets.Button, results: List[Dict[str, str]], layout: widgets.Box
) -> None:
    molecule_sheet = layout.children[LayoutPosition.SEARCH_OUTPUT.value]
    for i, button in enumerate(molecule_sheet.cells[0].value):
        if button == current:
            clear_search_output(layout)
            add_row_with_inchi(results[i]["name"], results[i]["inchi"])
            return
    row_display = widgets.HTML(value="Could not located clicked button!")
    layout.children = swap_layout(layout.children, LayoutPosition.SEARCH_OUTPUT.value, row_display)


def update_all_cell_values(sheet_key: str, value_list: List[List[Union[str, float]]]) -> None:
    current_sheet = ipysheet.sheet(sheet_key)
    current_sheet.rows = len(value_list)
    ipysheet.easy.cell_range(value_list)  # this appends to cells
    current_sheet.cells = (current_sheet.cells[-1],)  # so only keep the last one


def add_row_with_inchi(name: str, inchi: str):
    input_sheet = ipysheet.sheet("input")
    value_list = input_sheet.cells[0].value
    new_row_values = [name, inchi] + [""] * (input_sheet.columns - 2)
    for i, row in enumerate(value_list):
        if row == [""] * input_sheet.columns:
            value_list[i] = new_row_values
            break
    else:
        value_list.append(new_row_values)
    update_all_cell_values("input", value_list)


def float_list_to_str(floats: List[float]) -> str:
    return "\n".join([f"{x:.6f}" for x in floats])


def create_refs_sheet(refs_df: pd.DataFrame, num_mzs: int = 10) -> ipysheet.sheet:
    sheet = ipysheet.sheet(
        key="output",
        rows=len(refs_df),
        columns=len(OUTPUT_COLUMNS),
        column_headers=list(OUTPUT_COLUMNS.keys()),
        column_resizing=False,
        column_width=[3 if x == "Spectra" else 1 for x in OUTPUT_COLUMNS],
    )
    for i, ref_key in enumerate(OUTPUT_COLUMNS.values()):
        if ref_key == "spectrum":
            ipysheet.column(i, [x.widget() for x in refs_df[ref_key].to_list()], read_only=True)
        elif ref_key == "most_intense_mzs":
            ipysheet.column(
                i,
                [float_list_to_str(x.mz_by_intensity()[:num_mzs]) for x in refs_df["spectrum"]],
                read_only=True,
            )
        elif ref_key in ["exact_mass", "precursor_mz"]:
            ipysheet.column(i, refs_df[ref_key].to_list(), numeric_format="0.000000", read_only=True)
        else:
            ipysheet.column(i, refs_df[ref_key].to_list(), read_only=True)
    return sheet


def get_collision_energy(file_name: str) -> str:
    """Extract collision energy from file name or return default value"""
    all_fields = Path(file_name).stem.split("_")
    optional_fields = all_fields[H5_FILE_NAME_OPTIONAL_POS].split("-")
    for field in optional_fields:
        if field.startswith("CE"):
            return field
    return DEFAULT_COLLISION_ENERGY


def get_instrument(file_name: str) -> str:
    """Extract instrumnet from file name"""
    return Path(file_name).stem.split("_")[H5_FILE_NAME_SYSTEM_POS]
