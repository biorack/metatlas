""" For minipulating msms_refs files """
# pylint: disable=too-few-public-methods,missing-function-docstring,too-many-arguments

import json
import logging
import math
import os
import uuid

from typing import Any, cast, Dict, Optional, List, Mapping, Sequence, Tuple, TypedDict, Union

import ipysheet
import ipywidgets as widgets
import matchms
import numpy as np
import pandas as pd
import traitlets

from pandas.api.types import CategoricalDtype
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
from traitlets import Float, HasTraits, Instance, Int, TraitError, TraitType, Unicode, validate

from metatlas.datastructures import metatlas_objects as metob
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.plots import dill2plots as dp
from metatlas.tools import cheminfo

logger = logging.getLogger(__name__)

COLUMN_WIDTH = 5

NEW_REFS_DB_NAME = "NorthernLabAddition:NoDB"

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

HELP_TEXT = (
    "Compound Name must be in the synonym list for the corresponding database entry\n"
    "Inchi or Smiles cannot be an Inchi Key\n"
    "The supported adducts can be found at "
    "https://github.com/matchms/matchms/blob/master/matchms/data/known_adducts_table.csv\n"
    "Instrument should contain a model name\n"
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
        self.identifier = identifier
        self.label = label
        self.basic_type = basic_type
        self.validator = validator


INPUTS = [
    Input("name", "Compound Name", "text", lambda x: not is_bad_name(x["name"], x["molecule_id"])),
    Input(
        "molecule_id",
        "Inchi or Smiles",
        "text",
        lambda x: cheminfo.inchi_or_smiles_to_molecule(x["molecule_id"]) is not None,
    ),
    Input("adduct", "Adduct", "text", lambda x: cheminfo.valid_adduct(x["adduct"])),
    Input(
        "instrument",
        "Instrument",
        "text",
        lambda x: len(x["instrument"]) > 0,
    ),
    Input(
        "instrument_type",
        "Instrument Type",
        INSTRUMENT_TYPES,
        lambda x: x["instrument_type"] in INSTRUMENT_TYPES,
    ),
    Input(
        "fragmentation_method",
        "Frag. Method",
        FRAG_METHODS,
        lambda x: x["fragmentation_method"] in FRAG_METHODS,
    ),
    Input("mz_tolerance", "m/z Tol. [ppm]", "numeric", lambda x: is_pos_number(x["mz_tolerance"])),
    Input("rt_min", "Min RT [min.]", "numeric", lambda x: is_valid_rt_min(x["rt_min"], x["rt_max"])),
    Input("rt_max", "Max RT [min.]", "numeric", lambda x: is_valid_rt_max(x["rt_min"], x["rt_max"])),
    Input("h5_file_name", "File Name (.h5)", "text", lambda x: is_readable_file(x["h5_file_name"])),
]


class InputDict(TypedDict, total=False):
    """Type for holding one row from input sheet"""
    name: str
    molecule_id: str
    adduct: str
    instrument: str
    instrument_type: str
    fragmentation_method: str
    h5_file_name: str
    mz_tolerance: float
    rt_min: float
    rt_max: float


def to_input_dict(data: Mapping[str, Any]) -> InputDict:
    result = InputDict()
    for key, key_type in InputDict.__annotations__.items():  # pylint: disable=no-member
        if key in data:
            result[key] = key_type(data[key])  # type: ignore
    return result


def is_number(value: Any) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def to_float(value: str) -> Union[float, str]:
    try:
        return float(value)
    except ValueError:
        return value


def is_pos_number(value: Any) -> bool:
    return is_number(value) and float(value) >= 0


def get_compound(inchi_key: str) -> Optional[metob.Compound]:
    """
    Returns first compound from database matching inchi_key and with username pasteur
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
            self.intensities = list(intensities)
            self.mzs = list(mzs)

    def __repr__(self) -> str:
        """Return representation of data"""
        nested_list_form = [[f"{m:.5f}" for m in self.mzs], [f"{x:.3f}" for x in self.intensities]]
        return str(nested_list_form).replace("'", "")

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
    if spectrum_str is None or spectrum_str == "":
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
        results = metob.retrieve("compounds", username="*", inchi_key=self.inchi_key)
        if len(results) == 0:
            logger.warning(
                "Could not find inchi_key=%s in database (name=%s), so skipping some tests.",
                self.inchi_key,
                self.name,
            )
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
        if not cheminfo.is_synonym(self.name, ref_compound.synonyms):
            logger.error(
                "The entry with inchi_key=%s does not contain name '%s' in database.",
                self.inchi_key,
                self.name,
            )
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


def row_list_to_dict(values: List[str], input_defs: List[Input]) -> InputDict:
    return to_input_dict(
        {x.identifier: v if x.basic_type != "numeric" else to_float(v) for x, v in zip(input_defs, values)}
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
        f"{most_intense['collision_energy']:.1f}eV",
    )


def build_msms_ref(in_dict: InputDict) -> MsmsRef:
    """MsmsRef factory"""
    ref_keys = MsmsRef().class_trait_names()
    ref_dict = {k: v for k, v in in_dict.items() if k in ref_keys}
    try:
        (
            ref_dict["spectrum"],
            _,
            ref_dict["exact_mass"],
            ref_dict["precursor_mz"],
            ref_dict["collision_energy"],
        ) = extract_most_intense(
            in_dict["h5_file_name"],
            in_dict["molecule_id"],
            in_dict["adduct"],
            in_dict["rt_min"],
            in_dict["rt_max"],
            in_dict["mz_tolerance"],
        )
    except IndexError as err:
        logger.error(f"Matching spectrum not found for {in_dict['name']}")
        raise err
    ref_dict["polarity"] = (
        "positive" if cheminfo.is_positive_mode(str(ref_dict["adduct"])) else "negative"
    )
    mol = cheminfo.normalize_molecule(cheminfo.inchi_or_smiles_to_molecule(in_dict["molecule_id"]))
    ref_dict["inchi"] = Chem.inchi.MolToInchi(mol)
    ref_dict["smiles"] = Chem.MolToSmiles(mol)
    ref_dict["inchi_key"] = Chem.inchi.InchiToInchiKey(ref_dict["inchi"])
    ref_dict["database"] = NEW_REFS_DB_NAME
    ref_dict["id"] = str(uuid.uuid4())
    return MsmsRef(**ref_dict)


def generate_msms_refs(
    existing_refs_file_name: Optional[str],
    output_file_name: str,
    sheet: ipysheet.sheet,
    validate_existing: bool = False,
) -> None:
    """Create CSV file containing old and new MSMS refs"""
    refs_df = get_empty_refs() if existing_refs_file_name is None else read_msms_refs(existing_refs_file_name)
    logger.info(
        "Number of existing msms reference records not passing validation is %d", get_num_bad_refs(refs_df)
    )
    if validate_existing and get_num_bad_refs(refs_df) > 0:
        logger.error("All existing MSMS references must pass validation before spectrum extraction")
        return
    new_refs = [build_msms_ref(row_list_to_dict(row, INPUTS)) for row in sheet.cells[0].value]
    new_df = refs_list_to_df(new_refs)
    out_df = pd.concat([refs_df, new_df])
    out_df.to_csv(output_file_name, sep="\t", index=False)
    logger.info("New MSMS references file with %d records written to %s.", len(out_df), output_file_name)


def display_inputs_ui(
    existing_refs_file_name: Optional[str], output_file_name: str, validate_existing: bool, num_rows: int
) -> widgets.Box:
    """Display spreadsheet for entering input values"""
    if existing_refs_file_name is not None and not is_readable_file(existing_refs_file_name):
        logger.error("%s does not exist or is not readable.", existing_refs_file_name)
        return widgets.Box()
    col_headers = [x.label for x in INPUTS]
    sheet = ipysheet.sheet(
        rows=num_rows,
        columns=len(INPUTS),
        column_headers=col_headers,
        column_resizing=False,
        column_width=COLUMN_WIDTH,
    )
    ipysheet.easy.cell_range([[""] * len(INPUTS)] * num_rows)
    extract = widgets.Button(description="Execute")
    docs = widgets.Button(description="Help")
    log_box = widgets.Output()

    def on_extract_clicked(_):
        """launch msms refs extraction and export"""
        log_box.clear_output()
        with log_box:
            invalid = get_invalid_cells(sheet, INPUTS)
            for row_num, col_name in invalid:
                logger.error("In row %d, invalid value for '%s'.", row_num + 1, col_name)
            if len(invalid) > 0:
                logger.error("All inputs must pass validation before spectrum extraction")
                return
            generate_msms_refs(existing_refs_file_name, output_file_name, sheet, validate_existing)

    extract.on_click(on_extract_clicked)

    def on_docs_clicked(_):
        """show help information below the sheet widget"""
        log_box.clear_output()
        with log_box:
            print(HELP_TEXT)

    docs.on_click(on_docs_clicked)
    return widgets.VBox([sheet, widgets.HBox([extract, docs]), log_box])
