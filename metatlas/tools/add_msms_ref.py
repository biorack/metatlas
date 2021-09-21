""" For minipulating msms_refs files """
import json
import logging
import math
import uuid

from typing import cast, Optional, List, TypedDict

import ipysheet
import ipywidgets as widgets
import numpy as np
import pandas as pd
import traitlets

from ipysheet import sheet, row, to_dataframe
from pandas.api.types import CategoricalDtype
from rdkit import Chem
from traitlets import TraitError, default, validate
from traitlets import Float, HasTraits, Instance, Int, TraitType, Unicode

from metatlas.datastructures import metatlas_objects as metob

logger = logging.getLogger(__name__)

REFS_V3_FILE_NAME = "/global/project/projectdirs/metatlas/projects/spectral_libraries/msms_refs_v3.tab"

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
    "cid",
]

polarity_type = CategoricalDtype(categories=POLARITIES, ordered=True)
frag_method_type = CategoricalDtype(categories=FRAG_METHODS, ordered=False)
instrument_type_type = CategoricalDtype(categories=INSTRUMENT_TYPES, ordered=False)

REFS_TYPES = {
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


def get_compound(inchi_key: str) -> Optional[metob.Compound]:
    """
    Returns first compound from database matching inchi_key and with username pasteur
    or None if not found
    """
    try:
        return metob.retrieve("Compounds", inchi_key=inchi_key, username="pasteur")[0]
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

    def __init__(self, mzs: List[float], intensities: List[float], **kwargs) -> None:
        """required fields are inputs"""
        with self.hold_trait_notifications():
            super().__init__(**kwargs)
            self.intensities = intensities
            self.mzs = mzs

    @validate("intensities")
    def _valid_intensities(self, proposal: Proposal) -> List[float]:
        """validate positive values, not empty, and same length as mzs list"""
        value = cast(List[float], proposal["value"])
        if len(value) == 0:
            raise TraitError("length of intensities must be greater than 0")
        if len(value) != len(self.mzs):
            raise TraitError("length of intensities and mzs must be equal")
        if any(x <= 0 for x in value):
            raise TraitError("intensities must be positive")
        return value

    @validate("mzs")
    def _valid_mzs(self, proposal: Proposal) -> List[float]:
        """validate positive values, not empty, and same length as intensities list"""
        value = cast(List[float], proposal["value"])
        if len(value) == 0:
            raise TraitError("length of mzs must be greater than 0")
        if len(value) != len(self.intensities):
            raise TraitError("length of intensities and mzs must be equal")
        if value != sorted(value):
            raise TraitError("mzs values must be monotonically increasing")
        if any(x <= 0 for x in value):
            raise TraitError("mzs values must be positive")
        return value


def str_to_spectrum(spectrum_str: str) -> Spectrum:
    """Converts a spectrum string into a Spectrum class instance"""
    decoded = json.loads(spectrum_str)
    return Spectrum(mzs=decoded[0], intensities=decoded[1])


class MsmsRef(HasTraits):
    # pylint: disable=too-few-public-methods,too-many-instance-attributes
    """one line from msms_refs file"""
    database: str = Unicode()
    id: str = Unicode(default=uuid.uuid4())
    name: str = Unicode()
    spectrum: Spectrum = Instance(klass=Spectrum)
    decimal: np.ushort = Int(default_value=4)
    precursor_mz: np.float64 = Float()
    polarity: str = Unicode()
    adduct: str = Unicode()
    fragmentation_method: str = Unicode()
    collision_energy: str = Unicode()
    instrument: str = Unicode()
    instrument_type: str = Unicode()
    formula: str = Unicode()
    exact_mass: np.float64 = Float()
    inchi_key: str = Unicode()
    inchi: str = Unicode()
    smiles: str = Unicode()

    # pylint: disable=no-self-use,too-many-arguments
    def __init__(
        self,
        name: str,
        spectrum: Spectrum,
        precursor_mz: np.float64,
        polarity: str,
        adduct: str,
        fragmentation_method: str,
        collision_energy: str,
        instrument: str,
        instrument_type: str,
        formula: str,
        exact_mass: np.float64,
        inchi_key: str,
        **kwargs,
    ) -> None:
        """required fields are inputs"""
        with self.hold_trait_notifications():
            super().__init__(**kwargs)
            self.name = name
            self.spectrum = spectrum
            self.precursor_mz = precursor_mz
            self.polarity = polarity
            self.adduct = adduct
            self.fragmentation_method = fragmentation_method
            self.collision_energy = collision_energy
            self.instrument = instrument
            self.instrument_type = instrument_type
            self.formula = formula
            self.exact_mass = exact_mass
            self.inchi_key = inchi_key
        if self.is_bad():
            raise ValueError("MSMS Ref does not pass validation")

    def _valid_enum(self, proposal, name, values_list):
        """generic validation for enumerated type"""
        if proposal["value"] not in values_list:
            raise TraitError(f"{name} must be one of {', '.join(values_list)}")
        return proposal["value"]

    def _valid_not_len_zero(self, proposal, name):
        """generic validation for length greater than 0"""
        if len(proposal["value"]) == 0:
            raise TraitError(f"{name} cannot have a length of zero")
        return proposal["value"]

    def _valid_positive(self, proposal, name):
        """generic validation for positive value"""
        if proposal["value"] < 0:
            raise TraitError(f"{name} must be positive")
        return proposal["value"]

    @validate("database")
    def _valid_database(self, proposal):
        """valid if database string has positive length"""
        return self._valid_not_len_zero(proposal, "database")

    @validate("id")
    def _valid_id(self, proposal):
        """valid if id string has positive length"""
        return self._valid_not_len_zero(proposal, "id")

    @validate("name")
    def _valid_name(self, proposal):
        """valid if name string has positive length"""
        return self._valid_not_len_zero(proposal, "name")

    @validate("decimal")
    def _valid_decimal(self, proposal):
        """valid if decimal is positive"""
        return self._valid_positive(proposal, "decimal")

    @validate("precursor_mz")
    def _valid_precursor_mz(self, proposal):
        """valid if precursor_mz is positive"""
        return self._valid_positive(proposal, "precursor_mz")

    @validate("polarity")
    def _valid_polarity(self, proposal):
        """valid if polarity is in POLARITIES"""
        return self._valid_enum(proposal, "polarity", POLARITIES)

    @validate("adduct")
    def _valid_adduct(self, proposal):
        """valid if adduct string has positive length"""
        return self._valid_not_len_zero(proposal, "adduct")

    @validate("fragmentation_method")
    def _valid_fragmentation_method(self, proposal):
        """valid if fragmentation_method in FRAG_METHODS"""
        return self._valid_enum(proposal, "fragmentation_method", FRAG_METHODS)

    @validate("collision_energy")
    def _valid_collision_energy(self, proposal):
        """valid if collision_energy has positive length"""
        return self._valid_not_len_zero(proposal, "collision_energy")

    @validate("instrument")
    def _valid_instrument(self, proposal):
        """valid if instrument has positive length"""
        return self._valid_not_len_zero(proposal, "instrument")

    @validate("instrument_type")
    def _valid_instrument_type(self, proposal):
        """valid if instrument_type is in INSTRUMENT_TYPES"""
        return self._valid_enum(proposal, "instrument_type", INSTRUMENT_TYPES)

    @validate("formula")
    def _valid_formula(self, proposal):
        """valid if formula has positive length"""
        return self._valid_not_len_zero(proposal, "formula")

    @validate("exact_mass")
    def _valid_exact_mass(self, proposal):
        """valid if exact_mass is positive"""
        return self._valid_positive(proposal, "exact_mass")

    @validate("inchi_key")
    def _valid_inchi_key(self, proposal):
        """valid if inchi_key has positive length"""
        return self._valid_not_len_zero(proposal, "inchi_key")

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

    @default("smiles")
    def _get_default_smiles(self):
        """generate smiles from inchi"""
        return Chem.MolToSmiles(Chem.inchi.MolFromInchi(self.inchi))

    def is_bad(self):
        """
        If returns True, then the inputs are bad, but if returns False do not assume the inputs are good
        returning False only means that there is no evidence the inputs are bad. Conclusively saying
        the inputs are good for unusual chemicals that are not in databases is hard.
        """
        # pylint: disable=too-many-return-statements
        if self.fragmentation_method not in FRAG_METHODS:
            logger.error('Invalid fragmentation method "%s" for %s.', self.fragmentation_method, self.name)
            return True
        if not is_valid_inchi_pair(self.inchi, self.inchi_key):
            logger.error("Invalid inchi/inchi_key pair for %s.", self.name)
            return True
        results = metob.retrieve("compounds", username="*", inchi_key=self.inchi_key)
        if len(results) == 0:
            return False
        ref_compound = results[0]
        if self.formula != ref_compound.formula:
            logger.error(
                'Formula "%s" for %s does not match value "%s" in database.',
                self.formula,
                self.name,
                ref_compound.formula,
            )
            return True
        if not math.isclose(self.exact_mass, ref_compound.mono_isotopic_molecular_weight, rel_tol=1e-9):
            logger.error(
                "Exact mass %s for %s does not match value %s in database.",
                self.exact_mass,
                self.name,
                ref_compound.mono_isotopic_molecular_weight,
            )
            return True
        if not is_synonym(self.name, ref_compound.synonyms):
            logger.error("Inchi_key %s does not contain name %s in database.", self.inchi_key, self.name)
            return True
        return False


def read_msms_refs(file_name: str, sep="\t", **kwargs) -> pd.DataFrame:
    """Read in msms refs from file with correct types"""
    logger.info("Reading in existing references from %s", file_name)
    return pd.read_csv(file_name, sep=sep, dtype=REFS_TYPES, **kwargs)


def get_empty_refs() -> pd.DataFrame:
    """Returns an empty MSMS refs DataFrame with the correct columns and types"""
    return pd.DataFrame(data={k: [] for k, v in REFS_TYPES.items()}).astype(REFS_TYPES)


def valid_refs(refs_df: pd.DataFrame) -> int:
    """Return number of rows that fail validation in refs_df. Info on failures to logger"" """
    return sum([1 if MsmsRef(**row).is_bad() else 0 for row in refs_df.rows()])


def add_cells_for_last_row(refs_sheet: sheet) -> None:
    """Creates cells for the last row. Required to be able to export the row later"""
    row(refs_sheet.rows - 1, [None for _, _ in enumerate(REFS_TYPES.keys())])


def display_refs_edit_ui(
    input_file_name: Optional[str], output_file_name: str, validate_input_file: bool = False
) -> widgets.Box:
    """Create GUI spreadsheet for edited MSMS references"""
    if input_file_name is None:
        old_df = get_empty_refs()
    else:
        logger.info("Reading in existsing references from %s", input_file_name)
        old_df = read_msms_refs(input_file_name)
        if validate_input_file:
            old_pass = valid_refs(old_df)
        logger.info("%s of %s references passed validation.", old_pass, len(old_df))

    sheet1 = ipysheet.sheet(rows=1, columns=len(REFS_TYPES), column_headers=list(REFS_TYPES.keys()))
    add_cells_for_last_row(sheet1)
    auto_populate = widgets.Button(description="Auto-populate")
    add_row = widgets.Button(description="Add row")
    export = widgets.Button(description="Export")

    def on_auto_populate_clicked(_):
        logger.info("Auto-populate button clicked.")

    def on_add_row_clicked(_):
        sheet1.rows += 1
        add_cells_for_last_row(sheet1)

    def on_export_clicked(_):
        old_df = read_msms_refs(input_file_name)
        to_add_df = to_dataframe(sheet1).dropna()
        updated_df = pd.concat([old_df, to_add_df], ignore_index=True)
        updated_df.to_csv(output_file_name, sep="\t", index=False)

    auto_populate.on_click(on_auto_populate_clicked)
    add_row.on_click(on_add_row_clicked)
    export.on_click(on_export_clicked)
    return widgets.VBox([sheet1, widgets.HBox([add_row, auto_populate, export])])
