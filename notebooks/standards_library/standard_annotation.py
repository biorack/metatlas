import pandas as pd
import numpy as np
import os
import glob
from difflib import get_close_matches
from tqdm.notebook import tqdm
import re
import math
import pubchempy as pcp
import base64
import io
import pickle
import subprocess
import shutil
import itertools

import ipywidgets as widgets
from IPython.display import display, clear_output, HTML
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from PIL import Image, ImageDraw, ImageFont

from metatlas.datastructures.groups import group_name
from metatlas.io import feature_tools as ft
from metatlas.tools.cheminfo import inchi_or_smiles_to_molecule, get_precursor_mz
from metatlas.datastructures import metatlas_objects as metob
from metatlas.tools import cheminfo
from metatlas.tools import extract_msms as exms
from metatlas.datastructures.utils import get_atlas
from metatlas.plots.dill2plots import make_atlas_from_spreadsheet
from metatlas.io.metatlas_get_data_helper_fun import make_atlas_df

from matchms.filtering.filter_utils.load_known_adducts import load_known_adducts

from scipy.signal import find_peaks

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors

from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import AllChem, Draw, MolFromSmiles
from rdkit.Chem import rdMolDescriptors as Descriptors
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from pyteomics import mgf

from typing import Dict, List, Tuple, Union, Any, Optional, Set

#####################
### LCMSrun tools ###
#####################

def get_lcmsrun_params(lcmsrun_path: str) -> str:
    """
    Extracts the polarity from the LCMS run file name.

    Args:
        lcmsrun_path (str): Path to the LCMS run file.

    Returns:
        str: Polarity extracted from the file name.
    """
    lcmsrun_basename = os.path.basename(lcmsrun_path)
    lcmsrun_polarity = lcmsrun_basename.split('_')[-2]
    return lcmsrun_polarity

def get_run_num(lcmsrun_path: str) -> int:
    """
    Extracts the run number from the LCMS run file name.

    Args:
        lcmsrun_path (str): Path to the LCMS run file.

    Returns:
        int: Run number extracted from the file name.
    """
    lcmsrun_basename = os.path.basename(lcmsrun_path)
    lcmsrun_run = lcmsrun_basename.split('_')[-1][:-3]
    run_num = int(re.sub(r'\D', '', lcmsrun_run))
    return run_num

def get_closest_injbl(lcmsrun_path: str, injbl_pattern: str = '-InjBL-') -> Optional[str]:
    """
    Retrieve the closest injection blank file before the standard was injected.

    Args:
        lcmsrun_path (str): Path to the LCMS run file.
        injbl_pattern (str, optional): Pattern to identify injection blank files. Defaults to '-InjBL-'.

    Returns:
        Optional[str]: Path to the closest injection blank file, or None if no match is found.
    """
    raw_data_dir = os.path.dirname(lcmsrun_path)
    lcmsrun_num = get_run_num(lcmsrun_path)
    
    all_injbl_files = glob.glob(os.path.join(raw_data_dir, f"*{injbl_pattern}*.h5"))
    injbl_run_nums = {get_run_num(injbl_path): injbl_path for injbl_path in all_injbl_files}
    
    closest_run_num = max((run_num for run_num in injbl_run_nums if run_num <= lcmsrun_num), default=None)

    if closest_run_num is None:
        # No injection blank found before the LCMS run, use the closest one afterwards as a holder
        closest_run_num = min((run_num for run_num in injbl_run_nums if run_num > lcmsrun_num), default=None)

    return injbl_run_nums.get(closest_run_num)


def get_rt_range(lcmsrun_path: str, polarity: str) -> Tuple[float, float]:
    """
    Get the retention time range for a given LCMS run and polarity.

    Args:
        lcmsrun_path (str): Path to the LCMS run file.
        polarity (str): Polarity of the LCMS run (e.g., 'POS' or 'NEG').

    Returns:
        Tuple[float, float]: Minimum and maximum retention times.
    """
    lcmsrun_data = ft.df_container_from_metatlas_file(lcmsrun_path, desired_key=f"ms1_{polarity.lower()}")
    
    rt_min = 0.0
    rt_max = round(lcmsrun_data.rt.max(), 2)
    return rt_min, rt_max

def is_matching_group(lcmsrun_path: str, group: str) -> bool:
    """
    Check if the LCMS run file belongs to the specified group.

    Args:
        lcmsrun_path (str): Path to the LCMS run file.
        group (str): Group name to match.

    Returns:
        bool: True if the LCMS run belongs to the specified group, False otherwise.
    """
    lcmsrun_basename = os.path.basename(lcmsrun_path)
    lcmsrun_group = lcmsrun_basename.split('_')[12]
    return lcmsrun_group == group


def get_file_polarity(lcmsrun_path: str) -> str:
    """
    Extract the polarity from the LCMS run file name.

    Args:
        lcmsrun_path (str): Path to the LCMS run file.

    Returns:
        str: Polarity extracted from the file name.
    """
    lcmsrun_basename = os.path.basename(lcmsrun_path)
    lcmsrun_polarity = lcmsrun_basename.split('_')[9]
    return lcmsrun_polarity


def is_matching_polarity(lcmsrun_path: str, included_polarities: List[str]) -> bool:
    """
    Check if the LCMS run file has a polarity that matches the included polarities.

    Args:
        lcmsrun_path (str): Path to the LCMS run file.
        included_polarities (List[str]): List of polarities to match.

    Returns:
        bool: True if the LCMS run's polarity matches one of the included polarities, False otherwise.
    """
    lcmsrun_polarity = get_file_polarity(lcmsrun_path)
    return lcmsrun_polarity in included_polarities


def get_matching_lcmsruns(row: pd.Series, include_polarities: List[str], include_chromatographies: List[str], raw_data_dir: str) -> List[str]:
    """
    Retrieve LCMS run files that match the specified polarities and chromatographies.

    Args:
        row (pd.Series): Row of a DataFrame containing experiment and group information.
        include_polarities (List[str]): List of polarities to include (e.g., ['POS', 'NEG']).
        include_chromatographies (List[str]): List of chromatographies to include (e.g., ['C18', 'HILIC']).
        raw_data_dir (str): Path to the directory containing raw LCMS data.

    Returns:
        List[str]: List of matching LCMS run file paths.
    """
    all_standard_files = []
    for chrom in include_chromatographies:
        all_chrom_files = glob.glob(os.path.join(raw_data_dir, row[f"{chrom.lower()}_experiment"], '*.h5'))
        standard_files = [path for path in all_chrom_files if is_matching_group(path, row[f"{chrom.lower()}_group"])]
        standard_files = [path for path in standard_files if is_matching_polarity(path, include_polarities)]
        all_standard_files += standard_files
    return all_standard_files

def build_adduct_annotated_table(
    standard_lcmsruns_table: pd.DataFrame,
    include_adducts: List[str] = ['[M+H]+', '[M+Na]+', '[M-H2O+H]+', '[M+K]+', '[M+NH4]+', '[M]+', '[M+2H]2+', '[M-H]-', '[M+Cl]-', '[M-H2O-H]-', '[M]-', '[M-2H]2-']
) -> pd.DataFrame:
    """
    Build a table of LCMS runs annotated with adduct information.

    Args:
        standard_lcmsruns_table (pd.DataFrame): A DataFrame containing LCMS run data.
        include_adducts (List[str], optional): A list of adducts to include in the annotation. Defaults to a predefined list.

    Returns:
        pd.DataFrame: A DataFrame with annotated adduct information, including precursor m/z values.
    """
    standard_lcmsruns_table['polarity'] = standard_lcmsruns_table.apply(lambda row: get_file_polarity(row.standard_lcmsrun), axis=1)
    standard_lcmsruns_table['exact_mass'] = standard_lcmsruns_table.apply(lambda row: inchi_or_smiles_to_mass(row.smiles), axis=1)
    standard_lcmsruns_table['all_adducts'] = standard_lcmsruns_table[['exact_mass', 'polarity']].apply(lambda row: calc_all_adducts(row.exact_mass, row.polarity, include_adducts), axis=1)
    standard_lcmsruns_table = standard_lcmsruns_table.explode('all_adducts').reset_index(drop=True).rename(columns={'all_adducts': 'adduct_data'})
    standard_lcmsruns_table[['adduct', 'precursor_mz']] = pd.DataFrame(standard_lcmsruns_table['adduct_data'].tolist(), index=standard_lcmsruns_table.index)
    return standard_lcmsruns_table

def build_standard_lcmsrun_table(
    config: Dict[str, Any],
    raw_data_dir: str = '/global/cfs/cdirs/metatlas/raw_data/*/'
) -> pd.DataFrame:
    """
    Build a table of LCMS run files that match the specified polarities, chromatographies, and adducts.

    Args:
        config (Dict[str, Any]): Configuration dictionary containing metadata for ref std project
        raw_data_dir (str, optional): Path to the directory containing raw LCMS data. Defaults to '/global/cfs/cdirs/metatlas/raw_data/*/'.

    Returns:
        pd.DataFrame: A DataFrame containing the LCMS run files with annotated adducts.
    """
    standard_info = pd.read_csv(config['project']['standards_input_file'], keep_default_na=False, dtype=str)
    standard_info['standard_lcmsruns'] = standard_info.apply(
        lambda row: get_matching_lcmsruns(row, config['project']['include_polarities'], config['project']['include_chromatographies'], raw_data_dir), axis=1
    )
    standard_lcmsruns_table = standard_info.explode('standard_lcmsruns').reset_index(drop=True).rename(
        columns={'standard_lcmsruns': 'standard_lcmsrun'}
    )
    
    standard_lcmsruns_table_with_adducts = build_adduct_annotated_table(
        standard_lcmsruns_table, include_adducts=config['project']['include_adducts']
    )

    return standard_lcmsruns_table_with_adducts

###########################
### Tools for chemistry ###
###########################

def get_pubchem_synonyms_from_inchi_keys(inchi_keys: List[str]) -> Dict[str, str]:
    """
    Retrieve synonyms for multiple compounds from PubChem using their InChIKeys and return them as strings.

    Args:
        inchi_keys (List[str]): A list of InChIKeys for the compounds.

    Returns:
        Dict[str, str]: A dictionary where keys are InChIKeys and values are strings of synonyms separated by "///".
                        If no synonyms are found, the value will be an empty string.
    """
    # Use PubChemPy to search for compounds by InChIKeys
    results = pcp.get_compounds(inchi_keys, namespace='inchikey')
    synonyms_dict = {}

    for compound in results:
        if compound and compound.inchikey:
            synonyms = compound.synonyms if compound.synonyms else []
            synonyms_dict[compound.inchikey] = "///".join(synonyms)
        else:
            synonyms_dict[compound.inchikey] = ""  # Empty string if no synonyms are found

    return synonyms_dict

def get_pubchem_cids_from_inchi_keys(inchi_keys: List[str]) -> Dict[str, Optional[int]]:
    """
    Retrieve PubChem Compound IDs (CIDs) for multiple compounds using their InChIKeys.

    Args:
        inchi_keys (List[str]): A list of InChIKeys for the compounds.

    Returns:
        Dict[str, Optional[int]]: A dictionary where keys are InChIKeys and values are the corresponding PubChem CIDs.
                                   If no CID is found, the value will be None.
    """
    # Use PubChemPy to search for compounds by InChIKeys
    results = pcp.get_compounds(inchi_keys, namespace='inchikey')
    cid_dict = {}

    for compound in results:
        if compound and compound.inchikey:
            cid_dict[compound.inchikey] = compound.cid  # Assign the CID
        else:
            cid_dict[compound.inchikey] = None  # None if no CID is found

    return cid_dict

def inchi_or_smiles_to_mass(molecule_id: str) -> float:
    """
    Convert an InChI or SMILES string to its exact molecular mass.

    Args:
        molecule_id (str): The InChI or SMILES string of the molecule.

    Returns:
        float: The exact molecular mass of the molecule.
    """
    mol = inchi_or_smiles_to_molecule(molecule_id)
    molweight = ExactMolWt(mol)
    return molweight

def load_adducts_dict() -> Dict[str, Dict[str, Union[str, float]]]:
    """
    Load a dictionary of known adducts with their properties.

    Returns:
        Dict[str, Dict[str, Union[str, float]]]: A dictionary where keys are adduct names and values are dictionaries of adduct properties.
    """
    return load_known_adducts().set_index('adduct').to_dict(orient='index')

def load_selected_adducts(include_adducts: List[str]) -> Dict[str, Dict[str, Union[str, float]]]:
    """
    Filter the known adducts dictionary to include only the specified adducts.

    Args:
        include_adducts (List[str]): A list of adduct names to include.

    Returns:
        Dict[str, Dict[str, Union[str, float]]]: A dictionary of the selected adducts and their properties.
    """
    all_adducts = load_adducts_dict()
    filtered_adducts = {adduct: all_adducts[adduct] for adduct in include_adducts}
    return filtered_adducts

def separate_adducts_dict(all_adducts: Dict[str, Dict[str, Union[str, float]]]) -> Tuple[Dict[str, Dict[str, Union[str, float]]], Dict[str, Dict[str, Union[str, float]]]]:
    """
    Separate adducts into positive and negative ion modes.

    Args:
        all_adducts (Dict[str, Dict[str, Union[str, float]]]): A dictionary of adducts and their properties.

    Returns:
        Tuple[Dict[str, Dict[str, Union[str, float]]], Dict[str, Dict[str, Union[str, float]]]]:
            A tuple containing two dictionaries: one for positive ion mode adducts and one for negative ion mode adducts.
    """
    pos_adducts = {}
    neg_adducts = {}
    for adduct, data in all_adducts.items():
        if data['ionmode'] == 'positive':
            pos_adducts[adduct] = data
        elif data['ionmode'] == 'negative':
            neg_adducts[adduct] = data
    return pos_adducts, neg_adducts

def calc_all_adducts(exact_mass: float, polarity: str, include_adducts: List[str]) -> List[Tuple[str, float]]:
    """
    Calculate the precursor m/z values for all selected adducts based on the exact mass and polarity.

    Args:
        exact_mass (float): The exact molecular mass of the compound.
        polarity (str): The ionization polarity ('POS' or 'NEG').
        include_adducts (List[str]): A list of adducts to include in the calculation.

    Returns:
        List[Tuple[str, float]]: A list of tuples where each tuple contains an adduct name and its corresponding precursor m/z value.
    """
    assert polarity in ['POS', 'NEG'], "Polarity must be 'POS' or 'NEG'."
    
    all_adducts = load_selected_adducts(include_adducts)
    pos_adducts, neg_adducts = separate_adducts_dict(all_adducts)
    
    adduct_pmzs = []
    if polarity == 'POS':
        for adduct in pos_adducts.keys():
            adduct_pmzs.append((adduct, get_precursor_mz(exact_mass, adduct)))
    elif polarity == 'NEG':
        for adduct in neg_adducts.keys():
            adduct_pmzs.append((adduct, get_precursor_mz(exact_mass, adduct)))
    return adduct_pmzs

def inchi_to_inchikey(inchi: str) -> str:
    """
    Convert an InChI string to an InChIKey.

    Args:
        inchi (str): The InChI string of the molecule.

    Returns:
        str: The corresponding InChIKey.
    """
    return AllChem.InchiToInchiKey(inchi)

def neutralize_inchi(inchi: str) -> str:
    """
    Neutralize a molecule represented by an InChI string.

    Args:
        inchi (str): The InChI string of the molecule.

    Returns:
        str: The neutralized InChI string.
    """
    mol = AllChem.MolFromInchi(inchi)
    neutral_mol = cheminfo.normalize_molecule(mol)
    neutralized_inchi = AllChem.MolToInchi(neutral_mol)
    return neutralized_inchi

def charge_from_inchi(inchi: str) -> int:
    """
    Calculate the formal charge of a molecule from its InChI string.

    Args:
        inchi (str): The InChI string of the molecule.

    Returns:
        int: The formal charge of the molecule.
    """
    mol = AllChem.MolFromInchi(inchi)
    charge = AllChem.GetFormalCharge(mol)
    return charge

def formula_from_inchi(inchi: str) -> Optional[str]:
    """
    Get the molecular formula of a molecule from its InChI string.

    Args:
        inchi (str): The InChI string of the molecule.

    Returns:
        Optional[str]: The molecular formula of the molecule, or None if it cannot be determined.
    """
    mol = AllChem.MolFromInchi(inchi)
    formula = Descriptors.CalcMolFormula(mol)
    return formula

def monoisotopic_mass_from_inchi(inchi: str) -> float:
    """
    Calculate the monoisotopic mass of a molecule from its InChI string.

    Args:
        inchi (str): The InChI string of the molecule.

    Returns:
        float: The monoisotopic mass of the molecule.
    """
    mol = AllChem.MolFromInchi(inchi)
    monoisotopic_mass = Descriptors.CalcExactMolWt(mol)
    return monoisotopic_mass

def get_adduct(compound_name: str) -> str:
    """
    Extract the adduct from a compound name.

    Args:
        compound_name (str): The compound name in the format "compound_adduct".

    Returns:
        str: The adduct part of the compound name.
    """
    return compound_name.split('_')[1]

def extract_adducts(eics: Dict[str, pd.DataFrame]) -> Set[str]:
    """
    Extracts unique adducts from EIC data.

    Args:
        eics (Dict[str, pd.DataFrame]): Dictionary of EIC data.

    Returns:
        Set[str]: Set of unique adducts.
    """
    eic_adducts = []
    for eic in eics.values():
        eic['adduct'] = eic.label.apply(lambda x: x.split('_')[-1])
        eic_adducts += eic.adduct.tolist()
    return set(eic_adducts)

def get_compound_name(compound_name: str) -> str:
    """
    Extract the compound name from a compound name with an adduct.

    Args:
        compound_name (str): The compound name in the format "compound_adduct".

    Returns:
        str: The compound name without the adduct.
    """
    return compound_name.split('_')[0]


def smiles_from_inchi_key(inchi_key: str) -> Optional[str]:
    """
    Retrieve the SMILES string for a given InChIKey.

    Args:
        inchi_key (str): The InChIKey of the compound.

    Returns:
        Optional[str]: The SMILES string of the compound, or None if not found.
    """
    compound = pcp.get_compounds(inchi_key, 'inchikey')
    return compound[0].canonical_smiles if compound else None


def display_smiles(smiles: str) -> Optional[Image.Image]:
    """
    Generate an image of a molecule from its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        Optional[Image.Image]: An image of the molecule, or None if the SMILES string is invalid.
    """
    mol = AllChem.MolFromSmiles(smiles)
    if mol:
        return Draw.MolToImage(mol)
    else:
        print("Invalid SMILES string.")
        return None

def enrich_metadata(refs: List[Dict[str, Any]], msms_refs_metadata: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Enrich metadata for a list of references using MS/MS reference metadata.

    Args:
        refs (List[Dict[str, Any]]): A list of reference dictionaries to be enriched.
        msms_refs_metadata (Dict[str, Any]): A dictionary containing metadata fields such as 
                                             'frag_method', 'instrument_type', 'decimal', and 'msms_refs_prefix'.

    Returns:
        List[Dict[str, Any]]: The enriched list of references.
    """
    exms.enrich_metadata(refs, msms_refs_metadata['frag_method'],
                         msms_refs_metadata['instrument_type'], 
                         msms_refs_metadata['decimal'],
                         msms_refs_metadata['msms_refs_prefix'])
    return refs

def get_collision_energy(lcmsrun_path: str) -> str:
    """
    Extract the collision energy from the LCMS run file name.

    Args:
        lcmsrun_path (str): Path to the LCMS run file.

    Returns:
        str: The extracted collision energy as a string.
    """
    if not lcmsrun_path or (isinstance(lcmsrun_path, float) and np.isnan(lcmsrun_path)):
        return np.nan
    lcmsrun = os.path.basename(lcmsrun_path)
    collision_energy = lcmsrun.split('_')[-2].split('-')[1][2:]
    return collision_energy

def make_text_spectrum(spectrum: Any) -> str:
    """
    Convert a spectrum object into a text representation.

    Args:
        spectrum (Any): The spectrum object to be converted.

    Returns:
        str: The text representation of the spectrum.
    """
    return exms.make_text_spectrum(spectrum)

def get_chromatography(lcmsrun_path: str) -> str:
    """
    Extracts the chromatography type from the LCMS run file name.

    Args:
        lcmsrun_path (str): Path to the LCMS run file.

    Returns:
        str: Chromatography type extracted from the file name.
    """
    lcmsrun_basename = os.path.basename(lcmsrun_path)
    lcmsrun_chrom = lcmsrun_basename.split('_')[7]
    return lcmsrun_chrom


#########################
### Metatlas DB tools ###
#########################

def store_in_metatlas_db(cid_not_in_db: Any) -> None:
    """
    Store a compound identification object in the MetAtlas database.

    Args:
        cid_not_in_db (Any): The compound identification object to be stored.
    """
    print("Storing compound in the metatlas database compounds table...")
    metob.store(cid_not_in_db)

def check_db_deposit(new_entries_df: pd.DataFrame) -> None:
    """
    Checks if compounds in the given DataFrame exist in the "Compounds" table of the metatlas database.
    Prints missing entries and displays them as a DataFrame if any are not found.

    Args:
        new_entries_df (pd.DataFrame): DataFrame containing 'inchi_key' and 'label' columns for compounds to check.

    Returns:
        None
    """
    print("Running double check for compounds in metatlas db Compounds table...")
    missing = {}
    all_molecules_check = new_entries_df[new_entries_df['adduct'] != "Ambiguous"]

    for inchi_key, label in zip(all_molecules_check['inchi_key'], all_molecules_check['label']):
        entry = test_metatlas_db_insertion(inchi_key=inchi_key, table="Compounds")
        if not entry:
            #print(f"{label} entry not found in Compounds table for InChIKey {inchi_key}")
            missing[inchi_key] = label
    if missing:
        print("\tSome compounds still missing from database:")
        display(pd.DataFrame.from_dict(missing, orient='index', columns=['label']))
    else:
        print("\tAll new entries found in the database.\n")

def test_metatlas_db_insertion(inchi_key: str, table: str) -> list:
    """
    Test if a compound with the given InChIKey exists in the specified MetAtlas database table.

    Args:
        inchi_key (str): The InChIKey of the compound to search for.
        table (str): The name of the database table to search in.

    Returns:
        list: A list of matching entries from the database.
    """
    return metob.retrieve(table, inchi_key=inchi_key)


#######################
### Save/load tools ###
#######################

def handle_data(
    mode: str, 
    config: Dict[str, Any], 
    timestamp: Optional[str] = None, 
    data: Optional[Tuple[Any, ...]] = None, 
    file_suffix: str = "data"
) -> Optional[Union[Tuple[Any, ...], None]]:
    """
    Handles saving or loading data based on the mode and file suffix.

    Args:
        mode (str): Operation mode, either "save" to save data or "load" to load data.
        config (Dict[str, Any]): Configuration dictionary containing the save path.
        timestamp (Optional[str], optional): Timestamp for saving the file. Required for "save" mode. Defaults to None.
        data (Optional[Tuple[Any, ...]], optional): Data to save. Required for "save" mode. Defaults to None.
        file_suffix (str): Suffix for the file name (e.g., "filtered", "rt_correction", "selected"). Defaults to "data".

    Returns:
        Optional[Union[Tuple[Any, ...], None]]: Loaded data as a tuple if mode is "load", or None if mode is "save".
    """
    data_path = os.path.join(config['project']['standards_output_path'], "cache")
    if mode == "save":
        if not os.path.isdir(data_path):
            os.mkdir(data_path)
        if timestamp is None or data is None:
            raise ValueError("Timestamp and data are required for saving.")
        write_fname = os.path.join(data_path, f"{timestamp}_ref_stds_{file_suffix}.pkl")
        print(f"Saving data to: {write_fname}")
        with open(write_fname, 'wb') as f:
            pickle.dump(data, f)
    elif mode == "load":
        pkl_files = glob.glob(os.path.join(data_path, f"*_ref_stds_{file_suffix}.pkl"))
        try:
            most_recent_pkl = max(pkl_files, key=os.path.getmtime)
        except ValueError:
            print(f"No pkl files found in {data_path}/")
            return None
        print(f"Loading most recent pkl file: {most_recent_pkl}")
        with open(most_recent_pkl, 'rb') as f:
            return pickle.load(f)
    else:
        raise ValueError("Invalid mode. Use 'save' or 'load'.")
    
#############################
### Data extraction tools ###
#############################

def create_compound_atlas(
    group: pd.DataFrame,
    ppm_tolerance: float,
    rt_min: float,
    rt_max: float,
    polarity: str,
    extra_time: float = 0.0
) -> pd.DataFrame:
    """
    Create a compound atlas DataFrame with metadata for each compound.

    Args:
        group (pd.DataFrame): DataFrame containing compound information.
        ppm_tolerance (float): Tolerance for mass-to-charge ratio (ppm).
        rt_min (float): Minimum retention time.
        rt_max (float): Maximum retention time.
        polarity (str): Polarity of the LCMS run (e.g., 'POS' or 'NEG').
        extra_time (float, optional): Additional time to extend retention time range. Defaults to 0.0.

    Returns:
        pd.DataFrame: DataFrame representing the compound atlas.
    """
    atlas = group.rename(columns={'precursor_mz': 'mz'})
    atlas['label'] = group.apply(lambda row: f"{row.compound_name}_{row.adduct}", axis=1)
    atlas['polarity'] = polarity
    atlas['ppm_tolerance'] = ppm_tolerance
    atlas['extra_time'] = extra_time
    atlas_cols = ['label', 'mz', 'rt_min', 'rt_max', 'rt_peak', 'smiles', 'adduct', 'ppm_tolerance', 'polarity', 'extra_time']
    
    atlas['rt_min'] = rt_min
    atlas['rt_max'] = rt_max
    atlas['rt_peak'] = rt_max / 2  # Dummy value
    
    return atlas[atlas_cols]

def collect_eics_and_ms2(group: pd.DataFrame, ppm_tolerance: float) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Collect Extracted Ion Chromatograms (EICs) and MS2 data for a given group of compounds.

    Args:
        group (pd.DataFrame): A DataFrame containing compound information, including LCMS run details.
        ppm_tolerance (float): The mass accuracy tolerance in parts per million (ppm).

    Returns:
        Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], pd.DataFrame]:
            - eics: A dictionary where keys are LCMS run identifiers and values are DataFrames of EICs.
            - ms2_data: A dictionary where keys are LCMS run identifiers and values are DataFrames of MS2 summaries.
            - atlas: A DataFrame representing the compound atlas with metadata for each compound.
    """
    first_row = group.iloc[0]

    lcmsrun_polarity = first_row.polarity
    lcmsrun_params = get_lcmsrun_params(first_row.standard_lcmsrun)
    lcmsrun_chrom = get_chromatography(first_row.standard_lcmsrun)
    lcmsrun_metadata = "{}_{}_{}".format(lcmsrun_chrom, lcmsrun_polarity, lcmsrun_params)
    
    # Set min and max RT for search space
    rt_min, rt_max = get_rt_range(first_row.standard_lcmsrun, lcmsrun_polarity)
    # If provided, narrow the search space based on input csv
    def is_valid(val):
        return val not in [None, '']
    if lcmsrun_chrom.lower() == "hilicz":
        if is_valid(first_row.hilicz_rt_min) and is_valid(first_row.hilicz_rt_max):
            rt_min, rt_max = first_row.hilicz_rt_min, first_row.hilicz_rt_max
        elif not is_valid(first_row.hilicz_rt_min) and is_valid(first_row.hilicz_rt_max):
            rt_max = first_row.hilicz_rt_max
        elif is_valid(first_row.hilicz_rt_min) and not is_valid(first_row.hilicz_rt_max):
            rt_min = first_row.hilicz_rt_min
    elif lcmsrun_chrom.lower() == "c18":
        if is_valid(first_row.c18_rt_min) and is_valid(first_row.c18_rt_max):
            rt_min, rt_max = first_row.c18_rt_min, first_row.c18_rt_max
        elif not is_valid(first_row.c18_rt_min) and is_valid(first_row.c18_rt_max):
            rt_max = first_row.c18_rt_max
        elif is_valid(first_row.c18_rt_min) and not is_valid(first_row.c18_rt_max):
            rt_min = first_row.c18_rt_min

    atlas = create_compound_atlas(group, ppm_tolerance, float(rt_min), float(rt_max), lcmsrun_polarity, extra_time=0.0)
    
    # Get experimental input for reference standards
    files = group.standard_lcmsrun.tolist() + [first_row.closest_injbl]
    experiment_input = ft.setup_file_slicing_parameters(
        atlas, files, base_dir=os.getcwd(), ppm_tolerance=ppm_tolerance, extra_time=0.0, polarity=lcmsrun_polarity.lower()
    )

    eics = {}
    ms2_data = {}
    for file_input in experiment_input:
        data = ft.get_data(file_input, save_file=False, return_data=True, ms1_feature_filter=False)
        adduct_eics = ft.group_duplicates(data['ms1_data'], 'label', make_string=False)
        
        ms2_summary = ft.calculate_ms2_summary(data['ms2_data'])

        if not data['ms1_data'].empty:
            eics[file_input['lcmsrun']] = adduct_eics
        if not ms2_summary.empty:
            ms2_data[file_input['lcmsrun']] = ms2_summary
    
    return eics, ms2_data, atlas

def get_all_ms1_and_ms2(
    eics: Dict[str, pd.DataFrame],
    ms2_data: Dict[str, pd.DataFrame],
    standard_lcmsrun: str,
    theoretical_mzs: Dict[str, float],
    compound_smiles: str,
    adduct_to_polarity: Dict[str, str],
    prominence_percentage: float = 0.25,
    width_threshold: Optional[float] = None,
    distance_threshold: Optional[float] = None
) -> Tuple[List[pd.DataFrame], pd.DataFrame]:
    """
    Extract MS1 retention time (RT) peaks and match them with the closest MS2 spectra for a given LCMS run.

    Args:
        eics (Dict[str, pd.DataFrame]): Extracted ion chromatograms (EICs) for all LCMS runs.
        ms2_data (Dict[str, pd.DataFrame]): MS2 spectra data for all LCMS runs.
        standard_lcmsrun (str): The LCMS run identifier.
        theoretical_mzs (Dict[str, float]): Theoretical m/z values for each adduct.
        compound_smiles (str): SMILES string of the compound.
        adduct_to_polarity (Dict[str, str]): Mapping of adducts to their polarities.
        prominence_percentage (float, optional): Percentage of the maximum intensity to use as the prominence threshold. Defaults to 0.25.
        width_threshold (Optional[float], optional): Minimum width of peaks to consider. Defaults to None.
        distance_threshold (Optional[float], optional): Minimum distance between peaks to consider. Defaults to None.

    Returns:
        Tuple[List[pd.DataFrame], pd.DataFrame]: A list of DataFrames containing MS1 RT peaks and a DataFrame of matched MS2 spectra.
    """
    rt_peaks = []
    if standard_lcmsrun in eics.keys():
        predicted_peaks = []
        for _, eic_row in eics[standard_lcmsrun].iterrows():
            adduct = get_adduct(eic_row['label'])
            chromatography = get_chromatography(standard_lcmsrun)
            polarity = adduct_to_polarity[adduct]
            
            # Sort retention times and intensities
            rt_sort = np.argsort(eic_row['rt'])
            sorted_rt = eic_row['rt'][rt_sort]
            sorted_intensity = eic_row['i'][rt_sort]
            
            # Calculate dynamic prominence threshold
            max_intensity = np.max(sorted_intensity)
            prominence_threshold = max_intensity * prominence_percentage
            
            # Find peaks using scipy.signal.find_peaks with dynamic prominence
            peaks, _ = find_peaks(
                sorted_intensity, 
                prominence=prominence_threshold, 
                width=width_threshold, 
                distance=distance_threshold
            )
            # If no peaks are found, fall back to using intensity method
            if peaks.size == 0:
                peaks = [np.argmax(eic_row['i'][rt_sort])]

            # Filter peaks to retain only the one with the highest intensity if they are within 0.05 RT minutes
            filtered_peaks = []
            for i, peak_index in enumerate(peaks):
                if i == 0 or (sorted_rt[peak_index] - sorted_rt[filtered_peaks[-1]] > 0.05):
                    filtered_peaks.append(peak_index)
                else:
                    # Replace the last peak if the current one has higher intensity
                    if sorted_intensity[peak_index] > sorted_intensity[filtered_peaks[-1]]:
                        filtered_peaks[-1] = peak_index
            
            # Sort filtered peaks by intensity in descending order and keep only the top 5
            filtered_peaks = sorted(filtered_peaks, key=lambda idx: sorted_intensity[idx], reverse=True)[:5]

            # Create a DataFrame for the filtered peaks
            peak_data = []
            for i, peak_index in enumerate(filtered_peaks):
                peak_data.append({
                    'lcmsrun': standard_lcmsrun,
                    'chromatography': chromatography,
                    'compound_name': get_compound_name(eic_row['label']),
                    'adduct': adduct,
                    'polarity': polarity,
                    'rt_peak': sorted_rt[peak_index],
                    'intensity': sorted_intensity[peak_index],
                    'mz_observed': eic_row['mz'][rt_sort][peak_index],
                    'mz_theoretical': theoretical_mzs[adduct],
                    'ppm_error': ((theoretical_mzs[adduct] - eic_row['mz'][rt_sort][peak_index]) / theoretical_mzs[adduct]) * 1e6,
                    'smiles': compound_smiles,
                    'peak_index': f"peak{i+1}"  # Adding unique peak index
                })
            predicted_peaks.extend(peak_data)

    try:
        rt_peaks.append(pd.DataFrame(predicted_peaks))
    except:
        rt_peaks.append(pd.DataFrame())

    top_spectra = []
    if standard_lcmsrun in ms2_data.keys():
        standard_ms2_data = ms2_data[standard_lcmsrun]
        standard_ms2_data['adduct'] = standard_ms2_data['label'].apply(lambda x: x.split('_')[-1])
        standard_ms2_data['lcmsrun'] = standard_lcmsrun

        for rt_peak_df in rt_peaks:
            for _, rt_peak in rt_peak_df.iterrows():
                adduct = rt_peak['adduct']
                rt_value = rt_peak['rt_peak']
                
                # Filter MS2 data for the specific adduct
                standard_ms2_data.loc[:, 'total_intensity'] = standard_ms2_data['spectrum'].apply(lambda x: x[1].sum())
                adduct_ms2_data = standard_ms2_data[standard_ms2_data['adduct'] == adduct]
                
                if not adduct_ms2_data.empty:
                    try:
                        # Find the MS2 spectrum closest to the RT peak and add some info from rt_peak
                        closest_idx = (adduct_ms2_data['rt'] - rt_value).abs().argsort()[:1]
                        closest_spectrum = adduct_ms2_data.iloc[closest_idx].copy()  # Create explicit copy
                        
                        # Check if the closest spectrum is within 0.3 minutes of the RT value
                        if abs(closest_spectrum['rt'].values[0] - rt_value) > 0.3:
                            continue  # Skip this spectrum if it is too far away from peak rt

                        # Use .loc to safely modify the dataframe copy
                        closest_spectrum.loc[:, 'total_intensity_fraction'] = closest_spectrum['total_intensity'] / standard_ms2_data['total_intensity'].max()
                        closest_spectrum.loc[:, 'peak_index'] = rt_peak['peak_index']
                        closest_spectrum.loc[:, 'chromatography'] = rt_peak['chromatography']
                        closest_spectrum.loc[:, 'compound_name'] = rt_peak['compound_name']
                        closest_spectrum.loc[:, 'polarity'] = rt_peak['polarity']
                        
                        top_spectra.append(closest_spectrum)
                    except:
                        continue

    try:
        rt_matched_spectra = pd.concat(top_spectra).reset_index(drop=True)
    except:
        rt_matched_spectra = pd.DataFrame()

    return rt_peaks, rt_matched_spectra


def get_top_ms1_and_ms2(
    eics: Dict[str, pd.DataFrame],
    ms2_data: Dict[str, pd.DataFrame],
    standard_lcmsrun: str,
    theoretical_mzs: Dict[str, float],
    compound_smiles: str,
    adduct_to_polarity: Dict[str, str]
) -> Tuple[List[pd.DataFrame], pd.DataFrame]:
    """
    Extract the most intense MS1 RT peak and the top MS2 spectrum for each adduct in a given LCMS run.

    Args:
        eics (Dict[str, pd.DataFrame]): Extracted ion chromatograms (EICs) for all LCMS runs.
        ms2_data (Dict[str, pd.DataFrame]): MS2 spectra data for all LCMS runs.
        standard_lcmsrun (str): The LCMS run identifier.
        theoretical_mzs (Dict[str, float]): Theoretical m/z values for each adduct.
        compound_smiles (str): SMILES string of the compound.
        adduct_to_polarity (Dict[str, str]): Mapping of adducts to their polarities.

    Returns:
        Tuple[List[pd.DataFrame], pd.DataFrame]: A list of DataFrames containing MS1 RT peaks and a DataFrame of top MS2 spectra.
    """
    rt_peaks = []
    if standard_lcmsrun in eics.keys():
        peak_rows = []
        for _, eic_row in eics[standard_lcmsrun].iterrows():
            adduct = get_adduct(eic_row['label'])
            chromatography = get_chromatography(standard_lcmsrun)
            polarity = adduct_to_polarity[adduct]
            
            rt_peak = {'lcmsrun': standard_lcmsrun, 'chromatography': None, 'compound_name': None, 'adduct': None , 'polarity': None, 'rt_peak': None, 
                       'intensity': None, 'mz_observed': None, 'mz_theoretical': None, 'ppm_error': None, 'smiles': None}
            rt_sort = np.argsort(eic_row['rt'])
            peak = np.argmax(eic_row['i'][rt_sort])
            
            rt_peak['chromatography'] = chromatography
            rt_peak['compound_name'] = get_compound_name(eic_row['label'])
            rt_peak['adduct'] = adduct
            rt_peak['polarity'] = polarity
            rt_peak['rt_peak'] = round(eic_row['rt'][rt_sort][peak], 2)
            rt_peak['intensity'] = eic_row['i'][rt_sort][peak]
            rt_peak['mz_observed'] = eic_row['mz'][rt_sort][peak]
            rt_peak['mz_theoretical'] = theoretical_mzs[adduct]
            rt_peak['ppm_error'] = ((theoretical_mzs[adduct] - rt_peak['mz_observed']) / theoretical_mzs[adduct]) * 1e6
            rt_peak['smiles'] = compound_smiles
            rt_peak['peak_index'] = "peak1"
            
            peak_rows.append(rt_peak)
        if peak_rows:
            rt_peaks.append(pd.DataFrame(peak_rows))
        else:
            rt_peaks.append(pd.DataFrame())
    else:
        rt_peaks.append(pd.DataFrame())

    if standard_lcmsrun in ms2_data.keys():
        standard_ms2_data = ms2_data[standard_lcmsrun]
        standard_ms2_data['total_intensity'] = standard_ms2_data['spectrum'].apply(lambda x: x[1].sum())
        standard_ms2_data['adduct'] = standard_ms2_data['label'].apply(lambda x: x.split('_')[-1])
        standard_ms2_data['lcmsrun'] = standard_lcmsrun

        top_spectra = standard_ms2_data.sort_values('total_intensity', ascending=False).groupby('adduct').head(1).reset_index(drop=True)
        top_spectra['peak_index'] = "peak1"
    else:
        top_spectra = pd.DataFrame()
        
    return rt_peaks, top_spectra


def extract_data(
    lcmsruns_table: pd.DataFrame, 
    config: Dict[str, Any],
    method: str = "intensity"
) -> Tuple[
    List[Dict[str, pd.DataFrame]], 
    List[pd.DataFrame], 
    List[Tuple[str, str, str]], 
    List[pd.DataFrame], 
    List[pd.DataFrame], 
    Dict[int, str]
]:
    """
    Extracts data from LCMS runs, including EICs, MS2 spectra, retention time peaks, and molecular grid images.

    Args:
        lcmsruns_table (pd.DataFrame): A DataFrame containing LCMS run data with columns such as 'compound_name', 
                                       'standard_lcmsrun', 'adduct_data', 'smiles', 'adduct', and 'polarity'.
        method (str, optional): The method to use for peak extraction. Options are "intensity" or "find_peaks". 
                                Defaults to "intensity".
        config (Dict[str, Any]): Configuration dictionary containing parameters for data extraction.

    Returns:
        Tuple[
            List[Dict[str, pd.DataFrame]],  # List of EICs for each group
            List[pd.DataFrame],            # List of top MS2 spectra for each group
            List[Tuple[str, str, str]],    # List of group names (compound_name, standard_lcmsrun, smiles)
            List[pd.DataFrame],            # List of retention time peaks for all groups
            List[pd.DataFrame],            # List of atlas data for each group
            Dict[int, str]                 # Molecular grid images for each run number
        ]
    """
    molecular_grid_image = generate_gridded_molecular_images(lcmsruns_table)

    lcmsruns_table['closest_injbl'] = lcmsruns_table['standard_lcmsrun'].apply(get_closest_injbl)
    grouped_lcmsruns_table = lcmsruns_table.groupby(['compound_name', 'standard_lcmsrun'])
    
    eics_list = []
    top_spectra_list = []
    group_name_list = []
    rt_peak_list = []
    atlas_list = []

    for group_name, group in tqdm(grouped_lcmsruns_table, total=len(grouped_lcmsruns_table), desc=" Extracting data from lcmsruns", unit='compound group'):
        theoretical_mzs = dict(group['adduct_data'].tolist())
        compound_smiles = group['smiles'].iloc[0]
        adduct_to_polarity = dict(zip(group['adduct'].tolist(), group['polarity'].tolist()))
        group_name = (group_name[0], group_name[1], compound_smiles)

        eics, ms2_data, atlas = collect_eics_and_ms2(group, config['project']['ppm_tolerance'])

        if method == "intensity":
            rt_peaks, top_spectra = get_top_ms1_and_ms2(eics, ms2_data, group_name[1], theoretical_mzs, compound_smiles, adduct_to_polarity)
        elif method == "find_peaks":
            rt_peaks, top_spectra = get_all_ms1_and_ms2(eics, ms2_data, group_name[1], theoretical_mzs, compound_smiles, adduct_to_polarity)
        else:
            raise ValueError("Invalid method. Use 'intensity' or 'find_peaks'.")

        eics_list.append(eics)
        top_spectra_list.append(top_spectra)
        group_name_list.append(group_name)
        rt_peak_list.extend(rt_peaks)
        atlas_list.append(atlas)

    return eics_list, top_spectra_list, group_name_list, rt_peak_list, atlas_list, molecular_grid_image

##################################
### Existing data search tools ###
##################################

def search_for_matches_in_metatlas_db(
    all_molecules: pd.DataFrame, 
    check_by_flat: bool = True
) -> Tuple[pd.DataFrame, List]:
    """
    Search for matches of molecules in the MetAtlas database using their InChI keys.

    Args:
        all_molecules (pd.DataFrame): DataFrame containing molecule information with columns such as 'label', 'inchi', and 'inchi_key'.
        check_by_flat (bool, optional): Whether to perform a secondary search using a "flat" InChI key. Defaults to True.

    Returns:
        Tuple[pd.DataFrame, List]: 
            - A DataFrame summarizing matches found in the MetAtlas database.
            - A list of compounds not found in the database, formatted for storage.
    """
    matches_dict: Dict[str, List[Union[str, List[str]]]] = {}
    nonmatches_dict: Dict[str, pd.Series] = {}
    all_molecules_keep = all_molecules[all_molecules['adduct'] != "Ambiguous"]

    for _, molecule in tqdm(all_molecules_keep.iterrows(), total=len(all_molecules_keep), desc=" Searching for matches in metatlas db", unit=" compound"):
        
        molecule_subset = molecule[['label', 'inchi', 'inchi_key']]
        if pd.notna(molecule_subset['inchi_key']) is False:
            pass
        
        inchi_key_parts = molecule_subset.inchi_key.split('-')
        inchi_key_parts[1] = '%'
        flat_inchi_key = '-'.join(inchi_key_parts)
        
        db_entry = metob.retrieve('compound', inchi_key=molecule_subset.inchi_key)
        
        if db_entry == []:
            if check_by_flat is True:
                flat_entry = metob.retrieve('compound', inchi_key=flat_inchi_key)
                if flat_entry == []:
                    nonmatches_dict[molecule_subset.label] = molecule
                else:
                    matches_dict[molecule_subset.label] = ["flat_inchi_key", flat_inchi_key, list({entry.inchi_key for entry in flat_entry})]
            else:
                nonmatches_dict[molecule_subset.label] = molecule
        else:
            matches_dict[molecule_subset.label] = ["inchi_key", molecule_subset.inchi_key, list({entry.inchi_key for entry in db_entry})]

    # Convert matches dictionary to DataFrame
    if len(matches_dict) != 0:
        matches_df = pd.DataFrame(
            [(key, value[0], value[1], value[2]) for key, value in matches_dict.items()],
            columns=['query_label', 'query_matching_criterion', 'query_to_db', 'db_match']
        )
        print("\nSummary of compounds already in the metatlas database:\n")
        display(matches_df)
    else:
        print("\nNo compounds were found in the metatlas database.\n")
        matches_df = pd.DataFrame()
    if len(nonmatches_dict) != 0:
        nonmatches_df = pd.concat(nonmatches_dict.values(), axis=1).T.reset_index(drop=True)
        attributes_to_save = ['label', 'inchi', 'inchi_key', 'neutralized_inchi', 'neutralized_inchi_key', 'permanent_charge', 'formula', 'monoisotopic_mass']
        nonmatches_df = nonmatches_df[attributes_to_save].drop_duplicates()
        nonmatches_df.drop_duplicates(inplace=True)
        print("\nThese compounds are not yet in the metatlas database:\n")
        display(nonmatches_df)
        nonmatches_list = format_for_atlas_store(nonmatches_df)
    else:
        print("\nAll compounds are already in the metatlas database.\n")
        nonmatches_list = []

    return matches_df, nonmatches_list

def search_for_matches_in_msms_refs(
    all_molecules: pd.DataFrame, 
    msms_refs: pd.DataFrame, 
    check_by_flat: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Search for matches of molecules in MS/MS reference data using their InChI keys and metadata.

    Args:
        all_molecules (pd.DataFrame): DataFrame containing molecule information with columns such as 'name', 'adduct', 'inchi', and 'inchi_key'.
        msms_refs (pd.DataFrame): DataFrame containing MS/MS reference data.
        check_by_flat (bool, optional): Whether to perform a secondary search using a "flat" InChI key. Defaults to True.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: 
            - A DataFrame summarizing matches found in the MS/MS references.
            - A DataFrame summarizing compounds not found in the references.
    """
    matches_dict: Dict[str, List[Union[str, List[str]]]] = {}
    nonmatches_dict: Dict[str, pd.Series] = {}

    for _, molecule in tqdm(all_molecules.iterrows(), total=len(all_molecules), desc=" Searching for matches in MSMS refs", unit=" compound"):
        
        molecule_subset = molecule[['name', 'adduct', 'inchi', 'inchi_key', 
                                    'database', 'instrument', 'collision_energy', 'polarity', 
                                    'fragmentation_method', 'instrument_type']]
        if pd.notna(molecule_subset['inchi_key']) is False:
            pass
        
        inchi_key_parts = molecule_subset.inchi_key.split('-')
        inchi_key_parts[1] = '%'
        flat_inchi_key = '-'.join(inchi_key_parts)

        # Search for exact matches with all criteria
        matching_refs = msms_refs.loc[
            (msms_refs['inchi_key'] == molecule_subset.inchi_key) &
            (msms_refs['adduct'] == molecule_subset.adduct) &
            (msms_refs['database'] == molecule_subset.database) &
            (msms_refs['instrument'] == molecule_subset.instrument) &
            (msms_refs['collision_energy'] == molecule_subset.collision_energy) &
            (msms_refs['polarity'] == molecule_subset.polarity) &
            (msms_refs['fragmentation_method'] == molecule_subset.fragmentation_method) &
            (msms_refs['instrument_type'] == molecule_subset.instrument_type)
        ]
        
        if matching_refs.empty:
            if check_by_flat is True:
                # Search using flat InChI key
                flat_matching_refs = msms_refs.loc[
                    (msms_refs['inchi_key'] == flat_inchi_key) &
                    (msms_refs['adduct'] == molecule_subset.adduct) &
                    (msms_refs['database'] == molecule_subset.database) &
                    (msms_refs['instrument'] == molecule_subset.instrument) &
                    (msms_refs['collision_energy'] == molecule_subset.collision_energy) &
                    (msms_refs['polarity'] == molecule_subset.polarity) &
                    (msms_refs['fragmentation_method'] == molecule_subset.fragmentation_method) &
                    (msms_refs['instrument_type'] == molecule_subset.instrument_type)
                ]
                if flat_matching_refs.empty:
                    nonmatches_dict[molecule_subset.name] = molecule
                else:
                    matches_dict[molecule_subset.name] = [
                        molecule_subset.adduct, "flat_inchi_key", flat_inchi_key, flat_matching_refs.iloc[0].inchi_key
                    ]
            else:
                nonmatches_dict[molecule_subset.name] = molecule
        else:
            matches_dict[molecule_subset.name] = [
                molecule_subset.adduct, "inchi_key", molecule_subset.inchi_key, matching_refs.iloc[0].inchi_key
            ]

    # Convert matches dictionary to DataFrame
    matches_df = pd.DataFrame(
        [(key, value[0], value[1], value[2], value[3]) for key, value in matches_dict.items()],
        columns=['query_name', 'query_adduct', 'query_matching_criterion', 'query_to_refs', 'msms_match']
    )
    nonmatches_df = pd.concat(nonmatches_dict.values(), axis=1).T.reset_index(drop=True)
    nonmatches_df.drop_duplicates(inplace=True)

    # if not matches_df.empty:
    #     print("\nSummary of compounds+adducts already in MSMS refs:\n")
    #     display(matches_df)
    if not nonmatches_df.empty:
        print(f"\n{nonmatches_df.shape[0]} compounds+adducts are not yet in MSMS refs. Check notin_msms_refs to view.\n")

    return matches_df, nonmatches_df

def search_for_matches_in_atlases(
    query_entries: pd.DataFrame, 
    ema_atlases_dfs: Dict[str, Dict[str, Union[str, pd.DataFrame]]], 
    cutoff: float = 0.8
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Search for matches of query entries in existing atlases using exact and fuzzy matching.

    Args:
        query_entries (pd.DataFrame): DataFrame containing query entries with columns such as 'label', 'inchi', 'inchi_key', and 'adduct'.
        ema_atlases_dfs (Dict[str, Dict[str, Union[str, pd.DataFrame]]]): Nested dictionary with chromatography and polarity keys, 
            containing atlas DataFrames.
        cutoff (float, optional): Similarity cutoff for fuzzy matching. Defaults to 0.8.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: 
            - A DataFrame summarizing matches found in the atlases.
            - A DataFrame summarizing compounds not found in the atlases.
    """
    matches_dict: Dict[str, List[Union[str, List[str]]]] = {}
    nonmatches_dict: Dict[str, pd.Series] = {}

    for chrom, polarities in ema_atlases_dfs.items():
        for polarity, atlas_data in polarities.items():
            if isinstance(atlas_data, pd.DataFrame):  # Ensure we are working with a DataFrame
                atlas_data['label'] = atlas_data['label'].astype(str)
                if 'compound_name' in atlas_data.columns:
                    atlas_data['compound_name'] = atlas_data['compound_name'].astype(str)

                relevant_query_entries = query_entries[(query_entries['polarity'] == polarity) & (query_entries['chromatography'] == chrom)]
                if relevant_query_entries.empty:
                    print(f"No compounds need to be searched against the {chrom} {polarity} atlas. Skipping.")
                    continue
                
                for _, query_row in tqdm(relevant_query_entries.iterrows(), total=relevant_query_entries.shape[0], desc=f" Searching in {chrom} {polarity} atlas", unit=" compound"):
                    query_label = str(query_row['label'])
                    query_compound_name = str(query_row['compound_name'])
                    query_inchikey = str(query_row['inchi_key'])
                    query_inchi = str(query_row['inchi'])
                    query_adduct = str(query_row['adduct'])
                    query_polarity = str(query_row['polarity'])
                    query_chromatography = str(query_row['chromatography'])
                    query_unique_id = f"{query_label};;{query_adduct};;{query_polarity};;{query_chromatography}"

                    # Check for exact matches in inchi
                    matching_rows = atlas_data[(atlas_data['inchi'] == query_inchi) & (atlas_data['adduct'] == query_adduct)]
                    if not matching_rows.empty:            
                        source_atlases = matching_rows['source_atlas'].tolist()
                        atlas_entry = matching_rows['inchi'].tolist()
                        matches_dict[query_unique_id] = [query_inchi, atlas_entry, source_atlases]
                        continue

                    # Check for exact matches in inchi_key
                    matching_rows = atlas_data[(atlas_data['inchi_key'] == query_inchikey) & (atlas_data['adduct'] == query_adduct)]
                    if not matching_rows.empty:
                        source_atlases = matching_rows['source_atlas'].tolist()
                        atlas_entry = matching_rows['inchi_key'].tolist()
                        matches_dict[query_unique_id] = [query_inchi, atlas_entry, source_atlases]
                        continue

                    # Check for fuzzy matches in compound_name
                    if 'compound_name' in atlas_data.columns:
                        compound_name_matches = get_close_matches(query_compound_name, atlas_data['compound_name'].tolist(), n=1, cutoff=cutoff)
                        if compound_name_matches:
                            matching_rows = atlas_data[(atlas_data['compound_name'] == compound_name_matches[0]) & (atlas_data['adduct'] == query_adduct)]
                            if not matching_rows.empty:
                                source_atlases = matching_rows['source_atlas'].tolist()
                                atlas_entry = matching_rows['compound_name'].tolist()
                                matches_dict[query_unique_id] = [compound_name_matches[0], atlas_entry, source_atlases]
                                continue

                    # Check for fuzzy matches in label
                    label_matches = get_close_matches(query_label, atlas_data['label'].tolist(), n=1, cutoff=cutoff)
                    if label_matches:
                        matching_rows = atlas_data[(atlas_data['label'] == label_matches[0]) & (atlas_data['adduct'] == query_adduct)]
                        if not matching_rows.empty:
                            source_atlases = matching_rows['source_atlas'].tolist()
                            atlas_entry = matching_rows['label'].tolist()
                            matches_dict[query_unique_id] = [label_matches[0], atlas_entry, source_atlases]
                            continue

                    ## Fingerprint similarity?

                    # If no match is found, add to nonmatches_dict
                    nonmatches_dict[query_unique_id] = query_row

    # Convert matches dictionary to DataFrame
    if len(matches_dict) != 0:
        matches_df = pd.DataFrame(
            [(key, value[0], value[1], value[2]) for key, value in matches_dict.items()],
            columns=['query_unique_id', 'query_to_atlas', 'atlas_matches', 'atlas_source']
        )
        print("\nSummary of compounds+adducts already in the atlases:\n")
        display(matches_df)
    else:
        print("\nNone of the compounds+adducts searched were found in the atlases.")
        matches_df = pd.DataFrame()
    if len(nonmatches_dict) != 0:
        nonmatches_df = pd.concat(nonmatches_dict.values(), axis=1).T.reset_index(drop=True)
        print(f"\nThere are {len(nonmatches_df)} compounds+adducts are not yet in any atlases. View with 'nonmatches_to_atlases'.\n")
        #display(nonmatches_df)
    else:
        print("\nAll compounds+adducts are already in the atlases.\n")
        nonmatches_df = pd.DataFrame()

    return matches_df, nonmatches_df

#####################
### General tools ###
#####################

def filter_by_selected(
    eics_full: List[Dict[str, pd.DataFrame]],
    rt_peaks_full: List[pd.DataFrame],
    top_spectra_full: List[pd.DataFrame],
    selection_results_dict: Dict[str, Tuple[List[str], List[str], List[bool]]],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Filters EICs, RT peaks, and MS2 spectra based on selected adducts and compounds.

    Args:
        eics_full (List[Dict[str, pd.DataFrame]]): List of dictionaries containing EIC data.
        rt_peaks_full (List[pd.DataFrame]): List of DataFrames containing RT peak data.
        top_spectra_full (List[pd.DataFrame]): List of DataFrames containing MS2 spectra data.
        selection_results_dict Dict[str, Tuple[List[str], List[str], List[bool]]]: A dictionary mapping unique compound IDs to tuples
            containing a list of selected adduct-peak combinations and a list of booleans indicating whether each adduct is the best choice.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: Filtered EICs, RT peaks, and MS2 spectra DataFrames.
    """

    eics_selected = pd.concat(
        [df.assign(key=key) for d in eics_full for key, df in d.items()],
        ignore_index=True
    ).rename(columns={'key': 'standard_lcmsrun'})
    eics_selected['compound_name'] = eics_selected['label'].apply(lambda x: x.split('_')[0])
    eics_selected['adduct'] = eics_selected['label'].apply(lambda x: x.split('_')[1])
    eics_selected = select_compounds_from_gui(eics_selected, selection_results_dict)

    rt_peaks_selected = pd.concat(rt_peaks_full).rename(columns={'lcmsrun': 'standard_lcmsrun'})
    rt_peaks_selected = select_compounds_from_gui(rt_peaks_selected, selection_results_dict)

    top_spectra_selected = pd.concat(top_spectra_full, ignore_index=True).rename(columns={'lcmsrun': 'standard_lcmsrun'})
    top_spectra_selected['compound_name'] = top_spectra_selected['label'].apply(lambda x: x.split('_')[0])
    top_spectra_selected = select_compounds_from_gui(top_spectra_selected, selection_results_dict)

    all_data = {"eics": eics_selected,
                "rt_peaks": rt_peaks_selected,
                "top_spectra": top_spectra_selected}

    print(f"\nAll compounds selected: {all_data['eics']['compound_name'].nunique()}")
    print(f"All compound+adduct entries selected: {all_data['eics']['label'].nunique()}")
    print(f"All compound+adduct+peak entries selected: {all_data['rt_peaks'].shape[0]}")
    print(f"All compound+adduct+peak+spectra entries selected: {all_data['top_spectra'].shape[0]}")

    best_data = {"eics": eics_selected[eics_selected['best_adduct'] == True],
                "rt_peaks": rt_peaks_selected[rt_peaks_selected['best_adduct'] == True],
                "top_spectra": top_spectra_selected[top_spectra_selected['best_adduct'] == True]}

    print(f"\nBest compounds selected: {best_data['eics']['compound_name'].nunique()}")
    print(f"Best compound+adduct entries selected: {best_data['eics']['label'].nunique()}")
    print(f"Best compound+adduct+peak entries selected: {best_data['rt_peaks'].shape[0]}")
    print(f"Best compound+adduct+peak+spectra entries selected: {best_data['top_spectra'].shape[0]}")

    return all_data


def format_rt_peaks(
    rt_peaks: pd.DataFrame,
) -> pd.DataFrame:
    """
    Format and validate retention time (RT) peak data for downstream analysis.

    This function processes a dataset of RT peaks (`rt_peaks`) to:
    - Add unique labels for peaks with multiple indices.
    - Check for consistency in retention times (RTs) across collision energies (CEs) and polarities.
    - Identify and handle isomers based on monoisotopic mass.
    - Select the best collision energy row by intensity for the top adduct(s) per compound.
    - Retain the "best_adduct" column indicating best (True) or not best (False).

    Args:
        rt_peaks (pd.DataFrame): A DataFrame containing RT peaks with required columns:
            ['compound_name', 'standard_lcmsrun', 'adduct', 'rt_peak', 'collision_energy', 
             'polarity', 'monoisotopic_mass', 'intensity', 'best_adduct'].

    Returns:
        pd.DataFrame: The formatted RT peaks DataFrame, retaining the "best_adduct" column.
    """
    print("\nFormatting RT peaks dataset")

    rt_data_formatted = rt_peaks.copy()
    rt_data_formatted['label'] = rt_data_formatted['compound_name']

    # Put peak index into label if there are multiple peaks selected for a given adduct
    rt_data_formatted_labeled_list = []
    rt_data_formatted_grouped = rt_data_formatted.groupby(['standard_lcmsrun', 'label', 'adduct'])
    for _, group in rt_data_formatted_grouped:
        if group.shape[0] > 1:
            group = group.copy()
            group['label'] = group.apply(
                lambda row: f"{row['label']} ({row['peak_index']})", axis=1
            )
        rt_data_formatted_labeled_list.append(group)
    rt_data_formatted_labeled = pd.concat(rt_data_formatted_labeled_list, ignore_index=True)

    print("\tChecking for differing RTs between CEs and polarities, which are unexpected...")
    dataset_grouped = rt_data_formatted_labeled.groupby(['chromatography', 'label'])
    for group_name, group in dataset_grouped:
        rt_values = group['rt_peak']
        cutoff = 0.05
        ces = group['collision_energy'].unique()
        pols = group['polarity'].unique()
        diff = rt_values.max() - rt_values.min()
        if diff <= cutoff or np.isnan(diff):
            pass
        else:
            print(f"\t\tWarning! Group {group_name}: RT values for {ces} and {pols} are NOT within {cutoff} mins of each other ({round(rt_values.max(),4)} - {round(rt_values.min(),4)} = {round(diff, 4)}).")

    print("\n\tChecking monoisotopic mass of all BEST compounds+adducts to identify isomers...")
    rt_data_formatted_labeled_best = rt_data_formatted_labeled[rt_data_formatted_labeled['best_adduct'] == True].copy()
    dataset_grouped = rt_data_formatted_labeled_best.groupby(['monoisotopic_mass', 'polarity', 'chromatography'])
    grouped_compounds = dataset_grouped['compound_name'].nunique()
    multiple_compounds_per_mim = grouped_compounds[grouped_compounds > 1]
    if not multiple_compounds_per_mim.empty:
        for isomer_mim in multiple_compounds_per_mim.index:
            isomer_data = rt_data_formatted_labeled_best[
                (rt_data_formatted_labeled_best['monoisotopic_mass'] == isomer_mim[0]) &
                (rt_data_formatted_labeled_best['polarity'] == isomer_mim[1]) &
                (rt_data_formatted_labeled_best['chromatography'] == isomer_mim[2])
            ]
            unique_adducts = isomer_data['adduct'].unique()
            if len(unique_adducts) == 1:
                print(f"\t\tNote: Found isomers in {isomer_mim[2]} {isomer_mim[1]} mode at {isomer_mim[0]} ({list(isomer_data['label'])}) but they had matching selected adducts {unique_adducts[0]}.\n")
            else:
                print(f"\t\tWarning! Adducts for isomers do not agree. See data for monoisotopic mass {isomer_mim[0]}:\n")
                display(isomer_data[['label', 'adduct', 'inchi', 'monoisotopic_mass']])
                #print("\n\t\tPlease return to the GUI to select a matching adduct for isomers.\n")
    else:
        print(f"\t\tNo isomers found in RT peaks data.\n")

    print("Selecting best collision energy row by intensity for the best adduct(s) per compound...")
    def mark_best_adduct(group):
        # If all adducts are Ambiguous, just set the first as best
        if (group['adduct'] == "Ambiguous").all():
            group['best_adduct'] = False
            group.iloc[0, group.columns.get_loc('best_adduct')] = True
            return group

        # Otherwise, only consider rows where best_adduct is True (user-selected)
        mask = (group['best_adduct'] == True)
        if mask.any():
            intensities = group.loc[mask, 'intensity']
            if intensities.notna().any():
                idx_max = intensities.idxmax()
                group.loc[mask, 'best_adduct'] = False
                group.loc[idx_max, 'best_adduct'] = True
                for idx in group.loc[mask].index:
                    if idx != idx_max:
                        ce = group.loc[idx, 'collision_energy']
                        compound = group.loc[idx, 'compound_name']
                        adduct = group.loc[idx, 'adduct']
                        print(f"\tRemoving CE:{ce} from best adduct selection for {compound} ({adduct})")
            else:
                group.loc[mask, 'best_adduct'] = False
        return group

    group_cols = ['chromatography', 'polarity', 'label', 'adduct']
    rt_data_formatted_labeled_marked = (
        rt_data_formatted_labeled
        .groupby(group_cols, group_keys=False)
        .apply(mark_best_adduct)
        .reset_index(drop=True)
    )

    print(f"\nRT peaks dataset: {rt_data_formatted_labeled_marked.shape[0]} total compound peaks.\n")

    return rt_data_formatted_labeled_marked


def format_and_select_best_adducts(
    rt_peaks: pd.DataFrame, 
    best_adducts: Dict[str, List[str]]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Format and filter retention time (RT) peaks to select the top adducts for each compound.

    Args:
        rt_peaks (pd.DataFrame): A DataFrame containing retention time peaks with columns such as 
                                 'compound_name', 'standard_lcmsrun', 'adduct', 'rt_peak', etc.
        best_adducts (Dict[str, List[str]]): A dictionary where keys are compound identifiers 
                                            (formatted as "compound_name;;standard_lcmsrun") and 
                                            values are lists of top adducts to select.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: 
            - The unfiltered RT peaks DataFrame with additional formatting.
            - A filtered DataFrame containing only the top adducts for each compound.
    """
    # Create a copy of the original dataframe for processing
    unfiltered_rt_peaks = rt_peaks.copy()
    unfiltered_rt_peaks['label'] = unfiltered_rt_peaks['compound_name']

    # Initialize lists to store selected and unfiltered data
    selected_best_adducts = []
    unfiltered_best_adducts = []

    # Process rows matching the best_adducts dictionary
    for key, value in best_adducts.items():
        label = key.split(';;')[0]
        standard_lcmsrun = key.split(';;')[1]
        
        # Filter rows matching the current key and adducts
        selected_adduct = unfiltered_rt_peaks[
            (unfiltered_rt_peaks['standard_lcmsrun'] == standard_lcmsrun) &
            (unfiltered_rt_peaks['label'] == label) &
            (unfiltered_rt_peaks['adduct'].isin(value))
        ].copy()
        
        # Handle multiple peaks per top adduct for the filtered dataframe
        if selected_adduct.shape[0] > 1:
            for adduct in selected_adduct['adduct'].unique():
                adduct_subset = selected_adduct[selected_adduct['adduct'] == adduct].copy()
                if adduct_subset.shape[0] > 1:
                    display(adduct_subset)
                    adduct_subset['label'] = adduct_subset.apply(
                        lambda row: f"{row['label']} ({row['peak_index']})", axis=1
                    )
                selected_best_adducts.append(adduct_subset)
        else:
            selected_best_adducts.append(selected_adduct)

    # Handle multiple peaks per top adduct for the unfiltered dataframe
    unfiltered_rt_peaks_grouped = unfiltered_rt_peaks.groupby(['standard_lcmsrun', 'label', 'adduct'])
    for _, group in unfiltered_rt_peaks_grouped:
        if group.shape[0] > 1:
            group = group.copy()
            group['label'] = group.apply(
                lambda row: f"{row['label']} ({row['peak_index']})", axis=1
            )
        unfiltered_best_adducts.append(group)

    # Combine selected rows into a dataframe
    if selected_best_adducts:
        selected_best_adducts_df = pd.concat(selected_best_adducts, ignore_index=True)
    else:
        selected_best_adducts_df = pd.DataFrame()
    # Combine unfiltered rows into a dataframe
    if unfiltered_best_adducts:
        unfiltered_adducts_df = pd.concat(unfiltered_best_adducts, ignore_index=True)
    else:
        unfiltered_adducts_df = pd.DataFrame()

    datasets_dict = {'top': selected_best_adducts_df, 'all': unfiltered_adducts_df}
    
    for dataset_name, dataset in datasets_dict.items():
        print(f"\nWorking on dataset: {dataset_name}")
        print("\tChecking for differing RTs between CEs and polarities, which are unexpected...")
        dataset_grouped = dataset.groupby(['chromatography', 'label'])
        for group_name, group in dataset_grouped:
            rt_values = group['rt_peak']
            cutoff = 0.05
            ces = group['collision_energy'].unique()
            pols = group['polarity'].unique()
            if rt_values.max() - rt_values.min() <= cutoff:
                print(f"\t\tGroup {group_name}: All RT values for {ces} and {pols} are within {cutoff} mins of each other ({round(rt_values.max() - rt_values.min(), 4)}).")
            else:
                print(f"\t\tGroup {group_name}: RT values for {ces} and {pols} are NOT within {cutoff} mins of each other ({round(rt_values.max() - rt_values.min(), 4)}).")

        print("\n\tGrouping by monoisotopic_mass and identify isomers in the datasets...")
        dataset_grouped = dataset.groupby(['monoisotopic_mass', 'polarity', 'chromatography'])
        grouped_compounds = dataset_grouped['compound_name'].nunique()
        multiple_compounds_per_mim = grouped_compounds[grouped_compounds > 1]
        if not multiple_compounds_per_mim.empty:
            for isomer_mim in multiple_compounds_per_mim.index:
                isomer_data = unfiltered_rt_peaks[
                    (unfiltered_rt_peaks['monoisotopic_mass'] == isomer_mim[0]) &
                    (unfiltered_rt_peaks['polarity'] == isomer_mim[1]) &
                    (unfiltered_rt_peaks['chromatography'] == isomer_mim[2])
                ]
                unique_adducts = isomer_data['adduct'].unique()
                if len(unique_adducts) == 1:
                    print(f"\t\tNote: Found isomers in {isomer_mim[2]} {isomer_mim[1]} mode at {isomer_mim[0]} ({list(isomer_data['label'])}) but they had matching selected adducts {unique_adducts[0]}.")
                else:
                    print(f"\t\tWarning! Adducts for isomers do not agree. See data for monoisotopic mass {isomer_mim[0]}:\n")
                    display(isomer_data[['label', 'adduct', 'inchi', 'monoisotopic_mass']])
                    print("\n\t\tPlease return to the GUI to select a matching adduct for isomers.\n")
                    return None, None
        else:
            print(f"\t\tNo isomers found in {dataset_name} data.\n")

        if dataset_name == 'top':
            print("\tSelecting best collision energy row by intensity for the top adduct(s) per compound...")
            selected_best_adducts_df_grouped = selected_best_adducts_df.groupby(['chromatography', 'polarity', 'label', 'adduct'])
            selected_best_adducts_df_grouped_best_ces = []
            for group_name, group in selected_best_adducts_df_grouped:
                ces = group.sort_values(by='intensity', ascending=False)
                if len(ces) >= 1:
                    selected_best_adducts_df_grouped_best_ces.append(ces.iloc[0])
                    print(f"\t\tSelected 1 row and removed {ces.shape[0] - 1} row(s) for {group_name}.")
                else:
                    print(f"\t\tWarning! No collision energy found for {group_name}.")
            selected_best_adducts_df_best_ces = pd.DataFrame(selected_best_adducts_df_grouped_best_ces)

    print(f"'All' peaks dataset: {unfiltered_rt_peaks.shape[0]} total compound peaks.")
    print(f"'Best' peaks dataset: {selected_best_adducts_df_best_ces.shape[0]} best compound peaks.\n")

    # Return both the unfiltered dataframe and the subsetted dataframe
    return unfiltered_rt_peaks, selected_best_adducts_df_best_ces



#######################
### MSMS refs tools ###
#######################

def get_scan_number(file_path: str) -> int:
    """
    Extract the scan number from the last line of the file containing 'scan=' and 'spectrum id'.

    Args:
        file_path (str): Path to the file to process.

    Returns:
        int: The extracted scan number.
    """
    try:
        # Run the shell command
        command = f"cat {file_path} | grep 'scan=' | grep 'spectrum id' | tail -n 1 | sed -n 's/.*scan=\\([0-9]*\\).*/\\1/p'"
        result = subprocess.check_output(command, shell=True, text=True).strip()
        
        # Convert the result to an integer
        scan_number = int(result)
        return scan_number
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
    except ValueError:
        print("No valid scan number found.")
    return None

def get_mgf_last_index(mgf_refs_path: str) -> str:
    """
    Index an MGF file and read specific spectra efficiently.

    Args:
        mgf_refs_path (str): Path to the MGF file.

    Returns:
        IndexedMGF: An indexed MGF object for efficient access.
    """
    try:
        # Create an IndexedMGF object
        indexed_mgf = mgf.IndexedMGF(mgf_refs_path, index_by_scans=True)
        
        # Access the last spectrum in the file
        last_spectrum = indexed_mgf[-1]

        # Access a spectrum by its SPECTRUMID (if available)
        if 'spectrumid' in last_spectrum['params']:
            id = last_spectrum['params']['spectrumid']
            return id
        else:
            print("Warning: SPECTRUMID parameter not found in the last spectrum. Check MGF file/path/format.")
    
    except FileNotFoundError:
        print(f"Error: File not found at {mgf_refs_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

def write_mgf_from_top_spectra(top_spectra_all: pd.DataFrame, 
                               rt_peaks_all: pd.DataFrame, 
                               config: Dict[str, Any],
                               timestamp: str
) -> None:
    """
    Write an MGF file from the spectra data in top_spectra_all.

    Args:
        top_spectra_all (pd.DataFrame): DataFrame containing spectra data.
        rt_peaks_all (pd.DataFrame): DataFrame containing retention time peak data.
        config (Dict[str, Any]): Configuration dictionary containing metadata for the project, including output path.
        timestamp (str): Timestamp for naming the output file.

    Returns:
        None
    """
    mgf_entries = []

    rt_peaks_and_top_spectra = rt_peaks_all.merge(top_spectra_all[['compound_name', 'adduct', 'polarity', 'standard_lcmsrun', 'peak_index', 'rt', 'precursor_mz', 'precursor_peak_height', 'spectrum']],
                                                on=['compound_name', 'adduct', 'polarity', 'standard_lcmsrun', 'peak_index'], how='left')

    if config['msms_refs']['standalone_mgf'] is True:
        current_id = 0
    elif config['msms_refs']['standalone_mgf'] is False:
        starting_mgf_id = get_mgf_last_index(config['msms_refs']['current_mgf_path'])
        print(f"Last ID from existing file to increment: {starting_mgf_id}")
        if "CCMSLIB" in starting_mgf_id:
            current_id = int(starting_mgf_id.split("CCMSLIB")[1])

    for _, row in rt_peaks_and_top_spectra.iterrows():
        spectrum = row['spectrum']
        if spectrum is None or (isinstance(spectrum, float) and pd.isna(spectrum)):
            continue
        mz_values, intensity_values = spectrum
        
        # Increment the SPECTRUMID for each entry
        if config['msms_refs']['standalone_mgf'] is True:
            current_id += 1
            new_spectrum_id = f"{config['msms_refs']['msms_refs_metadata']['msms_refs_prefix'].upper()}{current_id:011d}"
        elif config['msms_refs']['standalone_mgf'] is False:
            current_id += 1
            new_spectrum_id = f"CCMSLIB{current_id:011d}"

        # Find scan number
        mzml_file = row['standard_lcmsrun'].replace(".h5", ".mzML")
        scan_num = get_scan_number(mzml_file)

        # Calculate charge from adduct
        charge = row['adduct'].split("]")[-1]
        if charge in ["-", "+"]:
            charge = str(1) + charge

        # Create an MGF entry
        mgf_entry = {
            'm/z array': mz_values,
            'intensity array': intensity_values,
            'params': {
                'TITLE': f"MS/MS scan for {row['compound_name']} at {round(row['rt_peak'], 3)} min with intensity {round(row['intensity'], 3)}",
                'PEPMASS': round(row['precursor_mz'], 3),
                'CHARGE': charge,
                'MSLEVEL': 2,
                'SOURCE_INSTRUMENT': 'LC-ESI-Orbitrap',
                'FILENAME': row['standard_lcmsrun'],
                'SEQ': "*..*",
                'IONMODE': "Positive" if row['polarity'] == "POS" else "Negative",
                'ORGANISM': 'BERKELEY-LAB',
                'NAME': f"{row['compound_name']} {row['collision_energy']} {row['adduct']}",
                'PI': 'Trent Northen',
                'DATACOLLECTOR': 'JGI',
                'SMILES': row['smiles'],
                'INCHI': row['inchi'],
                'FORMULA': row['formula'],
                'INCHIAUX': "N/A",
                'PUBMED': "N/A",
                'SUBMITUSER': 'bkieft',
                'LIBRARYQUALITY': 3,
                'SPECTRUMID': new_spectrum_id,
                'SCANS': scan_num
            }
        }
        
        mgf_entries.append(mgf_entry)

    # Write the MGF file
    updated_refs_dir = f"{config['project']['standards_output_path']}/updated_MSMS_refs"
    standalone_tag = "_standalone" if config['msms_refs']['standalone_mgf'] else ""
    if not os.path.exists(updated_refs_dir):
        os.makedirs(updated_refs_dir)
    fname = f"{updated_refs_dir}/berkeley_lab_refs_{timestamp}{standalone_tag}.mgf"
    mgf.write(mgf_entries, fname)

    if config['msms_refs']['standalone_mgf'] is True:
        print(f"Standalone MGF file created: {fname}")
    elif config['msms_refs']['standalone_mgf'] is False:
        temp_fname = f"{fname}.tmp"
        command = f"(cat {config['msms_refs']['current_mgf_path']} && echo && cat {fname}) > {temp_fname}"
        subprocess.run(command, shell=True, check=True)
        os.replace(temp_fname, fname)
        print(f"Updated MGF file created: {fname}")


def get_msms_refs(msms_refs_path: str) -> pd.DataFrame:
    """
    Load MS/MS reference data from a file.

    Args:
        msms_refs_path (str): Path to the MS/MS reference file.

    Returns:
        pd.DataFrame: A DataFrame containing the MS/MS reference data.
    """
    if not os.path.exists(msms_refs_path):
        raise FileNotFoundError(f"MSMS refs file not found at {msms_refs_path}. Please check the path.")
    msms_refs = pd.read_csv(msms_refs_path, sep='\t', index_col=0, low_memory=False)
    print(f"Loaded MSMS refs with {msms_refs.shape[0]} rows and {msms_refs.shape[1]} columns.")
    return msms_refs

def merge_selected_peaks_with_top_spectra(
    selected_peaks: pd.DataFrame, 
    selected_spectra: pd.DataFrame
) -> pd.DataFrame:
    """
    Merge selected peaks with their corresponding top spectra based on shared columns.

    Args:
        selected_peaks (pd.DataFrame): DataFrame containing peak information, including columns like 'label', 'adduct', 'standard_lcmsrun', and 'peak_index'.
        selected_spectra (pd.DataFrame): DataFrame containing spectra information, including columns like 'compound_name', 'adduct', 'standard_lcmsrun', 'peak_index', and 'spectrum'.

    Returns:
        pd.DataFrame: A merged DataFrame containing information from both peaks and spectra.
    """
    selected_peaks_compounded = selected_peaks.copy()
    selected_peaks_compounded = selected_peaks_compounded[selected_peaks_compounded['adduct'] != "Ambiguous"]
    selected_peaks_compounded['compound_name'] = selected_peaks_compounded['label'].str.replace(r" \(peak\d+\)", "", regex=True)
    selected_spectra_compounded = selected_spectra[selected_spectra['adduct'] != "Ambiguous"]
    merge_on = ['compound_name', 'adduct', 'standard_lcmsrun', 'peak_index']
    merge_on_with_spectrum = merge_on + ['spectrum']  # Create a new list with 'spectrum' added
    merged_peaks_spectra = pd.merge(
        selected_peaks_compounded,
        selected_spectra_compounded[merge_on_with_spectrum],  # Subset columns correctly
        on=merge_on,
        how='left'
    )
    return merged_peaks_spectra

def format_for_msms_refs(
    input_df: pd.DataFrame, 
    spectra_df: pd.DataFrame,
    msms_refs: pd.DataFrame, 
    config: Dict[str, Any],
) -> pd.DataFrame:
    """
    Format a DataFrame for MS/MS reference data by enriching metadata and ensuring required columns.

    Args:
        input_df (pd.DataFrame): Input DataFrame containing compound data, including a 'spectrum' column.
        msms_refs (pd.DataFrame): Existing MS/MS reference DataFrame to align columns with.
        config (Dict[str, Any]): Configuration dictionary containing metadata for MS/MS references.
                                 This contains a metadata dictionary with fields such as 'ce_type', 
                                             'frag_method', 'instrument_type', 'decimal', and 'msms_refs_prefix'.

    Returns:
        pd.DataFrame: A formatted DataFrame ready for MS/MS reference storage.
    """
    # Add spectra from top spectra to the rt peaks data
    input_df_with_spectra = merge_selected_peaks_with_top_spectra(input_df, spectra_df)

    # Remove rows with NaN in the 'spectrum' column and print a warning
    rows_with_nan = input_df_with_spectra[input_df_with_spectra['spectrum'].isna()]
    if not rows_with_nan.empty:
        print("Warning: The following rows were removed due to NaN in the 'spectrum' column:")
        display(rows_with_nan)
        input_df_with_spectra = input_df_with_spectra.dropna(subset=['spectrum'])

    # Add all required columns for MS/MS refs
    output_df = input_df_with_spectra.copy()
    output_df['ce_type'] = config['msms_refs']['msms_refs_metadata']['ce_type']
    output_df['ce'] = output_df['standard_lcmsrun'].apply(get_collision_energy)
    output_df['file'] = output_df['standard_lcmsrun'].apply(os.path.basename)
    output_df.rename(columns={'mz_theoretical': 'mz'}, inplace=True)
    output_df = enrich_metadata(output_df, config['msms_refs']['msms_refs_metadata'])
    output_df['spectrum'] = output_df['spectrum'].apply(make_text_spectrum)
    output_df = output_df[msms_refs.columns.intersection(output_df.columns)]
    output_df = output_df.reset_index(drop=True)
    output_df.index = range(
        msms_refs.index.max() + 1,
        msms_refs.index.max() + 1 + len(output_df)
    )

    return output_df

###################
### Atlas Tools ###
###################

def convert_rt_peaks_to_atlas_format(rt_peaks: pd.DataFrame) -> pd.DataFrame:
    """
    Converts a DataFrame of retention time (RT) peaks into a format compatible with the MetAtlas atlas convention.

    Args:
        rt_peaks (pd.DataFrame): A DataFrame containing RT peak data with columns such as 'label', 'rt_peak', 
                                 'smiles', 'mz_theoretical', 'monoisotopic_mass', 'polarity', and 'chromatography'.

    Returns:
        pd.DataFrame: A formatted DataFrame with enriched metadata and columns renamed/dropped to match 
                      the MetAtlas atlas convention. The output is sorted by chromatography and polarity.
    """
    rt_peaks_unformatted = rt_peaks.copy()
    rt_peaks_unformatted = rt_peaks_unformatted[rt_peaks_unformatted['adduct'] != "Ambiguous"]
    rt_peaks_unformatted['compound_name'] = rt_peaks_unformatted['label']

    # Enrich for atlas-related metadata
    rt_peaks_unformatted['rt_min'] = rt_peaks_unformatted['rt_peak'] - 0.5
    rt_peaks_unformatted['rt_max'] = rt_peaks_unformatted['rt_peak'] + 0.5
    rt_peaks_unformatted['mz_tolerance'] = 5
    rt_peaks_unformatted['mz_tolerance_units'] = "ppm"
    rt_peaks_unformatted['inchi'] = rt_peaks_unformatted['smiles'].apply(lambda row: AllChem.MolToInchi(AllChem.MolFromSmiles(row)))
    rt_peaks_unformatted['inchi_key'] = rt_peaks_unformatted['inchi'].apply(inchi_to_inchikey)
    rt_peaks_unformatted['in_metatlas'] = "True"

    # Rename and drop columns to match MetAtlas atlas convention
    rt_peaks_unformatted.rename(columns={'mz_theoretical': 'mz', 'monoisotopic_mass': 'mono_isotopic_molecular_weight'}, inplace=True)
    rt_peaks_unformatted['polarity'] = rt_peaks_unformatted['polarity'].apply(lambda pol: 'positive' if pol == 'POS' else 'negative')
    rt_peaks_unformatted.drop(columns=['intensity', 'mz_observed', 'ppm_error'], inplace=True)

    # Export chromatographic and polarity-specific atlases
    pos_annotations = rt_peaks_unformatted[rt_peaks_unformatted['polarity'] == 'positive']
    neg_annotations = rt_peaks_unformatted[rt_peaks_unformatted['polarity'] == 'negative']
    c18_pos_annotations = pos_annotations[pos_annotations['chromatography'] == 'C18'].sort_values('rt_peak')
    c18_neg_annotations = neg_annotations[neg_annotations['chromatography'] == 'C18'].sort_values('rt_peak')
    hilic_pos_annotations = pos_annotations[pos_annotations['chromatography'] == 'HILICZ'].sort_values('rt_peak')
    hilic_neg_annotations = neg_annotations[neg_annotations['chromatography'] == 'HILICZ'].sort_values('rt_peak')

    rt_peaks_formatted = pd.concat([c18_pos_annotations, c18_neg_annotations, hilic_pos_annotations, hilic_neg_annotations], ignore_index=True)
    
    return rt_peaks_formatted

def get_ema_atlas_data(ema_atlases_path: Dict[str, Dict[str, str]]) -> Dict[str, Dict[str, Union[str, pd.DataFrame]]]:
    """
    Load and combine atlas data from multiple files into a single DataFrame.

    Args:
        ema_atlases_path (Dict[str, Dict[str, str]]): A dictionary where keys are chromatography types,
            and values are dictionaries mapping polarities to file paths or metatlas UIDs.

    Returns:
        Dict[str, Dict[str, Union[str, pd.DataFrame]]]: Same format as ema_atlases_path dict, but with DataFrames instead of file paths or UUIDs.
        This would be a nested dictionary with keys of both chromatography and polarity, and values of DataFrames.
    """
    atlas_dfs = ema_atlases_path.copy()
    for chrom_type, polarities in atlas_dfs.items():
        for polarity, file_path in polarities.items():
            if file_path.endswith('.tsv'):
                if not os.path.exists(file_path):
                    print(f"Warning: file {file_path} does not exist.")
                    #atlas_dfs[chrom_type][polarity] = None
                    continue
                else:
                    df = pd.read_csv(file_path, sep='\t', index_col=None)
                    df['source_atlas'] = os.path.basename(file_path)
                    atlas_dfs[chrom_type][polarity] = df
            elif file_path.endswith('.csv'):
                if not os.path.exists(file_path):
                    print(f"Warning: file {file_path} does not exist.")
                    #atlas_dfs[chrom_type][polarity] = pd.DataFrame()
                    continue
                else:
                    df = pd.read_csv(file_path, index_col=None)
                    df['source_atlas'] = os.path.basename(file_path)
                    atlas_dfs[chrom_type][polarity] = df
            else:
                try:
                    df = atlas_id_to_df(file_path)
                    df['source_atlas'] = f"{file_path}_{chrom_type}_{polarity}"
                except:
                    df = pd.DataFrame()
                    print(f"Warning: Unable to retrieve atlas with ID {file_path}. Returning empty dataframe")
                atlas_dfs[chrom_type][polarity] = df

    return atlas_dfs

def format_for_atlas_store(input_compounds: pd.DataFrame) -> List[Any]:
    """
    Format a DataFrame of compounds for storage in the MetAtlas database.

    Args:
        input_compounds (pd.DataFrame): DataFrame containing compound information, including columns such as 
                                        'label', 'inchi', 'inchi_key', 'neutralized_inchi', 'neutralized_inchi_key', 
                                        'formula', 'monoisotopic_mass', and 'permanent_charge'.

    Returns:
        List[Any]: A list of MetAtlas Compound objects ready for database storage.
    """
    # Second check for non-neutralized InChI key
    any_diff_inchikeys = input_compounds[input_compounds['inchi_key'] != input_compounds['neutralized_inchi_key']]
    if not any_diff_inchikeys.empty:
        print("Warning! The InChIKey and neutralized InChIKey do not match for the following compounds:")
        print(any_diff_inchikeys[['label', 'inchi_key', 'neutralized_inchi_key']])

    input_compounds.reset_index(drop=True, inplace=True)
    input_compounds.rename(columns={'label': 'compound_name'}, inplace=True)

    metatlas_compounds = []
    for _, row in input_compounds.iterrows():
        cpd = metob.Compound()
        cpd.name = row.compound_name
        cpd.inchi = row.inchi
        cpd.inchi_key = row.inchi_key
        cpd.neutralized_inchi = row.neutralized_inchi
        cpd.neutralized_inchi_key = row.neutralized_inchi_key
        cpd.formula = row.formula
        cpd.mono_isotopic_molecular_weight = row.monoisotopic_mass
        cpd.permanent_charge = row.permanent_charge
        cpd.description = ''
        metatlas_compounds.append(cpd)

    return metatlas_compounds

def atlas_id_to_df(atlas_unique_id: str) -> pd.DataFrame:
    """
    Retrieve an atlas from the database using its unique ID and convert it into a DataFrame.

    Args:
        atlas_unique_id (str): The unique identifier of the atlas.

    Returns:
        pd.DataFrame: A DataFrame containing compound identification data from the atlas.
    """
    atlas = get_atlas(atlas_unique_id)
    atlas_df = make_atlas_df(atlas)
        
    return atlas_df

def get_qc_files(files_path: str, chromatography: str, include_istds: bool = False) -> List[str]:
    """
    Retrieve all QC files from a specified raw data path.

    Args:
        files_path (str): Path to the directory containing raw data files.
        chromatography (str): Chromatography type ('C18' or 'HILIC').
        include_istds (bool, optional): Whether to include ISTD files. Defaults to False.

    Returns:
        List[str]: A list of file paths for QC files.
    """
    all_files = glob.glob(os.path.join(files_path, f"*.h5"))

    polarity = "fps"
    if chromatography == 'C18':
        chromatography = 'C18'
    elif chromatography == 'HILIC':
        chromatography = 'HILICZ'

    if include_istds:
        qc_files = [
            file for file in all_files if os.path.basename(file).split('_')[9].lower() == polarity and 
            os.path.basename(file).split('_')[7].lower() == chromatography.lower() and
            ('QC_' in file or 'ISTD_' in file)
        ]
    else:
        qc_files = [
            file for file in all_files if os.path.basename(file).split('_')[9].lower() == polarity and 
            os.path.basename(file).split('_')[7].lower() == chromatography.lower() and
            'QC_' in file
        ]

    return qc_files


def collect_qc_ms1_data(qc_atlas: pd.DataFrame, qc_files: List[str], polarity: str) -> pd.DataFrame:
    """
    Collect MS1 data for QC compounds from a list of QC files.

    Args:
        qc_atlas (pd.DataFrame): A DataFrame containing QC atlas data.
        qc_files (List[str]): A list of QC file paths.
        polarity (str): The ionization polarity ('pos' or 'neg').

    Returns:
        pd.DataFrame: A DataFrame containing the collected MS1 data.
    """
    experiment_input = ft.setup_file_slicing_parameters(qc_atlas, qc_files, base_dir=os.getcwd(), ppm_tolerance=10, polarity=polarity)

    ms1_data = []
    for file_input in tqdm(experiment_input, desc=" Collecting MS1 data for QC compounds", unit=' file'):
        data = ft.get_data(file_input, save_file=False, return_data=True)
        data['ms1_summary']['lcmsrun_observed'] = file_input['lcmsrun']
        ms1_data.append(data['ms1_summary'])

    return pd.concat(ms1_data)

def get_qc_atlas_dataframe(atlas_identifier: str) -> pd.DataFrame:
    """
    Retrieve the atlas as a DataFrame. If the identifier ends with '.tsv' or '.csv', 
    read it as a file. Otherwise, use the atlas_id_to_df function to retrieve it 
    from the database.

    Args:
        atlas_identifier (str): The atlas identifier, either a file path (ending in '.tsv' or '.csv') 
                                or an atlas ID.

    Returns:
        pd.DataFrame: The atlas as a DataFrame.
    """
    if atlas_identifier.endswith('.tsv'):
        if not os.path.exists(atlas_identifier):
            print(f"Warning: file {atlas_identifier} does not exist.")
            return
        else:
            return pd.read_csv(atlas_identifier, sep='\t')
    elif atlas_identifier.endswith('.csv'):
        if not os.path.exists(atlas_identifier):
            print(f"Warning: file {atlas_identifier} does not exist.")
            return
        else:
            return pd.read_csv(atlas_identifier)
    else:
        return atlas_id_to_df(atlas_identifier)
    
def run_rt_correction(
    nonmatches_to_atlases: pd.DataFrame, 
    config: Dict[str, Any]
) -> None:
    """
    Run retention time (RT) correction for compounds not matched to existing atlases.
    This function collects experimental QC data, retrieves baseline atlases, and performs RT correction.
    Args:
        nonmatches_to_atlases (pd.DataFrame): DataFrame containing compounds not matched to existing atlases.
        config (Dict[str, Any]): Configuration dictionary containing parameters for mapping chromatography types to baseline atlas file paths or IDs.
    Returns:
        Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]: A tuple containing two dictionaries:
            - baseline_correction_inputs: Dictionary where keys are chromatography types and values are DataFrames with experimental and baseline QC data.
            - baseline_correction_outputs: Dictionary where keys are chromatography types and values are DataFrames with RT correction outputs.
    """

    _, _, baseline_to_experimental_qc = get_qc_experimental_atlas(nonmatches_to_atlases, config, include_istds=True)
    baseline_correction_inputs = create_baseline_correction_input(nonmatches_to_atlases, baseline_to_experimental_qc, config)
    baseline_correction_outputs = rt_correction_from_baseline(baseline_correction_inputs, config)

    return baseline_correction_inputs, baseline_correction_outputs

def get_qc_experimental_atlas(
    nonmatches_to_atlases: pd.DataFrame, 
    config: Dict[str, Any],
    include_istds: bool = False
) -> Dict[str, pd.DataFrame]:
    """
    Generate experimental QC atlases by collecting MS1 data and comparing it with baseline atlases.

    Args:
        nonmatches_to_atlases (pd.DataFrame): DataFrame containing compounds not matched to existing atlases.
        config (Dict[str, Any]): Configuration dictionary containing parameters for mapping chromatography types to baseline atlas file paths or IDs.
        include_istds (bool, optional): Whether to include internal standards (ISTDs) in the QC files. Defaults to False.

    Returns:
        Dict[str, pd.DataFrame]: A dictionary where keys are chromatography types and values are DataFrames containing 
                                 combined QC data with experimental and baseline retention times.
    """
    chromatographies = nonmatches_to_atlases['chromatography'].unique()

    experimental_qc = {}
    baseline_qc = {}
    combined_qc = {}

    for chrom in chromatographies:
        projects = nonmatches_to_atlases['standard_lcmsrun'].apply(os.path.dirname).unique()
        project = [p for p in projects if chrom.lower() in os.path.basename(p).split('_')[7].lower()]
        if len(project) == 1:
            project = project[0]
        elif len(project) == 0:
            print(f"Warning: no project found for {chrom}")
            continue
        elif len(project) > 1:
            print(f"Warning: more than one project found for {chrom}: {project}")
            continue
        print(f"\tGetting all QC files for project {project}\n")
        qc_files = get_qc_files(project, chrom, include_istds)
        
        print(f"\tRetrieving baseline {chrom} QC atlas: {config['atlases']['current_qc_atlases'][chrom]}\n")
        baseline_atlas_df = get_qc_atlas_dataframe(config["atlases"]["current_qc_atlases"][chrom])
        baseline_qc[chrom] = baseline_atlas_df

        print(f"\tCollecting QC MS1 data for {chrom}...\n")
        ms1_summary = collect_qc_ms1_data(baseline_atlas_df, qc_files=qc_files, polarity="pos")
        ms1_summary_median = ms1_summary.groupby('label', as_index=False).agg({'rt_peak': 'median'})
        ms1_summary_median['rt_min'] = ms1_summary_median['rt_peak'] - 0.5
        ms1_summary_median['rt_max'] = ms1_summary_median['rt_peak'] + 0.5
        experimental_qc[chrom] = ms1_summary_median

        formatted = ms1_summary_median.copy()
        formatted.rename(columns={'rt_peak': 'rt_peak_experimental', 'rt_min': 'rt_min_experimental', 'rt_max': 'rt_max_experimental'}, inplace=True)
        formatted = pd.merge(formatted, baseline_atlas_df[['label', 'rt_peak', 'rt_min', 'rt_max']], on='label', how='left')
        formatted.rename(columns={'rt_peak': 'rt_peak_baseline', 'rt_min': 'rt_min_baseline', 'rt_max': 'rt_max_baseline'}, inplace=True)
        formatted['polarity'] = 'QC'
        combined_qc[chrom] = formatted

    return baseline_qc, experimental_qc, combined_qc

def create_baseline_correction_input(
    compounds_to_correct: pd.DataFrame, 
    baseline_to_experimental_qc: Dict[str, pd.DataFrame],
    config: Dict[str, Any]
) -> Dict[str, pd.DataFrame]:
    """
    Create input data for retention time (RT) correction by combining experimental and baseline QC data.

    Args:
        compounds_to_correct (pd.DataFrame): DataFrame containing compounds to correct with columns 
                                             ['label', 'polarity', 'chromatography', 'rt_peak', 'rt_min', 'rt_max'].
        baseline_to_experimental_qc (Dict[str, pd.DataFrame]): Dictionary where keys are chromatography types and 
                                                               values are DataFrames with experimental and baseline QC data.
        config (Dict[str, Any]): Configuration dictionary containing parameters for RT correction.

    Returns:
        Dict[str, pd.DataFrame]: A dictionary where keys are chromatography types and values are DataFrames 
                                 containing combined data for RT correction.
    """
    uncorrected_atlas = compounds_to_correct[['label', 'adduct', 'polarity', 'chromatography', 'rt_peak', 'rt_min', 'rt_max']]

    baseline_correction_inputs = {}
    for chrom in baseline_to_experimental_qc.keys():
        uncorrected_atlas_chrom = uncorrected_atlas[uncorrected_atlas['chromatography'] == chrom]
        uncorrected_atlas_chrom = uncorrected_atlas_chrom.drop(columns=['chromatography'])
        uncorrected_atlas_chrom = uncorrected_atlas_chrom.rename(columns={'rt_peak': 'rt_peak_experimental', 'rt_min': 'rt_min_experimental', 'rt_max': 'rt_max_experimental'})
        uncorrected_atlas_chrom.loc[:, 'rt_peak_baseline'] = np.nan
        uncorrected_atlas_chrom.loc[:, 'rt_min_baseline'] = np.nan
        uncorrected_atlas_chrom.loc[:, 'rt_max_baseline'] = np.nan

        qc_atlas_chrom = baseline_to_experimental_qc[chrom]

        baseline_correction_input = pd.concat([uncorrected_atlas_chrom, qc_atlas_chrom], axis=0, ignore_index=True)
        baseline_correction_inputs[chrom] = baseline_correction_input

    return baseline_correction_inputs


def rt_correction_from_baseline(
    baseline_correction_dfs: Dict[str, pd.DataFrame], 
    config: Dict[str, Any],
) -> Dict[str, pd.DataFrame]:
    """
    Perform RT correction for each chromatography in the given dictionary of DataFrames.

    Args:
        baseline_correction_dfs (Dict[str, pd.DataFrame]): A dictionary where keys are chromatography types 
            (e.g., 'C18', 'HILIC') and values are DataFrames containing RT experimental and baseline data.
        config (Dict[str, Any]): Configuration dictionary containing parameters for RT correction.

    Returns:
        Dict[str, pd.DataFrame]: A dictionary where keys are chromatography types and values are corrected DataFrames.
    """
    corrected_dfs: Dict[str, pd.DataFrame] = {}

    print("\tPerforming RT correction...\n")
    for chromatography, df in tqdm(baseline_correction_dfs.items(), desc="Calculating RT correction model", unit=' chromatography'):
        
        if config['atlases']['rt_correction_models'][chromatography] is None:
            fit_data = df.dropna(subset=["rt_peak_baseline", "rt_peak_experimental"])
            fit_data = fit_data[fit_data["polarity"] == "QC"]
            fit_data["rt_peak_experimental"] = pd.to_numeric(fit_data["rt_peak_experimental"], errors="coerce")
            fit_data["rt_peak_baseline"] = pd.to_numeric(fit_data["rt_peak_baseline"], errors="coerce")
            coefficients = np.polyfit(fit_data["rt_peak_experimental"], fit_data["rt_peak_baseline"], 2)
            polynomial = np.poly1d(coefficients)
            type_of_fit = "auto-generated"
        elif config['atlases']['rt_correction_models'][chromatography] is not None:
            coefficients = config['atlases']['rt_correction_models'][chromatography]
            polynomial = np.poly1d(coefficients)
            type_of_fit = "user-defined"

        # Apply the polynomial to all rows where rt_experimental is available
        df["rt_peak_corrected"] = df["rt_peak_experimental"].apply(
            lambda x: polynomial(x) if not np.isnan(x) else np.nan
        )
        df["rt_min_corrected"] = df.apply(
            lambda row: row["rt_peak_corrected"] - row["rt_peak_experimental"] + row["rt_min_experimental"]
            if not np.isnan(row["rt_peak_corrected"])
            else np.nan,
            axis=1
        )
        df["rt_max_corrected"] = df.apply(
            lambda row: row["rt_peak_corrected"] + row["rt_max_experimental"] - row["rt_peak_experimental"]
            if not np.isnan(row["rt_peak_corrected"])
            else np.nan,
            axis=1
        )

        # Check difference between rt_peak_experimental and rt_peak_corrected
        df['rt_diff_experimental_vs_corrected'] = df["rt_peak_experimental"] - df["rt_peak_corrected"]
        large_diff_rows = df[df['rt_diff_experimental_vs_corrected'].abs() > 0.5]
        if not large_diff_rows.empty:
            print(f"Warning: Large differences in some experimental vs predicted RTs for {chromatography} chromatography:\n")
            print(large_diff_rows[['label', 'polarity', 'rt_peak_experimental', 'rt_peak_corrected', 'rt_diff_experimental_vs_corrected']])

        # Store the corrected DataFrame in the result dictionary
        corrected_dfs[chromatography] = df[
            ['label', 'adduct', 'polarity', 'rt_peak_baseline', 'rt_peak_experimental', 
             'rt_peak_corrected', 'rt_min_corrected', 'rt_max_corrected', 'rt_diff_experimental_vs_corrected']
        ]

        print(f"\t{chromatography} RT correction results using {type_of_fit} polynomial coefficients: {coefficients}")
        display(corrected_dfs[chromatography])

    return corrected_dfs


def substitute_corrected_rt_values(
    nonmatches_to_atlases: pd.DataFrame, 
    baseline_correction_outputs: Dict[str, pd.DataFrame]
) -> pd.DataFrame:
    """
    Substitutes rt_peak, rt_min, and rt_max values in nonmatches_to_atlases
    with corrected values from baseline_correction_outputs for all chromatography keys.

    Args:
        nonmatches_to_atlases (pd.DataFrame): DataFrame containing the original RT values.
        baseline_correction_outputs (Dict[str, pd.DataFrame]): Dictionary containing corrected RT values for each chromatography.

    Returns:
        pd.DataFrame: Updated nonmatches_to_atlases DataFrame with substituted RT values.
    """

    # Make rt_peak_experimental a merge key to ensure all rows are uniquely aligned
    updated_atlas = nonmatches_to_atlases.copy()
    updated_atlas.rename(columns={'rt_peak': 'rt_peak_experimental'}, inplace=True)
    
    for chromatography, df in baseline_correction_outputs.items():

        # Merge the two DataFrames on 'label' and 'polarity' to align rows
        baseline_df = df.copy()
        baseline_df.loc[:, 'chromatography'] = str(chromatography)
        merged_df = updated_atlas.merge(
            baseline_df[['label', 'adduct', 'polarity', 'chromatography', 'rt_peak_experimental','rt_peak_corrected', 'rt_min_corrected', 'rt_max_corrected']],
            on=['label', 'adduct', 'chromatography', 'polarity', 'rt_peak_experimental'],
            how='left'
        )

        # Check if merged_df has any non-QC values in the 'polarity' column
        if merged_df[merged_df['polarity'] != 'QC'].empty: # There are no compounds to correct for this chromatography
            print(f"No compounds to correct for {chromatography} chromatography.")
            continue

        # Substitute the RT values with the corrected ones
        merged_df.loc[:, 'rt_peak_experimental'] = merged_df.loc[:, 'rt_peak_corrected'].combine_first(merged_df.loc[:, 'rt_peak_experimental'])
        merged_df.loc[:, 'rt_min'] = merged_df.loc[:, 'rt_min_corrected'].combine_first(merged_df.loc[:, 'rt_min'])
        merged_df.loc[:, 'rt_max'] = merged_df.loc[:, 'rt_max_corrected'].combine_first(merged_df.loc[:, 'rt_max'])

        # Drop the corrected columns to clean up
        merged_df.drop(columns=['rt_peak_corrected', 'rt_min_corrected', 'rt_max_corrected'], inplace=True)

        # Update the main DataFrame for the next iteration
        updated_atlas = merged_df
        print(f"Formatted {updated_atlas.shape[0]} RT-corrected compounds for insertion into {chromatography} atlases.")

    # Rename the columns back to their original names
    updated_atlas.rename(columns={'rt_peak_experimental': 'rt_peak'}, inplace=True)

    return updated_atlas

def update_and_save_ema_atlases(
    nonmatches_to_atlases_rt_corrected: pd.DataFrame,
    ema_atlases: Dict[str, Dict[str, Union[str, pd.DataFrame]]],
    config: Dict[str, Any],
    current_time: str,
) -> Tuple[Dict[str, Dict[str, Optional[str]]], Dict[str, Dict[str, Optional[str]]]]:
    """
    Update and save EMA atlases by merging new compounds and splitting them back into their original files.

    Args:
        nonmatches_to_atlases_rt_corrected (pd.DataFrame): DataFrame containing new compounds to be added to the atlases.
            Must include columns such as 'chromatography', 'polarity', 'inchi_key', and 'rt_peak'.
        ema_atlases (Dict[str, Dict[str, Union[str, pd.DataFrame]]]): Nested dictionary with chromatography and polarity keys,
            containing atlas DataFrames or file paths.
        config (Dict[str, Any]): Configuration dictionary containing metadata for the reference standard project.
            Includes keys like 'standards_output_path', 'standalone_ema_atlas', and 'msms_refs_metadata'.
        current_time (str): Timestamp to append to the updated atlas filenames.

    Returns:
        Tuple[Dict[str, Dict[str, Optional[str]]], Dict[str, Dict[str, Optional[str]]]]:
            - A dictionary of new atlas IDs for each chromatography and polarity.
            - A dictionary of new atlas names or file paths for each chromatography and polarity.
    """
    new_atlas_ids = {}
    new_atlas_names = {}
    for chrom, polarities in ema_atlases.items():
        new_atlas_ids[chrom] = {}
        new_atlas_names[chrom] = {}
        for polarity, atlas_data in polarities.items():
            new_atlas_ids[chrom][polarity] = ""
            new_atlas_names[chrom][polarity] = ""
            atlas_name = atlas_data['source_atlas'].iloc[0]

            # Filter nonmatches_to_atlases_rt_corrected to match the current chrom and polarity
            compounds_to_add_to_atlas = nonmatches_to_atlases_rt_corrected[
                (nonmatches_to_atlases_rt_corrected['chromatography'] == chrom) &
                (nonmatches_to_atlases_rt_corrected['polarity'] == polarity)
            ]

            # Check if there are any compounds to add
            if compounds_to_add_to_atlas.empty:
                print(f"No compounds to add to the {chrom} {polarity} atlas. Skipping.")
                continue

            # Add some columns to atlases depending on chrom (they have different formats) to get as much info as possible
            if chrom == "HILICZ":
                compounds_to_add_to_atlas = compounds_to_add_to_atlas.copy()
                compounds_to_add_to_atlas.rename(columns={'polarity': 'detected_polarity'}, inplace=True)
                compounds_to_add_to_atlas['identification_notes'] = ""
                compounds_to_add_to_atlas['ms1_notes'] = ""
                compounds_to_add_to_atlas['ms2_notes'] = ""
                compounds_to_add_to_atlas['rt_units'] = "min"
                synonyms_dict = get_pubchem_synonyms_from_inchi_keys(compounds_to_add_to_atlas['inchi_key'].tolist())
                compounds_to_add_to_atlas['synonyms'] = compounds_to_add_to_atlas['inchi_key'].map(synonyms_dict)
                pubchem_id_dict = get_pubchem_cids_from_inchi_keys(compounds_to_add_to_atlas['inchi_key'].tolist())
                compounds_to_add_to_atlas['pubchem_compound_id'] = compounds_to_add_to_atlas['inchi_key'].map(pubchem_id_dict)
                compounds_to_add_to_atlas['file_name'] = compounds_to_add_to_atlas['standard_lcmsrun']
                compounds_to_add_to_atlas['file_type'] = "lcms_file"
                compounds_to_add_to_atlas['source_atlas'] = config['msms_refs']['msms_refs_metadata']['msms_refs_prefix']
            if chrom == "C18":
                compounds_to_add_to_atlas = compounds_to_add_to_atlas.copy()
                compounds_to_add_to_atlas['file_name'] = compounds_to_add_to_atlas['standard_lcmsrun']
                compounds_to_add_to_atlas['file_type'] = "lcms_file"
                compounds_to_add_to_atlas['library'] = config['msms_refs']['msms_refs_metadata']['msms_refs_prefix']
                compounds_to_add_to_atlas['code'] = compounds_to_add_to_atlas['library'].str.upper() + compounds_to_add_to_atlas.index.astype(str)
                compounds_to_add_to_atlas['name'] = compounds_to_add_to_atlas['compound_name']
                compounds_to_add_to_atlas['identification_notes'] = ""
                compounds_to_add_to_atlas['ms1_notes'] = ""
                compounds_to_add_to_atlas['ms2_notes'] = ""                
                compounds_to_add_to_atlas['source_atlas'] = config['msms_refs']['msms_refs_metadata']['msms_refs_prefix']
            for col in [
                'file_name', 'file_type', 'identification_notes',
                'ms1_notes', 'ms2_notes', 'compound_classes', 'compound_pathways'
            ]:
                if col not in atlas_data.columns:
                    atlas_data[col] = ""

            # Format the columns of compounds to be added to match the atlas format so no columns are missing
            missing_columns = atlas_data.columns.difference(compounds_to_add_to_atlas.columns)
            for col in missing_columns:
                compounds_to_add_to_atlas = compounds_to_add_to_atlas.copy()
                compounds_to_add_to_atlas.loc[:, col] = np.nan
            common_columns = compounds_to_add_to_atlas.columns.intersection(atlas_data.columns)
            compounds_to_add_to_atlas_formatted = compounds_to_add_to_atlas[common_columns]

            # Assign unique permanent_index to new rows
            if 'permanent_index' in atlas_data.columns:
                max_permanent_index = atlas_data['permanent_index'].max()
                new_indices = range(max_permanent_index + 1, max_permanent_index + 1 + len(compounds_to_add_to_atlas_formatted))
                compounds_to_add_to_atlas_formatted = compounds_to_add_to_atlas_formatted.copy()
                compounds_to_add_to_atlas_formatted['permanent_index'] = [str(int(index)) for index in new_indices]

            # Merge the compounds_to_add_to_atlas_formatted with the current atlas
            column_order = atlas_data.columns
            if config['atlases']['standalone_ema_atlas'] is True:
                new_atlas_data = compounds_to_add_to_atlas_formatted
            elif config['atlases']['standalone_ema_atlas'] is False:
                new_atlas_data = pd.concat([atlas_data, compounds_to_add_to_atlas_formatted], ignore_index=False)
            new_atlas_data.drop(columns=['source_atlas'], inplace=True)
            new_atlas_data = new_atlas_data.reindex(columns=[col for col in column_order if col in new_atlas_data.columns])

            # Sort the new atlas by rt_peak lowest to highest
            new_atlas_data.sort_values(by=["rt_peak", "mz"], ascending=[True, False], inplace=True)

            # Set up for printing/depositing
            updated_atlas_dir = f"{config['project']['standards_output_path']}/updated_EMA_atlases"
            if not os.path.exists(updated_atlas_dir):
                os.makedirs(updated_atlas_dir)
            standalone_tag = "_standalone" if config['atlases']['standalone_ema_atlas'] else ""
            if atlas_name.endswith('.tsv'):
                new_atlas_name = atlas_name.replace(".tsv", f"_{current_time}{standalone_tag}.tsv")
            elif atlas_name.endswith('.csv'):
                new_atlas_name = atlas_name.replace(".csv", f"_{current_time}{standalone_tag}.tsv")
            else: # If atlas input from yaml is a metatlas UID
                new_atlas_name = f"{atlas_name}_{current_time}{standalone_tag}.csv"
            fname = f"{updated_atlas_dir}/{new_atlas_name}"
            
            # Save the new atlas to the original file
            if fname.endswith('.tsv'):
                new_atlas_data.to_csv(fname, sep='\t', index=False)
            elif fname.endswith('.csv'):
                new_atlas_data.to_csv(fname, index=False)
            
            # Print summary of the update
            if config['atlases']['standalone_ema_atlas'] is True:
                print(f"\nStandalone atlas (not appended to existing data) generated for {chrom} {polarity} EMA atlas.")
                print(f"New atlas has {len(new_atlas_data)} compounds.")
                print(f"Saved to: {fname}\n")
            elif config['atlases']['standalone_ema_atlas'] is False:
                print(f"\nCurrent {chrom} {polarity} EMA atlas: {atlas_name}")
                print(f"{len(atlas_data)} current compounds updated with {len(compounds_to_add_to_atlas_formatted)} new compounds for a total of {len(new_atlas_data)} compounds.")
                print(f"Updated {chrom} {polarity} EMA atlas saved to: {fname}\n")

            # Deposit atlas directly to metatlas db and retrieve the new atlas ID
            if config['atlases']['direct_deposit_new_emas']:
                if fname.endswith('.tsv'): # Must convert to csv for atlas deposit
                    input_atlas_file_name = fname.replace('.tsv', '.csv')
                    new_atlas_data.to_csv(input_atlas_file_name, index=False)
                else:
                    input_atlas_file_name = fname
                new_atlas_db_name = os.path.basename(input_atlas_file_name).replace(".csv","")
                print("Depositing atlas to metatlas database Atlas table...")
                atlas_deposited = make_atlas_from_spreadsheet(input_atlas_file_name, new_atlas_db_name, filetype="csv", polarity=polarity, store=True, mz_tolerance=config['project']['ppm_tolerance'])
                atlas_df = make_atlas_df(atlas_deposited)
                print(f"\tUpdated EMA atlas deposited to metatlas db with unique_id: {atlas_deposited.unique_id}")
                print(f"\tUpdated EMA atlas deposited to metatlas db with name: {atlas_deposited.name}")
                new_atlas_ids[chrom][polarity] = atlas_deposited.unique_id
                new_atlas_names[chrom][polarity] = atlas_deposited.name
            else:
                new_atlas_ids[chrom][polarity] = ""
                new_atlas_names[chrom][polarity] = fname
    
    return new_atlas_ids, new_atlas_names


def update_and_save_msms_refs(
    msms_refs: pd.DataFrame,
    rt_peaks_with_spectra: pd.DataFrame,
    config: Dict[str, Any],
    timestamp: str,
) -> None:
    """
    Update and save MS/MS reference data by appending new entries and optionally saving to disk.

    Args:
        msms_refs (pd.DataFrame): Existing MS/MS reference data.
        rt_peaks_with_spectra (pd.DataFrame): New MS/MS reference data to be added.
        config (Dict[str, Any]): Configuration dictionary containing metadata for ref std project (used for setting write path and project name).
        timestamp (str): Timestamp to append to the updated MS/MS reference filename.

    Returns:
        None
    """
    # Adjust the index of rt_peaks_with_spectra to start where msms_refs leaves off
    start_index = msms_refs.index.max() + 1
    rt_peaks_with_spectra = rt_peaks_with_spectra.reset_index(drop=True)
    rt_peaks_with_spectra.index = range(start_index, start_index + len(rt_peaks_with_spectra))

    # Check if the columns of msms_refs and rt_peaks_with_spectra match
    if not all(msms_refs.columns == rt_peaks_with_spectra.columns):
        print("Warning! Column names don't match between existing and new MSMS refs.")
        print(f"Existing columns: {msms_refs.columns}")
        print(f"New columns: {rt_peaks_with_spectra.columns}")
        return
    
    # Combine existing and new MS/MS refs
    if config['msms_refs']['standalone_msms_refs'] is True:
        new_msms_refs = rt_peaks_with_spectra
        print(f"Standalone MSMS refs generated with {new_msms_refs.shape[0]} compounds.")
    elif config['msms_refs']['standalone_msms_refs'] is False:
        new_msms_refs = pd.concat([msms_refs, rt_peaks_with_spectra])
        print(f"Existing MSMS refs went from {msms_refs.shape[0]} to {new_msms_refs.shape[0]} compounds.")
        if new_msms_refs.shape[0] != msms_refs.shape[0] + rt_peaks_with_spectra.shape[0]:
            print("Warning! Some new MSMS refs may not have been added correctly.")
            return
        if new_msms_refs.shape[1] != msms_refs.shape[1]:
            print("Warning! Column numbers don't match between existing and new MSMS refs.")
            return

    updated_refs_dir = f"{config['project']['standards_output_path']}/updated_MSMS_refs"
    standalone_tag = "_standalone" if config['msms_refs']['standalone_msms_refs'] else ""
    if not os.path.exists(updated_refs_dir):
        os.makedirs(updated_refs_dir)
    fname = f"{updated_refs_dir}/msms_refs_{timestamp}{standalone_tag}.tab"
    print(f"\tNew MSMS refs: {fname}")
    new_msms_refs.to_csv(fname, sep='\t')

    return

####################
### Upload tools ###
####################

def upload_to_google_drive(
    project_folder: str,
    project_name: str,
    timestamp: str,
    overwrite: bool = False
) -> bool:
    """
    Upload all contents of project_folder except the 'cache' subdirectory
    to Google Drive Outputs_Unreviewed/{project_name} using rclone.
    """
    dest_name = os.path.basename(project_name.replace('.csv',''))
    dest_folder = f"Reference_Standards_Annotation:/Outputs_Unreviewed/{dest_name}/{timestamp}_analysis"
    orig_folder = os.path.join(project_folder, f"{timestamp}_analysis")

    # Copy csv and yaml from top level of project_folder to orig_folder
    if os.path.exists(orig_folder):
        command = f"cp {project_name} {orig_folder}/ && cp {project_name.replace('.csv','.yaml')} {orig_folder}"
        try:
            subprocess.check_output(command, shell=True)
        except Exception as e:
            print(f"Warning! Failed to copy CSV and YAML files to {orig_folder} with exception: {e}")
            return False

    # Use --update if overwrite is True, else skip existing files
    update_flag = "--update" if overwrite else "--ignore-existing"
    if overwrite is False:
        print("Warning! Overwrite is set to False, existing files will not be replaced.\n")

    upload_command = (
        f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone copy --exclude cache/** {update_flag} '
        f'"{orig_folder}/" "{dest_folder}"'
    )
    try:
        print(f"Uploading to Google Drive with command:\n\t{upload_command}")
        subprocess.check_output(upload_command, shell=True)
    except Exception as e:
        print(f"Warning! Google Drive upload failed with exception: {e}\nCommand: {upload_command}")
        return False

    # Check that upload worked
    check_upload_command = (
        f'/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone ls "{dest_folder}" --max-depth 1'
    )
    try:
        check_upload_out = subprocess.check_output(check_upload_command, shell=True)
        if check_upload_out.decode('utf-8').strip():
            print(f"\nGoogle Drive upload confirmed with overwrite set to {overwrite}!")
            return True
        else:
            print("Warning! Google Drive upload check failed. Upload may not have been successful.")
            return False
    except Exception as e:
        print(f"Warning! Google Drive upload failed on upload check with exception: {e}\nCommand: {check_upload_command}")
        return False


###########################
### Visualization tools ###
###########################

def generate_adduct_colors(include_adducts: List[str]) -> Dict[str, str]:
    """
    Generates a dictionary mapping adducts to unique colors.

    Args:
        include_adducts (List[str]): List of adducts to include.

    Returns:
        Dict[str, str]: Dictionary mapping adducts to hex color codes.
    """
    adduct_color = {}
    colors = cm.rainbow(np.linspace(0, 1, len(include_adducts)))  # Generate evenly spaced colors
    for i, adduct in enumerate(include_adducts):
        # Convert RGBA array to hex string
        rgba = colors[i]
        hex_color = mcolors.to_hex(rgba)
        adduct_color[adduct] = hex_color
    return adduct_color

def select_compounds_from_gui(full_dataset: pd.DataFrame, selected_adducts_dict: dict) -> pd.DataFrame:
    """
    Select compounds from a dataset based on user-specified criteria, including compound name, LCMS run, adducts, and peak indices.

    Args:
        full_dataset (pd.DataFrame): The full dataset containing compound information.
        selected_adducts_dict (dict): Dictionary of selected adducts for each compound.

    Returns:
        pd.DataFrame: A filtered dataset containing only the selected compounds.
    """
    if not selected_adducts_dict:
        print("No compounds selected for annotation.")
        return pd.DataFrame()

    if 'smiles' in full_dataset.columns:
        full_dataset['inchi'] = full_dataset['smiles'].apply(lambda row: AllChem.MolToInchi(AllChem.MolFromSmiles(row)))
        full_dataset['inchi_key'] = full_dataset['inchi'].apply(inchi_to_inchikey)
        full_dataset['neutralized_inchi'] = full_dataset['inchi'].apply(neutralize_inchi)
        full_dataset['neutralized_inchi_key'] = full_dataset['neutralized_inchi'].apply(inchi_to_inchikey)
        full_dataset['permanent_charge'] = full_dataset['neutralized_inchi'].apply(charge_from_inchi)
        full_dataset['formula'] = full_dataset['neutralized_inchi'].apply(formula_from_inchi)
        full_dataset['monoisotopic_mass'] = full_dataset['neutralized_inchi'].apply(monoisotopic_mass_from_inchi)
        full_dataset['collision_energy'] = full_dataset['standard_lcmsrun'].apply(get_collision_energy)
        if not full_dataset[full_dataset['inchi_key'] != full_dataset['neutralized_inchi_key']].empty:
            print(f"Warning! The InChIKey and neutralized InChIKey do not match for the following selected compounds:")
            print(full_dataset[full_dataset['inchi_key'] != full_dataset['neutralized_inchi_key']][['compound_name', 'inchi_key', 'neutralized_inchi_key']])

    selected_rows_list = []
    for key, (adducts, peaks, bests) in selected_adducts_dict.items():
        compound = key.split(";;")[0]
        lcmsrun = key.split(";;")[1]

        if adducts != ["Ambiguous"] and peaks != ["Ambiguous"]:
            for adduct, peak, best in zip(adducts, peaks, bests):
                mask = (
                    (full_dataset['compound_name'] == compound) &
                    (full_dataset['standard_lcmsrun'] == lcmsrun) &
                    (full_dataset['adduct'] == adduct)
                )
                if 'peak_index' in full_dataset.columns:
                    mask = mask & (full_dataset['peak_index'] == peak)
                keep_row = full_dataset[mask].copy()  # Ensure keep_row is a copy
                if not keep_row.empty:
                    keep_row.loc[:, 'best_adduct'] = best  # Use .loc to set the value
                    selected_rows_list.append(keep_row)

        if adducts == ["Ambiguous"] and peaks == ["Ambiguous"]:
            mask = (
                (full_dataset['compound_name'] == compound) &
                (full_dataset['standard_lcmsrun'] == lcmsrun)
            )
            selected_compounds = full_dataset[mask].copy()  # Ensure selected_compounds is a copy
            # Replace values for ambiguous selection
            ambiguous_values = {
                'adduct': "Ambiguous",
                'peak_index': "Ambiguous",
                'rt_peak': np.nan,
                'rt': np.nan,
                'mz': np.nan,
                'i': np.nan,
                'label': lambda df: df['compound_name'],
                'spectrum': np.nan,
                'intensity': np.nan,
                'total_intensity_fraction': np.nan,
                'total_intensity': np.nan,
                'precursor_peak_height': np.nan,
                'precursor_mz': np.nan,
                'mz_observed': np.nan,
                'mz_theoretical': np.nan,
                'ppm_error': np.nan,
            }
            for col, val in ambiguous_values.items():
                if col in selected_compounds:
                    if callable(val):
                        selected_compounds.loc[:, col] = val(selected_compounds)
                    else:
                        selected_compounds.loc[:, col] = val
            if not selected_compounds.empty:
                selected_rows_list.append(selected_compounds)

    # Filter out empty DataFrames
    selected_rows_list = [df for df in selected_rows_list if not df.empty]

    if not selected_rows_list:
        return pd.DataFrame()

    selected_rows_from_full = pd.concat(selected_rows_list, ignore_index=True)
    key_cols = [col for col in ['compound_name', 'standard_lcmsrun', 'adduct', 'peak_index'] if col in selected_rows_from_full.columns]
    selected_rows_from_full.drop_duplicates(subset=key_cols, keep='first', inplace=True)
    return selected_rows_from_full

def select_compounds_from_gui_old(full_dataset: pd.DataFrame, selected_compounds_table: pd.DataFrame, method: str = "all") -> pd.DataFrame:
    """
    Select compounds from a dataset based on user-specified criteria, including compound name, LCMS run, adducts, and peak indices.

    Args:
        full_dataset (pd.DataFrame): The full dataset containing compound information.
        selected_compounds_table (pd.DataFrame): A table specifying the compounds, adducts, and peak indices to select.
        method (str): The method of selection. Defaults to "all", which selects all adducts from the selected_compounds_table, or "best" for best adducts

    Returns:
        pd.DataFrame: A filtered dataset containing only the selected compounds.
    """
    if len(selected_compounds_table) == 0:
        print("No compounds selected for annotation.")
        return pd.DataFrame()

    if 'smiles' in full_dataset.columns:
        full_dataset['inchi'] = full_dataset['smiles'].apply(lambda row: AllChem.MolToInchi(AllChem.MolFromSmiles(row)))
        full_dataset['inchi_key'] = full_dataset['inchi'].apply(inchi_to_inchikey)
        full_dataset['neutralized_inchi'] = full_dataset['inchi'].apply(neutralize_inchi)
        full_dataset['neutralized_inchi_key'] = full_dataset['neutralized_inchi'].apply(inchi_to_inchikey)
        full_dataset['permanent_charge'] = full_dataset['neutralized_inchi'].apply(charge_from_inchi)
        full_dataset['formula'] = full_dataset['neutralized_inchi'].apply(formula_from_inchi)
        full_dataset['monoisotopic_mass'] = full_dataset['neutralized_inchi'].apply(monoisotopic_mass_from_inchi)
        full_dataset['collision_energy'] = full_dataset['standard_lcmsrun'].apply(get_collision_energy)
        if not full_dataset[full_dataset['inchi_key'] != full_dataset['neutralized_inchi_key']].empty:
            print(f"Warning! The InChIKey and neutralized InChIKey do not match for the following selected compounds:")
            print(full_dataset[full_dataset['inchi_key'] != full_dataset['neutralized_inchi_key']][['compound_name', 'inchi_key', 'neutralized_inchi_key']])

    select_dataset = pd.DataFrame()
    for _, row in selected_compounds_table.iterrows():
        compound_name = row['compound_name']
        standard_lcmsrun = row['standard_lcmsrun']
        selected_adducts = row['selected_adducts']
        selected_peak_indices = row['selected_peak_indices']
        best_adduct = row['best_adduct']
        
        # Filter by compound name and standard_lcmsrun
        if selected_adducts != ["Ambiguous"] and selected_peak_indices != ["Ambiguous"]:
            if 'peak_index' in full_dataset.columns:
                mask = (full_dataset['compound_name'] == compound_name) & \
                            (full_dataset['standard_lcmsrun'] == standard_lcmsrun) & \
                            (full_dataset['adduct'].isin(selected_adducts)) & \
                            (full_dataset['peak_index'].isin(selected_peak_indices))
            else: # EIC doesn't contain peak index
                mask = (full_dataset['compound_name'] == compound_name) & \
                            (full_dataset['standard_lcmsrun'] == standard_lcmsrun) & \
                            (full_dataset['adduct'].isin(selected_adducts))
            if method == "best":
                # Handle missing or NaN best_adduct
                if not best_adduct or (isinstance(best_adduct, float) and pd.isna(best_adduct)):
                    continue  # Skip this row if best_adduct is not set
                # If best_adduct is a string, convert to list
                if isinstance(best_adduct, str):
                    best_adduct = [best_adduct]
                mask = mask & (full_dataset['adduct'].isin(best_adduct))

            selected_compounds = full_dataset[mask].copy()
            if not selected_compounds.empty:
                select_dataset = pd.concat([select_dataset, selected_compounds], ignore_index=True)

        elif selected_adducts == ["Ambiguous"] and selected_peak_indices == ["Ambiguous"]:
            mask = (full_dataset['compound_name'] == compound_name) & \
                   (full_dataset['standard_lcmsrun'] == standard_lcmsrun)
            
            selected_compounds = full_dataset[mask].copy()
            # This takes care of ambiguous compounds
            ambiguous_values = {
                'adduct': "Ambiguous",
                'peak_index': "Ambiguous",
                'rt_peak': np.nan,
                'rt': np.nan,
                'mz': np.nan,
                'i': np.nan,
                'label': lambda df: df['compound_name'],
                'spectrum': np.nan,
                'intensity': np.nan,
                'total_intensity_fraction': np.nan,
                'total_intensity': np.nan,
                'precursor_peak_height': np.nan,
                'precursor_mz': np.nan,
                'mz_observed': np.nan,
                'mz_theoretical': np.nan,
                'ppm_error': np.nan,
            }
            for col, val in ambiguous_values.items():
                if col in selected_compounds:
                    if callable(val):
                        selected_compounds[col] = val(selected_compounds)
                    else:
                        selected_compounds[col] = val
            if not selected_compounds.empty:
                select_dataset = pd.concat([select_dataset, selected_compounds], ignore_index=True)

    key_cols = [col for col in ['compound_name', 'standard_lcmsrun', 'adduct', 'peak_index'] if col in select_dataset.columns]
    select_dataset.drop_duplicates(subset=key_cols, keep='first', inplace=True)
    return select_dataset


def generate_molecular_image(
    smiles_string: str, 
    annotation: Optional[str] = None
) -> Optional[str]:
    """
    Generate an SVG representation of a molecule from its SMILES string, with optional annotation.

    Args:
        smiles_string (str): The SMILES string of the molecule.
        annotation (Optional[str]): Optional annotation text to include in the SVG. Defaults to None.

    Returns:
        Optional[str]: An SVG string representing the molecule, or None if the SMILES string is invalid.
    """
    mol = MolFromSmiles(smiles_string)
    if mol:
        inchi = Chem.MolToInchi(mol)  # Generate the InChI string
        if annotation is None:
            annotation_text = inchi_to_inchikey(inchi)  # Generate the InChIKey
        else:
            annotation_text = annotation
        svg = mol_to_svg(mol, molSize=(500, 500), kekulize=True, annotation=annotation_text)  # Pass the InChI as annotation
        return svg
    return None

def mol_to_svg(
    mol: Mol, 
    molSize: tuple[int, int] = (300, 300), 
    kekulize: bool = True, 
    annotation: Optional[str] = None
) -> str:
    """
    Generate an SVG representation of a molecule with optional annotation text (e.g., InChI string).

    Args:
        mol (Mol): RDKit molecule object.
        molSize (tuple[int, int]): Size of the molecule image in pixels (width, height). Defaults to (300, 300).
        kekulize (bool): Whether to kekulize the molecule structure. Defaults to True.
        annotation (Optional[str]): Optional annotation text to include below the molecule. Defaults to None.

    Returns:
        str: An SVG string representing the molecule.
    """
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1] + 50)  # Increase height for annotation
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # Add annotation (e.g., InChI string) below the molecule
    if annotation:
        svg_lines = svg.splitlines()
        
        # Define font size
        font_size = 25
        
        # Function to estimate text width
        def estimate_text_width(text: str, font_size: int = font_size) -> float:
            # Approximate width calculation (assuming average character width)
            avg_char_width = font_size / 2  # pixels per character
            return len(text) * avg_char_width
        
        max_width = molSize[0] - 20  # Leave some margin
        current_line = ""
        lines = []
        
        for char in annotation:
            test_line = current_line + char
            if estimate_text_width(test_line, font_size=font_size) > max_width:
                lines.append(current_line)
                current_line = char
            else:
                current_line = test_line
        
        if current_line:
            lines.append(current_line)
        
        # Insert each line into the SVG
        y_position = molSize[1] - 30
        line_height = font_size * 1.2  # Line height is 1.2 times the font size
        
        for line in lines:
            svg_lines.insert(-1, f'<text x="10" y="{y_position}" font-family="sans-serif" font-size="{font_size}">{line}</text>')
            y_position += line_height  # Move to the next line
        
        svg = '\n'.join(svg_lines)

    return svg.replace('svg:', '')


def generate_gridded_molecular_images(lcmsruns_table: pd.DataFrame) -> dict:
    """
    Generate molecular structure image grids for unique SMILES strings in each 'run_num' group,
    with the compound_name overlaid on each image.

    Args:
        lcmsruns_table (pd.DataFrame): DataFrame containing 'smiles' and 'compound_name' columns.

    Returns:
        dict: A dictionary where keys are 'run_num' and values are base64-encoded grids of molecular structure images.
    """
    lcmsruns_table['run_num'] = lcmsruns_table['standard_lcmsrun'].apply(get_run_num)
    grouped_lcmsruns_table = lcmsruns_table.groupby(['run_num'])

    runnum_to_structure_image_grid = {}
    unique_runs = lcmsruns_table['run_num'].dropna().unique()

    for runnum in unique_runs:
        group = grouped_lcmsruns_table.get_group((runnum,))
        unique_smiles = group['smiles'].unique()
        images = []

        for smiles in unique_smiles:
            try:
                mol = MolFromSmiles(smiles)
                if mol:
                    # Generate the image for the molecule
                    img = Draw.MolToImage(mol, size=(200, 200))

                    # Overlay the compound_name on the image
                    draw = ImageDraw.Draw(img)
                    compound_name = group.loc[group['smiles'] == smiles, 'compound_name'].iloc[0]
                    precursor_mz = round(group.loc[group['smiles'] == smiles, 'precursor_mz'].iloc[0], 4)
                    mono_mass = round(group.loc[group['smiles'] == smiles, 'exact_mass'].iloc[0], 4)
                    image_title = f"{compound_name}\n(MIM: {mono_mass})"

                    # Dynamically adjust font size
                    font_size = 20
                    min_font_size = 5  # Set a minimum font size
                    while font_size >= min_font_size:
                        try:
                            font = ImageFont.truetype("arial.ttf", font_size)
                        except IOError:
                            font = ImageFont.load_default()  # Fallback to default font if truetype font is unavailable

                        text_bbox = draw.textbbox((0, 0), image_title, font=font)
                        text_width = text_bbox[2] - text_bbox[0]
                        text_height = text_bbox[3] - text_bbox[1]

                        # Check if the text fits within the image width
                        if text_width <= img.width - 10:  # Leave some padding
                            break
                        font_size -= 1  # Reduce font size if it doesn't fit

                    if font_size < min_font_size:
                        print(f"Warning: Text '{image_title}' could not fit within the image.")

                    text_position = ((img.width - text_width) // 2, img.height - text_height - 5)
                    draw.text(text_position, image_title, fill="black", font=font)

                    images.append(img)
            except Exception as e:
                print(f"Error generating image for SMILES: {smiles}, Error: {e}")

        if images:
            # Create a grid of images
            grid_size = int(len(images) ** 0.5) + (1 if len(images) ** 0.5 % 1 > 0 else 0)
            grid_width = grid_size * 200
            grid_height = grid_size * 200
            grid_img = Image.new('RGB', (grid_width, grid_height), (255, 255, 255))

            for idx, img in enumerate(images):
                x_offset = (idx % grid_size) * 200
                y_offset = (idx // grid_size) * 200
                grid_img.paste(img, (x_offset, y_offset))

            # Convert the grid image to base64
            buffered = io.BytesIO()
            grid_img.save(buffered, format="PNG")
            img_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")

            # Store the base64 string in the dictionary
            runnum_to_structure_image_grid[runnum] = img_base64

    return runnum_to_structure_image_grid


def process_data_for_plotting(
    eics_list: List[Dict[str, pd.DataFrame]],
    top_spectra_list: List[pd.DataFrame],
    group_name_list: List[Tuple[str, str, str]],
    rt_peak_list: List[pd.DataFrame],
    config: Dict[str, Any]
) -> List[Dict[str, Any]]:
    """
    Processes data for plotting by combining EICs, spectra, and metadata for each group.

    Args:
        eics_list (List[Dict[str, pd.DataFrame]]): List of EIC data for each group.
        top_spectra_list (List[pd.DataFrame]): List of top spectra data for each group.
        group_name_list (List[Tuple[str, str, str]]): List of group metadata (compound name, file path, SMILES).
        rt_peak_list (List[pd.DataFrame]): List of retention time peak data for each group.
        config (Dict[str, Any]): Configuration dictionary containing metadata for ref std project. Uses:
                                analysis_sort (str): Sorting criteria, either "run" or "specs". Defaults to "run".
                                analysis_compound_subset (str): Compound name to filter the data. Defaults to None.
                                analysis_run_subset (int): Run number to filter the data. Defaults to None.

    Returns:
        List[Dict[str, Any]]: Processed data for each group, including metadata and EICs.
    """
    processed_data = []
    adduct_color = generate_adduct_colors(config['project']['include_adducts'])
    obj_lengths = [len(eics_list), len(top_spectra_list), len(group_name_list), len(rt_peak_list)]

    if len(set(obj_lengths)) != 1:
        print(f"Warning: Lists have inconsistent lengths: {obj_lengths}")
    number_of_groups = obj_lengths[0]

    for i in range(number_of_groups):
        eics = eics_list[i]
        top_spectra = top_spectra_list[i]
        group_name = group_name_list[i]
        rt_peaks = rt_peak_list[i]

        # Extract group-specific information
        compound_name = group_name[0]
        group_file = group_name[1]
        compound_smiles = group_name[2]
        group_run_number = get_run_num(group_file)
        group_chrom = get_chromatography(group_file)
        group_pol = get_file_polarity(group_file)
        group_params = get_lcmsrun_params(group_file)

        # Unique identifier for the group
        group_id = f"{compound_name}_{group_chrom}_{group_pol}_{group_params}_{group_run_number}"
        unique_id = f"{compound_name};;{group_file}"

        # Append processed data for further use
        processed_data.append({
            "group_id": group_id,
            "eics": eics,
            "top_spectra": top_spectra,
            "rt_peaks": rt_peaks,
            "compound_name": compound_name,
            "compound_smiles": compound_smiles,
            "group_file": group_file,
            "unique_id": unique_id,
            "group_run_number": group_run_number,
            "group_chrom": group_chrom,
            "group_pol": group_pol,
            "group_params": group_params,
            "adduct_color": adduct_color
        })

    if config['analysis']['analysis_sort'] == "run":
        processed_data.sort(key=lambda x: (x['group_run_number']))
    elif config['analysis']['analysis_sort'] == "specs":
        processed_data.sort(key=lambda x: (x['compound_name'], x['group_chrom'], x['group_params'], x['group_pol']))

    if config['analysis']['analysis_compound_subset'] is not None:
        processed_data = [
            entry for entry in processed_data
            if entry.get('compound_name') == config['analysis']['analysis_compound_subset']
        ]
    if config['analysis']['analysis_run_subset'] is not None:
        processed_data = [
            entry for entry in processed_data
            if entry.get('group_run_number') == config['analysis']['analysis_run_subset']
        ]

    return processed_data


def extract_selected_compounds(
    selected_dict: Dict[str, List[str]],
    best_adducts_dict: Dict[str, List[str]]
) -> pd.DataFrame:
    """
    Extracts selected compounds and merges them with top adducts.

    Args:
        selected_dict (Dict[str, List[str]]): Dictionary of selected compounds with adduct peaks.
        best_adducts_dict (Dict[str, List[str]]): Dictionary of top adducts for each compound.

    Returns:
        pd.DataFrame: Merged table of selected compounds and top adducts.
    """
    if len(selected_dict) == 0:
        print("Warning: No selected adducts for any compounds!")
        return
    
    else:
        selected_compounds_table = pd.DataFrame({
            'index': selected_dict.keys(),
            'selected_adduct_peaks': selected_dict.values()
        }).reset_index(drop=True)
        selected_compounds_table[['compound_name', 'standard_lcmsrun']] = selected_compounds_table['index'].str.split(';;', expand=True)
        selected_compounds_table['selected_adducts'] = selected_compounds_table['selected_adduct_peaks'].apply(
            lambda x: [item.split('||')[0] for item in x]
        )
        selected_compounds_table['selected_peak_indices'] = selected_compounds_table['selected_adduct_peaks'].apply(
            lambda x: [item.split('||')[1] for item in x]
        )
        selected_compounds_table = selected_compounds_table.drop(columns=['index', 'selected_adduct_peaks'])

        best_adducts_table = pd.DataFrame({
            'index': best_adducts_dict.keys(),
            'best_adduct': best_adducts_dict.values()
        }).reset_index(drop=True)
        best_adducts_table[['compound_name', 'standard_lcmsrun']] = best_adducts_table['index'].str.split(';;', expand=True)
        best_adducts_table = best_adducts_table.drop(columns=['index'])

        merged_table = pd.merge(selected_compounds_table, best_adducts_table, on=['compound_name', 'standard_lcmsrun'], how='left')

    return merged_table


def extract_ambiguous_compounds(
    ambiguous_dict: Dict[str, Any]
) -> pd.DataFrame:
    """
    Extracts ambiguous compounds from a dictionary and formats them into a DataFrame.

    Args:
        ambiguous_dict (Dict[str, Any]): Dictionary of ambiguous compounds.

    Returns:
        pd.DataFrame: DataFrame containing ambiguous compounds with their metadata.
    """
    if len(ambiguous_dict) == 0:
        ambiguous_adducts_table = pd.DataFrame()
    else:
        ambiguous_adducts_table = pd.DataFrame({
            'unique_id': ambiguous_dict.keys(),
            'combined': ambiguous_dict.values()
        }).reset_index(drop=True)

        ambiguous_adducts_table[['compound_name', 'standard_lcmsrun']] = ambiguous_adducts_table['unique_id'].str.split(';;', expand=True)

        ambiguous_adducts_table = ambiguous_adducts_table.drop(columns=['unique_id', 'combined'])

    return ambiguous_adducts_table


def generate_selection_summary_table(
    rt_peaks_data: pd.DataFrame,
    running_notes_dict: Dict[str, str],
    config: Dict[str, Any], 
    timestamp: str
) -> None:
    """
    Create a summary table by combining rt_peaks_all DataFrame with running_notes_dict.

    Parameters:
        rt_peaks_data (pd.DataFrame): DataFrame containing RT peak information.
        running_notes_dict (dict): Dictionary with keys as "compound_name;;standard_lcmsrun"
                                   and values as notes.
        ambiguous_adducts (Dict[str, List[str]]): Dictionary of ambiguous adducts for each compound.
        config (Dict[str, Any]): Configuration dictionary containing metadata used for writing summary to disk.
        timestamp (str): Timestamp to append to the summary filename.
    """

    standard_info = pd.read_csv(config['project']['standards_input_file'], keep_default_na=False, dtype=str)

    for chrom in tqdm(config['project']['include_chromatographies'], desc="Generating summary tables", unit="chromatography"):

        # Create summary table from input template
        compound_list = list(itertools.product(standard_info['compound_name'].unique(), config['project']['include_polarities']))
        compound_table = pd.DataFrame(compound_list, columns=['compound_name', 'polarity'])
        compound_table['chromatography'] = chrom

        # Create a new column in rt_peaks_all to store the notes
        rt_peaks_table = rt_peaks_data[rt_peaks_data['chromatography'] == chrom].copy()
        rt_peaks_table['notes'] = ""

        # Merge rt_peaks_table with summary_table on compound_name and polarity
        summary_table = pd.merge(compound_table, rt_peaks_table, on=['compound_name', 'chromatography', 'polarity'], how='left')
        summary_table['collision_energy'] = summary_table['standard_lcmsrun'].apply(get_collision_energy)
        summary_table['post_review_notes'] = ""

        # Iterate through the running_notes_dict and map notes to the DataFrame
        for key, notes in running_notes_dict.items():
            compound_name, standard_lcmsrun = key.split(";;")
            mask = (summary_table['compound_name'] == compound_name) & \
                (summary_table['standard_lcmsrun'] == standard_lcmsrun)
            summary_table.loc[mask, 'notes'] = notes

        # Format summary table
        summary_table.drop_duplicates(inplace=True)
        summary_table.sort_values(['compound_name', 'polarity'], ascending=[True, True], inplace=True)
        summary_table['ppm_error'] = abs(summary_table['mz_theoretical']-summary_table['mz_observed'])/summary_table['mz_theoretical']*1e6

        # Separate all and best
        all_summary_table = summary_table
        best_summary_table = summary_table[(summary_table['best_adduct'] == True) | (summary_table['best_adduct'].isna())]
        best_summary_table.loc[:, 'adduct'] = best_summary_table['adduct'].fillna("UNIDENTIFIED")
        summary_tables_dict = {'all': all_summary_table, 'best': best_summary_table}

        # Check in best summary table that each compound has at least one adduct in both polarities
        for compound_name in best_summary_table['compound_name'].unique():
            compound_data = best_summary_table[best_summary_table['compound_name'] == compound_name]
            polarities = compound_data['polarity'].unique()
            if len(polarities) < 2:
                print(f"Warning: Compound '{compound_name}' in chromatography '{chrom}' does not have at least one adduct in both polarities (only has {polarities}).")

        for summary_data_name, summary_data in summary_tables_dict.items():
            column_order = ['compound_name', 'adduct', 'best_adduct', 'polarity', 'collision_energy', 'formula', 'neutralized_inchi', 'monoisotopic_mass', 'mz_theoretical', 'mz_observed', 'ppm_error', 'standard_lcmsrun', 'rt_peak', 'intensity', 'notes', 'post_review_notes']
            summary_output = summary_data[column_order]
            summary_output = summary_output.fillna("-")

            # Save the summary table to a CSV file
            summaries_path = os.path.join(config['project']['standards_output_path'], f"{timestamp}_analysis", "summary_tables")
            os.makedirs(summaries_path, exist_ok=True)
            fname = f"{summaries_path}/{summary_data_name}_{chrom}_standards_summary_table_{timestamp}.csv"
            summary_output.to_csv(fname, index=False)
            #print(f"\n{summary_data_name} data for {chrom} summary table saved to: {fname}")

    return


def generate_static_summary_plots(
    processed_data: List[Dict[str, Any]],
    selection_results_dict: Dict[str, Tuple[List[str], List[str], List[bool]]],
    config: Dict[str, Any],
    timestamp: str
) -> None:
    """
    Generate static summary plots for selected compounds, including EICs, spectra, and molecular images.

    Args:
        processed_data (List[Dict[str, Any]]): A list of dictionaries containing processed data for each compound group.
            Each dictionary should include keys such as 'eics', 'rt_peaks', 'top_spectra', 'adduct_color', 'group_id',
            'unique_id', 'group_file', and 'compound_name'.
        selection_results_dict Dict[str, Tuple[List[str], List[str], List[bool]]]: A dictionary mapping unique compound IDs to tuples
            containing a list of selected adduct-peak combinations and a list of booleans indicating whether each adduct is the best choice.
        config (Dict[str, Any]): Configuration dictionary containing metadata for ref std project.
        timestamp (str): Timestamp to append to the summary filename.

    Returns:
        None
    """
    
    print("Writing summary plots for selected compounds...")
    for i in tqdm(range(len(processed_data)), desc="Summary plots", unit=" compound groups"):
        data = processed_data[i]
        eics = data['eics']
        rt_peaks = data['rt_peaks']
        top_spectra = data['top_spectra']
        adduct_color = data['adduct_color']
        group_id = data['group_id']
        unique_id = data['unique_id']
        group_file = data['group_file']
        compound_name = data['compound_name']
        group_display_xmin, group_display_xmax = get_rt_range(data['group_file'], data['group_pol'])

        if rt_peaks.empty:
            continue

        rt_peaks['inchi'] = rt_peaks['smiles'].apply(lambda row: AllChem.MolToInchi(AllChem.MolFromSmiles(row)))
        rt_peaks['inchi_key'] = rt_peaks['inchi'].apply(inchi_to_inchikey)
        rt_peaks['neutralized_inchi'] = rt_peaks['inchi'].apply(neutralize_inchi)
        rt_peaks['neutralized_inchi_key'] = rt_peaks['neutralized_inchi'].apply(inchi_to_inchikey)
        rt_peaks['permanent_charge'] = rt_peaks['neutralized_inchi'].apply(charge_from_inchi)
        rt_peaks['formula'] = rt_peaks['neutralized_inchi'].apply(formula_from_inchi)
        rt_peaks['monoisotopic_mass'] = rt_peaks['neutralized_inchi'].apply(monoisotopic_mass_from_inchi)

        # Subset top_spectra by just the selected adducts and peak indices
        if unique_id in selection_results_dict:
            adducts, peaks, bests = selection_results_dict[unique_id]
            if not adducts[0] == "Ambiguous":
                selected_adduct_peak_pairs = list(zip(adducts, peaks))
                top_spectra_subset = top_spectra[
                    top_spectra.apply(
                        lambda row: (row['adduct'], row['peak_index']) in selected_adduct_peak_pairs,
                        axis=1
                    )
                ]
            elif adducts[0] == "Ambiguous":
                top_spectra_subset = pd.DataFrame()
        else:
            top_spectra_subset = pd.DataFrame()

        # Subset rt_peaks by just the selected adducts and peak indices
        if unique_id in selection_results_dict:
            adducts, peaks, bests = selection_results_dict[unique_id]
            if not (adducts and adducts[0] == "Ambiguous"):
                selected_adduct_peak_pairs = list(zip(adducts, peaks))
                rt_peaks_subset = rt_peaks[
                    rt_peaks.apply(
                        lambda row: (row['adduct'], row['peak_index']) in selected_adduct_peak_pairs,
                        axis=1
                    )
                ]
            else:
                rt_peaks_subset = pd.DataFrame()
        else:
            continue

        num_columns = 3
        num_spectra = len(top_spectra_subset)
        num_spectra_rows = math.ceil(num_spectra / 2)  # Two spectra per row
        if num_spectra_rows == 0:
            num_spectra_rows = 1

        # Generate updated subplot titles
        eic_titles = ["EIC (Sample)", "EIC (Log Scale)", "Chosen Adducts Table"]
        zoomed_eic_titles = ["Zoomed EIC (Sample)", "Zoomed EIC (Log Scale)"]
        spectra_titles = [
            f"Spectrum: {row['adduct']} @ {round(row['rt'], 2)} min" for _, row in top_spectra_subset.iterrows()
        ]
        spectra_titles += [""] * (6 - len(spectra_titles))  # Pad to ensure 6 spectra slots

        # Add molecular structure title to the third row
        molecular_image_title = compound_name
        subplot_titles = eic_titles + zoomed_eic_titles + spectra_titles[:2] + [molecular_image_title] + spectra_titles[2:]

        # Create specs for subplots
        specs = [
            [{"type": "scatter"}, {"type": "scatter"}, {"type": "table", "rowspan": 2}],  # Row 1
            [{"type": "scatter"}, {"type": "scatter"}, None],  # Row 2 (merged with row 1 in column 3)
            [{"type": "scatter"}, {"type": "scatter"}, {"type": "image"}],  # Row 3 (Molecular Image with spectra)
        ]
        for _ in range(num_spectra_rows - 1):  # Add rows for additional spectra
            specs.append([{"type": "scatter"}, {"type": "scatter"}, None])

        # Ensure specs dimensions match the number of rows and columns
        while len(specs) < 3 + num_spectra_rows:
            specs.append([{"type": "scatter"}] * num_columns)

        fig = make_subplots(
            rows=3 + num_spectra_rows,
            cols=num_columns,
            shared_xaxes=False,
            vertical_spacing=0.3 / (3 + num_spectra_rows),
            horizontal_spacing=0.1,
            subplot_titles=subplot_titles,
            specs=specs
        )

        # Add EIC traces for "Sample" and "Sample (Log)"
        if rt_peaks_subset.empty:
            max_rt_zoom = max(rt_peaks['rt_peak']) + 2
            min_rt_zoom = min(rt_peaks['rt_peak']) - 2            
        else:
            max_rt_zoom = max(rt_peaks_subset['rt_peak']) + 0.5
            min_rt_zoom = min(rt_peaks_subset['rt_peak']) - 0.5       
        for lcmsrun_path, eic in eics.items():
            if "injbl" not in lcmsrun_path.lower(): # This removes blank
                for _, eic_row in eic.iterrows():
                    rt_sort = np.argsort(eic_row['rt'])
                    adduct = get_adduct(eic_row['label'])
                    color = adduct_color.get(adduct, 'gray')

                    # Add "Sample" trace
                    fig.add_trace(
                        go.Scatter(
                            x=eic_row['rt'][rt_sort],
                            y=eic_row['i'][rt_sort],
                            mode='lines',
                            name=f"{adduct} Sample",
                            line=dict(color=color)
                        ),
                        row=1,
                        col=1
                    )

                    # Add ZOOMED "Sample" trace
                    fig.add_trace(
                        go.Scatter(
                            x=eic_row['rt'][rt_sort],
                            y=eic_row['i'][rt_sort],
                            mode='lines',
                            name=f"{adduct} Sample",
                            line=dict(color=color)
                        ),
                        row=2,
                        col=1
                    )

                    # Add "Sample (Log)" trace
                    fig.add_trace(
                        go.Scatter(
                            x=eic_row['rt'][rt_sort],
                            y=np.log10(eic_row['i'][rt_sort].astype(float)),
                            mode='lines',
                            name=f"{adduct} Sample (Log)",
                            line=dict(color=color)
                        ),
                        row=1,
                        col=2
                    )

                    # Add ZOOMED "Sample (Log)" trace
                    fig.add_trace(
                        go.Scatter(
                            x=eic_row['rt'][rt_sort],
                            y=np.log10(eic_row['i'][rt_sort].astype(float)),
                            mode='lines',
                            name=f"{adduct} Sample (Log)",
                            line=dict(color=color)
                        ),
                        row=2,
                        col=2
                    )

                    # Add rt_peak markers
                    if not rt_peaks.empty:
                        adduct_peaks = rt_peaks[rt_peaks['adduct'] == adduct]
                        for _, peak_info in adduct_peaks.iterrows():
                            peak_rt = peak_info['rt_peak']
                            peak_intensity = peak_info['intensity']

                            # Add peak marker for "Sample"
                            fig.add_trace(
                                go.Scatter(
                                    x=[peak_rt],
                                    y=[peak_intensity],
                                    mode='markers',
                                    marker=dict(color=color, size=8),
                                    name=f"{adduct} RT Peak",
                                    showlegend=False
                                ),
                                row=1,
                                col=1
                            )

                            # Add ZOOMED peak marker for "Sample"
                            fig.add_trace(
                                go.Scatter(
                                    x=[peak_rt],
                                    y=[peak_intensity],
                                    mode='markers',
                                    marker=dict(color=color, size=8),
                                    name=f"{adduct} RT Peak",
                                    showlegend=False
                                ),
                                row=2,
                                col=1
                            )

                            # Add peak marker for "Sample (Log)"
                            fig.add_trace(
                                go.Scatter(
                                    x=[peak_rt],
                                    y=[np.log10(peak_intensity)],
                                    mode='markers',
                                    marker=dict(color=color, size=8),
                                    name=f"{adduct} RT Peak",
                                    showlegend=False
                                ),
                                row=1,
                                col=2
                            )

                            # Add ZOOMED peak marker for "Sample (Log)"
                            fig.add_trace(
                                go.Scatter(
                                    x=[peak_rt],
                                    y=[np.log10(peak_intensity)],
                                    mode='markers',
                                    marker=dict(color=color, size=8),
                                    name=f"{adduct} RT Peak",
                                    showlegend=False
                                ),
                                row=2,
                                col=2
                            )

                            if unique_id in selection_results_dict:
                                adducts, peaks, bests = selection_results_dict[unique_id]
                                if not (adducts and adducts[0] == "Ambiguous"):
                                    for selected_adduct, peak_index in zip(adducts, peaks):
                                        if adduct == selected_adduct and peak_index == peak_info['peak_index']:
                                            # Add star-shaped marker for "Sample"
                                            fig.add_trace(
                                                go.Scatter(
                                                    x=[peak_rt],
                                                    y=[peak_intensity],
                                                    mode='markers',
                                                    marker=dict(color=color, size=14, symbol='star'),
                                                    name=f"{adduct} Selected RT Peak",
                                                    showlegend=False
                                                ),
                                                row=1,
                                                col=1
                                            )

                                            # Add ZOOMED star-shaped marker for "Sample"
                                            fig.add_trace(
                                                go.Scatter(
                                                    x=[peak_rt],
                                                    y=[peak_intensity],
                                                    mode='markers',
                                                    marker=dict(color=color, size=14, symbol='star'),
                                                    name=f"{adduct} Selected RT Peak",
                                                    showlegend=False
                                                ),
                                                row=2,
                                                col=1
                                            )

                                            # Add star-shaped marker for "Sample (Log)"
                                            fig.add_trace(
                                                go.Scatter(
                                                    x=[peak_rt],
                                                    y=[np.log10(peak_intensity)],
                                                    mode='markers',
                                                    marker=dict(color=color, size=14, symbol='star'),
                                                    name=f"{adduct} Selected RT Peak (Log)",
                                                    showlegend=False
                                                ),
                                                row=1,
                                                col=2
                                            )

                                            # Add ZOOMED star-shaped marker for "Sample (Log)"
                                            fig.add_trace(
                                                go.Scatter(
                                                    x=[peak_rt],
                                                    y=[np.log10(peak_intensity)],
                                                    mode='markers',
                                                    marker=dict(color=color, size=14, symbol='star'),
                                                    name=f"{adduct} Selected RT Peak (Log)",
                                                    showlegend=False
                                                ),
                                                row=2,
                                                col=2
                                            )

                    # Add top_spectra markers with "X"
                    if not top_spectra_subset.empty:
                        adduct_spectra = top_spectra_subset[top_spectra_subset['adduct'] == adduct]
                        for _, spectrum_row in adduct_spectra.iterrows():
                            spectrum_rt = spectrum_row['rt']
                            spectrum_adduct = spectrum_row['adduct']
                            marker_color = adduct_color.get(spectrum_adduct, 'gray')

                            # Find closest point in the current EIC
                            sorted_rt = eic_row['rt'][rt_sort]
                            sorted_intensity = eic_row['i'][rt_sort]
                            
                            # Skip if no intensity data available
                            if len(sorted_rt) == 0 or len(sorted_intensity) == 0:
                                continue
                                
                            # Find the closest RT point in the EIC
                            closest_idx = np.argmin(np.abs(sorted_rt - spectrum_rt))
                            
                            if closest_idx >= len(sorted_rt):
                                # If the index is out of bounds, skip this spectrum
                                print(f"Warning: RT {spectrum_row['rt']} is out of bounds for EIC RT range.")
                                continue
                                
                            raw_intensity = sorted_intensity[closest_idx]

                            # Add "X" marker for "Sample"
                            fig.add_trace(
                                go.Scatter(
                                    x=[spectrum_rt],
                                    y=[raw_intensity],
                                    mode='markers',
                                    marker=dict(color=marker_color, size=15, symbol='x'),
                                    name=f"{adduct} Top Spectra",
                                    showlegend=False
                                ),
                                row=1,
                                col=1
                            )

                            # Add ZOOMED "X" marker for "Sample"
                            fig.add_trace(
                                go.Scatter(
                                    x=[spectrum_rt],
                                    y=[raw_intensity],
                                    mode='markers',
                                    marker=dict(color=marker_color, size=15, symbol='x'),
                                    name=f"{adduct} Top Spectra",
                                    showlegend=False
                                ),
                                row=2,
                                col=1
                            )

                            # Add "X" marker for "Sample (Log)"
                            fig.add_trace(
                                go.Scatter(
                                    x=[spectrum_rt],
                                    y=[np.log10(raw_intensity)],
                                    mode='markers',
                                    marker=dict(color=marker_color, size=15, symbol='x'),
                                    name=f"{adduct} Top Spectra (Log)",
                                    showlegend=False
                                ),
                                row=1,
                                col=2
                            )

                            # Add ZOOMED "X" marker for "Sample (Log)"
                            fig.add_trace(
                                go.Scatter(
                                    x=[spectrum_rt],
                                    y=[np.log10(raw_intensity)],
                                    mode='markers',
                                    marker=dict(color=marker_color, size=15, symbol='x'),
                                    name=f"{adduct} Top Spectra (Log)",
                                    showlegend=False
                                ),
                                row=2,
                                col=2
                            )

        # Constrain the x-axis for the two ZOOMED subplots
        fig.update_xaxes(range=[group_display_xmin, group_display_xmax], row=1, col=1)
        fig.update_xaxes(range=[group_display_xmin, group_display_xmax], row=1, col=2)
        fig.update_xaxes(range=[min_rt_zoom, max_rt_zoom], row=2, col=1)
        fig.update_xaxes(range=[min_rt_zoom, max_rt_zoom], row=2, col=2)

        # Add spectra plots for selected adducts and peak indices
        for i, (_, spectrum_row) in enumerate(top_spectra_subset.iterrows()):
            mz_values = spectrum_row['spectrum'][0]
            i_values = spectrum_row['spectrum'][1]
            adduct = spectrum_row['adduct']
            color = adduct_color.get(adduct, 'gray')

            spectrum_row_idx = 3 + (i // 2)  # Start from row 2
            spectrum_col = (i % 2) + 1  # Alternate between columns 1 and 2

            for mz, intensity in zip(mz_values, i_values):
                fig.add_trace(
                    go.Scatter(
                        x=[mz, mz],
                        y=[0, intensity],
                        mode='lines',
                        line=dict(color=color),
                        showlegend=False
                    ),
                    row=spectrum_row_idx,
                    col=spectrum_col
                )

            fig.add_trace(
                go.Scatter(
                    x=mz_values,
                    y=i_values,
                    mode='markers',
                    marker=dict(color=color, size=6),
                    name=f"Spectrum {i+1}: {adduct}",
                    showlegend=False
                ),
                row=spectrum_row_idx,
                col=spectrum_col
            )

        # Add table with adduct information
        if not rt_peaks_subset.empty and \
            unique_id in selection_results_dict and \
            not (selection_results_dict[unique_id][0] and selection_results_dict[unique_id][0][0] == "Ambiguous"):

            # Unpack the tuple
            adducts, peaks, bests = selection_results_dict[unique_id]

            # Create a dictionary to hold data for each peak
            peak_data = {}
            for idx, peak in rt_peaks_subset.iterrows():
                # Find the index in the selected adduct/peak pairs
                is_best = False
                for i, (a, p, s) in enumerate(zip(adducts, peaks, bests)):
                    if peak['adduct'] == a and peak['peak_index'] == p:
                        is_best = s
                        break
                peak['top'] = "Yes" if is_best else "No"
                peak_id = f"{peak['adduct']} ({peak['peak_index']})"
                peak_data[peak_id] = {
                    "RT Peak": round(peak['rt_peak'], 4),
                    "Peak Idx.": peak['peak_index'],
                    "Pol.": peak['polarity'],
                    "MZ Obs.": round(peak['mz_observed'], 4),
                    "MZ Theo.": round(peak['mz_theoretical'], 4),
                    "MI Mass": round(peak['monoisotopic_mass'], 4),
                    "Chrom.": peak['chromatography'] if 'chromatography' in peak else "N/A",
                    "Formula": peak['formula'],
                    "Best": peak['top'],
                    # "Inchi": peak['inchi'],
                }

            # Prepare table data with properties as rows and peaks as columns
            properties = ["RT Peak", "Peak Idx.", "Pol.", "MZ Obs.", "MZ Theo.", "MI Mass", "Chrom.", "Formula", "Best"]
            table_data = {
                " ": properties,
                **{peak_id: [peak_data[peak_id][prop] for prop in properties] for peak_id in peak_data.keys()}
            }

            fig.add_trace(
                go.Table(
                    columnwidth=[2.75] + [2.5] * len(peak_data.keys()),  # Adjust column widths (1 for Property (" "), 3 for others)
                    header=dict(
                        values=[" "] + [key.split(" (peak")[0] for key in peak_data.keys()],
                        fill_color='white',
                        align='left',
                        font=dict(size=12, color='black', weight='bold'),  # Make header text bold
                        line=dict(color='black', width=1)  # Add black outlines to header
                    ),
                    cells=dict(
                        values=[table_data[" "]] + [table_data[peak_id] for peak_id in peak_data.keys()],
                        fill_color='white',
                        align='left',
                        height=25,
                        font=dict(
                            size=[12] + [12] * len(peak_data.keys()),  # Larger font size for the first column
                            color=['black'] + ['black'] * len(peak_data.keys()),
                            weight=['bold'] + ['normal'] * len(peak_data.keys())  # Bold for the first column
                        ),
                        line=dict(color='black', width=1)  # Add black outlines to cells
                    )
                ),
                row=1,
                col=3
            )

        # Add molecular image
        if rt_peaks_subset.empty:
            smiles = rt_peaks['smiles'].unique()[0] if len(rt_peaks['smiles'].unique()) == 1 else ""
        else:
            smiles = rt_peaks_subset['smiles'].unique()[0] if len(rt_peaks_subset['smiles'].unique()) == 1 else ""
        if smiles:
            # Generate the high-resolution SVG molecular image
            svg_image = generate_molecular_image(smiles)
            
            if svg_image:
                # Convert the SVG content to base64
                encoded_svg = base64.b64encode(svg_image.encode('utf-8')).decode('utf-8')
                
                # Add the molecular image to the figure
                fig.add_layout_image(
                    dict(
                        source=f"data:image/svg+xml;base64,{encoded_svg}",
                        xref="x domain",  # Reference the x domain of the subplot in row 2, col 3
                        yref="y domain",  # Reference the y domain of the subplot in row 2, col 3
                        x=3.25,  # Adjust x position for column 3
                        y=-2.5,  # Adjust y position for row 2
                        sizex=2,  # Adjust size as necessary
                        sizey=2,  # Adjust size as necessary
                        xanchor="center",
                        yanchor="middle"
                    )
                )
                # Hide x and y axes for the image subplot
                fig.update_xaxes(visible=False, row=3, col=3)
                fig.update_yaxes(visible=False, row=3, col=3)

        # Update layout
        fig.update_layout(
            title=dict(
                text=f"{os.path.dirname(group_file)}<br>{os.path.basename(group_file)}<br>{group_id.replace('_', '  |  ')}",
                font=dict(size=14),
                y=0.98,
            ),
            height=1500 + 300 * num_spectra_rows,  # Adjust base height and per-row increment
            width=1300,
            plot_bgcolor="white",
            paper_bgcolor="white",
            legend=dict(
                orientation="h",  # Horizontal legend
                yanchor="bottom",  # Align to the bottom
                y=0,  # Position below the plot
                xanchor="center",  # Center the legend
                x=0.5  # Center horizontally
            )
        )
        

        export_dirname = os.path.join(config['project']["standards_output_path"], f"{timestamp}_analysis", "summary_figures")
        os.makedirs(export_dirname, exist_ok=True)
        fname = f"{export_dirname}/{group_id}_summary_figure_{timestamp}.pdf"
        fig.write_image(
            fname,
            engine="kaleido",
            width=1500,
            height=1000,
            scale=2
        )

def create_interactive_plots(
    processed_data: List[Dict[str, Any]],
    selection_results_dict: Dict[str, Tuple[List[str], List[str], List[bool]]],
    runnum_to_structure_image_grid: Dict[int, str],
    running_notes: Dict[str, str]
) -> None:
    """
    Create interactive plots for visualizing and annotating processed data.

    Args:
        processed_data (List[Dict[str, Any]]): A list of dictionaries containing processed data for each compound group.
            Each dictionary should include keys such as 'eics', 'rt_peaks', 'top_spectra', 'adduct_color', 'group_id',
            'unique_id', 'group_file', and 'compound_name'.
        selection_results_dict (Dict[str, Tuple[List[str], List[str], List[bool]]]): A dictionary mapping unique compound IDs to tuples
            containing a list of selected adduct-peak combinations and a list of booleans indicating whether each adduct is the best choice.
        runnum_to_structure_image_grid (Dict[int, str]): A dictionary mapping run numbers to base64-encoded molecular structure images.
        running_notes (Dict[str, str]): A dictionary mapping unique compound IDs to user notes.
        
    Returns:
        None
    """
    
    # Create a state object to track current values
    class PlotState:
        def __init__(self):
            self.current_index = 0
            self.current_unique_id = processed_data[0]['unique_id'] if processed_data else None
            # Add debug counters
            self.update_count = 0
    
    state = PlotState()

    # Create persistent widgets (those that don't need to be recreated)
    image_toggle = widgets.ToggleButton(
        value=False,
        description='Show Structures',
        tooltip='Toggle to show/hide the compound structure image',
        layout=widgets.Layout(width='150px', margin='5px 0 0 0')
    )
    yaxis_toggle = widgets.ToggleButton(
        value=False,
        description='Shared Y-Axis',
        tooltip='Toggle between unique and shared y-axes for non-log EIC plots',
        layout=widgets.Layout(width='150px', margin='30px 0 0 0')
    )
    next_button = widgets.Button(
        description="Next Group"
    )
    previous_button = widgets.Button(
        description="Previous Group"
    )
    progress_label = widgets.Label(
        value=f"1/{len(processed_data)} Groups Completed"
    )
    navigate_textbox = widgets.Text(
        placeholder='Search...',
        layout=widgets.Layout(width='50px')
    )
    navigate_button = widgets.Button(
        description="Go",
        layout=widgets.Layout(width='50px')
    )
    compound_image_widget = widgets.Image(
        format='png',
        layout=widgets.Layout(
            width='400px',
            height='400px',
            margin='0 0 0 50px',
            display='none'  # Initially hidden
        )
    )
    compound_image_label = widgets.Label(
        value="Compound Structures",
        layout=widgets.Layout(margin='0 0 10px 0')
    )
    notes_textarea = widgets.Textarea(
        placeholder='Enter notes about this compound...',
        disabled=False,
        layout=widgets.Layout(width='750px', height='35px')
    )
    completion_label = widgets.Label(value="", layout=widgets.Layout(margin="0 0 0 0"))
    output_container = widgets.Output()
    
    # Main output area that will be refreshed
    main_output = widgets.Output()
    
    # Define event handlers for persistent widgets
    def on_image_toggle_change(change):
        if image_toggle.value:
            image_toggle.description = 'Hide Structures'
            compound_image_widget.layout.display = 'block'
        else:
            image_toggle.description = 'Show Structures'
            compound_image_widget.layout.display = 'none'
    
    def update_progress_text():
        progress_label.value = f"{state.current_index + 1}/{len(processed_data)} Groups Completed"
    
    def on_toggle_change(change):
        yaxis_toggle.description = 'Shared Y-Axis' if not yaxis_toggle.value else 'Unique Y-Axis'
        update_display()
    
    def next_group(b):
        if state.current_index < len(processed_data) - 1:
            if state.current_unique_id in selection_results_dict:
                adducts, peaks, bests = selection_results_dict[state.current_unique_id]
                if not (adducts and adducts[0] == "Ambiguous") and not any(bests):
                    output_container.clear_output(wait=True)
                    with output_container:
                        completion_label.value = "Please select at least one best adduct for the current compound because Ambiguous is not selected."
                    return
            state.current_index += 1
            state.current_unique_id = processed_data[state.current_index]['unique_id']
            image_toggle.value = False
            completion_label.value = ""
            update_display()
        else:
            completion_label.value = "Analysis completed!"
    
    def previous_group(b):
        if state.current_index > 0:
            state.current_index -= 1
            state.current_unique_id = processed_data[state.current_index]['unique_id']
            image_toggle.value = False
            completion_label.value = ""
            update_display()
        else:
            output_container.clear_output(wait=True)
            with output_container:
                print("Already at the first group.")
    
    def navigate_to_group(b):
        try:
            target_index = None  # Initialize to None
            if navigate_textbox.value.isdigit():  # Index
                target_index = int(navigate_textbox.value) - 1
            else:  # Compound name
                for idx, entry in enumerate(processed_data):
                    if entry.get('unique_id', '').split(';;')[0] == navigate_textbox.value:
                        target_index = idx
                        break
            if target_index is not None and 0 <= target_index < len(processed_data):
                state.current_index = target_index
                state.current_unique_id = processed_data[state.current_index]['unique_id']
                completion_label.value = ""
                update_display()
            else:
                output_container.clear_output(wait=True)
                with output_container:
                   completion_label.value = f"Invalid index or compound name given. Please enter a valid compound name or number between 1 and {len(processed_data)}."
        except ValueError:
            output_container.clear_output(wait=True)
            with output_container:
                completion_label.value = "Invalid input. Please enter a valid integer for index or compound name."

    def on_notes_change(change):
        if state.current_unique_id:
            running_notes[state.current_unique_id] = notes_textarea.value
            
    # Bind event handlers
    image_toggle.observe(on_image_toggle_change, names='value')
    yaxis_toggle.observe(on_toggle_change, names='value')
    next_button.on_click(next_group)
    previous_button.on_click(previous_group)
    navigate_button.on_click(navigate_to_group)
    notes_textarea.observe(on_notes_change, names='value')
    
    # Helper function to directly update dictionaries when checkboxes change
    def create_checkbox_handlers(unique_id, all_adducts_checkboxes, adduct_peak_combinations, ambiguous_checkbox, best_adducts_checkboxes):
        def on_all_adducts_change(change):
            if ambiguous_checkbox.value:
                selection_results_dict[unique_id] = (["Ambiguous"], ["Ambiguous"], [False])
                for checkbox in all_adducts_checkboxes[:-1]:  # Exclude the ambiguous checkbox
                    checkbox.value = False
                    checkbox.disabled = True  # Disable all other checkboxes
            else:
                adducts = []
                peaks = []
                bests = []
                checkbox_dict = {c.description: c for c in all_adducts_checkboxes[:-1]}
                for combo in adduct_peak_combinations:
                    desc = combo['description']
                    if desc in checkbox_dict and checkbox_dict[desc].value:
                        adducts.append(combo['adduct'])
                        peaks.append(combo['peak_index'])
                        bests.append(False)  # Will update below if best
                # Update states for best adducts
                for i, adduct in enumerate(adducts):
                    for checkbox in best_adducts_checkboxes.children:
                        if checkbox.description == adduct and checkbox.value:
                            bests[i] = True
                selection_results_dict[unique_id] = (adducts, peaks, bests)
                # Enable/disable best checkboxes based on current selection
                for checkbox in best_adducts_checkboxes.children:
                    checkbox.disabled = checkbox.description not in adducts
                # Enable all other checkboxes
                for checkbox in all_adducts_checkboxes[:-1]:  # Exclude the ambiguous checkbox
                    checkbox.disabled = False

        def on_best_adducts_change(change):
            adducts, peaks, bests = selection_results_dict.get(unique_id, ([], [], []))
            for i, adduct in enumerate(adducts):
                for checkbox in best_adducts_checkboxes.children:
                    if checkbox.description == adduct:
                        bests[i] = checkbox.value
            selection_results_dict[unique_id] = (adducts, peaks, bests)

        def on_ambiguous_change(change):
            if ambiguous_checkbox.value:
                selection_results_dict[unique_id] = (["Ambiguous"], ["Ambiguous"], [False])
                for checkbox in all_adducts_checkboxes[:-1]:  # Exclude the ambiguous checkbox
                    checkbox.value = False
                    checkbox.disabled = True  # Disable all other checkboxes
            else:
                # Re-evaluate the state of other checkboxes
                on_all_adducts_change(None)
                for checkbox in all_adducts_checkboxes[:-1]:  # Exclude the ambiguous checkbox
                    checkbox.disabled = False  # Enable all other checkboxes

        return on_all_adducts_change, on_best_adducts_change, on_ambiguous_change

    # Layout Definitions
    def create_layout(all_adducts_checkboxes, best_adducts_checkboxes):
        checkbox_layout = widgets.VBox(
            children=[
                widgets.Label(value="Select good adduct(s):"),
                *all_adducts_checkboxes
            ],
            layout=widgets.Layout(
                border='1px solid black',
                padding='10px',
                margin='10px',
                width='300px',
                align_items='flex-start'  # Align items to the start (left)
            )
        )
        
        best_adducts_checkboxes_layout = widgets.VBox(
            children=[
                widgets.Label(value="Select best adduct(s):"),
                best_adducts_checkboxes
            ],
            layout=widgets.Layout(
                border='1px solid black',
                padding='10px',
                margin='10px',
                width='300px',
                align_items='flex-start'
            )
        )
        
        go_to_layout = widgets.HBox(
            [navigate_textbox, navigate_button],
            layout=widgets.Layout(
                justify_content='flex-start',  # Align to the far left
                width='200px',
                spacing='5px',
                margin='20px 0 0 0'  # Add space above the widget
            )
        )
        
        navigate_textbox.description = ""
        navigate_textbox.layout = widgets.Layout(width='100px')  # size of the search box
    
        compound_image_widget.layout.display = 'none'
        image_toggle.layout.margin = '0 0 0 50px'
    
        navigation_buttons_layout = widgets.HBox(
            [
                widgets.VBox([next_button, previous_button]),  # Stack Previous and Next buttons vertically
                image_toggle  # Place the Image Toggle button to the right
            ],
            layout=widgets.Layout(
                justify_content='flex-start',  # Align items to the left
                spacing='10px',  # Add spacing between elements
                margin='0 0 0 0'  # No margin for the navigation buttons
            )
        )
        
        button_layout = widgets.VBox(
            [navigation_buttons_layout, progress_label, go_to_layout, yaxis_toggle],
            layout=widgets.Layout(
                align_items='flex-start',
                spacing='5px'
            )
        )

        compound_image_container = widgets.VBox(
            [compound_image_label, compound_image_widget],
            layout=widgets.Layout(
                align_items='center',  # Center-align the content
                padding='10px',  # Add padding around the container
                border='1px solid lightgray',  # Optional: Add a border for better visibility
                width='500px'  # Ensure the container is slightly wider than the image
            )
        )
        
        notes_container = widgets.VBox(
            [
                widgets.Label(value="Annotation Notes:"),
                notes_textarea
            ],
            layout=widgets.Layout(
                border='1px solid black',
                padding='10px',
                margin='10px',
                width='780px',
                height='100px',
                align_items='flex-start'
            )
        )
        
        top_layout = widgets.HBox(
            [checkbox_layout, best_adducts_checkboxes_layout, button_layout, compound_image_container],
            layout=widgets.Layout(
                align_items='flex-start',
                justify_content='flex-start',
                spacing='10px'
            )
        )
        
        notes_row = widgets.HBox(
            [notes_container],
            layout=widgets.Layout(
                align_items='flex-start',
                justify_content='flex-start',
                margin='10px 0 0 0'
            )
        )
        
        final_layout = widgets.VBox(
            [completion_label, top_layout, notes_row],
            layout=widgets.Layout(
                align_items='flex-start',
                padding='0px',
            )
        )
        return final_layout
    
    # Main display update function that recreates the plots and widgets
    def update_display():
        state.update_count += 1
        main_output.clear_output(wait=True)
        
        # Get current data
        data = processed_data[state.current_index]
        unique_id = data['unique_id']
        state.current_unique_id = unique_id

        # Get data from current group
        eics = data['eics']
        top_spectra = data['top_spectra']
        rt_peaks = data['rt_peaks']
        adduct_color = data['adduct_color']
        group_id = data['group_id']
        unique_id = data['unique_id']
        group_run_number = data['group_run_number']
        group_display_xmin, group_display_xmax = get_rt_range(data['group_file'], data['group_pol'])

        # Update the notes text area with any existing notes
        if unique_id in running_notes:
            notes_textarea.value = running_notes[unique_id]
        else:
            notes_textarea.value = ""
            
        # Extract adduct-peak combinations from rt_peaks and top_spectra
        adduct_to_peaks = {}
        if not rt_peaks.empty:
            # Create a mapping from adducts to peak indices and intensities
            for _, peak_row in rt_peaks.iterrows():
                adduct = peak_row['adduct'] if 'adduct' in peak_row else None
                if adduct:
                    if adduct not in adduct_to_peaks:
                        adduct_to_peaks[adduct] = []
                    adduct_to_peaks[adduct].append({
                        'peak_index': peak_row['peak_index'],
                        'intensity': peak_row['intensity']
                    })

        # Create unique identifiers for each adduct-peak combination
        adduct_peak_combinations = []
        for adduct, peaks in adduct_to_peaks.items():
            max_intensity = max(peak['intensity'] for peak in peaks)
            for peak in peaks:
                # Check if there is an MS2 spectrum for this adduct+peak_index
                if top_spectra.empty:
                    has_ms2 = False
                else:
                    has_ms2 = not top_spectra[
                        (top_spectra['adduct'] == adduct) & 
                        (top_spectra['peak_index'] == peak['peak_index'])
                    ].empty

                # Add a star to the description if MS2 spectrum exists
                description = f"{adduct} ({peak['peak_index']}){' *' if has_ms2 else ''}"
                adduct_peak_combinations.append({
                    'adduct': adduct,
                    'peak_index': peak['peak_index'],
                    'description': description,
                    'max_intensity': max_intensity
                })

        # Sort adduct_peak_combinations by max_intensity in descending order
        adduct_peak_combinations.sort(key=lambda x: x['max_intensity'], reverse=True)

        # Create the summary EIC plot data
        group_run_eics = [
            eic for pdata in processed_data if pdata['group_run_number'] == group_run_number
            for eic in pdata['eics'].values()
        ]
        summary_traces = []
        for eic in group_run_eics:
            # Loop through each row in the eic DataFrame
            for _, eic_row in eic.iterrows():
                # Filter data where intensity is above a threshold
                valid_indices = eic_row['i'] > 1e4
                filtered_rt = eic_row['rt'][valid_indices]
                filtered_i = eic_row['i'][valid_indices]

                if len(filtered_rt) > 0:  # Ensure there are valid points
                    # Sort retention times
                    rt_sort = np.argsort(filtered_rt)
                    adduct = get_adduct(eic_row['label'])  # Extract adduct from the label
                    color = adduct_color.get(adduct, 'gray')  # Default to gray if adduct color is missing
                    label = eic_row['label']

                    # Add a trace for the current adduct
                    summary_traces.append(
                        go.Scatter(
                            x=filtered_rt[rt_sort],
                            y=filtered_i[rt_sort],
                            mode='lines',
                            name=f"{label}",
                            line=dict(color=color),
                            showlegend=False
                        )
                    )

        # Create the figure with subplots
        num_spectra = len(top_spectra)
        if num_spectra == 0:
            num_spectra = 1  # Ensure at least one row for empty top_spectra
        num_columns = 4
        num_spectra_rows = math.ceil(num_spectra / num_columns)

        # Adjust subplot titles and specifications
        subplot_titles = [
            "Sample",
            "Blank",
            f"EIC Summary for Run{group_run_number}",
            "Sample (Log)",
            "Blank (Log)",
            *(f"" if top_spectra.empty else f"{row['adduct']} @ {round(row['rt'], 2)} mins" for _, row in top_spectra.iterrows())
        ]

        specs = [
            [{"type": "scatter"}, {"type": "scatter"}, {"type": "scatter", "rowspan": 2, "colspan": 2}, None],
            [{"type": "scatter"}, {"type": "scatter"}, None, None],
            *[[{"type": "scatter"} for _ in range(4)] for _ in range(num_spectra_rows)]
        ]

        # Ensure there is at least one row for empty top_spectra
        if top_spectra.empty:
            subplot_titles.extend([""] * (num_spectra_rows * num_columns - len(subplot_titles) + 5))
            specs.extend([[{"type": "scatter"} for _ in range(4)] for _ in range(num_spectra_rows - 1)])

        fig = make_subplots(
            rows=2 + num_spectra_rows,
            cols=4,
            shared_xaxes=False,
            shared_yaxes=yaxis_toggle.value,
            vertical_spacing=0.3 / (2 + num_spectra_rows),
            horizontal_spacing=0.1,
            subplot_titles=subplot_titles,
            specs=specs
        )

        # Add fallback traces if top_spectra is empty
        if top_spectra.empty:
            fig.add_trace(
                go.Scatter(
                    x=[],
                    y=[],
                    mode='lines',
                    name="No Spectra Available",
                    line=dict(color='gray'),
                    showlegend=False
                ),
                row=3,
                col=1
            )

        # Add the summary traces to the spanning subplot
        for trace in summary_traces:
            fig.add_trace(trace, row=1, col=3)  # Add to row 1, col 3

        # Add EIC traces for each adduct/peak
        for idx, (lcmsrun_path, eic) in enumerate(eics.items()):
            for i, eic_row in eic.iterrows():
                rt_sort = np.argsort(eic_row['rt'])
                adduct = get_adduct(eic_row['label'])
                color = adduct_color[adduct]
                
                # Determine row and column for the current trace
                row = 1 if idx < 2 else 2
                col = (idx % 2) + 1
                
                # Dynamic facet_name determination 
                if row == 1 and col == 1:
                    facet_name = "Sample"
                elif row == 1 and col == 2:
                    facet_name = "Blank"
                elif row == 2 and col == 1:
                    facet_name = "Sample (Log)"
                elif row == 2 and col == 2:
                    facet_name = "Blank (Log)"

                # Add line traces for raw intensity
                trace_index = len(fig.data)
                fig.add_trace(
                    go.Scatter(
                        x=eic_row['rt'][rt_sort],
                        y=eic_row['i'][rt_sort],
                        mode='lines',
                        name=f"{adduct} {facet_name}",  # Include facet_name in legend
                        line=dict(color=color),
                        showlegend=True
                    ),
                    row=row,
                    col=col
                )

                # Recalculate facet_name for log-transformed traces
                if row + 1 == 2 and col == 1:
                    facet_name = "Sample (Log)"
                elif row + 1 == 2 and col == 2:
                    facet_name = "Blank (Log)"

                # Add line traces for log-transformed intensity
                trace_index = len(fig.data)
                fig.add_trace(
                    go.Scatter(
                        x=eic_row['rt'][rt_sort],
                        y=np.log10(eic_row['i'][rt_sort].astype(float)),
                        mode='lines',
                        name=f"{adduct} {facet_name}",  # Include updated facet_name in legend
                        line=dict(color=color),
                        showlegend=True
                    ),
                    row=row + 1,  # Log traces go to the next row
                    col=col
                )

                # Add peak markers for each peak associated with this adduct
                if not rt_peaks.empty:
                    if facet_name == "Sample" or facet_name == "Sample (Log)":
                        adduct_peaks = rt_peaks[rt_peaks['adduct'] == adduct]
                        for _, peak_info in adduct_peaks.iterrows():
                            peak_rt = peak_info['rt_peak']
                            peak_index = peak_info['peak_index']
                            peak_intensity = peak_info['intensity']

                            # Add marker for raw intensity
                            fig.add_trace(
                                go.Scatter(
                                    x=[peak_rt],
                                    y=[peak_intensity],
                                    mode='markers',
                                    marker=dict(color=color, size=10),
                                    name=f"{adduct} RT {peak_index}",
                                    showlegend=False
                                ),
                                row=row,
                                col=col
                            )

                            # Add marker for log-transformed intensity
                            fig.add_trace(
                                go.Scatter(
                                    x=[peak_rt],
                                    y=[np.log10(peak_intensity)],
                                    mode='markers',
                                    marker=dict(color=color, size=10),
                                    name=f"{adduct} RT {peak_index}",
                                    showlegend=False
                                ),
                                row=row + 1,  # Log traces go to the next row
                                col=col
                            )

                # Add MS2 spectra markers
                if not top_spectra.empty:
                    if facet_name == "Sample" or facet_name == "Sample (Log)":
                        adduct_spectra = top_spectra[top_spectra['adduct'] == adduct]
                        # Remove adduct filtering to show all MS2 spectra
                        for _, spectrum_row in adduct_spectra.iterrows():
                            spectrum_adduct = spectrum_row['adduct']
                            spectrum_peak_index = spectrum_row['peak_index']
                            rounded_rt = round(spectrum_row['rt'], 2)
                            marker_color = adduct_color.get(spectrum_adduct, 'gray')

                            # Find closest point in the current EIC
                            sorted_rt = eic_row['rt'][rt_sort]
                            sorted_intensity = eic_row['i'][rt_sort]
                            
                            # Skip if no intensity data available
                            if len(sorted_rt) == 0 or len(sorted_intensity) == 0:
                                continue
                                
                            # Find the closest RT point in the EIC
                            closest_idx = np.argmin(np.abs(sorted_rt - spectrum_row['rt']))
                            
                            if closest_idx >= len(sorted_rt):
                                # If the index is out of bounds, skip this spectrum
                                print(f"Warning: RT {spectrum_row['rt']} is out of bounds for EIC RT range.")
                                continue
                                
                            raw_intensity = sorted_intensity[closest_idx]
                            log_intensity = np.log10(raw_intensity)
                            
                            # Display marker for all spectra on the raw intensity plot
                            fig.add_trace(
                                go.Scatter(
                                    x=[spectrum_row['rt']],
                                    y=[raw_intensity],
                                    mode='markers',
                                    marker=dict(color=marker_color, symbol='x', size=10),
                                    name=f"MS2: {spectrum_adduct} ({spectrum_peak_index}) @ {rounded_rt}",
                                    showlegend=False
                                ),
                                row=row,
                                col=col
                            )
                            
                            # Display marker for all spectra on the log-transformed plot
                            fig.add_trace(
                                go.Scatter(
                                    x=[spectrum_row['rt']],
                                    y=[log_intensity],
                                    mode='markers',
                                    marker=dict(color=marker_color, symbol='x', size=10),
                                    name=f"MS2: {spectrum_adduct} ({spectrum_peak_index}) @ {rounded_rt}",
                                    showlegend=False
                                ),
                                row=row + 1,  # Log traces go to the next row
                                col=col
                            )

        # Update summary and Sample/Sample(Log) plots with correct axis bounds
        fig.update_xaxes(range=[group_display_xmin, group_display_xmax], row=1, col=3)
        fig.update_xaxes(range=[group_display_xmin, group_display_xmax], row=1, col=1)
        fig.update_xaxes(range=[group_display_xmin, group_display_xmax], row=2, col=1)

        # Add traces for Spectra plots
        if not top_spectra.empty:
            top_spectra_sorted = top_spectra.sort_values(['adduct', 'peak_index'])

            mz_list = [lst[0] for lst in top_spectra_sorted['spectrum'] if isinstance(lst, (list, np.ndarray)) and len(lst) > 0]
            mz_list_flattened = np.concatenate([np.ravel(arr) if isinstance(arr, np.ndarray) else np.array([arr]) for arr in mz_list])
            lowest_mz = np.min(mz_list_flattened)*0.9
            highest_mz = np.max(mz_list_flattened)*1.1

            for i, spectrum_row in enumerate(top_spectra_sorted.iterrows()):
                mz_values = spectrum_row[1]['spectrum'][0]
                i_values = spectrum_row[1]['spectrum'][1]
                adduct = spectrum_row[1]['adduct']
                color = adduct_color[adduct]
                precursor_mz = spectrum_row[1]['precursor_mz']
                peak_index = spectrum_row[1]['peak_index']
                spectrum_title = f"{adduct} ({peak_index}) @ {round(spectrum_row[1]['rt'], 2)} mins"

                # Determine the row and column for this spectrum
                spectrum_row_idx = 3 + (i // num_columns)  # Start after EIC rows
                spectrum_col = (i % num_columns) + 1

                # Update the x-axis range for the current subplot
                fig.update_xaxes(
                    range=[lowest_mz, highest_mz],  # Set x-axis limits
                    row=spectrum_row_idx,
                    col=spectrum_col
                )

                # Add vertical lines for each point
                for mz, intensity in zip(mz_values, i_values):
                    fig.add_trace(
                        go.Scatter(
                            x=[mz, mz],
                            y=[0, intensity],
                            mode='lines',
                            line=dict(color=color),
                            showlegend=False
                        ),
                        row=spectrum_row_idx,
                        col=spectrum_col
                    )

                # Add markers for each point
                fig.add_trace(
                    go.Scatter(
                        x=mz_values,
                        y=i_values,
                        mode='markers',
                        marker=dict(color=color, size=6),
                        name=f"Spectrum {i+1}: {adduct}",
                        showlegend=False
                    ),
                    row=spectrum_row_idx,
                    col=spectrum_col
                )

                # Add a black circle at precursor_mz (y=0)
                fig.add_trace(
                    go.Scatter(
                        x=[precursor_mz],
                        y=[0],
                        mode='markers',
                        marker=dict(color='black', symbol='circle', size=20),
                        name=f"Precursor MZ: {precursor_mz}",
                        showlegend=False
                    ),
                    row=spectrum_row_idx,
                    col=spectrum_col
                )
                
                fig.layout.annotations[5 + i].text = spectrum_title

        # Update layout
        fig_title = (group_id.replace('_', '  |  '))
        fig.update_layout(
            hoverlabel=dict(
                font_size=11,  # Increase font size for better readability
                namelength=-1  # Show the full name without truncation
            ),
            title=dict(text=fig_title, font=dict(size=14), x=0.5, xanchor="center"),
            height=700 + 300 * num_spectra_rows,
            width=1500,
            plot_bgcolor="white",
            paper_bgcolor="white",
            legend=dict(
                orientation="h",
                xanchor="center",
                yanchor="top",
                x=0.5,
                y=-0.2
            )
        )

        # Add black borders and keep gridlines
        fig.update_xaxes(
            showline=True,  # Show axis line
            linewidth=1,  # Set line width
            linecolor="black",  # Set line color to black
            showgrid=True,  # Keep gridlines
            gridcolor="lightgray"  # Set gridline color
        )
        fig.update_yaxes(
            showline=True,  # Show axis line
            linewidth=1,  # Set line width
            linecolor="black",  # Set line color to black
            showgrid=True,  # Keep gridlines
            gridcolor="lightgray"  # Set gridline color
        )

        # Update the compound image based on group_run_number
        if group_run_number in runnum_to_structure_image_grid:
            compound_image_widget.value = base64.b64decode(runnum_to_structure_image_grid[group_run_number])
            compound_image_label.value = f"Structures in Run{group_run_number}"  # Update the label text
        else:
            compound_image_widget.value = b''  # Clear the image if not found
            compound_image_label.value = "No Structures Found"  # Update the label text
        
        # Create new checkboxes for this unique_id
        all_adducts_checkboxes = []
        for combo in adduct_peak_combinations:
            is_selected = False
            if unique_id in selection_results_dict:
                adducts, peaks, bests = selection_results_dict[unique_id]
                # Check for both adduct and peak string match
                for a, p in zip(adducts, peaks):
                    if a == combo['adduct'] and p == combo['peak_index']:
                        is_selected = True
                        break
            checkbox = widgets.Checkbox(
                value=is_selected,
                description=combo['description'],
                disabled=False,
                layout=widgets.Layout(width='300px', margin='0 0 0 -75px')
            )
            all_adducts_checkboxes.append(checkbox)
        
        # Create ambiguous checkbox
        ambiguous_is_selected = unique_id in selection_results_dict and selection_results_dict[unique_id][0] == ["Ambiguous"]
        ambiguous_checkbox = widgets.Checkbox(
            value=ambiguous_is_selected,
            description="Ambiguous",
            disabled=False,
            layout=widgets.Layout(width='300px', margin='0 0 0 -75px')
        )
        all_adducts_checkboxes.append(ambiguous_checkbox)
        
        # Create best_adducts checkboxes
        best_adducts_options = list(dict.fromkeys(combo['adduct'] for combo in adduct_peak_combinations))
        best_adducts_children = []

        if unique_id in selection_results_dict:
            adducts, peaks, bests = selection_results_dict[unique_id]
            enabled_adducts = set(adducts)
            bests_dict = {}
            for a, s in zip(adducts, bests):
                # If multiple peaks for same adduct, keep True if any are True
                if a in bests_dict:
                    bests_dict[a] = bests_dict[a] or s
                else:
                    bests_dict[a] = s
        else:
            enabled_adducts = set()
            bests_dict = {}

        for adduct in best_adducts_options:
            is_selected = bests_dict.get(adduct, False)
            enabled = adduct in enabled_adducts
            checkbox = widgets.Checkbox(
                value=is_selected,
                description=adduct,
                disabled=not enabled,
                layout=widgets.Layout(width='275px', margin='0 0 0 -75px')
            )
            best_adducts_children.append(checkbox)

        best_adducts_checkboxes = widgets.VBox(children=best_adducts_children)
        
        # Create and attach handlers that explicitly use the current unique_id
        on_all_adducts_change, on_best_adducts_change, on_ambiguous_change = create_checkbox_handlers(
            unique_id, all_adducts_checkboxes, adduct_peak_combinations, ambiguous_checkbox, best_adducts_checkboxes
        )
        
        for checkbox in all_adducts_checkboxes[:-1]:  # Exclude the ambiguous checkbox
            checkbox.observe(on_all_adducts_change, names='value')
        
        ambiguous_checkbox.observe(on_ambiguous_change, names='value')
        
        for checkbox in best_adducts_checkboxes.children:
            checkbox.observe(on_best_adducts_change, names='value')
        
        # Create layout with current widgets
        layout = create_layout(all_adducts_checkboxes, best_adducts_checkboxes)
        
        # Update image
        group_run_number = data['group_run_number']
        if group_run_number in runnum_to_structure_image_grid:
            compound_image_widget.value = base64.b64decode(runnum_to_structure_image_grid[group_run_number])
            compound_image_label.value = f"Structures in Run{group_run_number}"
        else:
            compound_image_widget.value = b''
            compound_image_label.value = "No Structures Found"
        
        # Display everything in the main output
        with main_output:
            display(layout)
            display(fig)  # Display the plotly figure
        
        update_progress_text()
    
    # Initial display
    update_display()
    
    # Build the main UI
    ui = widgets.VBox([
        main_output,
        output_container
    ])
    
    # Show the UI
    display(ui)