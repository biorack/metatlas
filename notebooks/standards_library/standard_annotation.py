import pandas as pd
import numpy as np
import os
import glob
import matplotlib.colors as mcolors
from difflib import get_close_matches
from tqdm.notebook import tqdm
import re
import math
import pubchempy as pcp
import base64
import io
import pickle

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

from matchms.filtering.filter_utils.load_known_adducts import load_known_adducts

from scipy.signal import find_peaks

import matplotlib.pyplot as plt
from matplotlib import cm

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, MolFromSmiles
from rdkit.Chem import rdMolDescriptors as Descriptors
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

from typing import Dict, List, Tuple, Union, Any, Optional

##################################
### Get all matching LCMS runs ###
##################################

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
        all_chrom_files = glob.glob(os.path.join(raw_data_dir, row[f'{chrom.lower()}_experiment'], '*.h5'))
        standard_files = [path for path in all_chrom_files if is_matching_group(path, row[f'{chrom.lower()}_group'])]
        standard_files = [path for path in standard_files if is_matching_polarity(path, include_polarities)]
        all_standard_files += standard_files
    return all_standard_files


def build_standard_lcmsrun_table(
    csv_standard_info_path: str,
    include_polarities: List[str] = ['POS', 'NEG'],
    include_chromatographies: List[str] = ['C18', 'HILIC'],
    include_adducts: Optional[List[str]] = None,
    raw_data_dir: str = '/global/cfs/cdirs/metatlas/raw_data/*/'
) -> pd.DataFrame:
    """
    Build a table of LCMS run files that match the specified polarities, chromatographies, and adducts.

    Args:
        csv_standard_info_path (str): Path to the CSV file containing standard information.
        include_polarities (List[str], optional): List of polarities to include (e.g., ['POS', 'NEG']). Defaults to ['POS', 'NEG'].
        include_chromatographies (List[str], optional): List of chromatographies to include (e.g., ['C18', 'HILIC']). Defaults to ['C18', 'HILIC'].
        include_adducts (Optional[List[str]], optional): List of adducts to include. Defaults to None.
        raw_data_dir (str, optional): Path to the directory containing raw LCMS data. Defaults to '/global/cfs/cdirs/metatlas/raw_data/*/'.

    Returns:
        pd.DataFrame: A DataFrame containing the LCMS run files with annotated adducts.
    """
    standard_info = pd.read_csv(csv_standard_info_path)
    standard_info['standard_lcmsruns'] = standard_info.apply(
        lambda row: get_matching_lcmsruns(row, include_polarities, include_chromatographies, raw_data_dir), axis=1
    )
    standard_lcmsruns_table = standard_info.explode('standard_lcmsruns').reset_index(drop=True).rename(
        columns={'standard_lcmsruns': 'standard_lcmsrun'}
    )
    
    standard_lcmsruns_table_with_adducts = build_adduct_annotated_table(
        standard_lcmsruns_table, include_adducts=include_adducts
    )

    return standard_lcmsruns_table_with_adducts

##############################
### Calculate Adduct Table ###
##############################

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


###########################
### Tools for chemistry ###
###########################

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

#############################
### Extract data and plot ###
#############################

def store_in_metatlas_db(cid_not_in_db: Any) -> None:
    """
    Store a compound identification object in the MetAtlas database.

    Args:
        cid_not_in_db (Any): The compound identification object to be stored.
    """
    metob.store(cid_not_in_db)

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

def handle_data(
    mode: str, 
    save_path: str, 
    timestamp: Optional[str] = None, 
    data: Optional[Tuple[Any, ...]] = None, 
    file_suffix: str = "data"
) -> Optional[Union[Tuple[Any, ...], None]]:
    """
    Handles saving or loading data based on the mode and file suffix.

    Args:
        mode (str): Operation mode, either "save" to save data or "load" to load data.
        save_path (str): Directory path to save or load the data.
        timestamp (Optional[str], optional): Timestamp for saving the file. Required for "save" mode. Defaults to None.
        data (Optional[Tuple[Any, ...]], optional): Data to save. Required for "save" mode. Defaults to None.
        file_suffix (str): Suffix for the file name (e.g., "filtered", "rt_correction", "selected"). Defaults to "data".

    Returns:
        Optional[Union[Tuple[Any, ...], None]]: Loaded data as a tuple if mode is "load", or None if mode is "save".
    """
    write_path = os.path.join(save_path, "cache")
    if mode == "save":
        if not os.path.isdir(write_path):
            os.mkdir(write_path)
        if timestamp is None or data is None:
            raise ValueError("Timestamp and data are required for saving.")
        write_fname = os.path.join(write_path, f"{timestamp}_ref_stds_{file_suffix}.pkl")
        print(f"Saving data to: {write_fname}")
        with open(write_fname, 'wb') as f:
            pickle.dump(data, f)
    elif mode == "load":
        pkl_files = glob.glob(os.path.join(write_path, f"*_ref_stds_{file_suffix}.pkl"))
        try:
            most_recent_pkl = max(pkl_files, key=os.path.getmtime)
        except ValueError:
            print(f"No pkl files found in {write_path}/")
            return None
        print(f"Loading most recent pkl file: {most_recent_pkl}")
        with open(most_recent_pkl, 'rb') as f:
            return pickle.load(f)
    else:
        raise ValueError("Invalid mode. Use 'save' or 'load'.")
    
def filter_by_selected(
    eics_full: List[Dict[str, pd.DataFrame]],
    rt_peaks_full: List[pd.DataFrame],
    top_spectra_full: List[pd.DataFrame],
    selected_adducts_dict: Dict[str, List[str]],
    top_adducts_dict: Dict[str, List[str]]
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Filters EICs, RT peaks, and MS2 spectra based on selected adducts and compounds.

    Args:
        eics_full (List[Dict[str, pd.DataFrame]]): List of dictionaries containing EIC data.
        rt_peaks_full (List[pd.DataFrame]): List of DataFrames containing RT peak data.
        top_spectra_full (List[pd.DataFrame]): List of DataFrames containing MS2 spectra data.
        selected_adducts_dict (Dict[str, List[str]]): Dictionary of selected adducts for each compound.
        top_adducts_dict (Dict[str, List[str]]): Dictionary of top adducts for each compound.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: Filtered EICs, RT peaks, and MS2 spectra DataFrames.
    """
    selected_compounds_table = extract_selected_compounds(selected_adducts_dict, top_adducts_dict)

    eics_selected = pd.concat(
        [df.assign(key=key) for d in eics_full for key, df in d.items()],
        ignore_index=True
    ).rename(columns={'key': 'standard_lcmsrun'})
    eics_selected['compound_name'] = eics_selected['label'].apply(lambda x: x.split('_')[0])
    eics_selected = select_compounds_from_gui(eics_selected, selected_compounds_table)
        
    rt_peaks_selected = pd.concat(rt_peaks_full).rename(columns={'lcmsrun': 'standard_lcmsrun'})
    rt_peaks_selected = select_compounds_from_gui(rt_peaks_selected, selected_compounds_table)

    top_spectra_selected = pd.concat(top_spectra_full, ignore_index=True).rename(columns={'lcmsrun': 'standard_lcmsrun'})
    top_spectra_selected['compound_name'] = top_spectra_selected['label'].apply(lambda x: x.split('_')[0])
    top_spectra_selected = select_compounds_from_gui(top_spectra_selected, selected_compounds_table)

    print(f"\nTotal unique compounds retained: {eics_selected['compound_name'].nunique()}")
    print(f"Total unique compound+adduct entries retained: {eics_selected['label'].nunique()}\n")
    print(f"Total EICs selected: {eics_selected.shape[0]}")
    print(f"Total RT peaks selected: {rt_peaks_selected.shape[0]}")
    print(f"Total MS2 spectra selected: {top_spectra_selected.shape[0]}")

    return eics_selected, rt_peaks_selected, top_spectra_selected


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
    
    all_injbl_files = glob.glob(os.path.join(raw_data_dir, f'*{injbl_pattern}*.h5'))
    injbl_run_nums = {get_run_num(injbl_path): injbl_path for injbl_path in all_injbl_files}
    
    closest_run_num = max((run_num for run_num in injbl_run_nums if run_num <= lcmsrun_num), default=None)
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
    lcmsrun_data = ft.df_container_from_metatlas_file(lcmsrun_path, desired_key=f'ms1_{polarity.lower()}')
    
    rt_min = 0.0
    rt_max = round(lcmsrun_data.rt.max(), 2)
    return rt_min, rt_max


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
    
    rt_min, rt_max = get_rt_range(first_row.standard_lcmsrun, lcmsrun_polarity)
    atlas = create_compound_atlas(group, ppm_tolerance, rt_min, rt_max, lcmsrun_polarity, extra_time=0.0)
    
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

            # Filter peaks to retain only the one with the highest intensity if they are within 0.1 RT minutes
            filtered_peaks = []
            for i, peak_index in enumerate(peaks):
                if i == 0 or (sorted_rt[peak_index] - sorted_rt[filtered_peaks[-1]] > 0.1):
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
        rt_peaks = pd.DataFrame()

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
) -> Tuple[List[Dict[str, Union[str, float]]], pd.DataFrame]:
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
        Tuple[List[Dict[str, Union[str, float]]], pd.DataFrame]: A list of dictionaries containing MS1 RT peaks and a DataFrame of top MS2 spectra.
    """
    rt_peaks = []
    if standard_lcmsrun in eics.keys():
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
            
            rt_peaks.append(rt_peak)

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
    method: str = "intensity", 
    ppm_tolerance: int = 5
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
        ppm_tolerance (int, optional): The mass accuracy tolerance in parts per million (ppm). Defaults to 5.

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

    for group_name, group in tqdm(grouped_lcmsruns_table, total=len(grouped_lcmsruns_table), desc=" Extracting data from lcmsruns", unit=' compound group'):
        theoretical_mzs = dict(group['adduct_data'].tolist())
        compound_smiles = group['smiles'].iloc[0]
        adduct_to_polarity = dict(zip(group['adduct'].tolist(), group['polarity'].tolist()))
        group_name = (group_name[0], group_name[1], compound_smiles)

        eics, ms2_data, atlas = collect_eics_and_ms2(group, ppm_tolerance)

        if method == "intensity":
            rt_peaks, top_spectra = get_top_ms1_and_ms2(eics, ms2_data, group_name[1], theoretical_mzs, compound_smiles, adduct_to_polarity)
        elif method == "find_peaks":
            rt_peaks, top_spectra = get_all_ms1_and_ms2(eics, ms2_data, group_name[1], theoretical_mzs, compound_smiles, adduct_to_polarity)
        else:
            raise ValueError("Invalid method. Use 'intensity' or 'find_peaks'.")

        eics_list.append(eics)
        top_spectra_list.append(top_spectra)
        group_name_list.append(group_name)
        rt_peak_list += rt_peaks
        atlas_list.append(atlas)

    return eics_list, top_spectra_list, group_name_list, rt_peak_list, atlas_list, molecular_grid_image


# def extract_and_plot(standard_lcmsruns_table, ppm_tolerance=5, plot_output_dir='./annotation_plots', data_output_dir='./annotation_data', generate_plots=True):
#     standard_lcmsruns_table['closest_injbl'] = standard_lcmsruns_table['standard_lcmsrun'].apply(get_closest_injbl)
#     grouped_lcmsruns_table = standard_lcmsruns_table.groupby(['compound_name', 'standard_lcmsrun'])
    
#     if not os.path.isdir(plot_output_dir):
#         os.mkdir(plot_output_dir)
        
#     if not os.path.isdir(data_output_dir):
#         os.mkdir(data_output_dir)
    
#     all_top_spectra = []
#     all_rt_peaks = []

#     for group_name, group in tqdm(grouped_lcmsruns_table, unit=' Compound LCMSRun Group'):
#         theoretical_mzs = dict(group['adduct_data'].tolist())
#         compound_smiles = group['smiles'].iloc[0]
#         adduct_to_polarity = dict(zip(group['adduct'].tolist(), group['polarity'].tolist()))
        
#         eics, ms2_data = collect_eics_and_ms2(group, ppm_tolerance)
#         rt_peaks, top_spectra = get_top_ms1_and_ms2(eics, ms2_data, group_name[1], theoretical_mzs, compound_smiles, adduct_to_polarity)
        
#         if generate_plots:
#             save_annotation_fig(eics, top_spectra, group_name, plot_output_dir)
        
#         all_top_spectra.append(top_spectra)
#         all_rt_peaks += rt_peaks
    
#     if all_rt_peaks:
#         print("Writing RT peak data")
#         all_rt_peaks = pd.DataFrame(all_rt_peaks)
#         all_rt_peaks.to_csv(os.path.join(data_output_dir, 'rt_peak_annotations.csv'), index=False)
#     if all_top_spectra:
#         print("Writing top spectra data")
#         all_top_spectra = pd.concat(all_top_spectra).reset_index(drop=True)
#         all_top_spectra.to_json(os.path.join(data_output_dir, 'top_intensity_spectra.json'))

#     return all_rt_peaks, all_top_spectra

# def save_annotation_fig(eics, top_spectra, group_name, output_dir):
#     group_chrom = get_chromatography(group_name[1])
#     group_pol = get_file_polarity(group_name[1])
#     group_params = get_lcmsrun_params(group_name[1])
    
#     fig_title = f"{group_name[0]} {group_chrom} {group_pol} {group_params}".replace('/', 'or')
    
#     eic_adducts = []
#     for eic in eics.values():
#         eic['adduct'] = eic.label.apply(lambda x: x.split('_')[-1])
#         eic_adducts += eic.adduct.tolist()

#     eic_adducts = set(eic_adducts)

#     adduct_color = dict()
#     color = iter(cm.rainbow(np.linspace(0, 1, len(eic_adducts))))
#     for adduct in eic_adducts:
#         adduct_color[adduct] = next(color)
        
#     # Pass top_spectra to plot_adduct_eics to add the "X" marker
#     plot_adduct_eics(eics, adduct_color, fig_title, top_spectra)
#     plt.savefig(os.path.join(output_dir, f"{fig_title}_eics.pdf"))
#     plt.close()
        
#     plot_top_spectra(top_spectra, adduct_color, fig_title)
#     plt.savefig(os.path.join(output_dir, f"{fig_title}_top_spectra.pdf"))
#     plt.close()


# def plot_adduct_eics(eics, adduct_color, fig_title, top_spectra):
#     fig, axs = plt.subplots(2, 2, figsize=(20, 9), sharex=True, sharey=False)

#     for idx, (lcmsrun_path, eic) in enumerate(eics.items()):
#         total_eics = eic.shape[0]
#         for eic_idx, eic_row in eic.iterrows():

#             rt_sort = np.argsort(eic_row['rt'])
#             adduct = get_adduct(eic_row['label'])

#             # Plot on the left column for the first lcmsrun, right column for the second
#             ax_raw = axs[0, idx]  # Top row for raw values
#             ax_log = axs[1, idx]  # Bottom row for log-scale values

#             # Plot raw values in the top row
#             line_raw, = ax_raw.plot(eic_row['rt'][rt_sort], eic_row['i'][rt_sort], alpha=0.8, label=adduct, color=adduct_color[adduct])
#             line_log, = ax_log.plot(eic_row['rt'][rt_sort], np.log10(eic_row['i'][rt_sort].astype(float)), alpha=0.8, label=adduct, color=adduct_color[adduct])

#             peak = np.argmax(eic_row['i'][rt_sort])
#             ax_raw.scatter(eic_row['rt'][rt_sort][peak], eic_row['i'][rt_sort][peak], color=adduct_color[adduct])
#             ax_log.scatter(eic_row['rt'][rt_sort][peak], np.log10(eic_row['i'][rt_sort].astype(float))[peak], color=adduct_color[adduct])

#             # Add "X" marker for the retention time in top_spectra
#             if not top_spectra.empty:
#                 for _, spectrum_row in top_spectra.iterrows():
#                     if spectrum_row['adduct'] == adduct:
#                         rounded_rt = round(spectrum_row['rt'], 2)  # Round RT to 2 decimal places

#                         # Get the y-axis intensity value at the given RT
#                         raw_intensity = eic_row['i'][rt_sort][np.searchsorted(eic_row['rt'][rt_sort], spectrum_row['rt'])]
#                         log_intensity = np.log10(raw_intensity)

#                         # Plot the "X" marker at the intensity value
#                         x_marker = ax_raw.scatter(spectrum_row['rt'], raw_intensity, color=adduct_color[adduct], marker='x', s=100, label=f"Spectra RT: {rounded_rt}")
#                         ax_log.scatter(spectrum_row['rt'], log_intensity, color=adduct_color[adduct], marker='x', s=100)

#             lcms_params = get_lcmsrun_params(lcmsrun_path)
#             run_num = get_run_num(lcmsrun_path)

#             # Set subplot titles to the filenames
#             ax_raw.set_title(f"{lcms_params}_{run_num}")

#     # Set shared X axis labels for both rows
#     for ax in axs[1, :]:  # Bottom row
#         ax.set_xlabel('Retention Time (RT)')
#         ax.grid()

#     for ax in axs[0, :]:
#         ax.grid()

#     # Set Y axis labels for both rows
#     axs[0, 0].set_ylabel('Intensity (i)')  # Left, top row
#     axs[1, 0].set_ylabel('Log10 Intensity')  # Left, bottom row

#     axs[1, 0].sharey(axs[1, 1])  
#     axs[0, 1].sharey(axs[0, 0])  

#     # Collect legend handles and labels
#     handles_left, labels_left = axs[0, 0].get_legend_handles_labels()
#     handles_right, labels_right = axs[0, 1].get_legend_handles_labels()

#     # Ensure there are handles and labels before creating legends
#     if handles_left and labels_left:
#         leg1 = axs[0, 0].legend(handles=handles_left, labels=labels_left, loc='upper left', bbox_to_anchor=(-0.4, 1), fontsize=13, markerscale=0.2, handletextpad=0.5)
#         for legobj in leg1.legend_handles:
#             legobj.set_linewidth(8.0)

#     if handles_right and labels_right:
#         leg2 = axs[0, 1].legend(handles=handles_right, labels=labels_right, loc='upper right', bbox_to_anchor=(1.4, 1), fontsize=13, markerscale=0.2, handletextpad=0.5)
#         for legobj in leg2.legend_handles:
#             legobj.set_linewidth(8.0)

#     fig.suptitle(fig_title)
#     plt.tight_layout()


# def plot_top_spectra(top_spectra, adduct_color, fig_title):
#     if top_spectra.empty:
#         plt.figure()  
#         return
    
#     num_columns = 3
#     num_rows = math.ceil(len(top_spectra) / num_columns)
    
#     fig, axes = plt.subplots(nrows=num_rows, ncols=num_columns, figsize=(19.845, 4.5 * num_rows))
#     axes = axes.flatten()

#     for i, row in top_spectra.iterrows():
#         mz_values = row['spectrum'][0]
#         i_values = row['spectrum'][1]

#         markerline, stemlines, baseline = axes[i].stem(mz_values, i_values, basefmt=" ", markerfmt=" ")
#         plt.setp(stemlines, 'color', adduct_color[row.adduct])
#         axes[i].set_title(f"{row['adduct']} RT: {round(row['rt'], 2)}")

#     plt.tight_layout()

# def filter_by_selected_top_adduct(rt_peaks, top_adducts, rt_check_cutoff=0.05):
#     """
#     rt_peaks_filtered_top_adduct = sta.filter_by_selected_top_adduct(rt_peaks_filtered, top_adducts_dict)
#     """
#     unfiltered_rt_peaks = rt_peaks.copy()
#     unfiltered_rt_peaks['label'] = unfiltered_rt_peaks['compound_name']

#     selected_top_adducts = []
#     # Find the row that matches the top_adducts dict keys
#     for key, value in top_adducts.items():
#         label = key.split(';;')[0]
#         standard_lcmsrun = key.split(';;')[1]
        
#         # Filter rows matching the current key and adducts
#         selected_adduct = unfiltered_rt_peaks[
#             (unfiltered_rt_peaks['standard_lcmsrun'] == standard_lcmsrun) &
#             (unfiltered_rt_peaks['label'] == label) &
#             (unfiltered_rt_peaks['adduct'].isin(value))
#         ].copy()
        
#         # Handle multiple peaks per top adduct
#         if selected_adduct.shape[0] > 1:
#             for adduct in selected_adduct['adduct'].unique():
#                 adduct_subset = selected_adduct[selected_adduct['adduct'] == adduct].copy()
#                 if adduct_subset.shape[0] > 1:
#                     display(adduct_subset)
#                     adduct_subset['label'] = adduct_subset.apply(
#                         lambda row: f"{row['label']} ({row['peak_index']})", axis=1
#                     )
#                 selected_top_adducts.append(adduct_subset)
#         else:
#             selected_top_adducts.append(selected_adduct)

#     if selected_top_adducts:
#         selected_top_adducts_df = pd.concat(selected_top_adducts, ignore_index=True)
#     else:
#         selected_top_adducts_df = pd.DataFrame()  # Return an empty DataFrame if no matches are found

#     # Check that polarities and CEs for a given compound+adduct have the same RTs
#     selected_top_adducts_df_grouped = selected_top_adducts_df.groupby(['chromatography', 'label'])
#     for group_name, group in selected_top_adducts_df_grouped:
#         rt_values = group['rt_peak']
#         cutoff = rt_check_cutoff
#         ces = group['collision_energy'].unique()
#         pols = group['polarity'].unique()
#         if rt_values.max() - rt_values.min() <= cutoff:
#             print(f"Group {group_name}: All RT values for {ces} and {pols} are within {cutoff} mins of each other ({round(rt_values.max() - rt_values.min(), 4)}).\n")
#         else:
#             print(f"Group {group_name}: RT values for {ces} and {pols} are NOT within {cutoff} mins of each other ({round(rt_values.max() - rt_values.min(), 4)}).\n")

#     # Group by monoisotopic_mass and identify isomers if present, then check they have the same adduct
#     selected_top_adducts_grouped = selected_top_adducts_df.groupby(['monoisotopic_mass','polarity','chromatography'])
#     grouped_compounds = selected_top_adducts_grouped['compound_name'].nunique()
#     multiple_compounds_per_mim = grouped_compounds[grouped_compounds > 1]
#     if not multiple_compounds_per_mim.empty:
#         # Iterate over each monoisotopic mass with multiple compounds
#         for isomer_mim in multiple_compounds_per_mim.index:
#             isomer_data = selected_top_adducts_df[
#                 (selected_top_adducts_df['monoisotopic_mass'] == isomer_mim[0]) &
#                 (selected_top_adducts_df['polarity'] == isomer_mim[1]) &
#                 (selected_top_adducts_df['chromatography'] == isomer_mim[2])
#             ]
#             # Check if all adducts are the same
#             unique_adducts = isomer_data['adduct'].unique()
#             if len(unique_adducts) == 1:
#                 # All adducts are the same, do nothing
#                 print(f"Note: Found isomers in {isomer_mim[2]} {isomer_mim[1]} mode at {isomer_mim[0]} ({list(isomer_data['label'])}) but they had matching selected adducts {unique_adducts[0]}.")
#                 continue
#             else: # Adducts for isomers do not agree
#                 print(f"\tWarning! Adducts for isomers do not agree. Exiting. See data for monoisotopic mass {isomer_mim[0]}:\n")
#                 display(isomer_data[['label', 'adduct', 'inchi', 'monoisotopic_mass']])
#                 print("\nPlease return to the GUI to select a matching adduct for isomers.")
#                 return
#     else:
#         print("No isomers found in selected top adducts.")
    
#     # Select best collision energy row by intensity for each compound+adduct in selected_top_adducts_df
#     selected_top_adducts_df_grouped = selected_top_adducts_df.groupby(['chromatography', 'polarity', 'label', 'adduct'])
#     selected_top_adducts_df_grouped_best_ces = []
#     for group_name, group in selected_top_adducts_df_grouped:
#         ces = group.sort_values(by='intensity', ascending=False)
#         if len(ces) == 2:
#             selected_top_adducts_df_grouped_best_ces.append(ces.iloc[0])
#         elif len(ces) == 1:
#             selected_top_adducts_df_grouped_best_ces.append(ces.iloc[0])
#         else:
#             print("Unexpected number of collision energies found.")
#     selected_top_adducts_df_best_ces = pd.DataFrame(selected_top_adducts_df_grouped_best_ces)

#     print(f"\nFiltered {unfiltered_rt_peaks.shape[0]} compound peaks to {selected_top_adducts_df_best_ces.shape[0]} peaks by best adduct. Here are the compounds+adducts retained:\n")
#     display(selected_top_adducts_df_best_ces[['label', 'adduct', 'polarity', 'chromatography', 'inchi_key', 'monoisotopic_mass', 'collision_energy']].sort_values(by=['label','adduct']))

#     return selected_top_adducts_df_best_ces

# def filter_by_top_adduct(rt_peaks):
    
#     unfiltered_rt_peaks = rt_peaks.copy()

#     # Find the row with the highest intensity for each group
#     unfiltered_rt_peaks['label'] = unfiltered_rt_peaks['compound_name']
#     group_list = ['chromatography', 'polarity', 'label']
#     idx_max_intensity = unfiltered_rt_peaks.groupby(group_list)['intensity'].idxmax()
#     highest_intensity_row = unfiltered_rt_peaks.loc[idx_max_intensity]

#     # # Filter rows to keep only those with the same adduct as the highest intensity row
#     group_list.extend(['adduct', 'collision_energy'])
#     top_adducts_per_pol = unfiltered_rt_peaks.merge(
#         highest_intensity_row[group_list],
#         on=group_list,
#         how='inner'
#     )

#     # Find all other peaks for the selected adduct
#     top_adducts_per_pol_grouped = top_adducts_per_pol.groupby(group_list)
#     unfiltered_rt_peaks_grouped = unfiltered_rt_peaks.groupby(group_list)
#     all_peaks = []

#     for group_key, _ in top_adducts_per_pol_grouped:
#         # Check if the group_key exists in rt_peaks_grouped
#         if group_key in unfiltered_rt_peaks_grouped.groups:
#             # Retrieve all rows for the matching group
#             matching_rows = unfiltered_rt_peaks_grouped.get_group(group_key)
#             if matching_rows.shape[0] > 1: # Are there multiple peaks per chrom+polarity+compound+adduct+collision_energy?
#                 matching_rows.loc[:,'label'] = matching_rows.apply(lambda row: f"{row['label']} ({row['peak_index']})", axis=1)
#             all_peaks.append(matching_rows)

#     top_adducts_per_pol_allpeaks = pd.concat(all_peaks, ignore_index=True) if all_peaks else pd.DataFrame()

#     # Group by monoisotopic_mass and identify isomers if present
#     top_adducts_per_pol_allpeaks_isomer_grouping = top_adducts_per_pol_allpeaks.groupby(['monoisotopic_mass','polarity','chromatography'])
#     grouped_compounds = top_adducts_per_pol_allpeaks_isomer_grouping['compound_name'].nunique()
#     multiple_compounds_per_mim = grouped_compounds[grouped_compounds > 1]

#     if not multiple_compounds_per_mim.empty:
#         # Iterate over each monoisotopic mass with multiple compounds
#         for isomer_mim in multiple_compounds_per_mim.index:
#             isomer_data = top_adducts_per_pol_allpeaks[
#                 (top_adducts_per_pol_allpeaks['monoisotopic_mass'] == isomer_mim[0]) &
#                 (top_adducts_per_pol_allpeaks['polarity'] == isomer_mim[1]) &
#                 (top_adducts_per_pol_allpeaks['chromatography'] == isomer_mim[2])
#             ]
            
#             # Check if all adducts are the same
#             unique_adducts = isomer_data['adduct'].unique()
#             if len(unique_adducts) == 1:
#                 # All adducts are the same, do nothing
#                 print(f"Note: Found isomers in {isomer_mim[2]} {isomer_mim[1]} mode at {isomer_mim[0]} ({list(isomer_data['label'])}) but they had matching selected adducts {unique_adducts[0]}.")
#                 continue
#             else: # Adducts for isomers do not agree
#                 print(f"Warning! Adducts for isomers do not agree. See data for monoisotopic mass {isomer_mim[0]}:\n")
#                 display(isomer_data[['label', 'adduct', 'inchi', 'monoisotopic_mass']])
#                 print("\nPlease return to the GUI to select a matching adduct for isomers.")
#                 return
    
#     print(f"\nFiltered {unfiltered_rt_peaks.shape[0]} compound peaks to {top_adducts_per_pol_allpeaks.shape[0]} peaks by best adduct. Here are the compounds+adducts retained:\n")
#     display(top_adducts_per_pol_allpeaks[['label', 'adduct', 'polarity', 'chromatography', 'inchi_key', 'monoisotopic_mass']].sort_values(by=['label','adduct']))

#     return top_adducts_per_pol_allpeaks


# def filter_by_top_adduct(rt_peaks, top_adducts_per_pol_allpeaks=None, method="manual"):
#     import pandas as pd
#     from IPython.display import display, clear_output
#     import ipywidgets as widgets

#     unfiltered_rt_peaks = rt_peaks.copy()
#     unfiltered_rt_peaks['label'] = unfiltered_rt_peaks['compound_name']
#     all_peaks = []

#     # Get unique groups
#     unique_groups = unfiltered_rt_peaks[['label', 'polarity', 'chromatography']].drop_duplicates()
#     total_groups = len(unique_groups)
#     current_index = 0

#     # Widgets
#     progress_label = widgets.Label(value=f"Group {current_index + 1} of {total_groups}")
#     next_button = widgets.Button(description="Next Group")
#     output_container = widgets.Output()
#     table_container = widgets.Output()

#     selected_row = None  # Store the selected row index
#     def update_progress():
#         progress_label.value = f"Group {current_index + 1} of {total_groups}"

#     def process_group(index):
#         nonlocal current_index, selected_row
#         selected_row = None  # Reset the selected row for the new group
#         group = unique_groups.iloc[index]
#         label, polarity, chromatography = group['label'], group['polarity'], group['chromatography']

#         # Subset unfiltered_rt_peaks for the current group
#         group_peaks = unfiltered_rt_peaks[
#             (unfiltered_rt_peaks['label'] == label) &
#             (unfiltered_rt_peaks['polarity'] == polarity) &
#             (unfiltered_rt_peaks['chromatography'] == chromatography)
#         ]

#         if group_peaks.empty:
#             print("No peaks found for this group.")
#             return

#         if method == "intensity":
#             # Automatically select the row with the highest intensity
#             idx_max_intensity = group_peaks['intensity'].idxmax()
#             highest_intensity_row = group_peaks.loc[[idx_max_intensity]]
#             all_peaks.append(highest_intensity_row)
#         elif method == "manual":
#             # Manual selection
#             group_peaks = group_peaks.sort_values(by='intensity', ascending=False)
#             with table_container:
#                 clear_output(wait=True)
#                 print(f"\nManual selection for group: {label}, {polarity}, {chromatography}")
#                 display(group_peaks.style.set_table_attributes("style='display:inline'").set_table_styles(
#                     [{'selector': 'table', 'props': [('overflow-x', 'auto'), ('display', 'block')]}]
#                 ))

#             dropdown = widgets.Dropdown(
#                 options=group_peaks.index.tolist(),
#                 description='Select Row:',
#                 style={'description_width': 'initial'}
#             )
#             confirm_button = widgets.Button(description="Confirm Selection")

#             def on_confirm(b):
#                 nonlocal selected_row
#                 selected_row = dropdown.value
#                 with table_container:
#                     clear_output(wait=True)
#                     print(f"Selected row index: {selected_row}")
#                     display(group_peaks.loc[[selected_row]])

#                 # Add the selected row to all_peaks
#                 highest_intensity_row = group_peaks.loc[[selected_row]]
#                 all_peaks.append(highest_intensity_row)

#                 # Debug: Print all_peaks after appending
#                 print("Updated all_peaks:")
#                 print(pd.concat(all_peaks, ignore_index=True))

#             confirm_button.on_click(on_confirm)
#             with table_container:
#                 display(dropdown, confirm_button)

#     def on_next(b):
#         nonlocal current_index
#         if selected_row is not None or method == "intensity":  # Ensure a selection is made before proceeding
#             if current_index < total_groups - 1:
#                 current_index += 1
#                 update_progress()
#                 process_group(current_index)
#             else:
#                 with output_container:
#                     clear_output(wait=True)
#                     print("All groups processed!")
#                     # Combine all processed peaks
#                     if all_peaks:
#                         top_adducts_per_pol_allpeaks[0] = pd.concat(all_peaks, ignore_index=True)
#                     else:
#                         top_adducts_per_pol_allpeaks[0] = pd.DataFrame()  # Return an empty DataFrame if no peaks are selected

#                     # Display final results
#                     with output_container:
#                         clear_output(wait=True)
#                         if top_adducts_per_pol_allpeaks[0].empty:
#                             print("\nNo peaks were retained. Please check your input data or selection process.")
#                         else:
#                             print("Analysis complete!")
#                             #print(f"\nFiltered {unfiltered_rt_peaks.shape[0]} compound peaks to {top_adducts_per_pol_allpeaks[0].shape[0]} peaks by best adduct. Here are the compounds+adducts retained:\n")
#                             #display(top_adducts_per_pol_allpeaks[0][['label', 'adduct', 'polarity', 'chromatography', 'inchi_key', 'monoisotopic_mass']].sort_values(by=['label', 'adduct']))

#     # Attach event handler
#     next_button.on_click(on_next)

#     # Initial display
#     display(progress_label, next_button, output_container, table_container)
#     process_group(current_index)


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

##############################
### Existing data searches ###
##############################

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

    for _, molecule in tqdm(all_molecules.iterrows(), total=len(all_molecules), desc=" Searching for matches in metatlas db", unit=" compound"):
        
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

    if not matches_df.empty:
        print("\nSummary of compounds+adducts already in MSMS refs:\n")
        display(matches_df)
    if not nonmatches_df.empty:
        print("\nThese compounds+adducts are not yet in MSMS refs:\n")
        display(nonmatches_df)

    return matches_df, nonmatches_df

def search_for_matches_in_atlases(
    query_entries: pd.DataFrame, 
    atlases: pd.DataFrame, 
    cutoff: float = 0.8
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Search for matches of query entries in existing atlases using exact and fuzzy matching.

    Args:
        query_entries (pd.DataFrame): DataFrame containing query entries with columns such as 'label', 'inchi', 'inchi_key', and 'adduct'.
        atlases (pd.DataFrame): DataFrame containing atlas data.
        cutoff (float, optional): Similarity cutoff for fuzzy matching. Defaults to 0.8.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: 
            - A DataFrame summarizing matches found in the atlases.
            - A DataFrame summarizing compounds not found in the atlases.
    """
    matches_dict: Dict[str, List[Union[str, List[str]]]] = {}
    nonmatches_dict: Dict[str, pd.Series] = {}

    atlases['compound_name'] = atlases['compound_name'].astype(str)
    atlases['label'] = atlases['label'].astype(str)

    for _, query_row in tqdm(query_entries.iterrows(), total=query_entries.shape[0], desc=" Searching for matches in existing atlases", unit=" compound"):
        query_label = str(query_row['label'])
        query_compound_name = str(query_row['compound_name'])
        query_inchikey = str(query_row['inchi_key'])
        query_inchi = str(query_row['inchi'])
        query_adduct = str(query_row['adduct'])
        query_polarity = str(query_row['polarity'])
        query_chromatography = str(query_row['chromatography'])
        query_unique_id = f"{query_label};;{query_adduct};;{query_polarity};;{query_chromatography}"

        # Check for exact matches in inchi
        matching_rows = atlases[(atlases['inchi'] == query_inchi) & (atlases['adduct'] == query_adduct)]
        if not matching_rows.empty:            
            source_files = matching_rows['source_file'].tolist()
            atlas_entry = matching_rows['inchi'].tolist()
            matches_dict[query_unique_id] = [query_inchi, atlas_entry, source_files]
            continue

        # Check for exact matches in inchi_key
        matching_rows = atlases[(atlases['inchi_key'] == query_inchikey) & (atlases['adduct'] == query_adduct)]
        if not matching_rows.empty:
            source_files = matching_rows['source_file'].tolist()
            atlas_entry = matching_rows['inchi_key'].tolist()
            matches_dict[query_unique_id] = [query_inchi, atlas_entry, source_files]
            continue

        # Check for fuzzy matches in compound_name
        compound_name_matches = get_close_matches(query_compound_name, atlases['compound_name'].tolist(), n=1, cutoff=cutoff)
        if compound_name_matches:
            matching_rows = atlases[(atlases['compound_name'] == compound_name_matches[0]) & (atlases['adduct'] == query_adduct)]
            if not matching_rows.empty:
                source_files = matching_rows['source_file'].tolist()
                atlas_entry = matching_rows['compound_name'].tolist()
                matches_dict[query_unique_id] = [compound_name_matches[0], atlas_entry, source_files]
                continue

        # Check for fuzzy matches in label
        label_matches = get_close_matches(query_label, atlases['label'].tolist(), n=1, cutoff=cutoff)
        if label_matches:
            matching_rows = atlases[(atlases['label'] == label_matches[0]) & (atlases['adduct'] == query_adduct)]
            if not matching_rows.empty:
                source_files = matching_rows['source_file'].tolist()
                atlas_entry = matching_rows['label'].tolist()
                matches_dict[query_unique_id] = [label_matches[0], atlas_entry, source_files]
                continue

        ## Fingerprint similarity?

        # If no match is found, add to nonmatches_dict
        nonmatches_dict[query_unique_id] = query_row

    # Convert matches dictionary to DataFrame
    if len(matches_dict) != 0:
        matches_df = pd.DataFrame(
            [(key, value[0], value[1], value[2]) for key, value in matches_dict.items()],
            columns=['query_unique_id', 'query_to_atlas', 'atlas_matches', 'atlas_source_files']
        )
        print("\nSummary of compounds+adducts already in the atlases:\n")
        display(matches_df)
    else:
        print("\nNone of these compounds+adducts were found in the atlases.\n")
        matches_df = pd.DataFrame()
    if len(nonmatches_dict) != 0:
        nonmatches_df = pd.concat(nonmatches_dict.values(), axis=1).T.reset_index(drop=True)
        print("\nThese compounds+adducts are not yet in any atlases:\n")
        display(nonmatches_df)
    else:
        print("\nAll compounds+adducts are already in the atlases.\n")
        nonmatches_df = pd.DataFrame()

    return matches_df, nonmatches_df

############################
### Atlas/refs retrieval ###
############################

def get_ema_atlas_data(existing_atlases_path: Dict[str, Dict[str, str]]) -> pd.DataFrame:
    """
    Load and combine atlas data from multiple files into a single DataFrame.

    Args:
        existing_atlases_path (Dict[str, Dict[str, str]]): A dictionary where keys are chromatography types,
            and values are dictionaries mapping polarities to file paths.

    Returns:
        pd.DataFrame: A combined DataFrame containing all atlas data, with an additional 'source_file' column
            indicating the source file for each row.
    """
    atlas_dfs = []
    for chrom_type, polarities in existing_atlases_path.items():
        for polarity, file_path in polarities.items():
            if os.path.exists(file_path):  # Ensure the file exists
                df = pd.read_csv(file_path, sep='\t')
                df['source_file'] = os.path.basename(file_path)  # Add the file name as a new column
                atlas_dfs.append(df)

    atlases = pd.concat(atlas_dfs, ignore_index=True)
    return atlases

# def get_existing_atlases(existing_atlases_path: str) -> pd.DataFrame:
#     """
#     Load and combine all atlas files matching a given path pattern into a single DataFrame.
#
#     Args:
#         existing_atlases_path (str): A glob pattern to locate atlas files.
#
#     Returns:
#         pd.DataFrame: A combined DataFrame containing all atlas data, with an additional 'source_file' column
#             indicating the source file for each row.
#     """
#     all_atlases_paths = glob.glob(existing_atlases_path)
#     atlas_dfs = []
#     for df_path in all_atlases_paths:
#         df = pd.read_csv(df_path, sep='\t')
#         df['source_file'] = os.path.basename(df_path)  # Add the file name as a new column
#         atlas_dfs.append(df)
#
#     atlases = pd.concat(atlas_dfs)
#     return atlases

def get_qc_atlas_data(qc_atlases: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """
    Organize QC atlas data by source file.

    Args:
        qc_atlases (pd.DataFrame): A DataFrame containing QC atlas data, with a 'source_file' column
            indicating the source file for each row.

    Returns:
        Dict[str, pd.DataFrame]: A dictionary where keys are source file names and values are DataFrames
            containing the corresponding QC atlas data.
    """
    qc_atlas_data = {}
    for atlas_source in qc_atlases['source_file'].unique():
        qc_atlas = qc_atlases[qc_atlases['source_file'] == atlas_source]
        qc_atlas_data[atlas_source] = qc_atlas

    return qc_atlas_data

def get_msms_refs(msms_refs_path: str) -> pd.DataFrame:
    """
    Load MS/MS reference data from a file.

    Args:
        msms_refs_path (str): Path to the MS/MS reference file.

    Returns:
        pd.DataFrame: A DataFrame containing the MS/MS reference data.
    """
    msms_refs = pd.read_csv(msms_refs_path, sep='\t', index_col=0, low_memory=False)
    return msms_refs

##################
### Formatters ###
##################

def format_for_msms_refs(
    input_df: pd.DataFrame, 
    msms_refs: pd.DataFrame, 
    msms_refs_metadata: Dict[str, Any]
) -> pd.DataFrame:
    """
    Format a DataFrame for MS/MS reference data by enriching metadata and ensuring required columns.

    Args:
        input_df (pd.DataFrame): Input DataFrame containing compound data, including a 'spectrum' column.
        msms_refs (pd.DataFrame): Existing MS/MS reference DataFrame to align columns with.
        msms_refs_metadata (Dict[str, Any]): Metadata dictionary containing fields such as 'ce_type', 
                                             'frag_method', 'instrument_type', 'decimal', and 'msms_refs_prefix'.

    Returns:
        pd.DataFrame: A formatted DataFrame ready for MS/MS reference storage.
    """
    # Remove rows with NaN in the 'spectrum' column and print a warning
    rows_with_nan = input_df[input_df['spectrum'].isna()]
    if not rows_with_nan.empty:
        print("Warning: The following rows were removed due to NaN in the 'spectrum' column:")
        display(rows_with_nan)
        input_df = input_df.dropna(subset=['spectrum'])

    # Add all required columns for MS/MS refs
    output_df = input_df.copy()
    output_df['ce_type'] = msms_refs_metadata['ce_type']
    output_df['ce'] = output_df['standard_lcmsrun'].apply(get_collision_energy)
    output_df['file'] = output_df['standard_lcmsrun'].apply(os.path.basename)
    output_df.rename(columns={'mz_theoretical': 'mz'}, inplace=True)
    output_df = enrich_metadata(output_df, msms_refs_metadata)
    output_df['spectrum'] = output_df['spectrum'].apply(make_text_spectrum)
    output_df = output_df[msms_refs.columns.intersection(output_df.columns)]
    output_df = output_df.reset_index(drop=True)
    output_df.index = range(
        msms_refs.index.max() + 1,
        msms_refs.index.max() + 1 + len(output_df)
    )

    return output_df

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

def format_and_select_top_adducts(
    rt_peaks: pd.DataFrame, 
    top_adducts: Dict[str, List[str]]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Format and filter retention time (RT) peaks to select the top adducts for each compound.

    Args:
        rt_peaks (pd.DataFrame): A DataFrame containing retention time peaks with columns such as 
                                 'compound_name', 'standard_lcmsrun', 'adduct', 'rt_peak', etc.
        top_adducts (Dict[str, List[str]]): A dictionary where keys are compound identifiers 
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
    selected_top_adducts = []
    unfiltered_top_adducts = []

    # Process rows matching the top_adducts dictionary
    for key, value in top_adducts.items():
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
                selected_top_adducts.append(adduct_subset)
        else:
            selected_top_adducts.append(selected_adduct)

    # Handle multiple peaks per top adduct for the unfiltered dataframe
    unfiltered_rt_peaks_grouped = unfiltered_rt_peaks.groupby(['standard_lcmsrun', 'label', 'adduct'])
    for _, group in unfiltered_rt_peaks_grouped:
        if group.shape[0] > 1:
            group = group.copy()
            group['label'] = group.apply(
                lambda row: f"{row['label']} ({row['peak_index']})", axis=1
            )
        unfiltered_top_adducts.append(group)

    # Combine selected rows into a dataframe
    if selected_top_adducts:
        selected_top_adducts_df = pd.concat(selected_top_adducts, ignore_index=True)
    else:
        selected_top_adducts_df = pd.DataFrame()
    # Combine unfiltered rows into a dataframe
    if unfiltered_top_adducts:
        unfiltered_adducts_df = pd.concat(unfiltered_top_adducts, ignore_index=True)
    else:
        unfiltered_adducts_df = pd.DataFrame()

    datasets_dict = {'top': selected_top_adducts_df, 'all': unfiltered_adducts_df}
    
    for dataset_name, dataset in datasets_dict.items():
        print(f"Working on dataset: {dataset_name}")
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
            print("\t\tNo isomers found in unfiltered data.\n")

        if dataset_name == 'top':
            print("\tSelecting best collision energy row by intensity for the top adduct(s) per compound...\n")
            selected_top_adducts_df_grouped = selected_top_adducts_df.groupby(['chromatography', 'polarity', 'label', 'adduct'])
            selected_top_adducts_df_grouped_best_ces = []
            for group_name, group in selected_top_adducts_df_grouped:
                ces = group.sort_values(by='intensity', ascending=False)
                if len(ces) >= 1:
                    selected_top_adducts_df_grouped_best_ces.append(ces.iloc[0])
                else:
                    print(f"\t\tWarning! No collision energy found for {group_name}.")
            selected_top_adducts_df_best_ces = pd.DataFrame(selected_top_adducts_df_grouped_best_ces)

    print(f"Formatted {unfiltered_rt_peaks.shape[0]} total compound peaks.")
    print(f"Filtered {unfiltered_rt_peaks.shape[0]} total compound peaks to {selected_top_adducts_df_best_ces.shape[0]} peaks by best adduct.")

    # Return both the unfiltered dataframe and the subsetted dataframe
    return unfiltered_rt_peaks, selected_top_adducts_df_best_ces


############################
### Graphical interfaces ###
############################

from typing import List, Dict, Any
import pandas as pd

def select_compounds_from_gui(full_dataset: pd.DataFrame, selected_compounds_table: pd.DataFrame) -> pd.DataFrame:
    """
    Select compounds from a dataset based on user-specified criteria, including compound name, LCMS run, adducts, and peak indices.

    Args:
        full_dataset (pd.DataFrame): The full dataset containing compound information.
        selected_compounds_table (pd.DataFrame): A table specifying the compounds, adducts, and peak indices to select.

    Returns:
        pd.DataFrame: A filtered dataset containing only the selected compounds.
    """
    select_dataset = pd.DataFrame()
    
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

    for _, row in selected_compounds_table.iterrows():
        compound_name = row['compound_name']
        standard_lcmsrun = row['standard_lcmsrun']
        selected_adducts = row['selected_adducts']
        selected_peak_indices = row['selected_peak_indices']
        
        # Filter by compound name and standard_lcmsrun
        mask = (full_dataset['compound_name'] == compound_name) & (full_dataset['standard_lcmsrun'] == standard_lcmsrun)
        
        # Further filter by adduct and peak_index if they exist in the dataframe
        if 'adduct' in full_dataset.columns and len(selected_adducts) > 0:
            adduct_mask = full_dataset['adduct'].isin(selected_adducts)
            mask = mask & adduct_mask
            
        if 'peak_index' in full_dataset.columns and len(selected_peak_indices) > 0:
            peak_mask = full_dataset['peak_index'].isin(selected_peak_indices)
            mask = mask & peak_mask
        
        selected_compounds = full_dataset[mask].copy()
        select_dataset = pd.concat([select_dataset, selected_compounds], ignore_index=True)
    
    return select_dataset


def atlas_id_to_df(atlas_unique_id: str) -> pd.DataFrame:
    """
    Retrieve an atlas from the database using its unique ID and convert it into a DataFrame.

    Args:
        atlas_unique_id (str): The unique identifier of the atlas.

    Returns:
        pd.DataFrame: A DataFrame containing compound identification data from the atlas.
    """
    atlas = get_atlas(atlas_unique_id)
    
    atlas_df = []
    for cid in atlas.compound_identifications:
        row = {
            'label': cid.name,
            'adduct': cid.mz_references[0].adduct,
            'mz': cid.mz_references[0].mz,
            'rt_peak': cid.rt_references[0].rt_peak,
            'rt_min': cid.rt_references[0].rt_min,
            'rt_max': cid.rt_references[0].rt_max,
            'inchi': cid.compound[0].inchi,
            'inchi_key': cid.compound[0].inchi_key,
            'mz_tolerance': cid.mz_references[0].mz_tolerance,
            'polarity': cid.mz_references[0].detected_polarity
        }
        atlas_df.append(row)
        
    return pd.DataFrame(atlas_df)


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

def get_atlas_dataframe(atlas_identifier):
    """
    Retrieve the atlas as a DataFrame. If the identifier ends with '.tsv', read it as a file.
    Otherwise, use the atlas_id_to_df function to retrieve it from the database.

    Args:
        atlas_identifier (str): The atlas identifier, either a file path or an atlas ID.

    Returns:
        pd.DataFrame: The atlas as a DataFrame.
    """
    if atlas_identifier.endswith('.tsv'):
        # Read the file as a DataFrame
        return pd.read_csv(atlas_identifier, sep='\t')
    else:
        # Use the atlas_id_to_df function to retrieve the atlas
        return atlas_id_to_df(atlas_identifier)
    
def get_qc_experimental_atlas(
    nonmatches_to_atlases: pd.DataFrame, 
    current_atlases: Dict[str, str], 
    include_istds: bool = False
) -> Dict[str, pd.DataFrame]:
    """
    Generate experimental QC atlases by collecting MS1 data and comparing it with baseline atlases.

    Args:
        nonmatches_to_atlases (pd.DataFrame): DataFrame containing compounds not matched to existing atlases.
        current_atlases (Dict[str, str]): Dictionary mapping chromatography types to baseline atlas file paths or IDs.
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
        print(f"Getting raw QC files for {project}\n")
        qc_files = get_qc_files(project, chrom, include_istds)
        
        print(f"Retrieving baseline {chrom} QC atlas...\n")
        baseline_atlas_df = get_atlas_dataframe(current_atlases[chrom.lower()])
        baseline_qc[chrom] = baseline_atlas_df

        print(f"Collecting QC MS1 data for {chrom}...\n")
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

    return combined_qc

def create_baseline_correction_input(
    compounds_to_correct: pd.DataFrame, 
    baseline_to_experimental_qc: Dict[str, pd.DataFrame]
) -> Dict[str, pd.DataFrame]:
    """
    Create input data for retention time (RT) correction by combining experimental and baseline QC data.

    Args:
        compounds_to_correct (pd.DataFrame): DataFrame containing compounds to correct with columns 
                                             ['label', 'polarity', 'chromatography', 'rt_peak', 'rt_min', 'rt_max'].
        baseline_to_experimental_qc (Dict[str, pd.DataFrame]): Dictionary where keys are chromatography types and 
                                                               values are DataFrames with experimental and baseline QC data.

    Returns:
        Dict[str, pd.DataFrame]: A dictionary where keys are chromatography types and values are DataFrames 
                                 containing combined data for RT correction.
    """
    uncorrected_atlas = compounds_to_correct[['label', 'polarity', 'chromatography', 'rt_peak', 'rt_min', 'rt_max']]

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
    include_chromatographies: List[str]
) -> Dict[str, pd.DataFrame]:
    """
    Perform RT correction for each chromatography in the given dictionary of DataFrames.

    Args:
        baseline_correction_dfs (Dict[str, pd.DataFrame]): A dictionary where keys are chromatography types 
            (e.g., 'C18', 'HILIC') and values are DataFrames containing RT experimental and baseline data.
        include_chromatographies (List[str]): List of chromatography types to include in the correction process.

    Returns:
        Dict[str, pd.DataFrame]: A dictionary where keys are chromatography types and values are corrected DataFrames.
    """
    corrected_dfs: Dict[str, pd.DataFrame] = {}

    for chromatography, df in baseline_correction_dfs.items():
        # Step 1: Filter rows with both rt_peak_baseline and rt_peak_experimental
        fit_data = df.dropna(subset=["rt_peak_baseline", "rt_peak_experimental"])
        fit_data = fit_data[fit_data["polarity"] == "QC"]
        
        # Ensure numeric data types
        fit_data["rt_peak_experimental"] = pd.to_numeric(fit_data["rt_peak_experimental"], errors="coerce")
        fit_data["rt_peak_baseline"] = pd.to_numeric(fit_data["rt_peak_baseline"], errors="coerce")
        
        # Step 2: Generate 2nd order polynomial from valid rows
        coefficients = np.polyfit(fit_data["rt_peak_experimental"], fit_data["rt_peak_baseline"], 2)
        polynomial = np.poly1d(coefficients)

        # Step 3: Apply the polynomial to all rows where rt_experimental is available
        df["rt_peak_corrected"] = df["rt_peak_experimental"].apply(
            lambda x: polynomial(x) if not np.isnan(x) else np.nan
        )

        # Step 4: Compute rt_min_corrected
        df["rt_min_corrected"] = df.apply(
            lambda row: row["rt_peak_corrected"] - row["rt_peak_experimental"] + row["rt_min_experimental"]
            if not np.isnan(row["rt_peak_corrected"])
            else np.nan,
            axis=1
        )

        # Step 5: Compute rt_max_corrected
        df["rt_max_corrected"] = df.apply(
            lambda row: row["rt_peak_corrected"] + row["rt_max_experimental"] - row["rt_peak_experimental"]
            if not np.isnan(row["rt_peak_corrected"])
            else np.nan,
            axis=1
        )

        # Step 6: Check difference between rt_peak_experimental and rt_peak_corrected
        df['rt_diff_experimental_vs_corrected'] = df["rt_peak_experimental"] - df["rt_peak_corrected"]
        large_diff_rows = df[df['rt_diff_experimental_vs_corrected'].abs() > 0.5]
        if not large_diff_rows.empty:
            print(f"Warning: Large differences in some experimental vs predicted RTs for {chromatography} chromatography:\n")
            print(large_diff_rows[['label', 'polarity', 'rt_peak_experimental', 'rt_peak_corrected', 'rt_diff_experimental_vs_corrected']])

        # Store the corrected DataFrame in the result dictionary
        corrected_dfs[chromatography] = df[
            ['label', 'polarity', 'rt_peak_baseline', 'rt_peak_experimental', 
             'rt_peak_corrected', 'rt_min_corrected', 'rt_max_corrected', 'rt_diff_experimental_vs_corrected']
        ]

    for chrom in include_chromatographies:
        chrom = chrom if 'HILIC' not in chrom else 'HILICZ'
        print(f"Backward RT correction data for {chrom}:")
        display(corrected_dfs[chrom])

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
    updated_atlas = nonmatches_to_atlases.copy()

    for chromatography, df in baseline_correction_outputs.items():
        # Merge the two DataFrames on 'label' and 'polarity' to align rows
        df = df.copy()
        df.loc[:, 'chromatography'] = str(chromatography)
        merged_df = updated_atlas.merge(
            df[['label', 'polarity', 'chromatography', 'rt_peak_corrected', 'rt_min_corrected', 'rt_max_corrected']],
            on=['label', 'chromatography', 'polarity'],
            how='left'
        )

        # Substitute the RT values with the corrected ones
        merged_df.loc[:, 'rt_peak'] = merged_df.loc[:, 'rt_peak_corrected'].combine_first(merged_df.loc[:, 'rt_peak'])
        merged_df.loc[:, 'rt_min'] = merged_df.loc[:, 'rt_min_corrected'].combine_first(merged_df.loc[:, 'rt_min'])
        merged_df.loc[:, 'rt_max'] = merged_df.loc[:, 'rt_max_corrected'].combine_first(merged_df.loc[:, 'rt_max'])

        # Drop the corrected columns to clean up
        merged_df.drop(columns=['rt_peak_corrected', 'rt_min_corrected', 'rt_max_corrected'], inplace=True)

        updated_atlas = merged_df
        print(f"Formatted {updated_atlas.shape[0]} RT-corrected compounds for insertion into {chromatography} atlases.")

    return updated_atlas

def update_and_save_atlases(
    ema_atlases: pd.DataFrame,
    nonmatches_to_atlases_rt_corrected: pd.DataFrame,
    current_time: str,
    atlas_save_path: str,
    save_atlas: bool = False
) -> None:
    """
    Update and save EMA atlases by merging new compounds and splitting them back into their original files.

    Args:
        ema_atlases (pd.DataFrame): Existing EMA atlases with compound data.
        nonmatches_to_atlases_rt_corrected (pd.DataFrame): New compounds to be added to the atlases.
        current_time (str): Timestamp to append to the updated atlas filenames.
        atlas_save_path (str): Directory path to save the updated atlases.
        save_atlas (bool, optional): Whether to save the updated atlases to disk. Defaults to False.

    Returns:
        None
    """
    # Add source file information to nonmatches_to_atlases_rt_corrected
    atlas_names = ema_atlases['source_file'].unique()
    nonmatches_to_atlases_rt_corrected_sourced = nonmatches_to_atlases_rt_corrected.copy()
    for atlas_name in atlas_names:
        for index, row in nonmatches_to_atlases_rt_corrected_sourced.iterrows():
            pol = row['polarity'].lower()
            chrom = row['chromatography'].lower().replace('hilicz', 'hilic')
            if pol in atlas_name.lower() and chrom in atlas_name.lower():
                nonmatches_to_atlases_rt_corrected_sourced.at[index, 'source_file'] = atlas_name

    # Format the columns of missing compounds to match the atlas format
    missing_columns = ema_atlases.columns.difference(nonmatches_to_atlases_rt_corrected_sourced.columns)
    for col in missing_columns:
        nonmatches_to_atlases_rt_corrected_sourced[col] = np.nan
    common_columns = nonmatches_to_atlases_rt_corrected_sourced.columns.intersection(ema_atlases.columns)
    nonmatches_to_atlases_rt_corrected_sourced_formatted = nonmatches_to_atlases_rt_corrected_sourced[common_columns]

    # Check if source_file is missing from any compound lines before merging and splitting
    missed_source_files = nonmatches_to_atlases_rt_corrected_sourced_formatted[
        ~nonmatches_to_atlases_rt_corrected_sourced_formatted['source_file'].isin(atlas_names)
    ]
    if not missed_source_files.empty:
        print("Warning: Some rows have a source_file that does not match any atlas_name.")
        print(missed_source_files)
        return

    # Merge the nonmatches_to_atlases_rt_corrected_sourced_formatted with the ema_atlases
    new_ema_atlases = pd.concat([ema_atlases, nonmatches_to_atlases_rt_corrected_sourced_formatted], ignore_index=False).sort_values('rt_peak')

    # Split EMA atlases back to original files based on source_file
    old_ema_atlases_split = {}
    new_ema_atlases_split = {}
    for atlas_name in atlas_names:
        new_atlas_name = atlas_name.replace(".tsv", f"_{current_time}.tsv")
        old_ema_atlases_split[atlas_name] = ema_atlases[ema_atlases.loc[:, 'source_file'] == atlas_name].drop(columns=['source_file'])
        new_ema_atlases_split[atlas_name] = new_ema_atlases[new_ema_atlases.loc[:, 'source_file'] == atlas_name].drop(columns=['source_file'])

        print(f"\nOld atlas: {atlas_name}. Shape {old_ema_atlases_split[atlas_name].shape}.")
        print(f"New atlas: {new_atlas_name}. Shape {new_ema_atlases_split[atlas_name].shape}.")

        # Save the new EMA atlases to the original files
        if save_atlas:
            updated_atlas_dir = f'{atlas_save_path}/updated_EMA_atlases'
            if not os.path.exists(updated_atlas_dir):
                os.makedirs(updated_atlas_dir)
            fname = f'{updated_atlas_dir}/{new_atlas_name}'
            print(f"\tSaving new EMA atlas to disk at {fname}")
            new_ema_atlases_split[atlas_name].to_csv(fname, sep='\t', index=False)


def update_and_save_msms_refs(
    msms_refs: pd.DataFrame,
    rt_peaks_with_spectra: pd.DataFrame,
    msms_refs_save_path: str,
    timestamp: str,
    save_refs: bool = False
) -> None:
    """
    Update and save MS/MS reference data by appending new entries and optionally saving to disk.

    Args:
        msms_refs (pd.DataFrame): Existing MS/MS reference data.
        rt_peaks_with_spectra (pd.DataFrame): New MS/MS reference data to be added.
        msms_refs_save_path (str): Directory path to save the updated MS/MS references.
        timestamp (str): Timestamp to append to the updated MS/MS reference filename.
        save_refs (bool, optional): Whether to save the updated MS/MS references to disk. Defaults to False.

    Returns:
        None
    """
    # Adjust the index of rt_peaks_with_spectra to start where msms_refs leaves off
    start_index = msms_refs.index.max() + 1
    rt_peaks_with_spectra = rt_peaks_with_spectra.reset_index(drop=True)
    rt_peaks_with_spectra.index = range(start_index, start_index + len(rt_peaks_with_spectra))

    # Combine existing and new MS/MS refs
    new_msms_refs = pd.concat([msms_refs, rt_peaks_with_spectra])
    print(f"Existing MSMS refs: {msms_refs.shape}")
    print(f"New MSMS refs: {new_msms_refs.shape}")
    if new_msms_refs.shape[0] != msms_refs.shape[0] + rt_peaks_with_spectra.shape[0]:
        print("Warning! Some new MSMS refs may not have been added correctly.")
        return
    if new_msms_refs.shape[1] != msms_refs.shape[1]:
        print("Warning! Column numbers don't match between existing and new MSMS refs.")
        return

    if save_refs:
        updated_refs_dir = f'{msms_refs_save_path}/updated_MSMS_refs'
        if not os.path.exists(updated_refs_dir):
            os.makedirs(updated_refs_dir)
        fname = f"{updated_refs_dir}/msms_refs_{timestamp}.tab"
        print(f"\nSaving new MSMS refs to disk at {fname}")
        new_msms_refs.to_csv(fname, sep='\t', index=False)
    else:
        print("\nNew MSMS refs not saved to disk. Here is what would be added with the save_refs=True option:")
        print(rt_peaks_with_spectra)

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
    selected_peaks_compounded['compound_name'] = selected_peaks_compounded['label'].str.replace(r" \(peak\d+\)", "", regex=True)
    merge_on = ['compound_name', 'adduct', 'standard_lcmsrun', 'peak_index']
    merge_on_with_spectrum = merge_on + ['spectrum']  # Create a new list with 'spectrum' added
    merged_peaks_spectra = pd.merge(
        selected_peaks_compounded,
        selected_spectra[merge_on_with_spectrum],  # Subset columns correctly
        on=merge_on,
        how='left'
    )
    return merged_peaks_spectra

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
                    image_title = f"{compound_name}\n({precursor_mz})"

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
    include_adducts: Optional[List[str]] = None,
    sort_by: str = "run",
    subset_by_compound: Optional[str] = None,
    subset_by_run: Optional[int] = None
) -> List[Dict[str, Any]]:
    """
    Processes data for plotting by combining EICs, spectra, and metadata for each group.

    Args:
        eics_list (List[Dict[str, pd.DataFrame]]): List of EIC data for each group.
        top_spectra_list (List[pd.DataFrame]): List of top spectra data for each group.
        group_name_list (List[Tuple[str, str, str]]): List of group metadata (compound name, file path, SMILES).
        rt_peak_list (List[pd.DataFrame]): List of retention time peak data for each group.
        include_adducts (Optional[List[str]]): List of adducts to include in the color mapping. Defaults to None.
        sort_by (str): Sorting criteria, either "run" or "specs". Defaults to "run".
        subset_by_compound (Optional[str]): Compound name to filter the data. Defaults to None.
        subset_by_run (Optional[int]): Run number to filter the data. Defaults to None.

    Returns:
        List[Dict[str, Any]]: Processed data for each group, including metadata and EICs.
    """
    processed_data = []
    adduct_color = generate_adduct_colors(include_adducts)
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

    if sort_by == "run":
        processed_data.sort(key=lambda x: (x['group_run_number']))
    elif sort_by == "specs":
        processed_data.sort(key=lambda x: (x['compound_name'], x['group_chrom'], x['group_pol']))

    if subset_by_compound:
        processed_data = [
            entry for entry in processed_data
            if entry.get('compound_name') == subset_by_compound
        ]
    if subset_by_run:
        processed_data = [
            entry for entry in processed_data
            if entry.get('group_run_number') == subset_by_run
        ]

    return processed_data


def extract_selected_compounds(
    selected_dict: Dict[str, List[str]],
    top_adducts_dict: Dict[str, List[str]]
) -> pd.DataFrame:
    """
    Extracts selected compounds and merges them with top adducts.

    Args:
        selected_dict (Dict[str, List[str]]): Dictionary of selected compounds with adduct peaks.
        top_adducts_dict (Dict[str, List[str]]): Dictionary of top adducts for each compound.

    Returns:
        pd.DataFrame: Merged table of selected compounds and top adducts.
    """
    if len(selected_dict) == 0:
        selected_compounds_table = pd.DataFrame()
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

        top_adducts_table = pd.DataFrame({
            'index': top_adducts_dict.keys(),
            'top_adduct': top_adducts_dict.values()
        }).reset_index(drop=True)
        top_adducts_table[['compound_name', 'standard_lcmsrun']] = top_adducts_table['index'].str.split(';;', expand=True)
        top_adducts_table = top_adducts_table.drop(columns=['index'])

        merged_table = pd.merge(selected_compounds_table, top_adducts_table, on=['compound_name', 'standard_lcmsrun'], how='left')

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


def generate_static_summary_plots(
    processed_data: List[Dict[str, Any]],
    selected_good_adducts: Dict[str, List[str]],
    top_adducts: Dict[str, List[str]],
    export_dir: Optional[str]
) -> None:
    """
    Generate static summary plots for selected compounds, including EICs, spectra, and molecular images.

    Args:
        processed_data (List[Dict[str, Any]]): A list of dictionaries containing processed data for each compound group.
            Each dictionary should include keys such as 'eics', 'rt_peaks', 'top_spectra', 'adduct_color', 'group_id',
            'unique_id', 'group_file', and 'compound_name'.
        selected_good_adducts (Dict[str, List[str]]): A dictionary mapping unique compound IDs to lists of selected
            adduct-peak combinations in the format "adduct||peak_index".
        top_adducts (Dict[str, List[str]]): A dictionary mapping unique compound IDs to lists of top adducts.
        export_dir (Optional[str]): Directory path to save the generated plots. If None, plots are not saved.

    Returns:
        None
    """
    
    for i in tqdm(range(len(processed_data)), desc=" Writing summary plots for selected compounds", unit=" compound group"):
        data = processed_data[i]
        eics = data['eics']
        rt_peaks = data['rt_peaks']
        top_spectra = data['top_spectra']
        adduct_color = data['adduct_color']
        group_id = data['group_id']
        unique_id = data['unique_id']
        group_file = data['group_file']
        compound_name = data['compound_name']

        rt_peaks['inchi'] = rt_peaks['smiles'].apply(lambda row: AllChem.MolToInchi(AllChem.MolFromSmiles(row)))
        rt_peaks['inchi_key'] = rt_peaks['inchi'].apply(inchi_to_inchikey)
        rt_peaks['neutralized_inchi'] = rt_peaks['inchi'].apply(neutralize_inchi)
        rt_peaks['neutralized_inchi_key'] = rt_peaks['neutralized_inchi'].apply(inchi_to_inchikey)
        rt_peaks['permanent_charge'] = rt_peaks['neutralized_inchi'].apply(charge_from_inchi)
        rt_peaks['formula'] = rt_peaks['neutralized_inchi'].apply(formula_from_inchi)
        rt_peaks['monoisotopic_mass'] = rt_peaks['neutralized_inchi'].apply(monoisotopic_mass_from_inchi)

        # Subset top_spectra by just the selected adducts and peak indices
        if unique_id in selected_good_adducts:
            selected_peaks = selected_good_adducts[unique_id]
            selected_adduct_peak_pairs = [
                (selected_peak.split("||")[0], selected_peak.split("||")[1])
                for selected_peak in selected_peaks
            ]
            top_spectra_subset = top_spectra[
                top_spectra.apply(
                    lambda row: (row['adduct'], row['peak_index']) in selected_adduct_peak_pairs,
                    axis=1
                )
            ]
        else:
            top_spectra_subset = pd.DataFrame()  # Empty DataFrame if no adducts are selected

        # Subset rt_peaks by just the selected adducts and peak indices
        if unique_id in selected_good_adducts:
            selected_peaks = selected_good_adducts[unique_id]
            selected_adduct_peak_pairs = [
                (selected_peak.split("||")[0], selected_peak.split("||")[1])
                for selected_peak in selected_peaks
            ]
            rt_peaks_subset = rt_peaks[
                rt_peaks.apply(
                    lambda row: (row['adduct'], row['peak_index']) in selected_adduct_peak_pairs,
                    axis=1
                )
            ]
        else:
            continue # Skip if no adduct peaks are selected

        num_columns = 3
        num_spectra = len(top_spectra_subset)
        num_spectra_rows = math.ceil(num_spectra / 2)  # Two spectra per row
        if num_spectra_rows == 0:
            num_spectra_rows = 1

        # Generate updated subplot titles
        eic_titles = ["EIC (Sample)", "EIC (Log Scale)", "Adduct Details"]
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
        max_rt = max(rt_peaks_subset['rt_peak']) + 0.5
        min_rt = min(rt_peaks_subset['rt_peak']) - 0.5       
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

                            if unique_id in selected_good_adducts:
                                selected_peaks = selected_good_adducts[unique_id]
                                for selected_peak in selected_peaks:
                                    selected_adduct, peak_index = selected_peak.split("||")
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
        fig.update_xaxes(range=[min_rt, max_rt], row=2, col=1)
        fig.update_xaxes(range=[min_rt, max_rt], row=2, col=2)
        

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
            unique_id in selected_good_adducts and \
            unique_id in top_adducts:

            # Create a dictionary to hold data for each peak
            peak_data = {}
            top_adduct = top_adducts[unique_id]
            for idx, peak in rt_peaks_subset.iterrows():
                if peak['adduct'] in top_adduct:
                    peak['top'] = "Yes"
                else:
                    peak['top'] = "No"
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
                    #"Inchi": peak['inchi'],
                }

            # Prepare table data with properties as rows and peaks as columns
            properties = ["RT Peak", "Peak Idx.", "Pol.", "MZ Obs.", "MZ Theo.", "MI Mass", "Chrom.", "Formula", "Best"]
            table_data = {
                " ": properties,
                **{peak_id: [peak_data[peak_id][prop] for prop in properties] for peak_id in peak_data.keys()}
            }

            fig.add_trace(
                go.Table(
                    columnwidth=[2] + [2.75] * len(peak_data.keys()),  # Adjust column widths (1 for Property (" "), 3 for others)
                    header=dict(
                        values=[" "] + [key.split(" (peak")[0] for key in peak_data.keys()],
                        fill_color='white',
                        align='left',
                        font=dict(size=13, color='black', weight='bold'),  # Make header text bold
                        line=dict(color='black', width=1)  # Add black outlines to header
                    ),
                    cells=dict(
                        values=[table_data[" "]] + [table_data[peak_id] for peak_id in peak_data.keys()],
                        fill_color='white',
                        align='left',
                        height=25,
                        font=dict(
                            size=[12] + [16] * len(peak_data.keys()),  # Larger font size for the first column
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
        smiles = rt_peaks_subset['smiles'].unique()[0] if not rt_peaks_subset.empty and len(rt_peaks_subset['smiles'].unique()) == 1 else ""
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
        # Export the figure if export_dir is provided
        if export_dir:
            export_dirname = os.path.join(export_dir, "selection_summary_plots")
            os.makedirs(export_dirname, exist_ok=True)
            fname=f"{export_dirname}/{group_id}_summary_plot.png"
            fig.write_image(
                fname,
                engine="kaleido",
                width=1500,
                height=1000,
                scale=2
            )

def create_interactive_plots(
    processed_data: List[Dict[str, Any]],
    runnum_to_structure_image_grid: Dict[int, str],
    selected_good_adducts: Dict[str, List[str]],
    ambiguous_adducts: Dict[str, str],
    top_adducts: Dict[str, List[str]]
) -> None:
    """
    Create interactive plots for visualizing and annotating processed data.

    Args:
        processed_data (List[Dict[str, Any]]): A list of dictionaries containing processed data for each compound group.
            Each dictionary should include keys such as 'eics', 'rt_peaks', 'top_spectra', 'adduct_color', 'group_id',
            'unique_id', 'group_file', and 'compound_name'.
        runnum_to_structure_image_grid (Dict[int, str]): A dictionary mapping run numbers to base64-encoded molecular structure images.
        selected_good_adducts (Dict[str, List[str]]): A dictionary mapping unique compound IDs to lists of selected
            adduct-peak combinations in the format "adduct||peak_index".
        ambiguous_adducts (Dict[str, str]): A dictionary mapping unique compound IDs to ambiguous annotations.
        top_adducts (Dict[str, List[str]]): A dictionary mapping unique compound IDs to lists of top adducts.

    Returns:
        None
    """

    # Widget Creation
    image_toggle = widgets.ToggleButton(
        value=False,  # Default to hidden
        description='Show Structures',
        tooltip='Toggle to show/hide the compound structure image',
        layout=widgets.Layout(width='150px', margin='5px 0 0 0')
    )
    yaxis_toggle = widgets.ToggleButton(
        value=False,  # Default to unique y-axis
        description='Shared Y-Axis',  # Description when toggled to shared y-axis
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
        placeholder='Index...',
        description='Go to:',
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
        )
    )
    compound_image_label = widgets.Label(
        value="Compound Structures",
        layout=widgets.Layout(margin='0 0 10px 0')  # Add some margin below the label
    )
    completion_label = widgets.Label(value="", layout=widgets.Layout(margin="0 0 0 0"))
    output_container = widgets.Output()

    # Event Handlers
    def on_image_toggle_change(change):
        if image_toggle.value:
            image_toggle.description = 'Hide Structures'
            compound_image_widget.layout.display = 'block'
        else:
            image_toggle.description = 'Show Structures'
            compound_image_widget.layout.display = 'none'

    def update_progress_text():
        progress_label.value = f"{current_index + 1}/{len(processed_data)} Groups Completed"

    def on_toggle_change(change):
        yaxis_toggle.description = 'Shared Y-Axis' if not yaxis_toggle.value else 'Unique Y-Axis'
        update_plot(current_index)

    def navigate_to_group(b):
        nonlocal current_index
        try:
            target_index = int(navigate_textbox.value) - 1
            if 0 <= target_index < len(processed_data):
                current_index = target_index
                update_plot(current_index)
            else:
                output_container.clear_output(wait=True)
                with output_container:
                    print(f"Invalid index. Please enter a number between 1 and {len(processed_data)}.")
        except ValueError:
            output_container.clear_output(wait=True)
            with output_container:
                print("Invalid input. Please enter a valid integer.")

    def next_group(b):
        nonlocal current_index
        if current_index < len(processed_data) - 1:            
            current_index += 1
            update_plot(current_index)
            image_toggle.value = False
            image_toggle.description = 'Show Structures'
        else:
            completion_label.value = "Analysis completed!"

    def previous_group(b):
        nonlocal current_index
        if current_index > 0:
            current_index -= 1
            update_plot(current_index)
            image_toggle.value = False
            image_toggle.description = 'Show Structures'
            completion_label.value = ""
        else:
            output_container.clear_output(wait=True)
            with output_container:
                print("Already at the first group.")

    # Attach Event Handlers
    image_toggle.observe(on_image_toggle_change, names='value')
    yaxis_toggle.observe(on_toggle_change, names='value')
    next_button.on_click(next_group)
    previous_button.on_click(previous_group)
    navigate_button.on_click(navigate_to_group)

    # Layout Definitions
    def create_layout(all_adducts_checkboxes, top_adducts_checkboxes):
        checkbox_layout = widgets.VBox(
            children=[
                widgets.Label(value="Select all good adducts:"),
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
        top_adducts_checkboxes_layout = widgets.VBox(
            children=[
                widgets.Label(value="Select best adduct(s) - default intensity:"),
                top_adducts_checkboxes
            ],
            layout=widgets.Layout(
                border='1px solid black',
                padding='10px',
                margin='10px',
                width='300px',
                align_items='flex-start'
            )
        )
        go_to_label = widgets.Label(value="Go To:")
        go_to_layout = widgets.HBox(
            [go_to_label, navigate_textbox, navigate_button],
            layout=widgets.Layout(
                justify_content='flex-start',  # Align to the far left
                spacing='5px',
                margin='30px 0 0 0'  # Add space above the widget
            )
        )
        navigate_textbox.description = ""
        navigate_textbox.layout = widgets.Layout(width='150px')  # Decrease the size of the search box

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
        top_layout = widgets.HBox(
            [checkbox_layout, top_adducts_checkboxes_layout, button_layout, compound_image_container],
            layout=widgets.Layout(
                align_items='flex-start',
                justify_content='flex-start',
                spacing='10px'
            )
        )
        final_layout = widgets.VBox(
            [completion_label, top_layout],
            layout=widgets.Layout(
                align_items='flex-start',  # Center-align the content
                padding='0px',  # Add padding around the container
            )
        )
        return final_layout

    # Plot Update Logic
    def update_plot(index):
        nonlocal current_index
        data = processed_data[index]

        eics = data['eics']
        top_spectra = data['top_spectra']
        rt_peaks = data['rt_peaks']
        adduct_color = data['adduct_color']
        group_id = data['group_id']
        unique_id = data['unique_id']
        group_run_number = data['group_run_number']

        # Extract adduct-peak combinations from rt_peaks and top_spectra
        adduct_peak_combinations = []
        if isinstance(rt_peaks, pd.DataFrame) and not rt_peaks.empty:
            # Create a mapping from adducts to peak indices and intensities
            adduct_to_peaks = {}
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
        summary_xmin_list = []
        summary_xmax_list = []
        for eic in group_run_eics:
            # Loop through each row in the eic DataFrame
            for _, eic_row in eic.iterrows():
                # Filter data where intensity is above 1e5
                valid_indices = eic_row['i'] > 1e5
                filtered_rt = eic_row['rt'][valid_indices]
                filtered_i = eic_row['i'][valid_indices]

                if len(filtered_rt) > 0:  # Ensure there are valid points
                    # Sort retention times
                    rt_sort = np.argsort(filtered_rt)
                    adduct = get_adduct(eic_row['label'])  # Extract adduct from the label
                    color = adduct_color.get(adduct, 'gray')  # Default to gray if adduct color is missing
                    label = eic_row['label']

                    # Update x_min and x_max based on filtered data
                    summary_xmin_list.append(filtered_rt.min())
                    summary_xmax_list.append(filtered_rt.max())

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
        x_min = min(summary_xmin_list) if summary_xmin_list else None
        x_max = max(summary_xmax_list) if summary_xmax_list else None

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
            "EIC Summary",
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
        fig.update_xaxes(range=[x_min, x_max], row=1, col=3)  # Set x-axis bounds for the summary graph
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
                                row=row + 1,
                                col=col
                            )

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

        # Create checkboxes for selecting good adducts with peak indices
        all_adducts_checkboxes = [
            widgets.Checkbox(
                value=False,
                description=combo['description'],
                disabled=False,
                layout=widgets.Layout(width='300px', margin='0 0 0 -75px')
            )
            for combo in adduct_peak_combinations
        ]
        ambiguous_checkbox = widgets.Checkbox(
            value=False,
            description="Ambiguous",
            disabled=False,
             layout=widgets.Layout(width='300px', margin='0 0 0 -75px')
        )
        all_adducts_checkboxes.append(ambiguous_checkbox)
        checkbox_dict = {checkbox.description: checkbox for checkbox in all_adducts_checkboxes}

        # Initialize top_adducts with the default value if not already set
        if unique_id not in top_adducts and adduct_peak_combinations:
            top_adducts[unique_id] = [adduct_peak_combinations[0]['adduct']]

        # Create top_adducts_checkboxes for adduct selection
        top_adducts_checkboxes = widgets.VBox(
            children=[
                widgets.Checkbox(
                    value=(adduct in top_adducts.get(unique_id, [])),
                    description=adduct,
                    layout=widgets.Layout(width='275px', margin='0 0 0 -75px')
                )
                for adduct in dict.fromkeys(combo['adduct'] for combo in adduct_peak_combinations)
            ]
        )
        
        def on_all_adducts_checkboxes_change(change):
            if ambiguous_checkbox.value:
                ambiguous_adducts[unique_id] = unique_id
                selected_good_adducts.pop(unique_id, None)
                top_adducts.pop(unique_id, None)
            else:
                selected_good_adducts[unique_id] = [
                    combo['adduct'] + "||" + str(combo['peak_index'])
                    for combo in adduct_peak_combinations 
                    if checkbox_dict[combo['description']].value
                ]
                ambiguous_adducts.pop(unique_id, None)

        def on_top_adducts_checkboxes_change(change):
            selected_adducts = [
                checkbox.description
                for checkbox in top_adducts_checkboxes.children
                if checkbox.value
            ]
            top_adducts[unique_id] = selected_adducts

        for checkbox in all_adducts_checkboxes:
            checkbox.observe(on_all_adducts_checkboxes_change, names='value')

        for checkbox in top_adducts_checkboxes.children:
            checkbox.observe(on_top_adducts_checkboxes_change, names='value')

        # Set previously selected values
        if unique_id in selected_good_adducts:
            for selected_combo in selected_good_adducts[unique_id]:
                adduct, peak_index = selected_combo.split("||")
                for combo in adduct_peak_combinations:
                    if combo['adduct'] == adduct and str(combo['peak_index']) == peak_index:
                        combo_description = combo['description']
                        if combo_description in checkbox_dict:
                            checkbox_dict[combo_description].value = True
        elif unique_id in ambiguous_adducts:
            ambiguous_checkbox.value = True

        # Set previously selected value for top_adducts_checkboxes
        if unique_id in top_adducts:
            for checkbox in top_adducts_checkboxes.children:
                if checkbox.description in top_adducts[unique_id]:
                    checkbox.value = True
    
        # Create the layout with all_adducts_checkboxes and top_adducts_checkboxes
        top_layout = create_layout(all_adducts_checkboxes, top_adducts_checkboxes)
        plot_and_checkboxes = widgets.VBox(
            [top_layout, widgets.Output()],
            layout=widgets.Layout(align_items='flex-start')
        )

        clear_output(wait=True)
        update_progress_text()
        with plot_and_checkboxes.children[1]:
            display(fig)
        display(plot_and_checkboxes)
        display(output_container)

    # Initialize
    current_index = 0
    update_plot(current_index)