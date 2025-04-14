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
from io import BytesIO

from PIL import Image, ImageDraw, ImageFont

from metatlas.datastructures.groups import group_name
from metatlas.io import feature_tools as ft
from metatlas.tools.cheminfo import inchi_or_smiles_to_molecule, get_precursor_mz
from metatlas.datastructures import metatlas_objects as metob
from metatlas.tools import cheminfo
from metatlas.tools import extract_msms as exms

from matchms.filtering.filter_utils.load_known_adducts import load_known_adducts

from scipy.signal import find_peaks

import matplotlib.pyplot as plt
from matplotlib import cm

from rdkit.Chem import AllChem, Draw, MolFromSmiles
from rdkit.Chem import rdMolDescriptors as Descriptors
from rdkit.Chem.Descriptors import ExactMolWt


################################
# Get all matching LCMS runs ###
################################

def is_matching_group(lcmsrun_path, group):
    lcmsrun_basename = os.path.basename(lcmsrun_path)
    lcmsrun_group = lcmsrun_basename.split('_')[12]
    
    return lcmsrun_group == group


def get_file_polarity(lcmsrun_path):
    lcmsrun_basename = os.path.basename(lcmsrun_path)
    lcmsrun_polarity = lcmsrun_basename.split('_')[9]
    
    return lcmsrun_polarity


def is_matching_polarity(lcmsrun_path, included_polarities):
    lcmsrun_polarity = get_file_polarity(lcmsrun_path)
    
    return lcmsrun_polarity in included_polarities


def get_matching_lcmsruns(row, include_polarities, include_chromatographies, raw_data_dir):
    
    all_standard_files = []
    for chrom in include_chromatographies:
        all_chrom_files = glob.glob(os.path.join(raw_data_dir, row[f'{chrom.lower()}_experiment'], '*.h5'))
        standard_files = [path for path in all_chrom_files if is_matching_group(path, row[f'{chrom.lower()}_group'])]
        standard_files = [path for path in standard_files if is_matching_polarity(path, include_polarities)]
        
        all_standard_files += standard_files
                                                                              
    return all_standard_files


def build_standard_lcmsrun_table(csv_standard_info_path, include_polarities=['POS', 'NEG'], include_chromatographies=['C18', 'HILIC'], raw_data_dir='/global/cfs/cdirs/metatlas/raw_data/*/'):
    standard_info = pd.read_csv(csv_standard_info_path)
    standard_info['standard_lcmsruns'] = standard_info.apply(lambda row: get_matching_lcmsruns(row, include_polarities, include_chromatographies, raw_data_dir), axis=1)
    standard_lcmsruns_table = standard_info.explode('standard_lcmsruns').reset_index(drop=True).rename(columns={'standard_lcmsruns': 'standard_lcmsrun'})
    
    return standard_lcmsruns_table


##########################
# Calculate Adduct Table #
##########################

def inchi_or_smiles_to_mass(molecule_id):
    mol = inchi_or_smiles_to_molecule(molecule_id)
    
    return ExactMolWt(mol)

def load_adducts_dict():
    return load_known_adducts().set_index('adduct').to_dict(orient='index')


def load_filtered_adducts(include_adducts):
    all_adducts = load_adducts_dict()
    filtered_adducts = {adduct: all_adducts[adduct] for adduct in include_adducts}
    
    return filtered_adducts
    

def separate_adducts_dict(all_adducts):
    pos_adducts = {}
    neg_adducts = {}
    for adduct, data in all_adducts.items():
        if data['ionmode'] == 'positive':
            pos_adducts[adduct] = data
        elif data['ionmode'] == 'negative':
            neg_adducts[adduct] = data
            
    return pos_adducts, neg_adducts


def calc_all_adducts(exact_mass, polarity, include_adducts):
    assert polarity == 'POS' or polarity == 'NEG'
    
    all_adducts = load_filtered_adducts(include_adducts)
    pos_adducts, neg_adducts = separate_adducts_dict(all_adducts)   
            
    adduct_pmzs = []
    if polarity == 'POS':
        for adduct in pos_adducts.keys():
            adduct_pmzs.append((adduct, get_precursor_mz(exact_mass, adduct)))
    elif polarity == 'NEG':
        for adduct in neg_adducts.keys():
            adduct_pmzs.append((adduct, get_precursor_mz(exact_mass, adduct)))
    elif polarity == 'FPS':
        for adduct in all_adducts.keys():
            adduct_pmzs.append((adduct, get_precursor_mz(exact_mass, adduct)))
            
    return adduct_pmzs

def build_adduct_annotated_table(standard_lcmsruns_table, include_adducts=['[M+H]+', '[M+Na]+', '[M-H2O+H]+', '[M+K]+', '[M+NH4]+', '[M]+', '[M+2H]2+', '[M-H]-', '[M+Cl]-', '[M-H2O-H]-', '[M]-', '[M-2H]2-']):
    standard_lcmsruns_table['polarity'] = standard_lcmsruns_table.apply(lambda row: get_file_polarity(row.standard_lcmsrun), axis=1)
    standard_lcmsruns_table['exact_mass'] = standard_lcmsruns_table.apply(lambda row: inchi_or_smiles_to_mass(row.smiles), axis=1)
    standard_lcmsruns_table['all_adducts'] = standard_lcmsruns_table[['exact_mass', 'polarity']].apply(lambda row: calc_all_adducts(row.exact_mass, row.polarity, include_adducts), axis=1)
    standard_lcmsruns_table = standard_lcmsruns_table.explode('all_adducts').reset_index(drop=True).rename(columns={'all_adducts': 'adduct_data'})
    standard_lcmsruns_table[['adduct', 'precursor_mz']] = pd.DataFrame(standard_lcmsruns_table['adduct_data'].tolist(), index=standard_lcmsruns_table.index)

    return standard_lcmsruns_table

#######################
# Tools for chemistry #
#######################

def inchi_to_inchikey(inchi):
    return AllChem.InchiToInchiKey(inchi)

def neutralize_inchi(inchi):
    mol = AllChem.MolFromInchi(inchi)
    neutral_mol = cheminfo.normalize_molecule(mol)
    neutralized_inchi = AllChem.MolToInchi(neutral_mol)
    
    return neutralized_inchi


def charge_from_inchi(inchi):
    mol = AllChem.MolFromInchi(inchi)
    charge = AllChem.GetFormalCharge(mol)
    
    return charge


def formula_from_inchi(inchi):
    mol = AllChem.MolFromInchi(inchi)
    formula = Descriptors.CalcMolFormula(mol)
    
    return formula


def monoisotopic_mass_from_inchi(inchi):
    mol = AllChem.MolFromInchi(inchi)
    monoisotopic_mass = Descriptors.CalcExactMolWt(mol)
    
    return monoisotopic_mass


#########################
# Extract data and plot #
#########################

def extract_adducts(eics):
    eic_adducts = []
    for eic in eics.values():
        eic['adduct'] = eic.label.apply(lambda x: x.split('_')[-1])
        eic_adducts += eic.adduct.tolist()
    return set(eic_adducts)

def generate_adduct_colors(include_adducts):
    adduct_color = {}
    colors = cm.rainbow(np.linspace(0, 1, len(include_adducts)))  # Generate evenly spaced colors
    for i, adduct in enumerate(include_adducts):
        # Convert RGBA array to hex string
        rgba = colors[i]
        hex_color = mcolors.to_hex(rgba)
        adduct_color[adduct] = hex_color
    return adduct_color

def get_lcmsrun_params(lcmsrun_path):
    lcmsrun_basename = os.path.basename(lcmsrun_path)
    lcmsrun_polarity = lcmsrun_basename.split('_')[-2]
    
    return lcmsrun_polarity


def get_chromatography(lcmsrun_path):
    lcmsrun_basename = os.path.basename(lcmsrun_path)
    lcmsrun_chrom = lcmsrun_basename.split('_')[7]
    
    return lcmsrun_chrom

def get_run_num(lcmsrun_path):
    lcmsrun_basename = os.path.basename(lcmsrun_path)
    lcmsrun_run = lcmsrun_basename.split('_')[-1][:-3]
    
    run_num = int(re.sub(r'\D', '', lcmsrun_run))
    
    return run_num


def get_closest_injbl(lcmsrun_path, injbl_pattern='-InjBL-'):
    """Retrieve the closest injection blank before the standard was injected"""
    raw_data_dir = os.path.dirname(lcmsrun_path)
    
    lcmsrun_num = get_run_num(lcmsrun_path)
    
    all_injbl_files = glob.glob(os.path.join(raw_data_dir, '*{}*.h5'.format(injbl_pattern)))
    injbl_run_nums = {get_run_num(injbl_path): injbl_path for injbl_path in all_injbl_files}
    
    closest_run_num = max((run_num for run_num in injbl_run_nums if run_num <= lcmsrun_num), default=None)
        
    return injbl_run_nums[closest_run_num]


def get_rt_range(lcmsrun_path, polarity):
    lcmsrun_data = ft.df_container_from_metatlas_file(lcmsrun_path, desired_key='ms1_{}'.format(polarity.lower()))
    
    rt_min = 0
    rt_max = round(lcmsrun_data.rt.max(), 2)
    
    return rt_min, rt_max


def create_compound_atlas(group, ppm_tolerance, rt_min, rt_max, polarity):
    atlas = group.rename(columns={'precursor_mz': 'mz'})
    atlas['label'] = group.apply(lambda row: '{}_{}'.format(row.compound_name, row.adduct), axis=1)
    atlas_cols = ['label', 'mz', 'rt_min', 'rt_max', 'rt_peak', 'smiles', 'adduct']
    
    atlas['rt_min'] = rt_min
    atlas['rt_max'] = rt_max
    atlas['rt_peak'] = rt_max / 2 # dummy value
    
    return atlas[atlas_cols]
    

def collect_eics_and_ms2(group, ppm_tolerance):
    first_row = group.iloc[0]
    
    lcmsrun_polarity = first_row.polarity
    lcmsrun_params = get_lcmsrun_params(first_row.standard_lcmsrun)
    lcmsrun_chrom = get_chromatography(first_row.standard_lcmsrun)
    lcmsrun_metadata = "{}_{}_{}".format(lcmsrun_chrom, lcmsrun_polarity, lcmsrun_params)
    
    rt_min, rt_max = get_rt_range(first_row.standard_lcmsrun, lcmsrun_polarity)
    atlas = create_compound_atlas(group, ppm_tolerance, rt_min, rt_max, lcmsrun_polarity)
    
    files = group.standard_lcmsrun.tolist() + [first_row.closest_injbl]
    experiment_input = ft.setup_file_slicing_parameters(atlas, files, base_dir=os.getcwd(), ppm_tolerance=ppm_tolerance, extra_time=0.0, polarity=lcmsrun_polarity.lower())

    eics = {}
    ms2_data = {}
    for file_input in experiment_input:
        data = ft.get_data(file_input, save_file=False, return_data=True, ms1_feature_filter=False)
        adduct_eics = ft.group_duplicates(data['ms1_data'],'label', make_string=False)
        
        ms2_summary = ft.calculate_ms2_summary(data['ms2_data'])

        if not data['ms1_data'].empty:
            eics[file_input['lcmsrun']] = adduct_eics
        if not ms2_summary.empty:
            ms2_data[file_input['lcmsrun']] = ms2_summary
    
    return eics, ms2_data


def get_all_ms1_and_ms2(eics, ms2_data, standard_lcmsrun, theoretical_mzs, compound_smiles, adduct_to_polarity, prominence_percentage=0.25, \
                        width_threshold=None, distance_threshold=None):
    """Get MS1 RT peaks for an EIC and the top MS2 spectra per adduct."""
    
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

    rt_peaks.append(pd.DataFrame(predicted_peaks))

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
                        
                        # Use .loc to safely modify the dataframe copy
                        closest_spectrum.loc[:, 'total_intensity_fraction'] = closest_spectrum['total_intensity'] / standard_ms2_data['total_intensity'].max()
                        closest_spectrum.loc[:, 'peak_index'] = rt_peak['peak_index']
                        closest_spectrum.loc[:, 'chromatography'] = rt_peak['chromatography']
                        closest_spectrum.loc[:, 'compound_name'] = rt_peak['compound_name']
                        closest_spectrum.loc[:, 'polarity'] = rt_peak['polarity']
                        
                        top_spectra.append(closest_spectrum)
                    except: # There is no "closest" spectrum
                        continue

    rt_matched_spectra = pd.concat(top_spectra).reset_index(drop=True)

    return rt_peaks, rt_matched_spectra

    
def get_top_ms1_and_ms2(eics, ms2_data, standard_lcmsrun, theoretical_mzs, compound_smiles, adduct_to_polarity):
    """Get MS1 RT peaks for an EIC and the top MS2 spectra per adduct."""
        
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


def get_adduct(compound_name):
    return compound_name.split('_')[1]


def get_compound_name(compound_name):
    return compound_name.split('_')[0]

def smiles_from_inchi_key(inchi_key):
    """
    # Example usage
    inchikey = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
    name = "glucose"

    print(smiles_from_inchi_key(inchikey))
    """
    compound = pcp.get_compounds(inchi_key, 'inchikey')
    return compound[0].canonical_smiles if compound else None


def display_smiles(smiles):
    mol = AllChem.MolFromSmiles(smiles)
    if mol:
        return Draw.MolToImage(mol)
    else:
        print("Invalid SMILES string.")
        return None

def extract_data(lcmsruns_table, method="intensity", ppm_tolerance=5):
    lcmsruns_table['closest_injbl'] = lcmsruns_table['standard_lcmsrun'].apply(get_closest_injbl)
    grouped_lcmsruns_table = lcmsruns_table.groupby(['compound_name', 'standard_lcmsrun'])
    
    eics_list = []
    top_spectra_list = []
    group_name_list = []
    rt_peak_list = []

    for group_name, group in tqdm(grouped_lcmsruns_table, unit=' Compound LCMSRun Group'):
        theoretical_mzs = dict(group['adduct_data'].tolist())
        compound_smiles = group['smiles'].iloc[0]
        adduct_to_polarity = dict(zip(group['adduct'].tolist(), group['polarity'].tolist()))
        group_name = (group_name[0], group_name[1], compound_smiles)

        eics, ms2_data = collect_eics_and_ms2(group, ppm_tolerance)

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

    return eics_list, top_spectra_list, group_name_list, rt_peak_list


def extract_and_plot(standard_lcmsruns_table, ppm_tolerance=5, plot_output_dir='./annotation_plots', data_output_dir='./annotation_data', generate_plots=True):
    standard_lcmsruns_table['closest_injbl'] = standard_lcmsruns_table['standard_lcmsrun'].apply(get_closest_injbl)
    grouped_lcmsruns_table = standard_lcmsruns_table.groupby(['compound_name', 'standard_lcmsrun'])
    
    if not os.path.isdir(plot_output_dir):
        os.mkdir(plot_output_dir)
        
    if not os.path.isdir(data_output_dir):
        os.mkdir(data_output_dir)
    
    all_top_spectra = []
    all_rt_peaks = []

    for group_name, group in tqdm(grouped_lcmsruns_table, unit=' Compound LCMSRun Group'):
        theoretical_mzs = dict(group['adduct_data'].tolist())
        compound_smiles = group['smiles'].iloc[0]
        adduct_to_polarity = dict(zip(group['adduct'].tolist(), group['polarity'].tolist()))
        
        eics, ms2_data = collect_eics_and_ms2(group, ppm_tolerance)
        rt_peaks, top_spectra = get_top_ms1_and_ms2(eics, ms2_data, group_name[1], theoretical_mzs, compound_smiles, adduct_to_polarity)
        
        if generate_plots:
            save_annotation_fig(eics, top_spectra, group_name, plot_output_dir)
        
        all_top_spectra.append(top_spectra)
        all_rt_peaks += rt_peaks
    
    if all_rt_peaks:
        print("Writing RT peak data")
        all_rt_peaks = pd.DataFrame(all_rt_peaks)
        all_rt_peaks.to_csv(os.path.join(data_output_dir, 'rt_peak_annotations.csv'), index=False)
    if all_top_spectra:
        print("Writing top spectra data")
        all_top_spectra = pd.concat(all_top_spectra).reset_index(drop=True)
        all_top_spectra.to_json(os.path.join(data_output_dir, 'top_intensity_spectra.json'))

    return all_rt_peaks, all_top_spectra

def save_annotation_fig(eics, top_spectra, group_name, output_dir):
    group_chrom = get_chromatography(group_name[1])
    group_pol = get_file_polarity(group_name[1])
    group_params = get_lcmsrun_params(group_name[1])
    
    fig_title = f"{group_name[0]} {group_chrom} {group_pol} {group_params}".replace('/', 'or')
    
    eic_adducts = []
    for eic in eics.values():
        eic['adduct'] = eic.label.apply(lambda x: x.split('_')[-1])
        eic_adducts += eic.adduct.tolist()

    eic_adducts = set(eic_adducts)

    adduct_color = dict()
    color = iter(cm.rainbow(np.linspace(0, 1, len(eic_adducts))))
    for adduct in eic_adducts:
        adduct_color[adduct] = next(color)
        
    # Pass top_spectra to plot_adduct_eics to add the "X" marker
    plot_adduct_eics(eics, adduct_color, fig_title, top_spectra)
    plt.savefig(os.path.join(output_dir, f"{fig_title}_eics.pdf"))
    plt.close()
        
    plot_top_spectra(top_spectra, adduct_color, fig_title)
    plt.savefig(os.path.join(output_dir, f"{fig_title}_top_spectra.pdf"))
    plt.close()


def plot_adduct_eics(eics, adduct_color, fig_title, top_spectra):
    fig, axs = plt.subplots(2, 2, figsize=(20, 9), sharex=True, sharey=False)

    for idx, (lcmsrun_path, eic) in enumerate(eics.items()):
        total_eics = eic.shape[0]
        for eic_idx, eic_row in eic.iterrows():

            rt_sort = np.argsort(eic_row['rt'])
            adduct = get_adduct(eic_row['label'])

            # Plot on the left column for the first lcmsrun, right column for the second
            ax_raw = axs[0, idx]  # Top row for raw values
            ax_log = axs[1, idx]  # Bottom row for log-scale values

            # Plot raw values in the top row
            line_raw, = ax_raw.plot(eic_row['rt'][rt_sort], eic_row['i'][rt_sort], alpha=0.8, label=adduct, color=adduct_color[adduct])
            line_log, = ax_log.plot(eic_row['rt'][rt_sort], np.log10(eic_row['i'][rt_sort].astype(float)), alpha=0.8, label=adduct, color=adduct_color[adduct])

            peak = np.argmax(eic_row['i'][rt_sort])
            ax_raw.scatter(eic_row['rt'][rt_sort][peak], eic_row['i'][rt_sort][peak], color=adduct_color[adduct])
            ax_log.scatter(eic_row['rt'][rt_sort][peak], np.log10(eic_row['i'][rt_sort].astype(float))[peak], color=adduct_color[adduct])

            # Add "X" marker for the retention time in top_spectra
            if not top_spectra.empty:
                for _, spectrum_row in top_spectra.iterrows():
                    if spectrum_row['adduct'] == adduct:
                        rounded_rt = round(spectrum_row['rt'], 2)  # Round RT to 2 decimal places

                        # Get the y-axis intensity value at the given RT
                        raw_intensity = eic_row['i'][rt_sort][np.searchsorted(eic_row['rt'][rt_sort], spectrum_row['rt'])]
                        log_intensity = np.log10(raw_intensity)

                        # Plot the "X" marker at the intensity value
                        x_marker = ax_raw.scatter(spectrum_row['rt'], raw_intensity, color=adduct_color[adduct], marker='x', s=100, label=f"Spectra RT: {rounded_rt}")
                        ax_log.scatter(spectrum_row['rt'], log_intensity, color=adduct_color[adduct], marker='x', s=100)

            lcms_params = get_lcmsrun_params(lcmsrun_path)
            run_num = get_run_num(lcmsrun_path)

            # Set subplot titles to the filenames
            ax_raw.set_title(f"{lcms_params}_{run_num}")

    # Set shared X axis labels for both rows
    for ax in axs[1, :]:  # Bottom row
        ax.set_xlabel('Retention Time (RT)')
        ax.grid()

    for ax in axs[0, :]:
        ax.grid()

    # Set Y axis labels for both rows
    axs[0, 0].set_ylabel('Intensity (i)')  # Left, top row
    axs[1, 0].set_ylabel('Log10 Intensity')  # Left, bottom row

    axs[1, 0].sharey(axs[1, 1])  
    axs[0, 1].sharey(axs[0, 0])  

    # Collect legend handles and labels
    handles_left, labels_left = axs[0, 0].get_legend_handles_labels()
    handles_right, labels_right = axs[0, 1].get_legend_handles_labels()

    # Ensure there are handles and labels before creating legends
    if handles_left and labels_left:
        leg1 = axs[0, 0].legend(handles=handles_left, labels=labels_left, loc='upper left', bbox_to_anchor=(-0.4, 1), fontsize=13, markerscale=0.2, handletextpad=0.5)
        for legobj in leg1.legend_handles:
            legobj.set_linewidth(8.0)

    if handles_right and labels_right:
        leg2 = axs[0, 1].legend(handles=handles_right, labels=labels_right, loc='upper right', bbox_to_anchor=(1.4, 1), fontsize=13, markerscale=0.2, handletextpad=0.5)
        for legobj in leg2.legend_handles:
            legobj.set_linewidth(8.0)

    fig.suptitle(fig_title)
    plt.tight_layout()


def plot_top_spectra(top_spectra, adduct_color, fig_title):
    if top_spectra.empty:
        plt.figure()  
        return
    
    num_columns = 3
    num_rows = math.ceil(len(top_spectra) / num_columns)
    
    fig, axes = plt.subplots(nrows=num_rows, ncols=num_columns, figsize=(19.845, 4.5 * num_rows))
    axes = axes.flatten()

    for i, row in top_spectra.iterrows():
        mz_values = row['spectrum'][0]
        i_values = row['spectrum'][1]

        markerline, stemlines, baseline = axes[i].stem(mz_values, i_values, basefmt=" ", markerfmt=" ")
        plt.setp(stemlines, 'color', adduct_color[row.adduct])
        axes[i].set_title(f"{row['adduct']} RT: {round(row['rt'], 2)}")

    plt.tight_layout()


def save_rt_peaks_to_atlas_format(rt_peaks):

    # Find the row with the highest intensity for each group
    idx_max_intensity = rt_peaks.groupby(['chromatography', 'polarity', 'compound_name'])['intensity'].idxmax()
    highest_intensity_row = rt_peaks.loc[idx_max_intensity]

    # Filter rows to keep only those with the same adduct as the highest intensity row
    top_adducts_per_pol = rt_peaks.merge(
        highest_intensity_row[['chromatography', 'polarity', 'compound_name', 'adduct']],
        on=['chromatography', 'polarity', 'compound_name', 'adduct'],
        how='inner'
    )

    # Find other good peaks for the selected adduct
    top_adducts_per_pol_grouped = top_adducts_per_pol.groupby(['compound_name', 'chromatography', 'polarity', 'adduct'])
    for (compound, chromatography, polarity, adduct), group in top_adducts_per_pol_grouped:
        if len(group) > 1:  # Check if there are multiple rows in the group
            # Update the 'compound_name' column to include the peak_index
            top_adducts_per_pol.loc[group.index, 'compound_name'] = group['compound_name'] + ' ' + group['peak_index']

    # enrich for atlas related metadata
    top_adducts_per_pol['rt_min'] = top_adducts_per_pol['rt_peak'] - 0.5
    top_adducts_per_pol['rt_max'] = top_adducts_per_pol['rt_peak'] + 0.5
    top_adducts_per_pol['label'] = top_adducts_per_pol['compound_name']
    top_adducts_per_pol['mz_tolerance'] = 5
    top_adducts_per_pol['inchi'] = top_adducts_per_pol['smiles'].apply(lambda row: AllChem.MolToInchi(AllChem.MolFromSmiles(row)))
    top_adducts_per_pol['inchi_key'] = top_adducts_per_pol['inchi'].apply(inchi_to_inchikey)

    # rename and drop columns to match metatlas atlas convention
    top_adducts_per_pol.rename({'mz_theoretical': 'mz'})
    top_adducts_per_pol['polarity'] = top_adducts_per_pol['polarity'].apply(lambda pol: 'positive' if pol == 'POS' else 'negative')
    top_adducts_per_pol.drop(columns=['peak_index', 'standard_lcmsrun', 'intensity', 'mz_observed', 'ppm_error'], inplace=True)

    pos_annotations = top_adducts_per_pol[top_adducts_per_pol['polarity'] == 'positive']
    neg_annotations = top_adducts_per_pol[top_adducts_per_pol['polarity'] == 'negative']

    c18_pos_annotations = pos_annotations[pos_annotations['chromatography'] == 'C18'].sort_values('rt_peak')
    c18_neg_annotations = neg_annotations[neg_annotations['chromatography'] == 'C18'].sort_values('rt_peak')

    hilic_pos_annotations = pos_annotations[pos_annotations['chromatography'] == 'HILICZ'].sort_values('rt_peak')
    hilic_neg_annotations = neg_annotations[neg_annotations['chromatography'] == 'HILICZ'].sort_values('rt_peak')

    return c18_pos_annotations, c18_neg_annotations, hilic_pos_annotations, hilic_neg_annotations


def find_atlas_matches(atlases, new_entries, cutoff=0.8):
    """
    Compare columns between atlases and new_entries to find matches.
    Matches are checked in the order of inchi, inchi_key, compound_name, and label.
    If a match is found, it is added to a dictionary with the matching value and source_file(s).

    Parameters:
        atlases (pd.DataFrame): DataFrame containing atlas data.
        new_entries (pd.DataFrame): DataFrame containing new entries to compare.
        cutoff (float): Similarity cutoff for fuzzy matching (default: 0.8).
    """
    atlases['compound_name'] = atlases['compound_name'].astype(str)
    atlases['label'] = atlases['label'].astype(str)
    matches_dict = {}
    nonmatches_dict = {}

    for _, new_row in new_entries.iterrows():
        new_label = str(new_row['label'])
        new_compound_name = str(new_row['compound_name'])
        new_inchikey = str(new_row['inchi_key'])
        new_inchi = str(new_row['inchi'])

        # Check for exact matches in inchi
        if new_inchi in atlases['inchi'].values:
            matching_rows = atlases[atlases['inchi'] == new_inchi]
            source_files = matching_rows['source_file'].tolist()
            matches_dict[new_label] = [new_inchi, source_files]
            print(f"Exact match found for {new_label} in inchi: {new_inchi}")
            continue

        # Check for exact matches in inchi_key
        if new_inchikey in atlases['inchi_key'].values:
            matching_rows = atlases[atlases['inchi_key'] == new_inchikey]
            source_files = matching_rows['source_file'].tolist()
            matches_dict[new_label] = [new_inchikey, source_files]
            print(f"Exact match found for {new_label} in inchi_key: {new_inchikey}")
            continue

        # Check for fuzzy matches in compound_name
        compound_name_matches = get_close_matches(new_compound_name, atlases['compound_name'].tolist(), n=1, cutoff=cutoff)
        if compound_name_matches:
            matching_rows = atlases[atlases['compound_name'] == compound_name_matches[0]]
            source_files = matching_rows['source_file'].tolist()
            matches_dict[new_label] = [compound_name_matches[0], source_files]
            print(f"Fuzzy match found for {new_label} in compound_name: {compound_name_matches[0]}")
            continue

        # Check for fuzzy matches in label
        label_matches = get_close_matches(new_label, atlases['label'].tolist(), n=1, cutoff=cutoff)
        if label_matches:
            matching_rows = atlases[atlases['label'] == label_matches[0]]
            source_files = matching_rows['source_file'].tolist()
            matches_dict[new_label] = [label_matches[0], source_files]
            print(f"Fuzzy match found for {new_label} in label: {label_matches[0]}")
            continue

        ## Fingerprint similarity?

        # If no match is found, add to nonmatches_dict
        nonmatches_dict[new_label] = [[new_inchikey, new_inchi, new_label, new_compound_name], ""]
        print(f"No match found for {new_label} in any atlas field.")

    # Convert matches dictionary to DataFrame
    matches_df = pd.DataFrame(
        [(key, value[0], value[1]) for key, value in matches_dict.items()],
        columns=['new_label', 'matching_value', 'atlas_source_files']
    )
    nonmatches_df = pd.DataFrame(
        [(key, value[0], value[1]) for key, value in nonmatches_dict.items()],
        columns=['new_label', 'attempted_matching_values', 'atlas_source_files']
    )

    return matches_df, nonmatches_df


def search_for_matches_in_metatlas_db(all_molecules):
    in_db = {}
    notin_db = {}
    flat_in_db = {}
    for _, molecule in tqdm(all_molecules.iterrows()):
        inchi_key_parts = molecule.inchi_key.split('-')
        inchi_key_parts[1] = '%'
        flat_inchi_key = '-'.join(inchi_key_parts)
        print(f"\nSearching metatlas db for {molecule.label} ({molecule.inchi_key})")
        
        db_entry = metob.retrieve('compound', inchi_key=molecule.inchi_key)
        
        if db_entry == []:
            print(f"{molecule.label} ({molecule.inchi_key}) not found in metalas db. Trying flat inchi key ({flat_inchi_key})")
            flat_entry = metob.retrieve('compound', inchi_key=flat_inchi_key)
            if flat_entry == []:
                print(f"{molecule.label} ({molecule.inchi_key}) not found in metatlas db.")
                notin_db[molecule.label] = molecule.inchi_key
            else:
                print(f"Found {len(flat_entry)} entries in metatlas db with flat inchi search.")
                flat_in_db[molecule.label] = flat_entry[0].inchi_key
        else:
            print(f"Found {len(db_entry)} entries in metatlas db.")
            in_db[molecule.label] = molecule.inchi_key

    return in_db, notin_db, flat_in_db


def search_for_matches_in_msms_refs(all_molecules, msms_refs):
    in_msms_refs = {}
    notin_msms_refs = {}
    flatin_msms_refs = {}
    for _, molecule in tqdm(all_molecules.iterrows()):
        inchi_key_parts = molecule.inchi_key.split('-')
        inchi_key_parts[1] = '%'
        flat_inchi_key = '-'.join(inchi_key_parts)
        print(f"Searching MSMS refs for {molecule.label} ({molecule.inchi_key})")

        matching_refs = msms_refs.loc[msms_refs['inchi_key'] == molecule.inchi_key]
        
        if matching_refs.empty:
            print(f"{molecule.label} ({molecule.inchi_key}) not found in MSMS refs. Trying flat inchi key ({flat_inchi_key})")
            flat_matching_refs = msms_refs.loc[msms_refs['inchi_key'] == flat_inchi_key]
            if flat_matching_refs.empty:
                print(f"{molecule.label} ({molecule.inchi_key}) not found in MSMS refs.")
                notin_msms_refs[molecule.label] = molecule.inchi_key
            else:
                print(f"Found {len(flat_matching_refs)} molecules in MSMS refs with flat inchi search.")
                flatin_msms_refs[molecule.label] = flat_matching_refs.iloc[0].inchi_key
        else:
            print(f"Found {len(matching_refs)} molecules in MSMS refs.")
            in_msms_refs[molecule.label] = molecule.inchi_key

    return in_msms_refs, notin_msms_refs, flatin_msms_refs


def format_for_atlas_store(input_compounds):

    input_compounds.reset_index(drop=True, inplace=True)

    metatlas_compounds = []
    for idx, row in input_compounds.iterrows():
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


def enrich_metadata(refs):
    frag_method = 'HCD'
    instrument_type = 'Orbitrap'
    decimal = 4.0
    id_prefix = 'schellermetasci'
    exms.enrich_metadata(refs, frag_method, instrument_type, decimal, id_prefix)
    return refs

def get_collision_energy(lcmsrun_path):
    lcmsrun = os.path.basename(lcmsrun_path)
    collision_energy = lcmsrun.split('_')[-2].split('-')[1][2:]
    return collision_energy

def make_text_spectrum(spectrum):
    return exms.make_text_spectrum(spectrum)


def select_compounds_from_gui(unfiltered_df, selected_compounds_table):
    """
    Updated version to handle multiple peaks per adduct
    """
    result = pd.DataFrame()
    
    for _, row in selected_compounds_table.iterrows():
        compound_name = row['compound_name']
        standard_lcmsrun = row['standard_lcmsrun']
        selected_adducts = row['selected_adducts']
        selected_peak_indices = row['selected_peak_indices']
        
        # Filter by compound name and standard_lcmsrun
        mask = (unfiltered_df['compound_name'] == compound_name) & (unfiltered_df['standard_lcmsrun'] == standard_lcmsrun)
        
        # Further filter by adduct and peak_index if they exist in the dataframe
        if 'adduct' in unfiltered_df.columns and len(selected_adducts) > 0:
            adduct_mask = unfiltered_df['adduct'].isin(selected_adducts)
            mask = mask & adduct_mask
            
        if 'peak_index' in unfiltered_df.columns and len(selected_peak_indices) > 0:
            peak_mask = unfiltered_df['peak_index'].isin(selected_peak_indices)
            mask = mask & peak_mask
        
        subset = unfiltered_df[mask].copy()
        result = pd.concat([result, subset], ignore_index=True)
    
    return result


# def select_compounds_from_gui(full_table, select_table):
#     subset_table = full_table[
#         full_table.apply(
#             lambda row: any(
#                 (row['standard_lcmsrun'] == compound_row['standard_lcmsrun']) and
#                 (row['compound_name'] == compound_row['compound_name']) and
#                 (row['adduct'] in compound_row['selected_adducts'])
#                 for _, compound_row in select_table.iterrows()
#             ),
#             axis=1
#         )
#     ]
#     return subset_table

def generate_molecular_images(lcmsruns_table):
    """
    Generate molecular structure images for unique SMILES strings in the given DataFrame.

    Args:
        lcmsruns_table (pd.DataFrame): DataFrame containing a 'smiles' column.

    Returns:
        dict: A dictionary where keys are SMILES strings and values are base64-encoded molecular structure images.
    """
    smiles_to_image = {}
    unique_smiles = lcmsruns_table['smiles'].dropna().unique()

    for smiles in unique_smiles:
        try:
            mol = MolFromSmiles(smiles)
            if mol:
                # Generate the image
                img = Draw.MolToImage(mol, size=(200, 200))
                
                # Convert the image to base64
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                img_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
                
                # Store the base64 string in the dictionary
                smiles_to_image[smiles] = img_base64
        except Exception as e:
            print(f"Error generating image for SMILES: {smiles}, Error: {e}")

    return smiles_to_image

def generate_gridded_molecular_images(lcmsruns_table):
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
            buffered = BytesIO()
            grid_img.save(buffered, format="PNG")
            img_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")

            # Store the base64 string in the dictionary
            runnum_to_structure_image_grid[runnum] = img_base64

    return runnum_to_structure_image_grid