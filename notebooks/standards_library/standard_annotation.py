import pandas as pd
import numpy as np
import os
import glob

from metatlas.datastructures.groups import group_name
from metatlas.io import feature_tools as ft
from metatlas.tools.cheminfo import inchi_or_smiles_to_molecule, get_precursor_mz

from matchms.filtering.filter_utils.load_known_adducts import load_known_adducts
from rdkit.Chem.Descriptors import ExactMolWt

import matplotlib.pyplot as plt
from matplotlib import cm

import re
import math

from tqdm.notebook import tqdm


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


#########################
# Extract data and plot #
#########################


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
    
    
def get_top_ms1_and_ms2(eics, ms2_data, standard_lcmsrun, theoretical_mzs, compound_smiles, adduct_to_polarity):
    """Get MS1 RT peaks for an EIC and the top MS2 spectra per adduct."""
    
    if standard_lcmsrun in ms2_data.keys():
        standard_ms2_data = ms2_data[standard_lcmsrun]
        standard_ms2_data['total_intensity'] = standard_ms2_data['spectrum'].apply(lambda x: x[1].sum())
        standard_ms2_data['adduct'] = standard_ms2_data['label'].apply(lambda x: x.split('_')[-1])
        standard_ms2_data['lcmsrun'] = standard_lcmsrun

        top_spectra = standard_ms2_data.sort_values('total_intensity', ascending=False).groupby('adduct').head(1).reset_index(drop=True)
    else:
        top_spectra = pd.DataFrame()
        
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
            
            rt_peaks.append(rt_peak)
        
    return rt_peaks, top_spectra


def get_adduct(compound_name):
    return compound_name.split('_')[1]


def get_compound_name(compound_name):
    return compound_name.split('_')[0]


def plot_adduct_eics(eics, adduct_color, fig_title):
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
            ax_raw.plot(eic_row['rt'][rt_sort], eic_row['i'][rt_sort], alpha=0.8, label=adduct, color=adduct_color[adduct])
            ax_log.plot(eic_row['rt'][rt_sort], np.log10(eic_row['i'][rt_sort].astype(float)), alpha=0.8, label=adduct, color=adduct_color[adduct])

            peak = np.argmax(eic_row['i'][rt_sort])
            ax_raw.scatter(eic_row['rt'][rt_sort][peak], eic_row['i'][rt_sort][peak], color=adduct_color[adduct])
            ax_log.scatter(eic_row['rt'][rt_sort][peak], np.log10(eic_row['i'][rt_sort].astype(float))[peak], color=adduct_color[adduct])

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

    handles_left, labels_left = axs[0, 0].get_legend_handles_labels()
    handles_right, labels_right = axs[0, 1].get_legend_handles_labels()

    # Left legend for standard run
    leg1 = axs[0, 0].legend(handles=handles_left, labels=labels_left, loc='upper left', bbox_to_anchor=(-0.3, 1), fontsize=13)
    # Right legend for blank run
    leg2 = axs[0, 1].legend(handles=handles_right, labels=labels_right, loc='upper right', bbox_to_anchor=(1.3, 1), fontsize=13)

    # set the linewidth of each legend object
    for leg in [leg1, leg2]:
        for legobj in leg.legend_handles:
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
        
    plot_adduct_eics(eics, adduct_color, fig_title)
    plt.savefig(os.path.join(output_dir, f"{fig_title}_eics.pdf"))
    plt.close()
        
    plot_top_spectra(top_spectra, adduct_color, fig_title)
    plt.savefig(os.path.join(output_dir, f"{fig_title}_top_spectra.pdf"))
    plt.close()
    
    
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
    
    all_rt_peaks = pd.DataFrame(all_rt_peaks)
    all_top_spectra = pd.concat(all_top_spectra).reset_index(drop=True)
    
    all_rt_peaks.to_csv(os.path.join(data_output_dir, 'rt_peak_annotations.csv'))
    all_top_spectra.to_json(os.path.join(data_output_dir, 'top_intensity_spectra.json'))
        
    return all_rt_peaks, all_top_spectra