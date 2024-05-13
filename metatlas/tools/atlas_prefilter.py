import glob
import os
import pandas as pd
import numpy as np
import numpy.typing as npt
from typing import TypeAlias, Optional
from tqdm.notebook import tqdm

import matchms as mms
from matchms.similarity import CosineHungarian

from metatlas.io.metatlas_get_data_helper_fun import make_atlas_df
from metatlas.io import feature_tools as ft
import metatlas.plots.dill2plots as dp
from metatlas.datastructures.analysis_identifiers import AnalysisIdentifiers
from metatlas.datastructures.metatlas_dataset import MetatlasDataset
from metatlas.tools.notebook import in_papermill

# Typing:
MS2Spectrum: TypeAlias = npt.NDArray[npt.NDArray[float]]


def get_sample_file_paths(ids: AnalysisIdentifiers) -> list[str]:
    """Return lists of sample hdf5 file paths filtered by polarity."""
    sample_file_paths = [lcmsrun.hdf5_file for lcmsrun in ids.all_lcmsruns if 'QC' not in lcmsrun.hdf5_file]
    return sample_file_paths


def order_ms2_spectrum(spectrum: MS2Spectrum) -> MS2Spectrum:
    """Order spectrum by m/z from lowest to highest.

    Ordering spectrum by m/z prevents MatchMS errors during MS/MS scoring.
    """
    order_idx = np.argsort(spectrum[0])
    ordered_spec = np.array([spectrum[0][order_idx], spectrum[1][order_idx]])
    return ordered_spec


def get_sample_data(aligned_atlas_df: pd.DataFrame, sample_files: list[str], ppm_tolerance: int, extra_time: float, polarity: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collect MS1 and MS2 data from experimental sample data using aligned atlas."""

    experiment_input = ft.setup_file_slicing_parameters(aligned_atlas_df, sample_files, base_dir=os.getcwd(), ppm_tolerance=ppm_tolerance, extra_time=extra_time, polarity=polarity)

    ms1_data = []
    ms2_data = []

    for file_input in tqdm(experiment_input, unit="file", disable=in_papermill()):

        data = ft.get_data(file_input, save_file=False, return_data=True, ms1_feature_filter=False)
        data['ms1_summary']['lcmsrun_observed'] = file_input['lcmsrun']

        ms2_summary = ft.calculate_ms2_summary(data['ms2_data'])
        ms2_summary['lcmsrun_observed'] = file_input['lcmsrun']

        ms1_data.append(data['ms1_summary'])
        ms2_data.append(ms2_summary)

    ms1_data = pd.concat(ms1_data)
    ms2_data = pd.concat(ms2_data)

    df = pd.DataFrame(ms2_data['spectrum'])
    df.to_csv("/out/ms2_data_array.csv")

    ms2_data['spectrum'] = ms2_data['spectrum'].apply(order_ms2_spectrum)

    return (ms1_data, ms2_data)


def calculate_ms2_scores(ms2_data_refs_merge: pd.DataFrame, frag_tolerance: float) -> list[float]:
    """Calculate cosine similarity scores for each feature in the sample MS2 data.

    To maintain parity with the MSMS hits collection, scoring is performed using the Hungarian alignment algorithm.
    """

    cosine_hungarian = CosineHungarian(tolerance=frag_tolerance)
    scores = ms2_data_refs_merge.apply(lambda x: cosine_hungarian.pair(x.mms_spectrum_x, x.mms_spectrum_y), axis=1)

    return scores


def load_and_filter_msms_refs_file(msms_refs_path: str, polarity: str) -> pd.DataFrame:
    """Load and filter MSMS refs file.

    In addition to loading and filtering MSMS refs, spectral data is converted to Numpy array format.
    """

    ref_dtypes = {'database': str, 'id': str, 'name': str,
                  'spectrum': object, 'decimal': int, 'precursor_mz': float,
                  'polarity': str, 'adduct': str, 'fragmentation_method': str,
                  'collision_energy': str, 'instrument': str, 'instrument_type': str,
                  'formula': str, 'exact_mass': float,
                  'inchi_key': str, 'inchi': str, 'smiles': str}

    msms_refs_df = pd.read_csv(msms_refs_path, sep='\t', dtype=ref_dtypes)
    msms_refs_filtered = msms_refs_df[(msms_refs_df['database'] == 'metatlas') & (msms_refs_df['polarity'] == polarity)].copy()

    msms_refs_filtered['spectrum'] = msms_refs_filtered['spectrum'].apply(lambda x: np.asarray(eval(x)))
    msms_refs_filtered['spectrum'] = msms_refs_filtered['spectrum'].apply(order_ms2_spectrum)

    return msms_refs_filtered


def score_ms2_data(ms2_data: pd.DataFrame, aligned_atlas_df: pd.DataFrame,
                   polarity: str, msms_refs_path: str, frag_tolerance: float) -> pd.DataFrame:
    """Score collected MS2 data against reference spectra.

    This function merges the MSMS refs with the MS2 data collected from the samples for scoring
    and calculates scores and number of matching ions.
    """

    msms_refs = load_and_filter_msms_refs_file(msms_refs_path, polarity)

    ms2_data_refs_merge = pd.merge(ms2_data, aligned_atlas_df[['label', 'inchi_key']], on='label')
    ms2_data_refs_merge = pd.merge(ms2_data_refs_merge, msms_refs[['id', 'inchi_key', 'spectrum']], on='inchi_key')

    ms2_data_refs_merge['mms_spectrum_x'] = ms2_data_refs_merge.apply(lambda x: mms.Spectrum(x.spectrum_x[0], x.spectrum_x[1], metadata={'precursor_mz': x.precursor_mz}), axis=1)
    ms2_data_refs_merge['mms_spectrum_y'] = ms2_data_refs_merge.apply(lambda x: mms.Spectrum(x.spectrum_y[0], x.spectrum_y[1], metadata={'precursor_mz': x.precursor_mz}), axis=1)

    ms2_data_refs_merge['mms_out'] = calculate_ms2_scores(ms2_data_refs_merge, frag_tolerance)

    ms2_data_refs_merge['score'] = ms2_data_refs_merge['mms_out'].apply(lambda x: x['score'])
    ms2_data_refs_merge['matches'] = ms2_data_refs_merge['mms_out'].apply(lambda x: x['matches'])

    return ms2_data_refs_merge


def filter_atlas_labels(ms1_data: pd.DataFrame, ms2_data_scored: pd.DataFrame, peak_height: Optional[float], num_points: Optional[int],
                        msms_score: Optional[float], msms_matches: Optional[int]) -> set[str]:
    """Filter atlas labels to include only those that pass the MS1 and MS2 thresholds."""

    ms1_data_filtered = ms1_data[ms1_data['peak_height'] >= peak_height] if peak_height is not None else ms1_data
    ms1_data_filtered = ms1_data_filtered[ms1_data_filtered['num_datapoints'] >= num_points] if num_points is not None else ms1_data_filtered
    ms1_reduced_labels = set(ms1_data_filtered.label.tolist())

    if msms_score is not None or msms_matches is not None:
        ms2_data_filtered = ms2_data_scored[ms2_data_scored['score'] >= msms_score] if msms_score is not None else ms2_data_scored
        ms2_data_filtered = ms2_data_filtered[ms2_data_filtered['matches'] >= msms_matches] if msms_matches is not None else ms2_data_filtered
        ms2_reduced_labels = set(ms2_data_filtered.label.tolist())
    else:
        ms2_reduced_labels = ms1_reduced_labels

    reduced_labels = ms1_reduced_labels.intersection(ms2_reduced_labels)

    return reduced_labels


def filter_atlas(aligned_atlas, ids: AnalysisIdentifiers, analysis, data: MetatlasDataset):

    aligned_atlas_df = make_atlas_df(aligned_atlas)
    sample_file_paths = get_sample_file_paths(ids)

    if analysis.parameters.mz_tolerance_override is not None:
        mz_tolerance = analysis.parameters.mz_tolerance_override
    else:
        mz_tolerance = analysis.parameters.mz_tolerance_default

    ms1_data, ms2_data = get_sample_data(aligned_atlas_df, sample_file_paths, mz_tolerance, data.extra_time, ids.polarity)
    ms2_data_scored = score_ms2_data(ms2_data, aligned_atlas_df, ids.polarity,
                                     analysis.parameters.msms_refs, analysis.parameters.frag_mz_tolerance)

    reduced_labels = filter_atlas_labels(ms1_data, ms2_data_scored, analysis.parameters.peak_height, analysis.parameters.num_points,
                                         analysis.parameters.msms_score, analysis.parameters.msms_matches)

    aligned_filtered_atlas_df = aligned_atlas_df[aligned_atlas_df['label'].isin(reduced_labels)]
    aligned_filtered_atlas = dp.get_atlas(aligned_atlas.name, aligned_filtered_atlas_df, ids.polarity, mz_tolerance)

    return aligned_filtered_atlas


# Prototype Pre-Filter RT-Alignment Functions


def extract_file_polarity(file_path: str) -> str:
    """Extract file polarity from file path.

    Per file naming conventions, field 9 of each (underscore delimited) filename contains the polarity information for that particular file.
    """
    return os.path.basename(file_path).split('_')[9]


def subset_file_paths(raw_data_dir: str, experiment: str, polarity: str) -> tuple[list[str], list[str]]:
    """Return lists of QC file paths and sample file paths filtered by polarity.

    Filter polarity is used in this case to retrieve both fast polarity switching (FPS) files and files matching the defined polarity
    """

    if polarity == 'positive':
        file_polarity = 'POS'
        filter_polarity = 'NEG'
    else:
        file_polarity = 'NEG'
        filter_polarity = 'POS'

    all_files = glob.glob(os.path.join(raw_data_dir, experiment, '*.h5'))

    sample_files = [file for file in all_files if extract_file_polarity(file) == file_polarity and 'QC' not in file]
    qc_files = [file for file in all_files if extract_file_polarity(file) != filter_polarity and 'QC_' in file]

    return (sample_files, qc_files)


def get_rt_alignment_ms1_data(rt_alignment_atlas: pd.DataFrame, qc_files: list[str], ppm_tolerance: int,
                              extra_time: float, polarity: str) -> pd.DataFrame:
    """Collect all MS1 feature data for each entry in the retention time adjustment atlas."""

    experiment_input = ft.setup_file_slicing_parameters(rt_alignment_atlas, qc_files, base_dir=os.getcwd(), ppm_tolerance=ppm_tolerance, extra_time=extra_time, polarity=polarity)

    ms1_data = []
    for file_input in experiment_input:

        data = ft.get_data(file_input, save_file=False, return_data=True)
        data['ms1_summary']['lcmsrun_observed'] = file_input['lcmsrun']

        ms1_data.append(data['ms1_summary'])

    return pd.concat(ms1_data)


def align_rt_adjustment_peaks(ms1_data: pd.DataFrame, rt_alignment_atlas: pd.DataFrame) -> tuple[list[float], list[float]]:
    """align median experimental retention time peaks with rt adjustment atlas peaks."""

    median_experimental_rt_peaks = ms1_data[ms1_data['peak_height'] >= 1e4].groupby('label')['rt_peak'].median()
    rt_peaks_merged = pd.merge(rt_alignment_atlas[['label', 'rt_peak']], median_experimental_rt_peaks, on='label')

    original_rt_peaks = rt_peaks_merged['rt_peak_x'].tolist()
    experimental_rt_peaks = rt_peaks_merged['rt_peak_y'].tolist()

    return (original_rt_peaks, experimental_rt_peaks)


def adjust_template_atlas_rt_peaks(template_atlas: pd.DataFrame, original_rt_peaks: list[float], experimental_rt_peaks: list[float],
                                   rt_regression: bool, model_degree: int, rt_window: float) -> pd.DataFrame:
    """Build and use model to adjust template atlas retention time peaks to match experimental retention time space."""

    aligned_template_atlas = template_atlas.copy()

    if rt_regression:
        rt_alignment_model = np.polyfit(original_rt_peaks, experimental_rt_peaks, model_degree)
        aligned_template_atlas['rt_peak'] = aligned_template_atlas['rt_peak'].apply(lambda x: np.polyval(rt_alignment_model, x))

    else:
        median_offset = np.median(np.array(original_rt_peaks) - np.array(experimental_rt_peaks))
        aligned_template_atlas['rt_peak'] = aligned_template_atlas['rt_peak'] + median_offset

    aligned_template_atlas['rt_min'] = aligned_template_atlas['rt_peak'] - rt_window
    aligned_template_atlas['rt_max'] = aligned_template_atlas['rt_peak'] + rt_window

    return aligned_template_atlas
