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
from metatlas.tools.config import Analysis
from metatlas.tools.notebook import in_papermill
from metatlas.tools import spectralprocessing as sp

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

    ms2_data_refs_merge['match_reference_jaccard'] = ms2_data_refs_merge.apply(lambda row: sp.jaccard_set([round(x, 2) for x in row['spectrum_x'][0]], [round(y, 2) for y in row['spectrum_y'][0]]), axis=1)
    ms2_data_refs_merge['ref_frags'] = ms2_data_refs_merge['spectrum_y'].apply(lambda x: len(x[0]) if x.any() else 0)
    ms2_data_refs_merge['data_frags'] = ms2_data_refs_merge['spectrum_x'].apply(lambda x: len(x[0]) if x.any() else 0)
    ms2_data_refs_merge['reference_data_ratio'] = ms2_data_refs_merge['ref_frags'] / ms2_data_refs_merge['data_frags']
    ms2_data_refs_merge['match_reference_ratio'] = ms2_data_refs_merge['matches'] / ms2_data_refs_merge['ref_frags']

    ms2_data_refs_merge = ms2_data_refs_merge.drop(['ref_frags', 'data_frags'], axis=1)

    return ms2_data_refs_merge


def filter_atlas_labels(ms1_data: pd.DataFrame, ms2_data_scored: pd.DataFrame, peak_height: Optional[float], num_points: Optional[int],
                        msms_score: Optional[float], msms_matches: Optional[int], msms_frag_ratio: Optional[float], msms_frag_jaccard: Optional[float]) -> set[str]:
    """Filter atlas labels to include only those that pass the MS1 and MS2 thresholds."""

    ms1_data_filtered = ms1_data[ms1_data['peak_height'] >= peak_height] if peak_height is not None else ms1_data
    ms1_data_filtered = ms1_data_filtered[ms1_data_filtered['num_datapoints'] >= num_points] if num_points is not None else ms1_data_filtered
    ms1_reduced_labels = set(ms1_data_filtered.label.tolist())

    if msms_score is not None or msms_matches is not None or msms_frag_ratio is not None or msms_frag_jaccard is not None:
        ms2_data_filtered = ms2_data_scored[ms2_data_scored['score'] >= msms_score] if msms_score is not None else ms2_data_scored
        ms2_data_filtered = ms2_data_filtered[ms2_data_filtered['matches'] >= msms_matches] if msms_matches is not None else ms2_data_filtered
        ms2_data_filtered = ms2_data_filtered[ms2_data_filtered['match_reference_ratio'] >= msms_frag_ratio] if msms_frag_ratio is not None else ms2_data_filtered
        ms2_data_filtered = ms2_data_filtered[ms2_data_filtered['match_reference_jaccard'] >= msms_frag_jaccard] if msms_frag_jaccard is not None else ms2_data_filtered
        ms2_reduced_labels = set(ms2_data_filtered.label.tolist())
    else:
        ms2_reduced_labels = ms1_reduced_labels
    
    reduced_labels = list(ms1_reduced_labels.intersection(ms2_reduced_labels))

    # Put fragment ratio and jaccard similarity into the return to be appended to aligned filtered atlas
    if msms_frag_ratio is not None and msms_frag_jaccard is not None and ms2_data_filtered is not None:
        frag_info_df = pd.DataFrame({
            'label': list(reduced_labels),
            'max_match_reference_ratio': [ms2_data_filtered[ms2_data_filtered['label'] == label]['match_reference_ratio'].max() for label in reduced_labels],
            'max_jaccard_similarity': [ms2_data_filtered[ms2_data_filtered['label'] == label]['match_reference_jaccard'].max() for label in reduced_labels]
        })
        return reduced_labels, frag_info_df
    
    else:
        return reduced_labels, None


def filter_atlas(aligned_atlas, ids: AnalysisIdentifiers, analysis: Analysis, data: MetatlasDataset):

    aligned_atlas_df = make_atlas_df(aligned_atlas)
    sample_file_paths = get_sample_file_paths(ids)

    if analysis.parameters.mz_tolerance_override is not None:
        mz_tolerance = analysis.parameters.mz_tolerance_override
    else:
        mz_tolerance = analysis.parameters.mz_tolerance_default

    ms1_data, ms2_data = get_sample_data(aligned_atlas_df, sample_file_paths, mz_tolerance, data.extra_time, analysis.parameters.polarity)
    ms2_data_scored = score_ms2_data(ms2_data, aligned_atlas_df, analysis.parameters.polarity,
                                     analysis.parameters.msms_refs, analysis.parameters.frag_mz_tolerance)

    reduced_labels, msms_frag_info = filter_atlas_labels(ms1_data, ms2_data_scored, analysis.parameters.peak_height, analysis.parameters.num_points,
                                         analysis.parameters.msms_score, analysis.parameters.msms_matches, analysis.parameters.msms_frag_ratio,
                                         analysis.parameters.msms_frag_jaccard)

    aligned_filtered_atlas_df = aligned_atlas_df[aligned_atlas_df['label'].isin(reduced_labels)]
    
    if msms_frag_info is not None:
        aligned_filtered_atlas_df = aligned_filtered_atlas_df.merge(msms_frag_info, on='label', how='left')
    
    aligned_filtered_atlas = dp.get_atlas(aligned_atlas.name, aligned_filtered_atlas_df, analysis.parameters.polarity, mz_tolerance)

    return aligned_filtered_atlas
