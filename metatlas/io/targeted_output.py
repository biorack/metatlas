"""Generate standarized outputs for targeted analysis"""
# pylint: disable=too-many-arguments

import logging
import os
import tarfile

from collections import namedtuple

import numpy as np
import pandas as pd

from IPython.core.display import display, HTML

from metatlas.io import rclone
from metatlas.io.write_utils import export_dataframe_die_on_diff
from metatlas.plots import dill2plots as dp
from metatlas.tools import fastanalysis as fa

logger = logging.getLogger(__name__)

RCLONE_PATH = "/global/cfs/cdirs/m342/USA/shared-repos/rclone/bin/rclone"


def write_atlas_to_spreadsheet(metatlas_dataset, overwrite=False):
    """Save atlas as csv file. Will not overwrite existing file unless overwrite is True"""
    out_file_name = os.path.join(metatlas_dataset.ids.output_dir, f"{metatlas_dataset.atlas.name}_export.csv")
    out_df = dp.export_atlas_to_spreadsheet(metatlas_dataset.atlas)
    export_dataframe_die_on_diff(out_df, out_file_name, "atlas", overwrite=overwrite, float_format="%.6e")


def write_stats_table(
    metatlas_dataset,
    min_intensity=1e4,
    rt_tolerance=0.5,
    mz_tolerance=20,
    min_msms_score=0.6,
    min_num_frag_matches=1,
    min_relative_frag_intensity=0.001,
    allow_no_msms=True,
    overwrite=False,
):
    """
    inputs:
        metatlas_dataset: a MetatlasDataset instance
        min_intensity: intensity threshold; 1e5 is strict, 1e3 is loose
        rt_tolerance: RT tolerance threshold
                      shift of median RT across all files for given compound to reference
        mz_tolerance: MZ tolerance threshold
                      ppm of median mz across all files for given compound relative to reference
                      5 is strict, 25 is loose
        min_msms_score: score threshold
                        max dot-product score across all files for given compound relative to reference
                        Score values in [0-1]; 0.6 is strict, 0.3 is loose
        min_num_frag_matches: threshold of number of frag matches between sample and reference
        min_relative_frag_intensity: threshold ratio of second highest to first highest intensity
                                     of matching sample mzs
        allow_no_msms: if True evaluate only on MS1 thresholds if no MSMS data is found,
                       if False filter out row if MSMS thresholds are not passing
        overwrite: if True, will write over existing files
    """
    ids = metatlas_dataset.ids
    prefix = f"{ids.short_polarity}_"
    scores_path = os.path.join(ids.output_dir, f"{prefix}stats_tables", f"{prefix}compound_scores.csv")
    _ = metatlas_dataset.hits  # regenerate hits if needed before logging about scores
    logger.info("Calculating scores and exporting them to %s.", scores_path)
    scores_df = fa.make_scores_df(metatlas_dataset, metatlas_dataset.hits)
    scores_df["passing"] = fa.test_scores_df(
        scores_df,
        min_intensity,
        rt_tolerance,
        mz_tolerance,
        min_msms_score,
        allow_no_msms,
        min_num_frag_matches,
        min_relative_frag_intensity,
    )
    export_dataframe_die_on_diff(scores_df, scores_path, "scores", overwrite=overwrite, float_format="%.8e")
    fa.make_stats_table(
        input_dataset=metatlas_dataset,
        msms_hits=metatlas_dataset.hits,
        output_loc=ids.output_dir,
        output_sheetname=f"{ids.project}_{ids.output_type}_Identifications.xlsx",
        min_peak_height=1e5,
        use_labels=True,
        min_msms_score=0.01,
        min_num_frag_matches=1,
        include_lcmsruns=[],
        exclude_lcmsruns=["QC"],
        polarity=ids.short_polarity,
        overwrite=overwrite,
    )


def write_chromatograms(metatlas_dataset, overwrite=False, max_cpus=1):
    """
    inputs:
        metatlas_dataset: a MetatlasDataset instance
        group_by: 'index', 'page', or None for grouping of plots
        overwrite: if False raise error if file already exists
    """
    # overwrite checks done within dp.make_chromatograms
    logger.info("Exporting chromotograms to %s", metatlas_dataset.ids.output_dir)
    params = {
        "input_dataset": metatlas_dataset,
        "include_lcmsruns": [],
        "exclude_lcmsruns": ["InjBl", "QC", "Blank", "blank"],
        "share_y": True,
        "output_loc": metatlas_dataset.ids.output_dir,
        "polarity": metatlas_dataset.ids.short_polarity,
        "overwrite": overwrite,
        "max_cpus": max_cpus,
        "suffix": "_sharedY",
    }
    dp.make_chromatograms(**params)
    params["share_y"] = False
    params["suffix"] = "_independentY"
    dp.make_chromatograms(**params)


def write_identification_figure(metatlas_dataset, overwrite=False):
    """Save identification figure. Will not overwrite existing file unless overwrite is True"""
    # overwrite checks done within dp.make_identification_figure_v2
    logger.info("Exporting indentification figures to %s", metatlas_dataset.ids.output_dir)
    dp.make_identification_figure_v2(
        input_dataset=metatlas_dataset,
        msms_hits=metatlas_dataset.hits,
        use_labels=True,
        include_lcmsruns=[],
        exclude_lcmsruns=["InjBl", "QC", "Blank", "blank"],
        output_loc=metatlas_dataset.ids.output_dir,
        short_names_df=metatlas_dataset.ids.lcmsruns_short_names,
        polarity=metatlas_dataset.ids.short_polarity,
        overwrite=overwrite,
    )


def write_metrics_and_boxplots(metatlas_dataset, overwrite=False, max_cpus=1):
    """
    Save metrics dataframes as csv and boxplots as PDF.
    Will not overwrite existing file unless overwrite is True
    """
    config = [
        {"name": "peak_height", "label": "Peak Height"},
        {"name": "peak_area", "label": None},
        {"name": "mz_peak", "label": None},
        {"name": "rt_peak", "label": "RT Peak"},
        {"name": "mz_centroid", "label": "MZ Centroid"},
        {"name": "rt_centroid", "label": None},
    ]
    prefix = f"{metatlas_dataset.ids.short_polarity}_"
    for fields in config:
        df_dir = os.path.join(metatlas_dataset.ids.output_dir, f"{prefix}data_sheets")
        dataframe = dp.make_output_dataframe(
            fieldname=fields["name"],
            input_dataset=metatlas_dataset,
            output_loc=df_dir,
            short_names_df=metatlas_dataset.ids.lcmsruns_short_names,
            polarity=metatlas_dataset.ids.short_polarity,
            use_labels=True,
            overwrite=overwrite,
        )
        if fields["label"] is not None:
            for logy in [False, True]:
                plot_dir = os.path.join(
                    metatlas_dataset.ids.output_dir,
                    f"{prefix}boxplot_{fields['name']}" + '_log' if logy else '',
                )
                dp.make_boxplot_plots(
                    dataframe,
                    output_loc=plot_dir,
                    use_shortnames=True,
                    ylabel=fields["label"],
                    overwrite=overwrite,
                    max_cpus=max_cpus,
                    logy=logy
                )


Max = namedtuple("Max", ["file_idx", "pre_intensity_idx", "pre_intensity", "precursor_mz"])


def write_msms_fragment_ions(
    data, intensity_fraction=0.01, min_mz=450, max_mz_offset=5, scale_intensity=1e5, overwrite=False
):
    """
    inputs:
        data: metatlas_datset
        intensity_fraction: intensity threshold as fraction of max_msms_intensity (0-1]
        min_mz: minimum threshold MSMS mz value
        max_mz: maximum threshold MSMS mz value. Relative to precursor mz with highest intensity
        scale_intensity: If not None, normalize output intensity to maximum of scale_intensity
    """
    out = []
    for compound_idx, _ in enumerate(data[0]):
        max_vars = get_max_precursor_intensity(data, compound_idx)
        out.append(
            get_spectra_strings(
                data[max_vars.file_idx][compound_idx],
                max_vars.pre_intensity,
                min_mz,
                max_mz_offset + max_vars.precursor_mz,
                intensity_fraction,
                scale_intensity,
            )
        )
    out_df = pd.DataFrame(out)
    path = os.path.join(data.ids.output_dir, f"spectra_{intensity_fraction:.2f}pct_{int(min_mz)}cut.csv")
    export_dataframe_die_on_diff(out_df, path, "MSMS fragment ions", overwrite=overwrite, float_format="%.8e")
    return out_df


def get_max_precursor_intensity(data, compound_idx):
    """
    inputs:
        data: metatlas_dataset
        compound_idx: index of compound to search over
    returns Max object with file index of highest precursor intensity, associated intensity value, and mz
    """
    max_pre_intensity = max_precursor_mz = 0
    max_file_idx = max_pre_intensity_idx = None
    for file_idx, _ in enumerate(data):
        try:
            msms = data[file_idx][compound_idx]["data"]["msms"]["data"]
            if len(msms["precursor_intensity"]) == 0:
                continue
            pre_intensity_idx = msms["precursor_intensity"].argmax()
            pre_intensity = msms["precursor_intensity"][pre_intensity_idx]
            precursor_mz = msms["precursor_MZ"][pre_intensity_idx]
            rts = msms["rt"][pre_intensity_idx]
            rt_ref = data[file_idx][compound_idx]["identification"].rt_references[-1]
            if pre_intensity > max_pre_intensity and rt_ref.rt_min < rts < rt_ref.rt_max:
                max_file_idx = file_idx
                max_pre_intensity_idx = pre_intensity_idx
                max_pre_intensity = pre_intensity
                max_precursor_mz = precursor_mz
        except (AttributeError, IndexError):
            pass
    return Max(max_file_idx, max_pre_intensity_idx, max_pre_intensity, max_precursor_mz)


def get_spectra_strings(data, max_pre_intensity, min_mz, max_mz, intensity_fraction, scale_intensity):
    """
    inputs:
        data: metatlas_dataset[x][y]
        sample_idx: first index into data
        compound_idx: second into into data
        max_pre_intensity: highest msms precursor intensity for this compound across all samples
        min_mz: minimum threshold MSMS mz value
        max_mz: maximum threshold MSMS mz value
        intensity_fraction: intensity threshold as fraction of max_msms_intensity (0-1]
        scale_intensity: If not None, normalize output intensity to maximum of scale_intensity
    returns a dict containing compound name and string representations of the spectra
    """
    mz_list, intensity_list = get_spectra(
        data, max_pre_intensity, min_mz, max_mz, intensity_fraction, scale_intensity
    )
    mz_str = str([f"{x:.2f}" for x in mz_list]).replace("'", "")
    intensity_str = str([int(x) for x in intensity_list]).replace("'", "")
    spectra_str = str([mz_str, intensity_str]).replace("'", "")
    name = data["identification"].name
    return {"name": name, "spectrum": spectra_str, "mz": mz_str, "intensity": intensity_str}


def get_spectra(data, max_pre_intensity, min_mz, max_mz, intensity_fraction, scale_intensity):
    """
    inputs:
        data: metatlas_dataset[i][j]
        max_pre_intensity: highest msms precursor intensity for this compound across all samples
        min_mz: minimum threshold MSMS mz value
        max_mz: maximum threshold MSMS mz value
        intensity_fraction: intensity threshold as fraction of max_msms_intensity (0-1]
        scale_intensity: If not None, normalize output intensity to maximum of scale_intensity
    returns a tuple containing a list of mz values and a list intensity values that make a spectra
    returns None, None if no spectra meet the filtering thresholds
    """
    if max_pre_intensity != 0:
        msms = data["data"]["msms"]["data"]
        idx = np.argwhere(msms["precursor_intensity"] == max_pre_intensity).flatten()
        msms_mz = msms["mz"][idx]
        intensity = msms["i"][idx]
        max_msms_intensity = intensity.max()
        cutoff = intensity_fraction * max_msms_intensity
        conditions = (intensity > cutoff) & (min_mz < msms_mz) & (msms_mz < max_mz)
        if any(conditions):
            keep_idx = np.argwhere(conditions).flatten()
            if scale_intensity is not None:
                intensity = (intensity / max_msms_intensity * scale_intensity).astype(int)
            return msms_mz[keep_idx], intensity[keep_idx]
    return None, None


def archive_outputs(ids):
    """
    Creates a .tar.gz file containing all output files
    Inputs:
        ids: an AnalysisIds object
    """
    logger.info("Generating archive of output files.")
    suffix = "" if ids.output_type == "data_QC" else f"-{ids.short_polarity}"
    output_file = f"{ids.short_experiment_analysis}{suffix}.tar.gz"
    output_path = os.path.join(ids.project_directory, ids.experiment, output_file)
    with tarfile.open(output_path, "w:gz") as tar:
        tar.add(ids.output_dir, arcname=os.path.basename(ids.output_dir))
    logger.info("Generation of archive completed succesfully: %s", output_path)


def copy_outputs_to_google_drive(ids):
    """
    Recursively copy the output files to Google Drive using rclone
    Inputs:
        ids: an AnalysisIds object
    """
    logger.info("Copying output files to Google Drive")
    rci = rclone.RClone(RCLONE_PATH)
    fail_suffix = "not copying files to Google Drive"
    if rci.config_file() is None:
        logger.warning("RClone config file not found -- %s.", fail_suffix)
        return
    drive = rci.get_name_for_id(ids.google_folder)
    if drive is None:
        logger.warning("RClone config file missing JGI_Metabolomics_Projects -- %s.", fail_suffix)
        return
    folders = [ids.experiment, ids.analysis, ids.output_type]
    if ids.output_type != "data_QC":
        folders.append(ids.short_polarity)
    sub_folders_string = os.path.join("Analysis_uploads", *folders)
    rci.copy_to_drive(ids.output_dir, drive, sub_folders_string)
    logger.info("Done copying output files to Google Drive")
    path_string = f"{drive}:{sub_folders_string}"
    display(
        HTML(f'Data is now on Google Drive at <a href="{rci.path_to_url(path_string)}">{path_string}</a>')
    )
