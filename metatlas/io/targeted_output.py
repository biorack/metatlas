"""Generate standarized outputs for targeted analysis"""
# pylint: disable=too-many-arguments

import logging
import math

from collections import namedtuple
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

from matplotlib import gridspec
from matplotlib.axis import Axis
from tqdm.notebook import tqdm

from metatlas.datastructures.metatlas_dataset import MetatlasDataset
from metatlas.io.write_utils import export_dataframe_die_on_diff
from metatlas.plots import dill2plots as dp
from metatlas.tools import fastanalysis as fa
from metatlas.tools.config import Analysis, Workflow
from metatlas.tools.notebook import in_papermill
from metatlas.plots.tic import save_sample_tic_pdf

logger = logging.getLogger(__name__)


def write_atlas_to_csv(metatlas_dataset, overwrite=False):
    """Save atlas as csv file. Will not overwrite existing file unless overwrite is True"""
    prefix = "CompoundAtlas__"
    out_file_name = metatlas_dataset.ids.output_dir / f"{prefix}{metatlas_dataset.atlas.name}.csv"
    out_df = dp.export_atlas_to_spreadsheet(metatlas_dataset.atlas)
    export_dataframe_die_on_diff(out_df, out_file_name, "atlas", overwrite=overwrite, float_format="%.6e")


def write_identifications_spreadsheet(
    metatlas_dataset,
    analysis_parameters,
    metatlas_workflow,
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
    ids.set_output_state(analysis_parameters, "ids_spreadsheet")
    prefix = f"{ids.short_polarity}_"
    parent = ids.additional_output_dir
    scores_path = parent / f"{prefix}stats_tables" / f"{prefix}compound_scores.csv"
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
    output_sheetname = f"{ids.project}_{ids.workflow}_{ids.analysis}_Identifications.xlsx"
    fa.make_stats_table(
        workflow_name=metatlas_workflow.name,
        msms_sorting_method=analysis_parameters.msms_sorting_method,
        input_dataset=metatlas_dataset,
        msms_hits=metatlas_dataset.hits,
        output_loc=ids.additional_output_dir,
        output_sheetname=output_sheetname,
        min_peak_height=1e5,
        use_labels=True,
        min_msms_score=0.01,
        min_num_frag_matches=1,
        ppm_tolerance = mz_tolerance,
        polarity=ids.short_polarity,
        overwrite=overwrite,
        data_sheets=False,  # eliminates output of the legacy 'data_sheets' dir but not '{POL}_data_sheets'
    )
    xlsx = ids.additional_output_dir / output_sheetname
    xlsx.rename(ids.output_dir / xlsx.name)


def write_chromatograms(metatlas_dataset, analysis_parameters, overwrite=False, max_cpus=1):
    """
    inputs:
        metatlas_dataset: a MetatlasDataset instance
        group_by: 'index', 'page', or None for grouping of plots
        overwrite: if False raise error if file already exists
    """
    # overwrite checks done within dp.make_chromatograms
    metatlas_dataset.ids.set_output_state(analysis_parameters, "chromatograms")
    logger.info("Exporting chromatograms with shared Y-axis.")
    params = {
        "input_dataset": metatlas_dataset,
        "share_y": True,
        "output_loc": metatlas_dataset.ids.additional_output_dir,
        "polarity": metatlas_dataset.ids.short_polarity,
        "overwrite": overwrite,
        "max_cpus": max_cpus,
        "suffix": "_sharedY",
    }
    dp.make_chromatograms(**params)
    logger.info("Exporting chromatograms with independent Y-axis.")
    params["share_y"] = False
    params["suffix"] = "_independentY"
    dp.make_chromatograms(**params)


def write_tics(metatlas_dataset, x_min=None, x_max=None, y_min=0, overwrite=False):
    """
    Create PDF files with TIC plot for each sample
    One file with shared Y-axes, one file with independent Y-axes
    """
    parent = metatlas_dataset.ids.additional_output_dir
    prefix = f"{metatlas_dataset.ids.short_polarity}_"
    for suffix, sharey in [("_independentY", False), ("_sharedY", True)]:
        file_name = parent / f"{prefix}TICs{suffix}.pdf"
        save_sample_tic_pdf(
            metatlas_dataset,
            metatlas_dataset.ids.polarity,
            file_name,
            overwrite,
            x_min=x_min,
            x_max=x_max,
            y_min=y_min,
            sharey=sharey,
        )


def write_identification_figure(metatlas_dataset, analysis_parameters, overwrite=False):
    """Save identification figure. Will not overwrite existing file unless overwrite is True"""
    # overwrite checks done within dp.make_identification_figure_v2
    ids = metatlas_dataset.ids
    logger.info("Exporting indentification figures to %s", ids.output_dir)
    ids.set_output_state(analysis_parameters, "chromatograms")
    dp.make_identification_figure_v2(
        input_dataset=metatlas_dataset,
        msms_hits=metatlas_dataset.hits,
        use_labels=True,
        output_loc=metatlas_dataset.ids.output_dir,
        short_names_df=metatlas_dataset.ids.lcmsruns_short_names,
        polarity=metatlas_dataset.ids.short_polarity,
        overwrite=overwrite,
    )


def write_metrics_and_boxplots(metatlas_dataset, analysis_parameters, overwrite=False, max_cpus=1):
    """
    Save metrics dataframes as csv and boxplots as PDF.
    Will not overwrite existing file unless overwrite is True
    """
    config = [
        {"name": "peak_height", "label": "Peak Height", "additional": False},
        {"name": "peak_area", "label": None, "additional": True},
        {"name": "mz_peak", "label": None, "additional": True},
        {"name": "rt_peak", "label": "RT Peak", "additional": True},
        {"name": "mz_centroid", "label": "MZ Centroid", "additional": True},
        {"name": "rt_centroid", "label": None, "additional": True},
    ]
    ids = metatlas_dataset.ids
    prefix = f"{metatlas_dataset.ids.short_polarity}_"
    df_dir = ids.output_dir / f"{prefix}data_sheets"
    ids.set_output_state(analysis_parameters, "data_sheets")
    dfs = [
        dp.make_output_dataframe(
            fieldname=fields["name"],
            input_dataset=metatlas_dataset,
            output_loc=df_dir,
            short_names_df=ids.lcmsruns_short_names,
            polarity=ids.short_polarity,
            use_labels=True,
            overwrite=overwrite,
        )
        for fields in config
    ]
    ids.set_output_state(analysis_parameters, "box_plots")
    for dataframe, fields in zip(dfs, config):
        if fields["label"] is not None:
            for logy in [False, True]:
                parent = ids.additional_output_dir if fields["additional"] else ids.output_dir
                plot_dir = parent / f"{prefix}boxplot_{fields['name']}{'_log' if logy else ''}"
                dp.make_boxplot_plots(
                    dataframe,
                    output_loc=plot_dir,
                    use_shortnames=True,
                    ylabel=fields["label"],
                    overwrite=overwrite,
                    max_cpus=max_cpus,
                    logy=logy,
                )


Max = namedtuple("Max", ["file_idx", "pre_intensity_idx", "pre_intensity", "precursor_mz"])


def write_msms_fragment_ions(
    data, intensity_fraction=0.00001, min_mz=0, max_mz_offset=5, scale_intensity=1e5, overwrite=False
):
    """
    inputs:
        data: metatlas_datset
        intensity_fraction: intensity threshold as fraction of max_msms_intensity (0-1]
        min_mz: minimum threshold MSMS mz value
        max_mz: maximum threshold MSMS mz value. Relative to precursor mz with highest intensity
        scale_intensity: If not None, normalize output intensity to maximum of scale_intensity
    """
    full_out = []
    intensity_fraction_out = []
    for compound_idx, _ in enumerate(data[0]):
        max_vars = get_max_precursor_intensity(data, compound_idx)
        if max_vars.file_idx is not None:
            spectra_string_dict = get_spectra_strings(
                                        data,
                                        max_vars.file_idx,
                                        compound_idx,
                                        max_vars.pre_intensity,
                                        min_mz,
                                        max_mz_offset + max_vars.precursor_mz,
                                        intensity_fraction,
                                        scale_intensity,
                                    )
            for key, spectra_dict in spectra_string_dict.items():
                if key == 0:
                    full_out.append(spectra_dict)
                elif key == intensity_fraction:
                    intensity_fraction_out.append(spectra_dict)
                else:
                    raise ValueError(f"Unexpected key: {key}")
    full_out_df = pd.DataFrame(full_out)
    intensity_fraction_out_df = pd.DataFrame(intensity_fraction_out)
    full_out_path = data.ids.output_dir / f"spectra_top1000_{int(min_mz)}cut.csv"
    intensity_fraction_out_path = data.ids.output_dir / f"spectra_{intensity_fraction}pct_{int(min_mz)}cut.csv"
    export_dataframe_die_on_diff(full_out_df, full_out_path, "MSMS fragment ions, top 1000", overwrite=overwrite, float_format="%.8e", header=True)
    export_dataframe_die_on_diff(intensity_fraction_out_df, intensity_fraction_out_path, f"MSMS fragment ions, {intensity_fraction}pct cutoff", overwrite=overwrite, float_format="%.8e", header=True)
    return
    #return full_out_df

    # Write two dataframes with different mz significant figures
    # if "mz4" in out_df.columns:
    #     out_df_mz2 = out_df.drop(columns=["mz4"])
    #     out_df_mz2 = out_df_mz2.rename(columns={"mz2": "mz"})
    #     path = data.ids.output_dir / f"spectra_{intensity_fraction:.2f}pct_{int(min_mz)}cut_mz2.csv"
    #     export_dataframe_die_on_diff(out_df_mz2, path, "MSMS fragment ions", overwrite=overwrite, float_format="%.8e")
    # if "mz2" in out_df.columns:
    #     out_df_mz4 = out_df.drop(columns=["mz2"])
    #     out_df_mz4 = out_df_mz4.rename(columns={"mz4": "mz"})
    #     path = data.ids.output_dir / f"spectra_{intensity_fraction:.2f}pct_{int(min_mz)}cut_mz4.csv"
    #     export_dataframe_die_on_diff(out_df_mz4, path, "MSMS fragment ions", overwrite=overwrite, float_format="%.8e")

def get_max_precursor_intensity(data, compound_idx):
    """
    inputs:
        data: metatlas_dataset
        compound_idx: index of compound to search over
    returns Max object with file index of highest precursor intensity, associated intensity value, and mz
    """
    max_pre_intensity = max_precursor_mz = 0
    max_file_idx = max_pre_intensity_idx = None
    for file_idx, sample in enumerate(data):
        try:
            msms = sample[compound_idx]["data"]["msms"]["data"]
            if len(msms["precursor_intensity"]) == 0:
                continue
            pre_intensity_idx = msms["precursor_intensity"].argmax()
            pre_intensity = msms["precursor_intensity"][pre_intensity_idx]
            precursor_mz = msms["precursor_MZ"][pre_intensity_idx]
            rts = msms["rt"][pre_intensity_idx]
            rt_ref = sample[compound_idx]["identification"].rt_references[-1]
            if pre_intensity > max_pre_intensity and rt_ref.rt_min < rts < rt_ref.rt_max:
                max_file_idx = file_idx
                max_pre_intensity_idx = pre_intensity_idx
                max_pre_intensity = pre_intensity
                max_precursor_mz = precursor_mz
        except (AttributeError, IndexError, KeyError):
            pass
    return Max(max_file_idx, max_pre_intensity_idx, max_pre_intensity, max_precursor_mz)


def get_spectra_strings(
    data, sample_idx, compound_idx, max_pre_intensity, min_mz, max_mz, intensity_fraction, scale_intensity
):
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
    compound_data = data[sample_idx][compound_idx]
    exported_fragment_dict = get_spectra(
        compound_data, max_pre_intensity, min_mz, max_mz, intensity_fraction, scale_intensity
    )
    spectra_string_dict = {}
    for key, (mz_list, intensity_list) in exported_fragment_dict.items():
        mz_str_2 = str([f"{x:.2f}" for x in mz_list]).replace("'", "")
        mz_str_4 = str([f"{x:.4f}" for x in mz_list]).replace("'", "")
        intensity_str = str([int(x) for x in intensity_list]).replace("'", "")
        spectra_str = str([mz_str_2, mz_str_4, intensity_str]).replace("'", "")
        name = compound_data["identification"].name
        run = Path(compound_data["lcmsrun"].hdf5_file).stem
        mz_peak = compound_data["data"]["ms1_summary"]["mz_peak"]
        rt_peak = compound_data["data"]["ms1_summary"]["rt_peak"]
        spectra_string_dict[key] = {
            "compound_idx": compound_idx,
            "name": name,
            "run": run,
            "rt_peak": rt_peak,
            "mz_peak": mz_peak,
            "spectrum": spectra_str,
            "mz2": mz_str_2,
            "mz4": mz_str_4,
            "intensity": intensity_str,
        }
    return spectra_string_dict


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
        exported_fragments = {
            intensity_fraction: (None, None),
            0: (None, None)
        }
        for key in exported_fragments.keys():
            cutoff = key * max_msms_intensity
            conditions = (intensity > cutoff) & (min_mz < msms_mz) & (msms_mz < max_mz)
            if any(conditions):
                keep_idx = np.argwhere(conditions).flatten()
                msms_mz = msms_mz[keep_idx]
                intensity = intensity[keep_idx]
                if key == 0:
                    sorted_indices = np.argsort(intensity)[::-1][:1000]                    
                else:
                    sorted_indices = np.argsort(intensity)[::-1]
                msms_mz_keep = msms_mz[sorted_indices]
                intensity_keep = intensity[sorted_indices]
                if scale_intensity is not None:
                    intensity_keep = (intensity_keep / max_msms_intensity * scale_intensity).astype(int)
                exported_fragments[key] = (msms_mz_keep, intensity_keep)
            else:
                exported_fragments[key] = (None, None)
        return exported_fragments
    else:
        exported_fragments[key] = (None, None)
        return exported_fragments


def generate_standard_outputs(
    data: MetatlasDataset,
    workflow: Workflow,
    analysis: Analysis,
    overwrite: bool = False,
) -> None:
    """Generates the default set of outputs for a targeted experiment"""
    params = analysis.parameters
    data.ids.set_output_state(params, "gui")
    write_atlas_to_csv(data, overwrite=overwrite)
    write_tics(data, overwrite=overwrite, x_min=1.5)
    write_identifications_spreadsheet(data, params, workflow, overwrite=overwrite)
    write_chromatograms(data, params, overwrite=overwrite, max_cpus=data.max_cpus)
    write_identification_figure(data, params, overwrite=overwrite)
    write_metrics_and_boxplots(data, params, overwrite=overwrite, max_cpus=data.max_cpus)
    logger.info("Generation of standard output files completed sucessfully.")


def generate_qc_plots(data: MetatlasDataset) -> None:
    """Write plots that can be used to QC the experiment"""
    rts_df = get_rts(data)
    compound_atlas_rts_file_name = data.ids.output_dir / f"{data.ids.short_polarity}_Compound_Atlas_RTs.pdf"
    plot_compound_atlas_rts(len(data), rts_df, compound_atlas_rts_file_name)
    peak_heights_df = get_peak_heights(data)
    peak_heights_plot_file_name = (
        data.ids.output_dir / f"{data.ids.short_polarity}_Compound_Atlas_peak_heights.pdf"
    )
    plot_compound_atlas_peak_heights(len(data), peak_heights_df, peak_heights_plot_file_name)


def generate_qc_outputs(data: MetatlasDataset, analysis: Analysis) -> None:
    """Write outputs that can be used to QC the experiment"""
    ids = data.ids
    ids.set_output_state(analysis.parameters, "ids_spreadsheet")
    save_rt_peak(data, ids.output_dir / f"{ids.short_polarity}_rt_peak.tab")
    save_measured_rts(data, ids.output_dir / f"{ids.short_polarity}_QC_Measured_RTs.csv")
    generate_qc_plots(data)


def save_measured_rts(data: MetatlasDataset, file_name: Path) -> None:
    """Save RT values in csv format file"""
    rts_df = get_rts(data, include_atlas_rt_peak=False)
    export_dataframe_die_on_diff(rts_df, file_name, "measured RT values", float_format="%.6e")


def save_rt_peak(data: MetatlasDataset, file_name: Path) -> None:
    """Save peak RT values in tsv format file"""
    rts_df = dp.make_output_dataframe(input_dataset=data, fieldname="rt_peak", use_labels=True)
    export_dataframe_die_on_diff(rts_df, file_name, "peak RT values", sep="\t", float_format="%.6e")


def get_rts(data: MetatlasDataset, include_atlas_rt_peak: bool = True) -> pd.DataFrame:
    """Returns RT values in DataFrame format"""
    rts_df = dp.make_output_dataframe(
        input_dataset=data,
        fieldname="rt_peak",
        use_labels=True,
        summarize=True,
    )
    if include_atlas_rt_peak:
        rts_df["atlas RT peak"] = [
            compound["identification"].rt_references[0].rt_peak for compound in data[0]
        ]
    return order_df_columns_by_run(rts_df)


def get_peak_heights(data: MetatlasDataset) -> pd.DataFrame:
    """Returns peak heights in DataFrame format"""
    peak_height_df = dp.make_output_dataframe(
        input_dataset=data,
        fieldname="peak_height",
        use_labels=True,
        summarize=True,
    )
    return order_df_columns_by_run(peak_height_df)


def order_df_columns_by_run(dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a dataframe with re-ordered columns such that second column up to column 'mean'
    are ordered by run number from low to high
    """
    cols = dataframe.columns.tolist()
    stats_start_idx = cols.index("mean")
    to_sort = cols[:stats_start_idx]
    no_sort = cols[stats_start_idx:]
    to_sort.sort(
        key=lambda x: int(
            x.split(".")[0].split("_")[-1].lower().replace("run", "").replace("seq", "").replace("s", "")
        )
    )
    new_cols = to_sort + no_sort
    return dataframe[new_cols]


def plot_per_compound(
    field_name: str,
    num_files: int,
    data: pd.DataFrame,
    file_name: Path,
    fontsize: float = 2,
    pad: float = 0.1,
    cols: int = 8,
) -> None:
    """
    Writes plot of RT peak for vs file for each compound
    inputs:
        field_name: one of rt_peak or peak_height
        num_files: number of files in data set, ie len(data)
        data: Dataframe with RTs values
        file_name: where to save plot
        fontsize: size of text
        pad: padding size
        cols: number of columns in plot grid
    """
    logger.info("Plotting %s vs file for each compound", field_name)
    plot_df = (
        data.sort_values(by="standard deviation", ascending=False, na_position="last")
        .drop(["#NaNs"], axis=1)
        .dropna(axis=0, how="all", subset=data.columns[:num_files])
    )
    rows = int(math.ceil((data.shape[0] + 1) / cols))
    fig = plt.figure()
    grid = gridspec.GridSpec(rows, cols, figure=fig, wspace=0.2, hspace=0.4)
    for i, (_, row) in tqdm(
        enumerate(plot_df.iterrows()), total=len(plot_df), unit="plot", disable=in_papermill()
    ):
        a_x = fig.add_subplot(grid[i])
        range_columns = list(plot_df.columns[:num_files])
        file_vs_value_plot(a_x, field_name, row, range_columns, fontsize, pad)
    plt.savefig(file_name, bbox_inches="tight")
    plt.close()


def file_vs_value_plot(
    a_x: Axis, field_name: str, row: pd.DataFrame, range_columns: List[str], fontsize: float, pad: float
) -> None:
    """Create a dot plot with one point per file"""
    assert field_name in ["rt_peak", "peak_height"]
    a_x.tick_params(direction="in", length=1, pad=pad, width=0.1, labelsize=fontsize)
    num_files = len(range_columns)
    a_x.scatter(range(num_files), row[:num_files], s=0.2)
    if field_name == "rt_peak":
        a_x.axhline(y=row["atlas RT peak"], color="r", linestyle="-", linewidth=0.2)
        range_columns += ["atlas RT peak"]
        a_x.set_ylim(np.nanmin(row.loc[range_columns]) - 0.12, np.nanmax(row.loc[range_columns]) + 0.12)
    else:
        a_x.set_yscale("log")
        a_x.set_ylim(bottom=1e4, top=1e10)
    a_x.set_xlim(-0.5, num_files + 0.5)
    a_x.xaxis.set_major_locator(mticker.FixedLocator(np.arange(0, num_files, 1.0)))
    _ = [s.set_linewidth(0.1) for s in a_x.spines.values()]
    # truncate name so it fits above a single subplot
    a_x.set_title(row.name[:33], pad=pad, fontsize=fontsize)
    a_x.set_xlabel("Files", labelpad=pad, fontsize=fontsize)
    ylabel = "Actual RTs" if field_name == "rt_peak" else "Peak Height"
    a_x.set_ylabel(ylabel, labelpad=pad, fontsize=fontsize)


def plot_compound_atlas_rts(
    num_files: int,
    rts_df: pd.DataFrame,
    file_name: Path,
    fontsize: float = 2,
    pad: float = 0.1,
    cols: int = 8,
) -> None:
    """Plot filenames vs peak RT for each compound"""
    plot_per_compound("rt_peak", num_files, rts_df, file_name, fontsize, pad, cols)


def plot_compound_atlas_peak_heights(
    num_files: int,
    peak_heights_df: pd.DataFrame,
    file_name: Path,
    fontsize: float = 2,
    pad: float = 0.1,
    cols: int = 8,
) -> None:
    """Plot filenames vs peak height for each compound"""
    plot_per_compound("peak_height", num_files, peak_heights_df, file_name, fontsize, pad, cols)
