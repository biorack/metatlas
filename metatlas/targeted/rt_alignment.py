"""Generate Model for Retention Time Alignment"""
# pylint: disable=too-many-arguments

import logging
import math
import os

from copy import deepcopy
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple, Sequence

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import papermill

from matplotlib import gridspec
from matplotlib.axis import Axis
from sklearn.base import BaseEstimator
from sklearn.linear_model import LinearRegression, RANSACRegressor
from sklearn.preprocessing import PolynomialFeatures
from tqdm.notebook import tqdm

from metatlas.datastructures.analysis_identifiers import AnalysisIdentifiers, MSMS_REFS_PATH
from metatlas.datastructures import metatlas_objects as metob
from metatlas.datastructures.utils import get_atlas
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.io import targeted_output
from metatlas.io import write_utils
from metatlas.io.gdrive import copy_outputs_to_google_drive
from metatlas.plots import dill2plots as dp
from metatlas.tools import notebook
from metatlas.tools import parallel
from metatlas.tools.config import Config, Workflow
from metatlas.tools.util import or_default, repo_path

logger = logging.getLogger(__name__)

# metatlas_dataset type that isn't an instance of the MetatlasDataset class
SimpleMetatlasData = List[List[dict]]


class Model:
    """Encapsulate both linear and polynomial models in a consistent interface"""

    def __init__(self, sk_model: BaseEstimator, intercept: float, coefficents: np.ndarray):
        """
        inputs:
            sk_model: scikit-learn model object
            intercept: y-intercept value
            coefficents: a list of coefficents, with x^n coefficent at index n-1
        """
        self.sk_model = sk_model
        self.intercept = intercept
        if coefficents.shape == (1, 1):
            self.coefficents = [intercept, coefficents[0][0]]
        elif coefficents.shape == (1, 3):
            self.coefficents = coefficents[0].tolist()

    def __repr__(self) -> str:
        """Text description of the model function"""
        if self.order == 1:
            return f"Linear model with intercept={self.intercept:.3f} and slope={self.coefficents[1]:.5f}"
        coef_str = ", ".join([f"{c:.5f}" for c in self.coefficents])
        return f"Polynomial model with intercept={self.intercept:.3f} and coefficents=[{coef_str}]"

    @property
    def order(self) -> int:
        """Polynomial order of the model"""
        return len(self.coefficents) - 1

    @property
    def name(self) -> str:
        """Type of model as string"""
        return "linear" if self.order == 1 else "polynomial"

    def predict(self, x_values) -> List[float]:
        """Returns y values for input x"""
        x_transformed = x_values.reshape(-1, 1)
        if self.order > 1:
            poly_reg = PolynomialFeatures(degree=2)
            x_transformed = poly_reg.fit_transform(x_transformed)
        return self.sk_model.predict(x_transformed).flatten().tolist()


def generate_rt_alignment_models(
    ids: AnalysisIdentifiers,
    metatlas_dataset: SimpleMetatlasData,
    groups: Sequence[metob.Group],
    qc_atlas: metob.Atlas,
    qc_atlas_df: pd.DataFrame,
    workflow: Workflow,
) -> Tuple[Model, Model]:
    """
    Generate the RT correction models and model charaterization files
    inputs:
        ids: an AnalysisIds object matching the one selected_cold in the main notebook
        selected_col: name of column to use for model generation
        inchi_keys_not_in_model: InChi Keys that will be ignored when for model creation
    Returns a tuple with a linear and polynomial model
    """
    params = workflow.rt_alignment.parameters
    rts_df = get_rts(metatlas_dataset)
    actual, pred = subset_data_for_model_input(
        params.dependent_data_source, rts_df, qc_atlas_df, params.inchi_keys_not_in_model
    )
    linear, poly = generate_models(actual, pred)
    out_dir = Path(ids.output_dir)
    actual_rts, aligned_rts = actual_and_aligned_rts(rts_df, qc_atlas_df, params.inchi_keys_not_in_model)
    actual_vs_pred_file_name = out_dir / "Actual_vs_Aligned_RTs.pdf"
    plot_actual_vs_aligned_rts(aligned_rts, actual_rts, rts_df, str(actual_vs_pred_file_name), linear, poly)
    rt_comparison_file_name = out_dir / "RT_Alignment_Model_Comparison.csv"
    save_model_comparison(
        params.dependent_data_source, qc_atlas_df, rts_df, linear, poly, str(rt_comparison_file_name)
    )
    models_file_name = out_dir / "rt_alignment_model.txt"
    write_models(str(models_file_name), linear, poly, groups, qc_atlas)
    return (linear, poly)


def generate_outputs(ids: AnalysisIdentifiers, workflow: Workflow) -> None:
    """
    Generate the RT alignment models, associated atlases with relative RT values, follow up notebooks
    """
    # pylint: disable=too-many-locals
    params = workflow.rt_alignment.parameters
    assert params.stop_before in ["atlases", "notebook_generation", "notebook_execution", None]
    metatlas_dataset, groups, atlas, atlas_df = load_data(ids, workflow)
    linear, poly = generate_rt_alignment_models(ids, metatlas_dataset, groups, atlas, atlas_df, workflow)
    if params.stop_before in ["notebook_generation", "notebook_execution", None]:
        atlases = create_aligned_atlases(linear, poly, ids, workflow)
    if params.stop_before in ["notebook_execution", None]:
        notebook_file_names = write_notebooks(ids, atlases, workflow)
    if params.stop_before is None:
        for in_file_name in notebook_file_names:
            out_file_name = in_file_name.with_name(in_file_name.stem + "_SLURM.ipynb")
            papermill.execute_notebook(in_file_name, out_file_name, {}, kernel_name="papermill")
    copy_outputs_to_google_drive(ids)
    targeted_output.archive_outputs(ids)
    logger.info("RT_Alignment notebook complete. Switch to an analysis notebook to continue.")


def generate_qc_plots(metatlas_dataset: SimpleMetatlasData, ids: AnalysisIdentifiers) -> None:
    """Write plots that can be used to QC the experiment"""
    rts_df = get_rts(metatlas_dataset)
    compound_atlas_rts_file_name = os.path.join(
        ids.output_dir, f"{ids.short_polarity}_Compound_Atlas_RTs.pdf"
    )
    plot_compound_atlas_rts(len(metatlas_dataset), rts_df, compound_atlas_rts_file_name)
    peak_heights_df = get_peak_heights(metatlas_dataset)
    peak_heights_plot_file_name = os.path.join(
        ids.output_dir, f"{ids.short_polarity}_Compound_Atlas_peak_heights.pdf"
    )
    plot_compound_atlas_peak_heights(len(metatlas_dataset), peak_heights_df, peak_heights_plot_file_name)


def generate_qc_outputs(metatlas_dataset: SimpleMetatlasData, ids: AnalysisIdentifiers, cpus: int) -> None:
    """Write outputs that can be used to QC the experiment"""
    save_rt_peak(metatlas_dataset, os.path.join(ids.output_dir, f"{ids.short_polarity}_rt_peak.tab"))
    save_measured_rts(
        metatlas_dataset, os.path.join(ids.output_dir, f"{ids.short_polarity}_QC_Measured_RTs.csv")
    )
    generate_qc_plots(metatlas_dataset, ids)
    write_chromatograms(metatlas_dataset, ids.output_dir, max_cpus=cpus)
    hits = dp.get_msms_hits(metatlas_dataset, extra_time=0.2, ref_loc=MSMS_REFS_PATH)
    write_identification_figures(
        metatlas_dataset, hits, ids.output_dir, ids.lcmsruns_short_names, prefix=ids.short_polarity
    )


def load_data(
    ids: AnalysisIdentifiers,
    workflow: Workflow,
) -> Tuple[SimpleMetatlasData, List[metob.Group], metob.Atlas, pd.DataFrame]:
    """create metatlas_dataset, groups and atlas"""
    params = workflow.rt_alignment.parameters
    groups = get_groups(ids)
    files_df = get_files_df(groups)
    qc_atlas, qc_atlas_df = get_qc_atlas(workflow)
    # this metatlas_dataset is not a class instance. Only has metatlas_dataset[file_idx][compound_idx]...
    metatlas_dataset = load_runs(files_df, qc_atlas_df, qc_atlas, params.max_cpus)
    try:
        if len(metatlas_dataset) == 0:
            raise ValueError("No matching LCMS runs, terminating without generating outputs.")
    except ValueError as err:
        logger.exception(err)
        raise err
    return metatlas_dataset, groups, qc_atlas, qc_atlas_df


def write_identification_figures(
    data: SimpleMetatlasData,
    hits: pd.DataFrame,
    output_dir: str,
    run_short_names: pd.DataFrame,
    overwrite: bool = False,
    prefix: str = "",
) -> None:
    """
    inputs:
       data: a metatlas_dataset datastructures (does not need to be MetatlasDataset class)
       hits: msms hits
       output_dir: directory to write pdfs to
       run_short_names: short names for LCMS runs in a dataframe
       overwrite: if False, throw error if files already exist
    """
    dp.make_identification_figure_v2(
        input_dataset=data,
        msms_hits=hits,
        use_labels=True,
        include_lcmsruns=["QC"],
        exclude_lcmsruns=[],
        output_loc=output_dir,
        short_names_df=run_short_names,
        polarity=prefix,
        overwrite=overwrite,
    )


def get_groups(ids: AnalysisIdentifiers) -> List[metob.Group]:
    """
    Create all experiment groups if they don't already exist and return the subset matching include_list
    inputs:
        ids: instance of AnalysisIds
    """
    ordered_groups = sorted(ids.groups, key=lambda x: x.name)
    for grp in ordered_groups:
        logger.info("Selected group: %s, %s", grp.name, int_to_date_str(grp.last_modified))
    return ordered_groups


def int_to_date_str(i_time: int) -> str:
    """unix epoc time in seconds to YYYY-MM-DD hh:mm:ss"""
    return str(datetime.fromtimestamp(i_time))


def get_files_df(groups: Sequence[metob.Group]) -> pd.DataFrame:
    """Pandas Datafram with one row per file plus columns for accquistion_time and group name"""
    files_df = pd.DataFrame(columns=["file", "time", "group"])
    for group in groups:
        for lcms_run in group.items:
            try:
                time = lcms_run.accquistion_time
            except AttributeError:
                time = 0
            files_df = files_df.append({"file": lcms_run, "time": time, "group": group}, ignore_index=True)
    return files_df.sort_values(by=["time"])


def get_qc_atlas(workflow: Workflow) -> Tuple[metob.Atlas, pd.DataFrame]:
    """Retreives template QC atlas and return tuple (atlas, atlas_df)"""
    atlas = workflow.rt_alignment.atlas
    params = workflow.rt_alignment.parameters
    logger.info("Loading QC Atlas %s", atlas.name)
    original_atlas = metob.retrieve("Atlas", name=atlas.name, username=atlas.username)[0]
    atlas = adjust_atlas_rt_range(original_atlas, params.rt_min_delta, params.rt_max_delta)
    atlas_df = ma_data.make_atlas_df(atlas)
    atlas_df["label"] = [cid.name for cid in atlas.compound_identifications]
    return atlas, atlas_df


def adjust_atlas_rt_range(
    in_atlas: metob.Atlas, rt_min_delta: Optional[float], rt_max_delta: Optional[float]
) -> metob.Atlas:
    """Reset the rt_min and rt_max values by adding rt_min_delta or rt_max_delta to rt_peak"""
    if rt_min_delta is None and rt_max_delta is None:
        return in_atlas
    out_atlas = deepcopy(in_atlas)
    for cid in out_atlas.compound_identifications:
        rts = cid.rt_references[0]
        rts.rt_min = rts.rt_min if rt_min_delta is None else rts.rt_peak + rt_min_delta
        rts.rt_max = rts.rt_max if rt_max_delta is None else rts.rt_peak + rt_max_delta
    return out_atlas


def load_runs(
    files_df: pd.DataFrame, qc_atlas_df: pd.DataFrame, qc_atlas: metob.Atlas, cpus: int
) -> SimpleMetatlasData:
    """
    Loads MSMS data file files
    inputs:
        files_df: files to load
        qc_atlas_df: dataframe form of the QC atlas
        qc_atlas: atlas of QC compounds
        cpus: number of cpus to use
    """
    files = [(i[1].file, i[1].group, qc_atlas_df, qc_atlas) for i in files_df.iterrows()]
    logger.info("Loading LCMS data files")
    return parallel.parallel_process(ma_data.get_data_for_atlas_df_and_file, files, cpus, unit="sample")


def save_measured_rts(metatlas_dataset: SimpleMetatlasData, file_name: str) -> None:
    """Save RT values in csv format file"""
    rts_df = get_rts(metatlas_dataset, include_atlas_rt_peak=False)
    write_utils.export_dataframe_die_on_diff(rts_df, file_name, "measured RT values", float_format="%.6e")


def save_rt_peak(metatlas_dataset: SimpleMetatlasData, file_name: str) -> None:
    """Save peak RT values in tsv format file"""
    rts_df = dp.make_output_dataframe(input_dataset=metatlas_dataset, fieldname="rt_peak", use_labels=True)
    write_utils.export_dataframe_die_on_diff(
        rts_df, file_name, "peak RT values", sep="\t", float_format="%.6e"
    )


def get_rts(metatlas_dataset: SimpleMetatlasData, include_atlas_rt_peak: bool = True) -> pd.DataFrame:
    """Returns RT values in DataFrame format"""
    rts_df = dp.make_output_dataframe(
        input_dataset=metatlas_dataset,
        fieldname="rt_peak",
        use_labels=True,
        summarize=True,
    )
    if include_atlas_rt_peak:
        rts_df["atlas RT peak"] = [
            compound["identification"].rt_references[0].rt_peak for compound in metatlas_dataset[0]
        ]
    return order_df_columns_by_run(rts_df)


def get_peak_heights(metatlas_dataset: SimpleMetatlasData) -> pd.DataFrame:
    """Returns peak heights in DataFrame format"""
    peak_height_df = dp.make_output_dataframe(
        input_dataset=metatlas_dataset,
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
    file_name: str,
    fontsize: float = 2,
    pad: float = 0.1,
    cols: int = 8,
) -> None:
    """
    Writes plot of RT peak for vs file for each compound
    inputs:
        field_name: one of rt_peak or peak_height
        num_files: number of files in data set, ie len(metatlas_dataset)
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
    for i, (_, row) in tqdm(enumerate(plot_df.iterrows()), total=len(plot_df), unit="plot"):
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
    num_files: int, rts_df: pd.DataFrame, file_name: str, fontsize: float = 2, pad: float = 0.1, cols: int = 8
) -> None:
    """Plot filenames vs peak RT for each compound"""
    plot_per_compound("rt_peak", num_files, rts_df, file_name, fontsize, pad, cols)


def plot_compound_atlas_peak_heights(
    num_files: int,
    peak_heights_df: pd.DataFrame,
    file_name: str,
    fontsize: float = 2,
    pad: float = 0.1,
    cols: int = 8,
) -> None:
    """Plot filenames vs peak height for each compound"""
    plot_per_compound("peak_height", num_files, peak_heights_df, file_name, fontsize, pad, cols)


def generate_models(actual: List[float], pred: List[float]) -> Tuple[Model, Model]:
    """
    inputs:
        actual: experimental RTs
        pred: predicted RTs
    returns tuple containing two Model classes of order 1 and 2
    """
    transformed_actual = np.array(actual).reshape(-1, 1)
    transformed_pred = np.array(pred).reshape(-1, 1)

    ransac = RANSACRegressor(random_state=42)
    rt_model_linear = ransac.fit(transformed_pred, transformed_actual)
    linear = Model(
        rt_model_linear, rt_model_linear.estimator_.intercept_[0], rt_model_linear.estimator_.coef_
    )

    poly_reg = PolynomialFeatures(degree=2)
    x_poly = poly_reg.fit_transform(transformed_pred)
    rt_model_poly = LinearRegression().fit(x_poly, transformed_actual)
    poly = Model(rt_model_poly, rt_model_poly.intercept_[0], rt_model_poly.coef_)
    return linear, poly


def subset_data_for_model_input(
    selected_column: str,
    rts_df: pd.DataFrame,
    atlas_df: pd.DataFrame,
    inchi_keys_not_in_model: Optional[List[str]] = None,
) -> Tuple[List[float], List[float]]:
    """
    inputs:
        selected_column: column number in rts_df to use for actual values
        rts_df: dataframe of RT values
        atlas_df: QC atlas in dataframe format
        inchi_keys_not_in_model: InChi Keys that will be ignored for model creation
    return a tuple of (actual, pred)
    """
    keep_idxs = get_keep_idxs(selected_column, rts_df, atlas_df, inchi_keys_not_in_model)
    actual = rts_df.iloc[keep_idxs][selected_column].tolist()
    pred = atlas_df.iloc[keep_idxs]["rt_peak"].tolist()
    return actual, pred


def get_keep_idxs(
    selected_column: str,
    rts_df: pd.DataFrame,
    atlas_df: pd.DataFrame,
    inchi_keys_not_in_model: Optional[List[str]] = None,
) -> List[int]:
    """Indices in rts_df that should be used within the model"""
    keep_idxs = set(np.flatnonzero(~np.isnan(rts_df.loc[:, selected_column])))
    if inchi_keys_not_in_model:
        keep_idxs = keep_idxs.intersection(
            set(np.flatnonzero(~atlas_df["inchi_key"].isin(inchi_keys_not_in_model)))
        )
    return list(keep_idxs)


def actual_and_aligned_rts(
    rts_df: pd.DataFrame, atlas_df: pd.DataFrame, inchi_keys_not_in_model: Optional[List[str]] = None
) -> Tuple[List[List[float]], List[List[float]]]:
    """
    inputs:
        rts_df: dataframe of RT values
        atlas_df: QC atlas in dataframe format
        inchi_keys_not_in_model: InChi Keys that will be ignored for model creation
    return a tuple of lists of lists: (actual_rts, aligned_rts)
    """
    actual_rts = []
    aligned_rts = []
    for i in range(rts_df.shape[1] - 5):
        keep_idxs = get_keep_idxs(rts_df.columns[i], rts_df, atlas_df, inchi_keys_not_in_model)
        current_actual_df = rts_df.loc[:, rts_df.columns[i]]
        current_actual_df = current_actual_df.iloc[keep_idxs]
        current_pred_df = atlas_df.iloc[keep_idxs][["rt_peak"]]
        actual_rts.append(current_actual_df.values.tolist())
        aligned_rts.append(current_pred_df.values.tolist())
    return actual_rts, aligned_rts


def plot_actual_vs_aligned_rts(
    aligned_rts: Sequence[Sequence[float]],
    actual_rts: Sequence[Sequence[float]],
    rts_df: pd.DataFrame,
    file_name: str,
    linear: Model,
    poly: Model,
) -> None:
    """Write scatter plot showing linear vs polynomial fit"""
    # pylint: disable=too-many-locals
    rows = int(math.ceil((rts_df.shape[1] + 1) / 5))
    cols = 5
    fig = plt.figure(constrained_layout=False)
    grid = gridspec.GridSpec(rows, cols, figure=fig)
    plt.rc("font", size=6)
    plt.rc("axes", labelsize=6)
    plt.rc("xtick", labelsize=3)
    plt.rc("ytick", labelsize=3)
    for i in range(rts_df.shape[1] - 5):
        x_values = aligned_rts[i]
        y_values = actual_rts[i]
        if len(x_values) == 0 or len(y_values) == 0:
            continue
        sub = fig.add_subplot(grid[i])
        sub.scatter(x_values, y_values, s=2)
        spaced_x = np.linspace(0, max(x_values), 100)
        sub.plot(spaced_x, linear.predict(spaced_x), linewidth=0.5, color="red")
        sub.plot(spaced_x, poly.predict(spaced_x), linewidth=0.5, color="green")
        sub.set_title("File: " + str(i))
        sub.set_xlabel("relative RTs")
        sub.set_ylabel("actual RTs")
    fig_legend = [
        (
            "Red line: linear model;  Green curve: polynomial model. "
            "Default model is a polynomial model using the median data."
        ),
        "",
        "file_index  data_source",
    ] + [f"{i:2d}              {rts_df.columns[i]}" for i in range(rts_df.shape[1] - 5)]
    fig.tight_layout(pad=0.5)
    line_height = 0.03
    legend_offset = line_height * len(fig_legend)
    plt.text(0, -1 * legend_offset, "\n".join(fig_legend), transform=plt.gcf().transFigure)
    plt.savefig(file_name, bbox_inches="tight")


def save_model_comparison(
    selected_column: str,
    qc_atlas_df: pd.DataFrame,
    rts_df: pd.DataFrame,
    linear: Model,
    poly: Model,
    file_name: str,
) -> None:
    """
    Save csv format file with per-compound comparision of linear vs polynomial models
    inputs:
        selected_column: column number in rts_df to use for actual values
        qc_atlas_df: QC atlas in dataframe format
        rts_df: dataframe with RT values
        linear: instance of class Model with first order model
        poly: instance of class Model with second order model
        file_name: where to save the plot
    """
    qc_df = rts_df[[selected_column]].copy()
    qc_df.columns = ["RT Measured"]
    qc_df.loc[:, "RT Reference"] = qc_atlas_df["rt_peak"].to_numpy()
    qc_df.loc[:, "Relative RT Linear"] = pd.Series(
        linear.predict(qc_df["RT Reference"].to_numpy()), index=qc_df.index
    )
    qc_df.loc[:, "Relative RT Polynomial"] = pd.Series(
        poly.predict(qc_df["RT Reference"].to_numpy()), index=qc_df.index
    )
    qc_df["RT Diff Linear"] = qc_df["RT Measured"] - qc_df["Relative RT Linear"]
    qc_df["RT Diff Polynomial"] = qc_df["RT Measured"] - qc_df["Relative RT Polynomial"]
    write_utils.export_dataframe_die_on_diff(qc_df, file_name, "model comparision", float_format="%.6e")


def write_models(
    file_name: str, linear_model: Model, poly_model: Model, groups: Sequence[metob.Group], atlas: metob.Atlas
) -> None:
    """
    inputs:
        file_name: text file to save model information
        linear_model: instance of class Model with first order model
        poly_model: instance of class Model with second order model
        groups: list of groups used in model generation
        atlas: QC atlas
    """
    with open(file_name, "w", encoding="utf8") as out_fh:
        for model in [linear_model, poly_model]:
            out_fh.write(f"{model.sk_model.set_params()}\n")
            out_fh.write(f"{model}\n")
            group_names = ", ".join([g.name for g in groups])
            out_fh.write(f"groups = {group_names}\n")
            out_fh.write(f"atlas = {atlas.name}\n\n")


def get_atlas_name(template_name: str, model: Model, project_id: str, analysis_name: str) -> str:
    """
    input:
        template_name: name of template atlas
        ids: an AnalysisIds object matching the one used in the main notebook
        model: an instance of Model
    returns the name of the production atlas
    """
    prod_name = template_name.replace("TPL", "PRD")
    return f"{prod_name}_{model.name}_{project_id}_{analysis_name}"


def adjust_atlas(atlas: metob.Atlas, model: Model, ids: AnalysisIdentifiers) -> pd.DataFrame:
    """use model to adjust RTs within atlas"""
    atlas_df = ma_data.make_atlas_df(atlas)
    atlas_df["label"] = [cid.name for cid in atlas.compound_identifications]
    atlas_df["rt_peak"] = model.predict(atlas_df["rt_peak"].to_numpy())
    rt_offset = 0.2 if ids.chromatography == "C18" else 0.5
    atlas_df["rt_min"] = atlas_df["rt_peak"].apply(lambda rt: rt - rt_offset)
    atlas_df["rt_max"] = atlas_df["rt_peak"].apply(lambda rt: rt + rt_offset)
    return atlas_df


def create_aligned_atlases(
    linear: Model,
    poly: Model,
    ids: AnalysisIdentifiers,
    workflow: Workflow,
) -> List[str]:
    """
    input:
        linear_model: instance of class Model with first order model
        poly_model: instance of class Model with second order model
        ids: an AnalysisIdentifiers object
        workflow: a config Workflow object
    returns a list of the names of atlases
    """
    # pylint: disable=too-many-locals
    out_atlas_names = []
    model = poly if workflow.rt_alignment.parameters.use_poly_model else linear
    for analysis in tqdm(workflow.analyses, unit="atlas"):
        template_atlas = get_atlas(analysis.atlas.name, analysis.atlas.username)
        out_atlas_names.append(get_atlas_name(template_atlas.name, model, ids.project, analysis.name))
        logger.info("Creating atlas %s", out_atlas_names[-1])
        out_atlas_file_name = os.path.join(ids.output_dir, f"{out_atlas_names[-1]}.csv")
        out_atlas_df = adjust_atlas(template_atlas, model, ids)
        write_utils.export_dataframe_die_on_diff(
            out_atlas_df, out_atlas_file_name, "RT aligned atlas", index=False, float_format="%.6e"
        )
        dp.make_atlas_from_spreadsheet(
            out_atlas_df,
            out_atlas_names[-1],
            filetype="dataframe",
            sheetname="",
            polarity=analysis.parameters.polarity,
            store=True,
            mz_tolerance=10 if ids.chromatography == "C18" else 12,
        )
    return out_atlas_names


def write_notebooks(ids: AnalysisIdentifiers, atlases: Sequence[str], workflow: Workflow) -> List[Path]:
    """
    Creates Targeted analysis jupyter notebooks with pre-populated parameter sets
    Inputs:
        ids: an AnalysisIds object matching the one used in the main notebook
        workflow: a Workflow object
        atlases: list of atlas names to use as source atlases
    Returns a list of Paths to notebooks
    """
    out = []
    for atlas_name, analysis in zip(atlases, workflow.analyses):
        source = repo_path() / "notebooks" / "reference" / "Targeted.ipynb"
        dest = Path(ids.output_dir).resolve().parent / f"{ids.project}_{workflow.name}_{analysis.name}.ipynb"
        parameters = {
            "experiment": ids.experiment,
            "rt_alignment_number": ids.rt_alignment_number,
            "analysis_number": 0,
            "workflow_name": workflow.name,
            "analysis_name": analysis.name,
            "source_atlas": atlas_name,
            "copy_atlas": analysis.parameters.copy_atlas,
            "polarity": analysis.parameters.polarity,
            "include_groups": analysis.parameters.include_groups,
            "exclude_groups": analysis.parameters.exclude_groups,
            "groups_controlled_vocab": analysis.parameters.groups_controlled_vocab,
            "exclude_files": analysis.parameters.exclude_files,
            "generate_qc_outputs": analysis.parameters.generate_qc_outputs,
            "num_points": analysis.parameters.num_points,
            "peak_height": analysis.parameters.peak_height,
            "msms_score": analysis.parameters.msms_score,
            "filter_removed": analysis.parameters.filter_removed,
            "line_colors": analysis.parameters.line_colors,
            "require_all_evaluated": analysis.parameters.require_all_evaluated,
            "generate_post_annotation_outputs": analysis.parameters.generate_post_annotation_outputs,
            "exclude_groups_for_post_annotation_outputs": analysis.parameters.exclude_groups_for_post_annotation_outputs,
            "export_msms_fragment_ions": analysis.parameters.export_msms_fragment_ions,
            "clear_cache": analysis.parameters.clear_cache,
            "config_file_name": analysis.parameters.config_file_name,
            "source_code_version_id": analysis.parameters.source_code_version_id,
            "project_directory": or_default(analysis.parameters.project_directory, ids.project_directory),
            "google_folder": or_default(
                analysis.parameters.google_folder, workflow.rt_alignment.parameters.google_folder
            ),
            "max_cpus": analysis.parameters.max_cpus,
            "log_level": analysis.parameters.log_level,
        }
        notebook.create_notebook(source, dest, parameters)
        out.append(dest)
    return out


def write_chromatograms(
    metatlas_dataset: SimpleMetatlasData, output_dir: str, overwrite: bool = False, max_cpus: int = 1
) -> None:
    """
    inputs:
        metatlas_dataset: a metatlas_dataset datastructure
        output_dir: directory to save plots within
        overwrite: if False raise error if file already exists
        max_cpus: number of cpus to use
    """
    # overwrite checks done within dp.make_chromatograms
    logger.info("Exporting chromotograms to %s", output_dir)
    params = {
        "input_dataset": metatlas_dataset,
        "share_y": True,
        "output_loc": output_dir,
        "overwrite": overwrite,
        "max_cpus": max_cpus,
        "include_lcmsruns": ["QC"],
        "suffix": "_sharedY",
    }
    dp.make_chromatograms(**params)
    params["share_y"] = False
    params["suffix"] = "_independentY"
    dp.make_chromatograms(**params)


def run(
    experiment: str,
    rt_alignment_number: int,
    configuration: Config,
    workflow: Workflow,
) -> None:
    """Generates RT alignment model, applies to atlases, and generates all outputs"""
    # pylint: disable=too-many-locals
    params = workflow.rt_alignment.parameters
    ids = AnalysisIdentifiers(
        analysis_number=0,
        source_atlas=workflow.rt_alignment.atlas.name,
        source_atlas_username=workflow.rt_alignment.atlas.username,
        experiment=experiment,
        project_directory=params.project_directory,
        google_folder=params.google_folder,
        rt_alignment_number=rt_alignment_number,
        exclude_files=params.exclude_files,
        include_groups=params.include_groups,
        exclude_groups=params.exclude_groups,
        groups_controlled_vocab=params.groups_controlled_vocab,
        configuration=configuration,
        workflow=workflow.name,
    )
    generate_outputs(ids=ids, workflow=workflow)
