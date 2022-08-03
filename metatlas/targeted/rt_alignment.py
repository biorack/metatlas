"""Generate Model for Retention Time Alignment"""
# pylint: disable=too-many-arguments

import logging
import math
import os

from pathlib import Path
from typing import List, Optional, Tuple, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import papermill  # pylint: disable=wrong-import-order

from matplotlib import gridspec
from sklearn.base import BaseEstimator
from sklearn.linear_model import LinearRegression, RANSACRegressor
from sklearn.preprocessing import PolynomialFeatures
from tqdm.notebook import tqdm

from metatlas.datastructures.analysis_identifiers import AnalysisIdentifiers
from metatlas.datastructures.metatlas_dataset import MetatlasDataset
from metatlas.datastructures import metatlas_objects as metob
from metatlas.datastructures.utils import get_atlas
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.io import targeted_output
from metatlas.io import write_utils
from metatlas.io.gdrive import copy_outputs_to_google_drive
from metatlas.plots import dill2plots as dp
from metatlas.tools import notebook
from metatlas.tools.config import Config, Workflow, Analysis
from metatlas.tools.util import or_default, repo_path

logger = logging.getLogger(__name__)


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
    data: MetatlasDataset,
    workflow: Workflow,
) -> Tuple[Model, Model]:
    """
    Generate the RT correction models and model charaterization files
    Returns a tuple with a linear and polynomial model
    """
    params = workflow.rt_alignment.parameters
    rts_df = get_rts(data)
    actual, pred = subset_data_for_model_input(
        params.dependent_data_source, rts_df, data.atlas_df, params.inchi_keys_not_in_model
    )
    linear, poly = generate_models(actual, pred)
    out_dir = Path(data.ids.output_dir)
    actual_rts, aligned_rts = actual_and_aligned_rts(rts_df, data.atlas_df, params.inchi_keys_not_in_model)
    actual_vs_pred_file_name = out_dir / "Actual_vs_Aligned_RTs.pdf"
    plot_actual_vs_aligned_rts(aligned_rts, actual_rts, rts_df, str(actual_vs_pred_file_name), linear, poly)
    rt_comparison_file_name = out_dir / "RT_Alignment_Model_Comparison.csv"
    save_model_comparison(
        params.dependent_data_source, data.atlas_df, rts_df, linear, poly, str(rt_comparison_file_name)
    )
    models_file_name = out_dir / "rt_alignment_model.txt"
    write_models(str(models_file_name), linear, poly, data.ids.groups, data.atlas)
    return (linear, poly)


def generate_outputs(data: MetatlasDataset, workflow: Workflow, analysis: Analysis) -> None:
    """
    Generate the RT alignment models, associated atlases with relative RT values, follow up notebooks
    """
    # pylint: disable=too-many-locals
    params = workflow.rt_alignment.parameters
    ids = data.ids
    assert params.stop_before in ["atlases", "notebook_generation", "notebook_execution", None]
    linear, poly = generate_rt_alignment_models(data, workflow)
    if params.stop_before in ["notebook_generation", "notebook_execution", None]:
        atlases = create_aligned_atlases(linear, poly, ids, workflow)
    if params.stop_before in ["notebook_execution", None]:
        notebook_file_names = write_notebooks(ids, atlases, workflow)
    if params.stop_before is None:
        for in_file_name in notebook_file_names:
            out_file_name = in_file_name.with_name(in_file_name.stem + "_SLURM.ipynb")
            papermill.execute_notebook(in_file_name, out_file_name, {}, kernel_name="papermill")
    copy_outputs_to_google_drive(ids)
    logger.info("RT_Alignment notebook complete. Switch to an analysis notebook to continue.")


def get_rts(metatlas_dataset: MetatlasDataset, include_atlas_rt_peak: bool = True) -> pd.DataFrame:
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
    return targeted_output.order_df_columns_by_run(rts_df)


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


def get_atlas_name(ids: AnalysisIdentifiers, workflow: Workflow, analysis: Analysis, model: Model) -> str:
    """returns the name of the RT aligned atlas"""
    prod_name = analysis.atlas.name.replace("TPL", "PRD")
    return (
        f"{prod_name}_{model.name}_{ids.experiment_id}_"
        f"{workflow.name}_{analysis.name}_{ids.rt_alignment_number}"
    )


def align_atlas(atlas: metob.Atlas, model: Model, ids: AnalysisIdentifiers) -> pd.DataFrame:
    """use model to align RTs within atlas"""
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
) -> List[metob.Atlas]:
    """
    input:
        linear_model: instance of class Model with first order model
        poly_model: instance of class Model with second order model
        ids: an AnalysisIdentifiers object
        workflow: a config Workflow object
    returns a list of the names of atlases
    """
    # pylint: disable=too-many-locals
    out_atlases = []
    model = poly if workflow.rt_alignment.parameters.use_poly_model else linear
    for analysis in tqdm(workflow.analyses, unit="atlas"):
        template_atlas = get_atlas(analysis.atlas.unique_id, analysis.atlas.name)
        if analysis.atlas.do_alignment:
            name = get_atlas_name(ids, workflow, analysis, model)
            logger.info("Creating atlas %s", name)
            out_atlas_file_name = os.path.join(ids.output_dir, f"{name}.csv")
            out_atlas_df = align_atlas(template_atlas, model, ids)
            write_utils.export_dataframe_die_on_diff(
                out_atlas_df, out_atlas_file_name, "RT aligned atlas", index=False, float_format="%.6e"
            )
            atlas = dp.make_atlas_from_spreadsheet(
                out_atlas_df,
                name,
                filetype="dataframe",
                sheetname="",
                polarity=analysis.parameters.polarity,
                store=True,
                mz_tolerance=10 if ids.chromatography == "C18" else 12,
            )
            out_atlases.append(atlas)
        else:
            out_atlases.append(template_atlas)
    return out_atlases


def write_notebooks(
    ids: AnalysisIdentifiers, atlases: Sequence[metob.Atlas], workflow: Workflow
) -> List[Path]:
    """
    Creates Targeted analysis jupyter notebooks with pre-populated parameter sets
    Inputs:
        ids: an AnalysisIds object matching the one used in the main notebook
        workflow: a Workflow object
        atlases: list of atlases to use as source atlases
    Returns a list of Paths to notebooks
    """
    out = []
    for atlas, analysis in zip(atlases, workflow.analyses):
        source = repo_path() / "notebooks" / "reference" / "Targeted.ipynb"
        dest = (
            Path(ids.output_dir).resolve().parent
            / f"{ids.experiment_id}_{workflow.name}_{analysis.name}.ipynb"
        )
        parameters = {
            "experiment": ids.experiment,
            "rt_alignment_number": ids.rt_alignment_number,
            "analysis_number": 0,
            "workflow_name": workflow.name,
            "analysis_name": analysis.name,
            "source_atlas_unique_id": atlas.unique_id,
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
            "generate_analysis_outputs": analysis.parameters.generate_analysis_outputs,
            "exclude_groups_for_analysis_outputs": analysis.parameters.exclude_groups_for_analysis_outputs,
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


def run(
    experiment: str,
    rt_alignment_number: int,
    configuration: Config,
    workflow: Workflow,
    analysis: Analysis,
) -> MetatlasDataset:
    """Generates RT alignment model, applies to atlases, and generates all outputs"""
    params = workflow.rt_alignment.parameters
    ids = AnalysisIdentifiers(
        analysis_number=0,
        source_atlas_unique_id=workflow.rt_alignment.atlas.unique_id,
        copy_atlas=params.copy_atlas,
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
    metatlas_dataset = MetatlasDataset(ids=ids, max_cpus=params.max_cpus)
    generate_outputs(metatlas_dataset, workflow, analysis)
    return metatlas_dataset
