"""Generate Retention Time Correction Model"""
# pylint: disable=too-many-arguments

import itertools
import logging
import math
import os

from datetime import datetime
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

from matplotlib import gridspec
from sklearn.linear_model import LinearRegression, RANSACRegressor
from sklearn.preprocessing import PolynomialFeatures
from tqdm.notebook import tqdm

from metatlas.datastructures import metatlas_dataset as mads
from metatlas.datastructures import metatlas_objects as metob
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.io import targeted_output
from metatlas.io import write_utils
from metatlas.plots import dill2plots as dp
from metatlas.tools import notebook
from metatlas.tools import parallel

logger = logging.getLogger(__name__)

TEMPLATES = {
    "positive": {
        "HILICZ": [
            {"name": "HILICz150_ANT20190824_TPL_EMA_Unlab_POS", "username": "vrsingan"},
            {"name": "HILICz150_ANT20190824_TPL_QCv3_Unlab_POS", "username": "vrsingan"},
            {"name": "HILICz150_ANT20190824_TPL_ISv5_Unlab_POS", "username": "vrsingan"},
            {"name": "HILICz150_ANT20190824_TPL_ISv5_13C15N_POS", "username": "vrsingan"},
            {"name": "HILICz150_ANT20190824_TPL_IS_LabUnlab2_POS", "username": "vrsingan"},
        ],
        "C18": [
            {"name": "C18_20220208c_QC_POS", "username": "wjholtz"},
            {"name": "C18_20220118_TPL_POS", "username": "wjholtz"},
        ],
    },
    "negative": {
        "HILICZ": [
            {"name": "HILICz150_ANT20190824_TPL_EMA_Unlab_NEG", "username": "vrsingan"},
            {"name": "HILICz150_ANT20190824_TPL_QCv3_Unlab_NEG", "username": "vrsingan"},
            {"name": "HILICz150_ANT20190824_TPL_ISv5_Unlab_NEG", "username": "vrsingan"},
            {"name": "HILICz150_ANT20190824_TPL_ISv5_13C15N_NEG", "username": "vrsingan"},
            {"name": "HILICz150_ANT20190824_TPL_IS_LabUnlab2_NEG", "username": "vrsingan"},
        ],
        "C18": [
            {"name": "C18_20220208c_QC_NEG", "username": "wjholtz"},
            {"name": "C18_20220118_TPL_NEG", "username": "wjholtz"},
        ],
    },
}

QC_ATLASES = {
    "positive": {
        "HILICZ": {"name": "HILICz150_ANT20190824_TPL_QCv3_Unlab_POS", "username": "vrsingan"},
        "C18": {"name": "C18_20220208c_QC_POS", "username": "wjholtz"},
    },
    "negative": {
        "HILICZ": {"name": "HILICz150_ANT20190824_TPL_QCv3_Unlab_NEG", "username": "vrsingan"},
        "C18": {"name": "C18_20220208c_QC_NEG", "username": "wjholtz"},
    },
}


class Model:
    """Encapsulate both linear and polynomial models in a consistent interface"""

    def __init__(self, sk_model, intercept, coefficents):
        """
        inputs:
            sk_model: scikit-learn model object
            intercept: y-intercept value
            coefficents: a list of coefficents, with x^n coefficent at index n-1
        """
        self.sk_model = sk_model
        self.intercept = intercept
        if isinstance(coefficents, (list, np.ndarray)):
            self.coefficents = coefficents
        else:
            self.coefficents = [coefficents]

    def __repr__(self):
        """Text description of the model function"""
        if self.order == 1:
            return f"Linear model with intercept={self.intercept:.3f} and slope={self.coefficents[0]:.5f}"
        coef_str = ", ".join([f"{c:.5f}" for c in self.coefficents])
        return f"Polynomial model with intercept={self.intercept:.3f} and coefficents=[{coef_str}]"

    @property
    def order(self):
        """Polynomial order of the model"""
        return 1 if len(self.coefficents) == 1 else len(self.coefficents) - 1

    @property
    def name(self):
        """Type of model as string"""
        return "linear" if self.order == 1 else "polynomial"

    def predict(self, x_values):
        """Returns y values for input x"""
        x_transformed = x_values.reshape(-1, 1)
        if self.order > 1:
            x_transformed = np.array([[i[0] ** n for n in range(self.order + 1)] for i in x_transformed])
        return self.sk_model.predict(x_transformed)


def generate_rt_correction_models(
    ids: mads.AnalysisIdentifiers, cpus: int, selected_col
) -> Tuple[Model, Model]:
    """
    Generate the RT correction models and model charaterization files
    inputs:
        ids: an AnalysisIds object matching the one selected_cold in the main notebook
        cpus: max number of cpus to us
        selected_col: name of column to use for model generation
    Returns a tuple with a linear and polynomial model
    """
    # pylint: disable=too-many-locals
    groups = get_groups(ids)
    files_df = get_files_df(groups)
    qc_atlas, qc_atlas_df = get_qc_atlas(ids)
    # this metatlas_dataset is not a class instance. Only has metatlas_dataset[file_idx][compound_idx]...
    metatlas_dataset = load_runs(files_df, qc_atlas_df, qc_atlas, cpus)
    try:
        if len(metatlas_dataset) == 0:
            raise ValueError("No matching LCMS runs, terminating without generating outputs.")
    except ValueError as err:
        logger.exception(err)
        raise err
    save_rt_peak(metatlas_dataset, os.path.join(ids.output_dir, "rt_peak.tab"))
    save_measured_rts(metatlas_dataset, os.path.join(ids.output_dir, "QC_Measured_RTs.csv"))
    rts_df = get_rts(metatlas_dataset)
    compound_atlas_rts_file_name = os.path.join(ids.output_dir, "Compound_Atlas_RTs.pdf")
    plot_compound_atlas_rts(len(metatlas_dataset), rts_df, compound_atlas_rts_file_name)
    actual_df, pred_df = actual_and_predicted_df(selected_col, rts_df, qc_atlas_df)
    linear, poly = generate_models(actual_df, pred_df)
    actual_rts, pred_rts = actual_and_predicted_rts(rts_df, qc_atlas_df, actual_df, pred_df)
    actual_vs_pred_file_name = os.path.join(ids.output_dir, "Actual_vs_Predicted_RTs.pdf")
    plot_actual_vs_pred_rts(pred_rts, actual_rts, rts_df, actual_vs_pred_file_name, linear, poly)
    rt_comparison_file_name = os.path.join(ids.output_dir, "RT_Predicted_Model_Comparison.csv")
    save_model_comparison(selected_col, qc_atlas_df, rts_df, linear, poly, rt_comparison_file_name)
    models_file_name = os.path.join(ids.output_dir, "rt_model.txt")
    write_models(models_file_name, linear, poly, groups, qc_atlas)
    return (linear, poly)


def generate_outputs(
    ids, cpus, num_points=5, peak_height=5e5, use_poly_model=True, model_only=False, selected_col="median"
):
    """
    Generate the RT correction models, associated atlases with adjusted RT values, follow up notebooks,
    msms hits pickles
    inputs:
        ids: an AnalysisIds object matching the one used in the main notebook
        cpus: max number of cpus to use
        num_points: minimum number of data points in a peak
        peak_height: threshold intensity level for filtering
        use_poly_model: If True, use the polynomial model, else use linear model
                        Both types of models are always generated, this only determines which ones
                        are pre-populated into the generated notebooks
        model_only: If True, do not create atlases or notebooks, if False create them
        selected_col: name of column to use for model generation
    """
    linear, poly = generate_rt_correction_models(ids, cpus, selected_col)
    if not model_only:
        atlases = create_adjusted_atlases(linear, poly, ids)
        write_notebooks(ids, atlases, use_poly_model)
        pre_process_data_for_all_notebooks(ids, atlases, cpus, use_poly_model, num_points, peak_height)
    targeted_output.copy_outputs_to_google_drive(ids)
    targeted_output.archive_outputs(ids)
    logger.info("RT correction notebook complete. Switch to Targeted notebook to continue.")


def pre_process_data_for_all_notebooks(ids, atlases, cpus, use_poly_model, num_points, peak_height):
    """
    inputs:
        ids: an AnalysisIds object matching the one used in the main notebook
        atlases: list of atlas names to consider generating hits for
        cpus: max number of cpus to use
        use_poly_model: If True, use the polynomial model, else use linear model
                        Both types of models are always generated, this only determines which ones
                        are pre-populated into the generated notebooks
        num_points: minimum number of data points in a peak
        peak_height: threshold intensity level for filtering
    Calls MetatlasDataset().hits, which will create a hits cache file
    Filters compounds by signal strength to reduce atlas size
    """
    for atlas_name in atlases:
        if (use_poly_model and "linear" in atlas_name) or (not use_poly_model and "polynomial" in atlas_name):
            continue
        polarity = "positive" if "_POS_" in atlas_name else "negative"
        output_type = "FinalEMA-HILIC" if "EMA_Unlab" in atlas_name else "ISTDsEtc"
        current_ids = mads.AnalysisIdentifiers(
            source_atlas=atlas_name,
            experiment=ids.experiment,
            output_type=output_type,
            polarity=polarity,
            analysis_number=ids.analysis_number,
            project_directory=ids.project_directory,
            google_folder=ids.google_folder,
        )
        metatlas_dataset = mads.MetatlasDataset(ids=current_ids, max_cpus=cpus)
        _ = metatlas_dataset.hits
        if metatlas_dataset.ids.output_type in ["FinalEMA-HILIC"]:
            metatlas_dataset.filter_compounds_by_signal(num_points=num_points, peak_height=peak_height)


def get_groups(ids):
    """
    Create all experiment groups if they don't already exist and return the subset matching include_list
    inputs:
        ids: instance of AnalysisIds
    """
    ordered_groups = sorted(ids.groups, key=lambda x: x.name)
    for grp in ordered_groups:
        logger.info("Selected group: %s, %s", grp.name, int_to_date_str(grp.last_modified))
    return ordered_groups


def int_to_date_str(i_time):
    """unix epoc time in seconds to YYYY-MM-DD hh:mm:ss"""
    return str(datetime.fromtimestamp(i_time))


def get_files_df(groups):
    """Pandas Datafram with one row per file plus columns for accquistion_time and group name"""
    files_df = pd.DataFrame(columns=["file", "time", "group"])
    for group in groups:
        for run in group.items:
            try:
                time = run.accquistion_time
            except AttributeError:
                time = 0
            files_df = files_df.append({"file": run, "time": time, "group": group}, ignore_index=True)
    return files_df.sort_values(by=["time"])


def get_qc_atlas(ids):
    """Retreives template QC atlas and return tuple (atlas, atlas_df)"""
    qc_atlas_dict = QC_ATLASES[ids.polarity][ids.chromatography]
    qc_atlas_name = qc_atlas_dict["name"]
    username = qc_atlas_dict["username"]
    logger.info("Loading QC Atlas %s", qc_atlas_name)
    atlas = metob.retrieve("Atlas", name=qc_atlas_name, username=username)[0]
    atlas_df = ma_data.make_atlas_df(atlas)
    atlas_df["label"] = [cid.name for cid in atlas.compound_identifications]
    return atlas, atlas_df


def load_runs(files_df, qc_atlas_df, qc_atlas, cpus):
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


def save_measured_rts(metatlas_dataset, file_name):
    """Save RT values in csv format file"""
    rts_df = get_rts(metatlas_dataset, include_atlas_rt_peak=False)
    write_utils.export_dataframe_die_on_diff(rts_df, file_name, "measured RT values", float_format="%.6e")


def save_rt_peak(metatlas_dataset, file_name):
    """Save peak RT values in tsv format file"""
    rts_df = dp.make_output_dataframe(input_dataset=metatlas_dataset, fieldname="rt_peak", use_labels=True)
    write_utils.export_dataframe_die_on_diff(
        rts_df, file_name, "peak RT values", sep="\t", float_format="%.6e"
    )


def get_rts(metatlas_dataset, include_atlas_rt_peak=True):
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


def order_df_columns_by_run(dataframe):
    """
    Returns a dataframe with re-ordered columns such that second column up to column 'mean'
    are ordered by run number from low to high
    """
    cols = dataframe.columns.tolist()
    stats_start_idx = cols.index("mean")
    to_sort = cols[:stats_start_idx]
    no_sort = cols[stats_start_idx:]
    to_sort.sort(key=lambda x: int(x.split(".")[0].split("_Run")[1]))
    new_cols = to_sort + no_sort
    return dataframe[new_cols]


def plot_compound_atlas_rts(num_files, rts_df, file_name, fontsize=2, pad=0.1, cols=8):
    """
    Writes plot of RT peak for vs file for each compound
    inputs:
        num_files: number of files in data set, ie len(metatlas_dataset)
        rts_df: Dataframe with RTs values
        file_name: where to save plot
        fontsize: size of text
        pad: padding size
        cols: number of columns in plot grid
    """
    logger.info("Plotting RT Peak vs file for each compound")
    rts_df_plot = (
        rts_df.sort_values(by="standard deviation", ascending=False, na_position="last")
        .drop(["#NaNs"], axis=1)
        .dropna(axis=0, how="all", subset=rts_df.columns[:num_files])
    )
    rows = int(math.ceil((rts_df.shape[0] + 1) / 8))
    fig = plt.figure()
    grid = gridspec.GridSpec(rows, cols, figure=fig, wspace=0.2, hspace=0.4)
    for i, (_, row) in tqdm(enumerate(rts_df_plot.iterrows()), total=len(rts_df_plot), unit="plot"):
        a_x = fig.add_subplot(grid[i])
        a_x.tick_params(direction="in", length=1, pad=pad, width=0.1, labelsize=fontsize)
        a_x.scatter(range(num_files), row[:num_files], s=0.2)
        a_x.axhline(y=row["atlas RT peak"], color="r", linestyle="-", linewidth=0.2)
        a_x.set_xlim(-0.5, num_files + 0.5)
        a_x.xaxis.set_major_locator(mticker.FixedLocator(np.arange(0, num_files, 1.0)))
        range_columns = list(rts_df_plot.columns[:num_files]) + ["atlas RT peak"]
        a_x.set_ylim(np.nanmin(row.loc[range_columns]) - 0.12, np.nanmax(row.loc[range_columns]) + 0.12)
        _ = [s.set_linewidth(0.1) for s in a_x.spines.values()]
        # truncate name so it fits above a single subplot
        a_x.set_title(row.name[:33], pad=pad, fontsize=fontsize)
        a_x.set_xlabel("Files", labelpad=pad, fontsize=fontsize)
        a_x.set_ylabel("Actual RTs", labelpad=pad, fontsize=fontsize)
    plt.savefig(file_name, bbox_inches="tight")


def generate_models(actual_df, pred_df):
    """
    inputs:
        actual_df: dataframe with experimental RTs
        pred_df: dataframe with predicted RTs
    returns tuple containing two Model classes of order 1 and 2
    """
    ransac = RANSACRegressor(random_state=42)
    rt_model_linear = ransac.fit(pred_df, actual_df)
    linear = Model(
        rt_model_linear, rt_model_linear.estimator_.intercept_, rt_model_linear.estimator_.coef_[0]
    )

    poly_reg = PolynomialFeatures(degree=2)
    x_poly = poly_reg.fit_transform(pred_df)
    rt_model_poly = LinearRegression().fit(x_poly, actual_df)
    poly = Model(rt_model_poly, rt_model_poly.intercept_, rt_model_poly.coef_)
    return linear, poly


def actual_and_predicted_df(selected_column, rts_df, atlas_df):
    """
    inputs:
        selected_column: column number in rts_df to use for actual values
        rts_df: dataframe of RT values
        atlas_df: QC atlas in dataframe format
    return a tuple of (actual_df, pred_df)
    """
    actual_df = rts_df.loc[:, selected_column]
    bad_qc_compounds = np.where(~np.isnan(actual_df))
    actual_df = actual_df.iloc[bad_qc_compounds]
    pred_df = atlas_df.iloc[bad_qc_compounds][["rt_peak"]]
    return actual_df, pred_df


def actual_and_predicted_rts(rts_df, atlas_df, actual_df, pred_df):
    """
    inputs:
        rts_df: dataframe of RT values
        atlas_df: QC atlas in dataframe format
        acutal_df: dataframe of actual RT values
        pred_df: dataframe of predicted RT values
    return a tuple of lists of lists: (actual_rts, pred_rts)
    """
    actual_rts = [actual_df.values.tolist()]
    pred_rts = [pred_df.values.tolist()]
    for i in range(rts_df.shape[1] - 5):
        current_actual_df = rts_df.loc[:, rts_df.columns[i]]
        bad_qc_compounds = np.where(~np.isnan(current_actual_df))
        current_actual_df = current_actual_df.iloc[bad_qc_compounds]
        current_pred_df = atlas_df.iloc[bad_qc_compounds][["rt_peak"]]
        actual_rts.append(current_actual_df.values.tolist())
        pred_rts.append(current_pred_df.values.tolist())
    return actual_rts, pred_rts


def plot_actual_vs_pred_rts(pred_rts, actual_rts, rts_df, file_name, linear, poly):
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
        sub = fig.add_subplot(grid[i])
        x_values = list(itertools.chain(*pred_rts[i]))
        y_values = actual_rts[i]
        sub.scatter(x_values, y_values, s=2)
        spaced_x = np.linspace(0, max(x_values), 100)
        sub.plot(spaced_x, linear.predict(spaced_x), linewidth=0.5, color="red")
        sub.plot(spaced_x, poly.predict(spaced_x), linewidth=0.5, color="green")
        sub.set_title("File: " + str(i))
        sub.set_xlabel("predicted RTs")
        sub.set_ylabel("actual RTs")
    fig_legend = "FileIndex       FileName"
    for i in range(rts_df.shape[1] - 5):
        fig_legend = fig_legend + "\n" + str(i) + "        " + rts_df.columns[i]
    fig.tight_layout(pad=0.5)
    plt.text(0, -0.03 * rts_df.shape[1], fig_legend, transform=plt.gcf().transFigure)
    plt.savefig(file_name, bbox_inches="tight")


def save_model_comparison(selected_column, qc_atlas_df, rts_df, linear, poly, file_name):
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
    # qc_df["RT Reference"] = qc_atlas_df["rt_peak"]
    qc_df.loc[:, "RT Reference"] = qc_atlas_df["rt_peak"].to_numpy()
    qc_df.loc[:, "RT Linear Pred"] = pd.Series(
        linear.predict(qc_df["RT Reference"].to_numpy()), index=qc_df.index
    )
    qc_df.loc[:, "RT Polynomial Pred"] = pd.Series(
        poly.predict(qc_df["RT Reference"].to_numpy()), index=qc_df.index
    )
    # qc_df["RT Linear Pred"] = linear.predict(qc_df["RT Reference"].to_numpy())
    # qc_df["RT Polynomial Pred"] = poly.predict(qc_df["RT Reference"].to_numpy())
    qc_df["RT Diff Linear"] = qc_df["RT Measured"] - qc_df["RT Linear Pred"]
    qc_df["RT Diff Polynomial"] = qc_df["RT Measured"] - qc_df["RT Polynomial Pred"]
    write_utils.export_dataframe_die_on_diff(qc_df, file_name, "model comparision", float_format="%.6e")


def write_models(file_name, linear_model, poly_model, groups, atlas):
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


def get_atlas_name(template_name, ids, model, free_text):
    """
    input:
        template_name: name of template atlas
        ids: an AnalysisIds object matching the one used in the main notebook
        model: an instance of Model
        free_text: arbitrary string to append to atlas name
    returns the name of the production atlas
    """
    prod_name = template_name.replace("TPL", "PRD")
    prod_atlas_name = f"{prod_name}_{model.name}_{ids.project}_{ids.analysis}"
    if free_text != "":
        prod_atlas_name += f"_{free_text}"
    return prod_atlas_name


def adjust_atlas(atlas, model, ids):
    """use model to adjust RTs within atlas"""
    atlas_df = ma_data.make_atlas_df(atlas)
    atlas_df["label"] = [cid.name for cid in atlas.compound_identifications]
    atlas_df["rt_peak"] = model.predict(atlas_df["rt_peak"].to_numpy())
    rt_offset = 0.2 if ids.chromatography == "C18" else 0.5
    atlas_df["rt_min"] = atlas_df["rt_peak"].apply(lambda rt: rt - rt_offset)
    atlas_df["rt_max"] = atlas_df["rt_peak"].apply(lambda rt: rt + rt_offset)
    return atlas_df


def get_template_atlas(ids, polarity, idx):
    """Retreives a template atlas with the correct chromatorgraphy and polarity"""
    template = TEMPLATES[polarity][ids.chromatography][idx]
    return metob.retrieve("Atlas", **template)[-1]


def create_adjusted_atlases(linear, poly, ids, atlas_indices=None, free_text=""):
    """
    input:
        linear_model: instance of class Model with first order model
        poly_model: instance of class Model with second order model
        ids: an AnalysisIds object matching the one used in the main notebook
        atlas_indices: list of integers for which adjusted atlases to create
                        0: EMA_Unlab
                        1: QCv3_Unlab
                        2: ISv5_Unlab
                        3: ISv5_13C15N
                        4: IS_LabUnlab2
        free_text: arbitrary string to append to atlas name
    returns a list of the names of atlases
    """
    # pylint: disable=too-many-locals
    assert ids.chromatography in ["HILICZ", "C18"]
    default_atlas_indices = [0, 1] if ids.chromatography == "C18" else [0, 4]
    atlas_indices = default_atlas_indices if atlas_indices is None else atlas_indices
    plot_vars = [
        (polarity, idx, model)
        for polarity in ["positive", "negative"]
        for idx in atlas_indices
        for model in [linear, poly]
    ]
    out_atlas_names = []
    for polarity, idx, model in tqdm(plot_vars, unit="atlas"):
        template_atlas = get_template_atlas(ids, polarity, idx)
        out_atlas_names.append(get_atlas_name(template_atlas.name, ids, model, free_text))
        logger.info("Creating atlas %s", out_atlas_names[-1])
        out_atlas_file_name = os.path.join(ids.output_dir, f"{out_atlas_names[-1]}.csv")
        out_atlas_df = adjust_atlas(template_atlas, model, ids)
        write_utils.export_dataframe_die_on_diff(
            out_atlas_df, out_atlas_file_name, "predicted atlas", index=False, float_format="%.6e"
        )
        dp.make_atlas_from_spreadsheet(
            out_atlas_df,
            out_atlas_names[-1],
            filetype="dataframe",
            sheetname="",
            polarity=polarity,
            store=True,
            mz_tolerance=10 if ids.chromatography == "C18" else 12,
        )
    return out_atlas_names


def write_notebooks(ids, atlases, use_poly_model):
    """
    Creates Targeted analysis jupyter notebooks with pre-populated parameter sets
    Inputs:
        ids: an AnalysisIds object matching the one used in the main notebook
        atlases: list of atlas names to use as source atlases
        use_poly_model: if True use polynomial RT prediction model, else use linear model
                        this value is used to filter atlases from the input atlases list
    """
    for atlas_name in atlases:
        if (use_poly_model and "linear" in atlas_name) or (not use_poly_model and "polynomial" in atlas_name):
            continue
        polarity = "positive" if "_POS_" in atlas_name else "negative"
        short_polarity = "POS" if polarity == "positive" else "NEG"
        output_type = "FinalEMA-HILIC" if "EMA_Unlab" in atlas_name else "ISTDsEtc"
        repo_path = Path(__file__).resolve().parent.parent.parent
        source = repo_path / "notebooks" / "reference" / "Targeted.ipynb"
        dest = Path(ids.output_dir).resolve().parent / f"{ids.project}_{output_type}_{short_polarity}.ipynb"
        # include_groups and exclude_groups do not get passed to subsequent notebooks
        # as they need to be updated for each output type
        parameters = {
            "experiment": ids.experiment,
            "output_type": output_type,
            "polarity": polarity,
            "analysis_number": 0,
            "project_directory": ids.project_directory,
            "source_atlas": atlas_name,
            "exclude_files": ids.exclude_files,
            "groups_controlled_vocab": ids.groups_controlled_vocab,
        }
        notebook.create_notebook(source, dest, parameters)


def get_analysis_ids_for_rt_prediction(
    experiment,
    project_directory,
    google_folder,
    analysis_number=0,
    polarity="positive",
    exclude_files=None,
    include_groups=None,
    exclude_groups=None,
    groups_controlled_vocab=None,
):
    """
    Simplified interface for generating an AnalysisIds instance for use in rt prediction
    inputs:
        experiment: name of experiment as given in LCMS run names
        project_directory: directory where per-experiment output directory will be created
        google_folder: id from URL of base export folder on Google Drive
        analysis_number: integer, defaults to 0, increment if redoing analysis
        polarity: polarity to use for RT prediction, defaults to positive
        exclude_files: list of substrings that will be used to filter out lcmsruns
        include_groups: list of substrings that will used to filter groups
        exclude_groups list of substrings that will used to filter out groups
        groups_controlled_vocab: list of substrings that will group all matches into one group
    Returns an AnalysisIds instance
    """
    return mads.AnalysisIdentifiers(
        experiment=experiment,
        output_type="data_QC",
        analysis_number=analysis_number,
        project_directory=project_directory,
        polarity=polarity,
        google_folder=google_folder,
        exclude_files=exclude_files,
        include_groups=include_groups,
        exclude_groups=exclude_groups,
        groups_controlled_vocab=groups_controlled_vocab,
    )
