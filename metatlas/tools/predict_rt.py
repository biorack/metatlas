"""Generate Retention Time Correction Model"""
# pylint: disable=too-many-arguments

import itertools
import math
import os

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

from matplotlib import gridspec
from sklearn.linear_model import LinearRegression, RANSACRegressor
from sklearn.preprocessing import PolynomialFeatures

from metatlas.datastructures import metatlas_dataset as mads
from metatlas.datastructures import metatlas_objects as metob
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.plots import dill2plots as dp


TEMPLATES = {
    "postive": [
        "HILICz150_ANT20190824_TPL_EMA_Unlab_POS",
        "HILICz150_ANT20190824_TPL_QCv3_Unlab_POS",
        "HILICz150_ANT20190824_TPL_ISv5_Unlab_POS",
        "HILICz150_ANT20190824_TPL_ISv5_13C15N_POS",
        "HILICz150_ANT20190824_TPL_IS_LabUnlab2_POS",
    ],
    "negative": [
        "HILICz150_ANT20190824_TPL_EMA_Unlab_NEG",
        "HILICz150_ANT20190824_TPL_QCv3_Unlab_NEG",
        "HILICz150_ANT20190824_TPL_ISv5_Unlab_NEG",
        "HILICz150_ANT20190824_TPL_ISv5_13C15N_NEG",
        "HILICz150_ANT20190824_TPL_IS_LabUnlab2_NEG",
    ],
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
        if isinstance(coefficents, list):
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
        return len(self.coefficents)

    @property
    def name(self):
        """Type of model as string"""
        return "linear" if self.order == 1 else "polynomial"

    def predict(self, x_values):
        """Returns y values for input x"""
        return self.sk_model.predict(x_values)


def generate_rt_correction_models(
    ids, groups_controlled_vocab, exclude_files, include_groups, cpus, save_to_db=True
):
    """
    Generate the RT correction models and associated atlases with adjusted RT values
    inputs:
        ids: an AnalysisIds object matching the one used in the main notebook
        groups_controlled_vocab: list of strings that will group together when creating groups
                                 application of groups_controlled_vocab is case insensitive
        exclude_files: list of strings that will exclude files if they are substrings of the filename
        include_groups: group will only be used in correction if their name has a substring match
                        to this list of strings
        cpus: max number of cpus to use
        save_to_db: If True, save the new atlases to the database
    """
    # pylint: disable=too-many-locals
    metatlas_dataset = mads.MetatlasDataset(ids, groups_controlled_vocab, exclude_files, save_metadata=False)
    qc_dir = os.path.join(ids.output_dir, "data_QC")
    groups = get_groups(metatlas_dataset, include_groups)
    files_df = get_files_df(groups)
    qc_atlas, qc_atlas_df = get_qc_atlas(metatlas_dataset.ids)
    metatlas_dataset = load_runs(files_df, qc_atlas_df, qc_atlas, cpus)
    save_measured_rts(metatlas_dataset, os.path.join(qc_dir, "QC_Measured_RTs.csv"))
    rts_df = get_rts(metatlas_dataset)
    plot_compound_atlas_rts(metatlas_dataset, rts_df, os.path.join(qc_dir, "Compound_Atlas_RTs.pdf"))
    selected_column = 9  # need to deal with this parameter, index from rts_df.columns
    actual_df, pred_df = actual_and_predicted_df(selected_column, rts_df, qc_atlas_df)
    linear, poly = generate_models(actual_df, pred_df)
    actual_rts, pred_rts = actual_and_predicted_rts(rts_df, qc_atlas_df, actual_df, pred_df)
    actual_vs_pred_file_name = os.path.join(qc_dir, "Actual_vs_Predicted_RTs.pdf")
    plot_actual_vs_pred_rts(pred_rts, actual_rts, rts_df, actual_vs_pred_file_name, linear, poly)
    rt_comparison_file_name = os.path.join(qc_dir, "RT_Predicted_Model_Comparison.csv")
    save_model_comparison(selected_column, qc_atlas_df, rts_df, linear, poly, rt_comparison_file_name)
    models_file_name = os.path.join(qc_dir, "rt_model.txt")
    write_models(models_file_name, linear, poly, groups, qc_atlas)
    create_adjusted_atlases(linear, poly, qc_dir, save_to_db=save_to_db)


def get_groups(metatlas_dataset, include_groups):
    """
    Create all experiment groups if they don't already exist and return the subset matching include_list
    inputs:
        metatlas_datset: instance of MetatlasDataset
        include_groups: group will only be used in correction if their name has a substring match
                        to this list of strings
    """
    metatlas_dataset.store_groups(exist_ok=True)
    ids = metatlas_dataset.ids
    groups = dp.select_groups_for_analysis(
        name=f"{ids.experiment}_{ids.short_polarity}_%",
        most_recent=True,
        remove_empty=True,
        include_list=include_groups,
        exclude_list=ids.short_polarity_inverse,
    )
    return sorted(groups, key=lambda x: x.name)


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
    qc_atlas_name = f"HILICz150_ANT20190824_TPL_QCv3_Unlab_{ids.short_polarity}"
    atlas = metob.retrieve("Atlas", name=qc_atlas_name, username="vrsingan")[0]
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
    return mads.parallel_process(ma_data.get_data_for_atlas_df_and_file, files, cpus, unit="sample")


def save_measured_rts(metatlas_dataset, filename):
    """Save RT values in csv format file"""
    rts_df = dp.make_output_dataframe(
        input_dataset=metatlas_dataset,
        fieldname="rt_peak",
        use_labels=True,
        summarize=True,
    )
    rts_df.to_csv(filename)


def get_rts(metatlas_dataset):
    """Returns RT values in DataFrame format"""
    rts_df = dp.make_output_dataframe(
        input_dataset=metatlas_dataset,
        fieldname="rt_peak",
        use_labels=True,
        summarize=True,
    )
    rts_df["atlas RT peak"] = [
        compound["identification"].rt_references[0].rt_peak for compound in metatlas_dataset[0]
    ]
    return rts_df


def plot_compound_atlas_rts(num_files, rts_df, file_name):
    """
    Writes plot of RT peak for vs file for each compound
    inputs:
        num_files: number of files in data set, ie len(metatlas_dataset)
        rts_df: Dataframe with RTs values
        filename: where to save plot
    """
    # pylint: disable=too-many-locals
    # number of columns in rts_df that are not values from a specific input file
    num_not_files = len(rts_df.columns) - num_files
    rts_df_plot = (
        rts_df.sort_values(by="standard deviation", ascending=False, na_position="last")
        .drop(["#NaNs"], axis=1)
        .dropna(axis=0, how="all", subset=rts_df.columns[:-num_not_files])
    )
    fontsize = 2
    pad = 0.1
    cols = 8
    rows = int(math.ceil((rts_df.shape[0] + 1) / 8))
    fig = plt.figure()
    grid = gridspec.GridSpec(rows, cols, figure=fig, wspace=0.2, hspace=0.4)
    for i, (_, row) in enumerate(rts_df_plot.iterrows()):
        a_x = fig.add_subplot(grid[i])
        a_x.tick_params(direction="in", length=1, pad=pad, width=0.1, labelsize=fontsize)
        a_x.scatter(range(rts_df_plot.shape[1] - num_not_files), row[:-num_not_files], s=0.2)
        ticks_loc = np.arange(0, len(rts_df_plot.columns) - num_not_files, 1.0)
        a_x.a_xhline(y=row["atlas RT peak"], color="r", linestyle="-", linewidth=0.2)
        a_x.set_xlim(-0.5, len(rts_df_plot.columns) - num_not_files + 0.5)
        a_x.xa_xis.set_major_locator(mticker.FixedLocator(ticks_loc))
        range_columns = list(rts_df_plot.columns[:-num_not_files]) + ["atlas RT peak"]
        a_x.set_ylim(np.nanmin(row.loc[range_columns]) - 0.12, np.nanma_x(row.loc[range_columns]) + 0.12)
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
    actual_df = rts_df.loc[:, rts_df.columns[selected_column]]
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
        filename: where to save the plot
    """
    qc_df = rts_df[[rts_df.columns[selected_column]]].copy()
    qc_df.columns = ["RT Measured"]
    qc_df["RT Reference"] = qc_atlas_df["rt_peak"]
    qc_df["RT Linear Pred"] = linear.predict(qc_df["RT Reference"])
    qc_df["RT Polynomial Pred"] = poly.predict(qc_df["RT Reference"])
    qc_df["RT Diff Linear"] = qc_df["RT Measured"] - qc_df["RT Linear Pred"]
    qc_df["RT Diff Polynomial"] = qc_df["RT Measured"] - qc_df["RT Polynomial Pred"]
    qc_df.to_csv(file_name)


def write_models(file_name, linear_model, poly_model, groups, atlas):
    """
    inputs:
        filename: text file to save model information
        linear_model: instance of class Model with first order model
        poly_model: instance of class Model with second order model
        groups: list of groups used in model generation
        atlas: QC atlas
    """
    with open(file_name, "w") as out_fh:
        for model in [linear_model, poly_model]:
            out_fh.write(f"{model.sk_model.set_params()}\n")
            out_fh.write(f"{model}\n")
            group_names = ", ".join([g.name for g in groups])
            out_fh.write(f"groups = {group_names}\n")
            out_fh.write(f"atlas = {atlas.name}\n\n")


def create_adjusted_atlases(linear, poly, qc_dir, atlas_indices=None, free_text="", save_to_db=True):
    """
    input:
        linear_model: instance of class Model with first order model
        poly_model: instance of class Model with second order model
        qc_dir: directory to write csv files to
        atlas_indices: list of integers for which adjusted atlases to create
                        0: EMA_Unlab
                        1: QCv3_Unlab
                        2: ISv5_Unlab
                        3: ISv5_13C15N
                        4: IS_LabUnlab2
        free_text: arbitrary string to append to atlas name
        save_to_db: if True, save the atlases to the database
    """
    if atlas_indices is None:
        atlas_indices = [0, 4]
    for polarity in ["positive", "negative"]:
        for idx in atlas_indices:
            for model in [linear, poly]:
                template_name = TEMPLATES[polarity][idx]
                atlas = metob.retrieve("Atlas", name=template_name, username="vrsingan")[-1]
                prd_atlas_name = template_name.replace("TPL", "PRD") + f"_{model.name}"
                if free_text != "":
                    prd_atlas_name = prd_atlas_name + "_" + free_text
                prd_atlas_filename = prd_atlas_name + ".csv"
                prd_atlas_df = ma_data.make_atlas_df(atlas)
                prd_atlas_df["label"] = [cid.name for cid in atlas.compound_identifications]
                prd_atlas_df["rt_peak"] = model.predict(prd_atlas_df["rt_peak"])
                prd_atlas_df["rt_min"] = prd_atlas_df["rt_peak"].apply(lambda rt: rt - 0.5)
                prd_atlas_df["rt_max"] = prd_atlas_df["rt_peak"].apply(lambda rt: rt + 0.5)
                prd_atlas_df.to_csv(os.path.join(qc_dir, prd_atlas_filename), index=False)
                if save_to_db:
                    dp.make_atlas_from_spreadsheet(
                        prd_atlas_df,
                        prd_atlas_name,
                        filetype="dataframe",
                        sheetname="",
                        polarity=polarity,
                        store=True,
                        mz_tolerance=12,
                    )
                print(prd_atlas_name + " Created!")
