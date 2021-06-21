"""Generate Retention Time Correction Model"""

import os
import sys
import time
import multiprocessing
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from metatlas.datastructures import metatlas_dataset as mads
from metatlas.datastructures import metatlas_objects as metob
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.plots import dill2plots as dp


def generate_rt_correction_models(ids, groups_controlled_vocab, exclude_files, include_groups):
    metatlas_dataset = mads.MetatlasDataset(ids, groups_controlled_vocab, exclude_files, save_metadata=False)
    groups = get_groups(metatlas_dataset, include_groups)
    files_df = get_files_df(groups)
    qc_atlas, qc_atlas_df = get_qc_atlas(metatlas_dataset.ids)


def get_groups(metatlas_dataset, include_groups):
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
    files_df = pd.DataFrame(columns=["file", "time", "group"])
    for group in groups:
        for run in group.items:
            time = run.accquistion_time if hasattr(run, "acquisition_time") else 0
            files_df = files_df.append({"file": run, "time": time, "group": group}, ignore_index=True)
    return files_df.sort_values(by=["time"])


def get_qc_atlas(ids):
    qc_atlas_name = f"HILICz150_ANT20190824_TPL_QCv3_Unlab_{ids.short_polarity}"
    atlas = metob.retrieve('Atlas', name=qc_atlas_name, username='vrsingan')[0]
    atlas_df = ma_data.make_atlas_df(atlas)
    atlas_df['label'] = [cid.name for cid in atlas.compound_identifications]
    return atlas, atlas_df


def load_runs(metatlas_dataset, files_df, qc_atlas_df, qc_atlas):
    files = []
    for file_data in files_df.iterrows():
        files.append((file_data[1].file, file_data[1].group, qc_atlas_df, qc_atlas))
    if metatlas_dataset.max_cpus > 1 and len(files) > 1:
        with multiprocessing.Pool(processes=min(metatlas_dataset.max_cpus, len(files))) as pool:
            data = pool.map(ma_data.get_data_for_atlas_df_and_file, files)
    else:  # skip multiprocessing as this makes for easier debugging
        data = [ma_data.get_data_for_atlas_df_and_file(i) for i in files]
    return data


rts_df = dp.make_output_dataframe(
    input_dataset=metatlas_dataset,
    fieldname="rt_peak",
    use_labels=True,
    output_loc=output_data_qc,
    summarize=True,
)
rts_df.to_csv(os.path.join(output_data_qc, "QC_Measured_RTs.csv"))

import itertools
import math
from __future__ import division
from matplotlib import gridspec
import matplotlib.ticker as mticker

rts_df["atlas RT peak"] = [
    compound["identification"].rt_references[0].rt_peak for compound in metatlas_dataset[0]
]
# number of columns in rts_df that are not values from a specific input file
num_not_files = len(rts_df.columns) - len(metatlas_dataset)
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
gs = gridspec.GridSpec(rows, cols, figure=fig, wspace=0.2, hspace=0.4)
for i, (index, row) in enumerate(rts_df_plot.iterrows()):
    ax = fig.add_subplot(gs[i])
    ax.tick_params(direction="in", length=1, pad=pad, width=0.1, labelsize=fontsize)
    ax.scatter(range(rts_df_plot.shape[1] - num_not_files), row[:-num_not_files], s=0.2)
    ticks_loc = np.arange(0, len(rts_df_plot.columns) - num_not_files, 1.0)
    ax.axhline(y=row["atlas RT peak"], color="r", linestyle="-", linewidth=0.2)
    ax.set_xlim(-0.5, len(rts_df_plot.columns) - num_not_files + 0.5)
    ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    range_columns = list(rts_df_plot.columns[:-num_not_files]) + ["atlas RT peak"]
    ax.set_ylim(np.nanmin(row.loc[range_columns]) - 0.12, np.nanmax(row.loc[range_columns]) + 0.12)
    [s.set_linewidth(0.1) for s in ax.spines.values()]
    # truncate name so it fits above a single subplot
    ax.set_title(row.name[:33], pad=pad, fontsize=fontsize)
    ax.set_xlabel("Files", labelpad=pad, fontsize=fontsize)
    ax.set_ylabel("Actual RTs", labelpad=pad, fontsize=fontsize)

plt.savefig(os.path.join(output_data_qc, "Compound_Atlas_RTs.pdf"), bbox_inches="tight")

for i, a in enumerate(rts_df.columns):
    print(i, a)

selected_column = 9

from sklearn.linear_model import LinearRegression, RANSACRegressor
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_absolute_error as mae

actual_rts, pred_rts, polyfit_rts = [], [], []

current_actual_df = rts_df.loc[:, rts_df.columns[selected_column]]
bad_qc_compounds = np.where(~np.isnan(current_actual_df))
current_actual_df = current_actual_df.iloc[bad_qc_compounds]
current_pred_df = atlas_df.iloc[bad_qc_compounds][["rt_peak"]]
actual_rts.append(current_actual_df.values.tolist())
pred_rts.append(current_pred_df.values.tolist())

ransac = RANSACRegressor(random_state=42)
rt_model_linear = ransac.fit(current_pred_df, current_actual_df)
coef_linear = rt_model_linear.estimator_.coef_[0]
intercept_linear = rt_model_linear.estimator_.intercept_

poly_reg = PolynomialFeatures(degree=2)
X_poly = poly_reg.fit_transform(current_pred_df)
rt_model_poly = LinearRegression().fit(X_poly, current_actual_df)
coef_poly = rt_model_poly.coef_
intercept_poly = rt_model_poly.intercept_

for i in range(rts_df.shape[1] - 5):
    current_actual_df = rts_df.loc[:, rts_df.columns[i]]
    bad_qc_compounds = np.where(~np.isnan(current_actual_df))
    current_actual_df = current_actual_df.iloc[bad_qc_compounds]
    current_pred_df = atlas_df.iloc[bad_qc_compounds][["rt_peak"]]
    actual_rts.append(current_actual_df.values.tolist())
    pred_rts.append(current_pred_df.values.tolist())

# User can change to use particular qc file
import itertools
import math
from __future__ import division
from matplotlib import gridspec

x = list(itertools.chain(*pred_rts))
y = list(itertools.chain(*actual_rts))

rows = int(math.ceil((rts_df.shape[1] + 1) / 5))
cols = 5
fig = plt.figure(constrained_layout=False)

gs = gridspec.GridSpec(rows, cols, figure=fig)
plt.rc("font", size=6)
plt.rc("axes", labelsize=6)
plt.rc("xtick", labelsize=3)
plt.rc("ytick", labelsize=3)


for i in range(rts_df.shape[1] - 5):
    x = list(itertools.chain(*pred_rts[i]))
    y = actual_rts[i]

    ax = fig.add_subplot(gs[i])
    ax.scatter(x, y, s=2)
    ax.plot(
        np.linspace(0, max(x), 100),
        coef_linear * np.linspace(0, max(x), 100) + intercept_linear,
        linewidth=0.5,
        color="red",
    )
    ax.plot(
        np.linspace(0, max(x), 100),
        (coef_poly[1] * np.linspace(0, max(x), 100))
        + (coef_poly[2] * (np.linspace(0, max(x), 100) ** 2))
        + intercept_poly,
        linewidth=0.5,
        color="green",
    )
    ax.set_title("File: " + str(i))
    ax.set_xlabel("predicted RTs")
    ax.set_ylabel("actual RTs")

fig_legend = "FileIndex       FileName"
for i in range(rts_df.shape[1] - 5):
    fig_legend = fig_legend + "\n" + str(i) + "        " + rts_df.columns[i]

fig.tight_layout(pad=0.5)
plt.text(0, -0.03 * rts_df.shape[1], fig_legend, transform=plt.gcf().transFigure)
plt.savefig(os.path.join(output_data_qc, "Actual_vs_Predicted_RTs.pdf"), bbox_inches="tight")

qc_df = rts_df[[rts_df.columns[selected_column]]]
qc_df = qc_df.copy()
print("Linear Parameters :", coef_linear, intercept_linear)
print("Polynomial Parameters :", coef_poly, intercept_poly)

qc_df.columns = ["RT Measured"]
atlas_df.index = qc_df.index
qc_df["RT Reference"] = atlas_df["rt_peak"]
qc_df["RT Linear Pred"] = qc_df["RT Reference"].apply(lambda rt: coef_linear * rt + intercept_linear)
qc_df["RT Polynomial Pred"] = qc_df["RT Reference"].apply(
    lambda rt: (coef_poly[1] * rt) + (coef_poly[2] * (rt ** 2)) + intercept_poly
)
qc_df["RT Diff Linear"] = qc_df["RT Measured"] - qc_df["RT Linear Pred"]
qc_df["RT Diff Polynomial"] = qc_df["RT Measured"] - qc_df["RT Polynomial Pred"]
qc_df.to_csv(os.path.join(output_data_qc, "RT_Predicted_Model_Comparison.csv"))

qc_df

# CHOOSE YOUR MODEL HERE (linear / polynomial).
# model = 'linear'
model = "polynomial"

# Save model

with open(os.path.join(output_data_qc, "rt_model.txt"), "w") as f:
    if model == "linear":
        f.write(
            "coef = {}\nintercept = {}\nqc_actual_rts = {}\nqc_predicted_rts = {}".format(
                coef_linear, intercept_linear, ", ".join([g.name for g in groups]), myAtlas.name
            )
        )
        f.write("\n" + repr(rt_model_linear.set_params()))

    else:
        f.write(
            "coef = {}\nintercept = {}\nqc_actual_rts = {}\nqc_predicted_rts = {}".format(
                coef_poly, intercept_poly, ", ".join([g.name for g in groups]), myAtlas.name
            )
        )
        f.write("\n" + repr(rt_model_poly.set_params()))


pos_atlas_indices = [0, 1, 2, 3, 4]
neg_atlas_indices = [0, 1, 2, 3, 4]
free_text = ""  # this will be appended to the end of the csv filename exported
save_to_db = False

for ix in pos_atlas_indices:
    atlases = metob.retrieve("Atlas", name=pos_templates[ix], username="vrsingan")
    prd_atlas_name = pos_templates[ix].replace("TPL", "PRD")
    if free_text != "":
        prd_atlas_name = prd_atlas_name + "_" + free_text
    prd_atlas_filename = prd_atlas_name + ".csv"
    myAtlas = atlases[-1]
    PRD_atlas_df = ma_data.make_atlas_df(myAtlas)
    PRD_atlas_df["label"] = [cid.name for cid in myAtlas.compound_identifications]
    if model == "linear":
        PRD_atlas_df["rt_peak"] = PRD_atlas_df["rt_peak"].apply(
            lambda rt: coef_linear * rt + intercept_linear
        )
    else:
        PRD_atlas_df["rt_peak"] = PRD_atlas_df["rt_peak"].apply(
            lambda rt: (coef_poly[1] * rt) + (coef_poly[2] * (rt ** 2)) + intercept_poly
        )
    PRD_atlas_df["rt_min"] = PRD_atlas_df["rt_peak"].apply(lambda rt: rt - 0.5)
    PRD_atlas_df["rt_max"] = PRD_atlas_df["rt_peak"].apply(lambda rt: rt + 0.5)

    PRD_atlas_df.to_csv(os.path.join(output_data_qc, prd_atlas_filename), index=False)

    if save_to_db:
        dp.make_atlas_from_spreadsheet(
            PRD_atlas_df,
            prd_atlas_name,
            filetype="dataframe",
            sheetname="",
            polarity="positive",
            store=True,
            mz_tolerance=12,
        )
    print(prd_atlas_name + " Created!")

for ix in neg_atlas_indices:
    atlases = metob.retrieve("Atlas", name=neg_templates[ix], username="vrsingan")
    prd_atlas_name = neg_templates[ix].replace("TPL", "PRD")
    if free_text != "":
        prd_atlas_name = prd_atlas_name + "_" + free_text
    prd_atlas_filename = prd_atlas_name + ".csv"
    myAtlas = atlases[-1]
    PRD_atlas_df = ma_data.make_atlas_df(myAtlas)
    PRD_atlas_df["label"] = [cid.name for cid in myAtlas.compound_identifications]
    if model == "linear":
        PRD_atlas_df["rt_peak"] = PRD_atlas_df["rt_peak"].apply(
            lambda rt: coef_linear * rt + intercept_linear
        )
    else:
        PRD_atlas_df["rt_peak"] = PRD_atlas_df["rt_peak"].apply(
            lambda rt: (coef_poly[1] * rt) + (coef_poly[2] * (rt ** 2)) + intercept_poly
        )
    PRD_atlas_df["rt_min"] = PRD_atlas_df["rt_peak"].apply(lambda rt: rt - 0.5)
    PRD_atlas_df["rt_max"] = PRD_atlas_df["rt_peak"].apply(lambda rt: rt + 0.5)

    PRD_atlas_df.to_csv(os.path.join(output_data_qc, prd_atlas_filename), index=False)

    if save_to_db:
        dp.make_atlas_from_spreadsheet(
            PRD_atlas_df,
            prd_atlas_name,
            filetype="dataframe",
            sheetname="",
            polarity="negative",
            store=True,
            mz_tolerance=12,
        )

    print(prd_atlas_name + " Created!")
