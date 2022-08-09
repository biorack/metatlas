""" process a targeted experiment"""

import getpass
import logging
import shutil

from typing import Optional

from IPython.display import display

import metatlas.datastructures.analysis_identifiers as analysis_ids
import metatlas.plots.dill2plots as dp

from metatlas.datastructures.id_types import Experiment
from metatlas.datastructures.metatlas_dataset import MetatlasDataset
from metatlas.datastructures.utils import Username
from metatlas.io.gdrive import copy_outputs_to_google_drive
from metatlas.io.targeted_output import generate_all_outputs
from metatlas.tools.config import Config, Workflow, Analysis
from metatlas.tools.notebook import in_papermill
from metatlas.io.targeted_output import generate_qc_outputs

logger = logging.getLogger(__name__)


# pylint: disable=too-many-arguments,too-many-locals
def pre_annotation(
    experiment: Experiment,
    rt_alignment_number: int,
    analysis_number: int,
    source_atlas_unique_id: str,
    configuration: Config,
    workflow: Workflow,
    analysis: Analysis,
    clear_cache: bool = False,
    username: Optional[Username] = None,
) -> MetatlasDataset:
    """All data processing that needs to occur before the annotation GUI in Targeted notebook"""
    params = analysis.parameters
    ids = analysis_ids.AnalysisIdentifiers(
        workflow=workflow.name,
        analysis=analysis.name,
        source_atlas_unique_id=source_atlas_unique_id,
        copy_atlas=params.copy_atlas,
        polarity=params.polarity,
        analysis_number=analysis_number,
        experiment=experiment,
        include_groups=params.include_groups.gui,
        exclude_groups=params.exclude_groups.gui,
        groups_controlled_vocab=params.groups_controlled_vocab,
        include_lcmsruns=params.include_lcmsruns.gui,
        exclude_lcmsruns=params.exclude_lcmsruns.gui,
        rt_alignment_number=rt_alignment_number,
        project_directory=params.project_directory,
        google_folder=params.google_folder,
        username=getpass.getuser() if username is None else username,
        configuration=configuration,
    )
    if clear_cache:
        logger.info("Clearing cache.")
        shutil.rmtree(ids.cache_dir)
    metatlas_dataset = MetatlasDataset(ids=ids, max_cpus=params.max_cpus)
    metatlas_dataset.filter_compounds_by_signal(params.num_points, params.peak_height, params.msms_score)
    return metatlas_dataset


def annotation_gui(
    data: MetatlasDataset,
    compound_idx: int = 0,
    width: float = 15,
    height: float = 3,
    alpha: float = 0.5,
    colors="",
    adjustable_rt_peak=False,
) -> Optional[dp.adjust_rt_for_selected_compound]:
    """
    Opens the interactive GUI for setting RT bounds and annotating peaks
    inputs:
        compound_idx: number of compound-adduct pair to start at
        width: width of interface in inches
        height: height of each plot in inches
        alpha: (0-1] controls transparency of lines on EIC plot
        colors: list (color_id, search_string) for coloring lines on EIC plot
                based on search_string occuring in LCMS run filename
    """
    display(dp.LOGGING_WIDGET)  # surface event handler error messages in UI
    return dp.adjust_rt_for_selected_compound(
        data,
        msms_hits=data.hits,
        color_me=colors,
        compound_idx=compound_idx,
        alpha=alpha,
        width=width,
        height=height,
        adjustable_rt_peak=adjustable_rt_peak,
    )


# pylint:disable=unused-argument
def post_annotation(
    data: MetatlasDataset, configuration: Config, workflow: Workflow, analysis: Analysis
) -> None:
    """All data processing that needs to occur after the annotation GUI in Targeted notebook"""
    params = analysis.parameters
    if params.require_all_evaluated and not in_papermill():
        data.error_if_not_all_evaluated()
    if params.filter_removed:
        data.filter_compounds_ms1_notes_remove()
    data.extra_time = 0.5
    logger.info("extra_time set to 0.5 minutes for output generation.")
    data.update()  # update hits and data if they no longer are based on current rt bounds
    if params.generate_qc_outputs:
        data.ids.set_output_state(params, "qc_outputs")
        generate_qc_outputs(data)
    if params.generate_analysis_outputs:
        generate_all_outputs(data, workflow, analysis)
    if not in_papermill():
        copy_outputs_to_google_drive(data.ids)
    logger.info("DONE - execution of notebook%s is complete.", " in draft mode" if in_papermill() else "")
