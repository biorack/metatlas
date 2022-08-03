"""Manage configuration options using a YAML file"""
# pylint: disable=too-few-public-methods

import getpass
import os

from pathlib import Path
from string import ascii_letters, digits
from typing import Dict, List, Optional, Sequence, Tuple, TypeVar

import yaml

from pydantic import BaseModel, validator

import metatlas.datastructures.metatlas_objects as metob

from metatlas.datastructures.id_types import Polarity
from metatlas.tools.util import or_default

ALLOWED_NAME_CHARS = ascii_letters + digits + "-"


class BaseNotebookParameters(BaseModel):
    """Parameters common to both RT_Alignment and Targeted notebooks"""

    copy_atlas: bool = False
    source_atlas: Optional[str] = None
    include_groups: Optional[List[str]] = None
    exclude_groups: List[str] = []
    exclude_files: List[str] = []
    groups_controlled_vocab: List[str] = []
    rt_min_delta: Optional[float] = None
    rt_max_delta: Optional[float] = None
    config_file_name: Optional[str] = None
    source_code_version_id: Optional[str] = None
    project_directory: str = str(Path().home())
    google_folder: Optional[str] = None
    max_cpus: int = 4
    log_level: str = "INFO"


class AnalysisNotebookParameters(BaseNotebookParameters):
    """Parameters for Targeted notebooks"""

    polarity: Polarity = Polarity("positive")
    generate_qc_outputs: bool = False
    num_points: Optional[int] = None
    peak_height: Optional[float] = None
    msms_score: Optional[float] = None
    filter_removed: bool = False
    line_colors: Optional[List[Tuple[str, str]]] = None
    require_all_evaluated: bool = False
    generate_analysis_outputs: bool = False
    exclude_groups_for_analysis_outputs: List[str] = []
    export_msms_fragment_ions: bool = False
    clear_cache: bool = False


class RTAlignmentNotebookParameters(BaseNotebookParameters):
    """Parameters used in RT Alignment notebooks"""

    inchi_keys_not_in_model: List[str] = []
    dependent_data_source: str = "median"
    use_poly_model: bool = False
    stop_before: Optional[str] = None


class Atlas(BaseModel):
    """Atlas specification"""

    unique_id: str
    name: str
    do_alignment: bool = False

    @validator("unique_id")
    @classmethod
    def single_unique_id(cls, to_check):
        """Check that unique_id exists in the database"""
        num_found = len(metob.retrieve("Atlas", unique_id=to_check, username="*"))
        if num_found == 1:
            return to_check
        if num_found == 0:
            raise ValueError(f"Atlas with unique_id '{to_check}' not found in database.")
        raise ValueError(f"{num_found} copies of atlas with unique_id '{to_check}' in database.")

    @validator("name")
    @classmethod
    def altas_name_unique_it_match(cls, to_check, values):
        """test that atlas unique_id and name are consistent"""
        if "unique_id" not in values:
            raise ValueError("unique_id was not supplied for atlas")
        unique_id = values["unique_id"]
        if len(metob.retrieve("Atlas", unique_id=unique_id, name=to_check, username="*")) == 0:
            raise ValueError(f"Atlas with unique_id '{unique_id}' does not have name '{to_check}'.")
        return to_check


class RTAlignment(BaseModel):
    """Define an RT-Alingment workflow step"""

    name: str = "RT-Alignment"
    atlas: Atlas
    parameters: RTAlignmentNotebookParameters

    @validator("name")
    @classmethod
    def allowed_name_chars(cls, to_check):
        """Only contains letters, digits and dashes"""
        return validate_allowed_chars("RT alignment name", to_check)


class Analysis(BaseModel):
    """Define an analysis workflow step"""

    name: str
    atlas: Atlas
    parameters: AnalysisNotebookParameters

    @validator("name")
    @classmethod
    def allowed_name_chars(cls, to_check):
        """Only contains letters, digits and dashes"""
        return validate_allowed_chars("analysis names", to_check)


class Workflow(BaseModel):
    """Define a targeted analysis workflow"""

    name: str
    rt_alignment: RTAlignment
    analyses: List[Analysis]

    @validator("name")
    @classmethod
    def allowed_name_chars(cls, to_check):
        """Only contains letters, digits and dashes"""
        return validate_allowed_chars("workflow names", to_check)

    def get_analysis(self, analysis_name: str) -> Analysis:
        """Returns Analysis with analysis_name or ValueError"""
        for analysis in self.analyses:
            if analysis.name == analysis_name:
                return analysis
        raise ValueError(f"Analysis named '{analysis_name}' was not found within the workflow.")


class Chromatography(BaseModel):
    """A type of chromotography"""

    name: str
    aliases: List[str] = []

    @validator("name")
    @classmethod
    def allowed_name_chars(cls, to_check):
        """Only contains letters, digits and dashes"""
        return validate_allowed_chars("chromatography names", to_check)


class Config(BaseModel):
    """This class defines the config file"""

    chromatography_types: List[Chromatography] = []
    workflows: List[Workflow]

    @validator("workflows")
    @classmethod
    def unique_workflow_names(cls, to_check):
        """Do not allow duplicated workflow names"""
        dup_workflow_names = get_dups([w.name for w in to_check])
        if dup_workflow_names:
            raise ValueError(f"Workflow names were redefined: {', '.join(dup_workflow_names)}.")
        return to_check

    def get_workflow(self, workflow_name: str) -> Workflow:
        """Returns workflow with workflow_name or ValueError"""
        for workflow in self.workflows:
            if workflow.name == workflow_name:
                return workflow
        raise ValueError(f"Workflow named '{workflow_name}' was not found in the configuration file.")

    def update(self, override_parameters: Dict) -> None:
        """update all parameters within self.workflows with any non-None values in override_parameters"""
        for flow in self.workflows:
            for analysis in [flow.rt_alignment] + flow.analyses:
                if analysis.parameters.source_atlas is not None:
                    analysis.atlas.name = analysis.parameters.source_atlas
                    analysis.atlas.username = getpass.getuser()
                for name in analysis.parameters.__dict__.keys():
                    if name in override_parameters and override_parameters[name] is not None:
                        setattr(analysis.parameters, name, override_parameters[name])
                analysis.parameters.google_folder = or_default(
                    analysis.parameters.google_folder, flow.rt_alignment.parameters.google_folder
                )


def load_config(file_name: os.PathLike) -> Config:
    """Return configuration from file"""
    with open(file_name, "r", encoding="utf-8") as yaml_fh:
        return Config.parse_obj(yaml.safe_load(yaml_fh))


def get_config(override_parameters: Dict) -> Tuple[Config, Workflow, Analysis]:
    """loads configuration from file and updates workflow parameters with override_parameters"""
    config = load_config(override_parameters["config_file_name"])
    config.update(override_parameters)
    workflow = config.get_workflow(override_parameters["workflow_name"])
    if "analysis_name" in override_parameters:
        analysis = workflow.get_analysis(override_parameters["analysis_name"])
    else:
        analysis = workflow.rt_alignment
    return (config, workflow, analysis)


Generic = TypeVar("Generic")


def get_dups(seq: Sequence[Generic]) -> List[Generic]:
    """Returns each non-unique value from seq"""
    seen = set()
    dupes = []
    for value in seq:
        if value in seen:
            dupes.append(value)
        else:
            seen.add(value)
    return dupes


def validate_allowed_chars(variable_name, to_check):
    """returns to_check if valid, otherwise raises ValueError"""
    if any(c not in ALLOWED_NAME_CHARS for c in to_check):
        raise ValueError(f"Only letters, numbers, and '-' are allowed in {variable_name}.")
    return to_check
