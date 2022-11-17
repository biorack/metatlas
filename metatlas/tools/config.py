"""Manage configuration options using a YAML file"""
# pylint: disable=too-few-public-methods

import logging
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

logger = logging.getLogger(__name__)


class OutputLists(BaseModel):
    """Lists that are used to configure outputs"""

    always: List[str] = []
    rt_alignment: List[str] = []
    qc_outputs: List[str] = []
    gui: List[str] = []
    ids_spreadsheet: List[str] = []
    chromatograms: List[str] = []
    data_sheets: List[str] = []
    box_plots: List[str] = []

    def update(self, override_parameters: Dict) -> None:
        """
        update all parameters with any non-None values in override_parameters
        Note, the override_parameters input is not the full parameter set,
        it contains a dict that maps to an OutputList
        """
        for name in self.__dict__:
            if name in override_parameters:
                override_value = override_parameters[name]
                if override_value is not None:
                    setattr(self, name, override_value)

    def distribute_always_values(self) -> None:
        """Append self.always to all other attributes, set self.always to []"""
        if self.always:
            for attribute_name in self.__dict__:
                if attribute_name != "always":
                    getattr(self, attribute_name).extend(self.always)
            self.always = []


class BaseNotebookParameters(BaseModel):
    """Parameters common to both RT_Alignment and Targeted notebooks"""

    copy_atlas: bool = False
    source_atlas_unique_id: Optional[str] = None
    include_groups: OutputLists = OutputLists()
    exclude_groups: OutputLists = OutputLists()
    include_lcmsruns: OutputLists = OutputLists()
    exclude_lcmsruns: OutputLists = OutputLists()
    groups_controlled_vocab: List[str] = []
    rt_min_delta: Optional[float] = None
    rt_max_delta: Optional[float] = None
    config_file_name: Optional[str] = None
    source_code_version_id: Optional[str] = None
    project_directory: Path = Path().home() / "metabolomics_data"
    max_cpus: int = 4
    log_level: str = "INFO"

    def update(self, override_parameters: Dict) -> None:
        """update all parameters with any non-None values in override_parameters"""
        for name, config_value in self.__dict__.items():
            if name in override_parameters:
                override_value = override_parameters[name]
                if override_value is not None:
                    if isinstance(config_value, OutputLists):
                        config_value.update(override_value)
                    else:
                        setattr(self, name, override_value)

    def distribute_always_values(self) -> None:
        """distribute always values within each OutputLists attribute"""
        attribute_names = set(self.__dict__.keys())
        for name in attribute_names:
            attribute_object = getattr(self, name)
            if isinstance(attribute_object, OutputLists):
                attribute_object.distribute_always_values()


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
    export_msms_fragment_ions: bool = False
    slurm_execute: bool = False
    clear_cache: bool = False
    # these are populated from the workflow's RTAlignment if None
    google_folder: Optional[str] = None
    msms_refs: Optional[Path] = None


class RTAlignmentNotebookParameters(BaseNotebookParameters):
    """Parameters used in RT Alignment notebooks"""

    polarity: Polarity = Polarity("positive")
    inchi_keys_not_in_model: List[str] = []
    dependent_data_source: str = "median"
    use_poly_model: bool = False
    stop_before: Optional[str] = None
    google_folder: str
    msms_refs: Path


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
    def altas_name_unique_id_match(cls, to_check, values):
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

    def update(self, override_parameters: Dict) -> None:
        """update all parameters with any non-None values in override_parameters"""
        update_analysis(self, override_parameters)


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

    def update(self, override_parameters: Dict) -> None:
        """update all parameters with any non-None values in override_parameters"""
        update_analysis(self, override_parameters)


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
        if analysis_name == "RT_Alignment":
            return self.rt_alignment
        for analysis in self.analyses:
            if analysis.name == analysis_name:
                return analysis
        raise ValueError(f"Analysis named '{analysis_name}' was not found within the workflow.")

    def update(self, override_parameters: Dict) -> None:
        """update all parameters with any non-None values in override_parameters"""
        self.rt_alignment.update(override_parameters)
        for analysis in self.analyses:
            analysis.parameters.google_folder = or_default(
                analysis.parameters.google_folder, self.rt_alignment.parameters.google_folder
            )
            analysis.parameters.msms_refs = or_default(
                analysis.parameters.msms_refs, self.rt_alignment.parameters.msms_refs
            )
            analysis.update(override_parameters)

    def distribute_always_values(self) -> None:
        """
        distribute always values within each
        OutputLists attribute within parameters
        """
        self.rt_alignment.parameters.distribute_always_values()
        for analysis in self.analyses:
            analysis.parameters.distribute_always_values()


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
        logging.info("Applying override parameters to configuration")
        for flow in self.workflows:
            flow.update(override_parameters)

    def distribute_always_values(self) -> None:
        """
        distribute always values within each
        OutputLists attribute within parameters
        """
        logging.info("Distributing 'always' values across include/exclude dicts")
        for flow in self.workflows:
            flow.distribute_always_values()


def load_config(file_name: os.PathLike) -> Config:
    """Return configuration from file"""
    logging.info("Loading configuration from %s", file_name)
    with open(file_name, "r", encoding="utf-8") as yaml_fh:
        return Config.parse_obj(yaml.safe_load(yaml_fh))


def get_config(override_parameters: Dict) -> Tuple[Config, Workflow, Analysis]:
    """loads configuration from file and updates workflow parameters with override_parameters"""
    config = load_config(override_parameters["config_file_name"])
    config.update(override_parameters)
    config.distribute_always_values()
    workflow_name = override_parameters["workflow_name"]
    workflow = config.get_workflow(override_parameters["workflow_name"])
    if "analysis_name" in override_parameters:
        analysis_name = override_parameters["analysis_name"]
    else:
        analysis_name = "RT_Alignment"
    logging.info("Using workflow/analysis: %s/%s", workflow_name, analysis_name)
    analysis = workflow.get_analysis(analysis_name)
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


def update_analysis(analysis, override_parameters: Dict) -> None:
    """update all parameters with any non-None values in override_parameters"""
    if analysis.parameters.source_atlas_unique_id is not None:
        analysis.atlas.unique_id = analysis.parameters.source_atlas_unique_id
    analysis.parameters.update(override_parameters)
