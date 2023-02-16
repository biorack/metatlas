""" AnalysisIdentifiers object for use with MetatlasDataset """
# pylint: disable=too-many-lines

import getpass
import logging

from functools import lru_cache
from pathlib import Path
from typing import cast, Dict, List, Optional, Union

import pandas as pd
import traitlets

from traitlets import observe, validate, Callable, Bool, HasTraits, Int, Instance, TraitError, Unicode
from traitlets.traitlets import ObserveHandler

from IPython.display import display, HTML

from metatlas.datastructures.id_types import (
    IterationNumber,
    Experiment,
    FileMatchList,
    GroupList,
    GroupMatchList,
    LcmsRunDict,
    LcmsRunsList,
    Polarity,
    POLARITIES,
    Proposal,
    ShortPolarity,
    SHORT_POLARITIES,
)
import metatlas.datastructures.metatlas_objects as metob
import metatlas.plots.dill2plots as dp

from metatlas.datastructures.utils import AtlasName, get_atlas, Username
from metatlas.datastructures import groups
from metatlas.io import write_utils
from metatlas.tools.config import BaseNotebookParameters, Config, OutputLists
from metatlas.tools.util import or_default


logger = logging.getLogger(__name__)

MSMS_REFS_PATH = Path("/global/cfs/cdirs/metatlas/projects/spectral_libraries/msms_refs_v3.tab")


class AnalysisIdentifiers(HasTraits):
    """Names used in generating an analysis"""

    project_directory: Path = Instance(klass=Path, read_only=True)
    experiment: Experiment = Unicode(read_only=True)
    polarity: Polarity = Unicode(default_value="positive", read_only=True)
    analysis_number: IterationNumber = Int(default_value=0, read_only=True)
    rt_alignment_number: IterationNumber = Int(default_value=0, read_only=True)
    google_folder: str = Unicode(read_only=True)
    source_atlas_unique_id: str = Unicode(read_only=True)
    copy_atlas: bool = Bool(default_value=True, read_only=True)
    username: Username = Unicode(default_value=getpass.getuser(), read_only=True)
    include_lcmsruns: FileMatchList = traitlets.List(trait=Unicode(), default_value=[])
    exclude_lcmsruns: FileMatchList = traitlets.List(trait=Unicode(), default_value=[])
    include_groups: GroupMatchList = traitlets.List(trait=Unicode(), default_value=[])
    exclude_groups: GroupMatchList = traitlets.List(trait=Unicode(), default_value=[])
    groups_controlled_vocab: GroupMatchList = traitlets.List(
        trait=Unicode(), default_value=[], read_only=True
    )
    configuration: Config = Instance(klass=Config)
    workflow: str = Unicode(read_only=True)
    analysis: str = Unicode(read_only=True)
    all_lcmsruns: Optional[LcmsRunsList] = traitlets.List(allow_none=True, default_value=None, read_only=True)
    _lcmsruns: Optional[LcmsRunsList] = traitlets.List(allow_none=True, default_value=None, read_only=True)
    _all_groups: Optional[GroupList] = traitlets.List(allow_none=True, default_value=None, read_only=True)
    _groups: Optional[GroupList] = traitlets.List(allow_none=True, default_value=None, read_only=True)
    groups_invalidation_callbacks: List[Callable] = traitlets.List(trait=Callable(), default_value=[])

    # pylint: disable=too-many-arguments,too-many-locals
    def __init__(
        self,
        project_directory,
        experiment,
        configuration,
        workflow,
        analysis="RT_Alignment",
        rt_alignment_number=0,
        analysis_number=0,
        source_atlas_unique_id=None,
        username=None,
        lcmsruns=None,
        groups=None,
    ) -> None:
        super().__init__()
        analysis_obj = configuration.get_workflow(workflow).get_analysis(analysis)
        self.set_trait("project_directory", Path(project_directory))
        self.set_trait("experiment", experiment)
        self.set_trait("configuration", configuration)
        self.set_trait("workflow", workflow)
        self.set_trait("analysis", analysis)
        self.set_trait("rt_alignment_number", rt_alignment_number)
        self.set_trait("analysis_number", analysis_number)
        self.set_trait(
            "source_atlas_unique_id", or_default(source_atlas_unique_id, analysis_obj.atlas.unique_id)
        )
        self.set_trait("username", or_default(username, getpass.getuser()))
        self.set_trait("groups_controlled_vocab", analysis_obj.parameters.groups_controlled_vocab)
        self.set_trait("all_lcmsruns", self._get_lcmsruns(lcmsruns))
        self.set_trait("_all_groups", groups)
        self.set_trait("polarity", Polarity(analysis_obj.parameters.polarity))
        self.set_trait("google_folder", analysis_obj.parameters.google_folder)
        self.set_trait("copy_atlas", analysis_obj.parameters.copy_atlas)
        logger.info(
            "IDs: source_atlas_unique_id=%s, atlas=%s, output_dir=%s",
            self.source_atlas_unique_id,
            self.atlas,
            self.output_dir,
        )
        self.store_all_groups(exist_ok=True)
        self.set_output_state(analysis_obj.parameters, "gui")

    @validate("polarity")
    def _valid_polarity(self, proposal: Proposal) -> Polarity:
        if proposal["value"] not in POLARITIES:
            raise TraitError(f"Parameter polarity must be one of {', '.join(POLARITIES)}")
        return cast(Polarity, proposal["value"])

    @validate("source_atlas_unique_id")
    def _valid_source_atlas_unique_id(self, proposal: Proposal) -> Optional[str]:
        if proposal["value"] is not None:
            try:
                # raises error if not found or matches multiple
                get_atlas(str(proposal["value"]))
            except ValueError as err:
                raise TraitError(str(err)) from err
            return str(proposal["value"])
        return None

    @validate("analysis_number")
    def _valid_analysis_number(self, proposal: Proposal) -> IterationNumber:
        value = cast(IterationNumber, proposal["value"])
        if value < 0:
            raise TraitError("Parameter analysis_number cannot be negative.")
        return value

    @validate("rt_alignment_number")
    def _valid_rt_alignment_number(self, proposal: Proposal) -> IterationNumber:
        value = cast(IterationNumber, proposal["value"])
        if value < 0:
            raise TraitError("Parameter rt_alignment_number cannot be negative.")
        return value

    @validate("experiment")
    def _valid_experiment(self, proposal: Proposal) -> Experiment:
        value = cast(str, proposal["value"])
        num_fields = len(value.split("_"))
        error_msg = "Parameter 'experiment' should contains 9 fields when split on '_', but has %d."
        if num_fields > 9:
            logger.warning(error_msg, num_fields)
        if num_fields < 9:
            raise TraitError(error_msg, num_fields)
        return cast(Experiment, value)

    @property
    def _exp_tokens(self) -> List[str]:
        """Returns list of strings from the experiment name"""
        return self.experiment.split("_")

    @property
    def project(self) -> str:
        """
        Returns project identifier (proposal id)
        This is an integer for JGI, but a string for EGSB
        """
        return self._exp_tokens[3]

    @property
    def exp(self) -> str:
        """Returns the exp field of the experiment"""
        return self._exp_tokens[4]

    @property
    def sample_set(self) -> str:
        """Returns the sample set field of the experiment"""
        return self._exp_tokens[5]

    @property
    def experiment_id(self) -> str:
        """Returns a unique ID for an experiment"""
        return f"{self.project}_{self.exp}_{self.sample_set}"

    @property
    @lru_cache
    def source_atlas(self) -> AtlasName:
        """Returns name of atlas"""
        atlas = get_atlas(self.source_atlas_unique_id)
        return atlas.name

    @property
    def atlas(self) -> AtlasName:
        """Atlas identifier (name)"""
        if self.copy_atlas:
            return AtlasName(f"{self.experiment_id}_{self.source_atlas}_{self.execution}")
        return self.source_atlas

    @property
    def execution(self) -> str:
        """execution identifier"""
        return f"{self.username}_{self.rt_alignment_number}_{self.analysis_number}"

    @property
    def short_polarity(self) -> ShortPolarity:
        """Short polarity identifier: 3 letters, upper case"""
        return SHORT_POLARITIES[self.polarity]

    @property
    def short_polarity_inverse(self) -> List[ShortPolarity]:
        """Returns the short_polarity values not used in this analysis"""
        return list(set(SHORT_POLARITIES.values()) - {self.short_polarity})

    @property
    def output_dir(self) -> Path:
        """Creates the output directory and returns the path as a string"""
        sub_dirs = [
            self.experiment_id,
            f"{self.username}_{self.workflow}_{self.rt_alignment_number}_{self.analysis_number}",
            "Targeted",
            f"{self.workflow}_{self.experiment}",
            self.analysis,
        ]
        out = self.project_directory.joinpath(*sub_dirs)
        out.mkdir(parents=True, exist_ok=True)
        return out

    @property
    def additional_output_dir(self) -> Path:
        """Creates the "Additional_Info" output directory"""
        out = self.output_dir / "Additional_Info"
        out.mkdir(parents=True, exist_ok=True)
        return out

    @property
    def notebook_dir(self) -> Path:
        """Directoy where notebooks are saved"""
        return self.output_dir.resolve().parent.parent.parent.parent

    @property
    def cache_dir(self) -> Path:
        """Creates directory for storing cache files and returns the path"""
        out = self.project_directory / self.experiment_id / "cache"
        out.mkdir(parents=True, exist_ok=True)
        return out

    @property
    def lcmsruns(self) -> List[metob.LcmsRun]:
        """Get LCMS runs from DB matching experiment"""
        return groups.get_groups_and_runs(
            self.execution,
            self.groups_controlled_vocab,
            or_default(self.include_lcmsruns, []),
            or_default(self.exclude_lcmsruns, []),
            or_default(self.include_groups, []),
            or_default(self.exclude_groups, []),
            self.all_lcmsruns,
            self._all_groups,
        )[1]

    def _get_lcmsruns(self, all_lcmsruns: Optional[List[metob.LcmsRun]] = None) -> List[metob.LcmsRun]:
        """Get the set of lcmsruns that are currently selected"""
        return or_default(all_lcmsruns, dp.get_metatlas_files(experiment=self.experiment, name="%"))

    @property
    def lcmsruns_dataframe(self) -> pd.DataFrame:
        """Returns a pandas DataFrame with lcmsrun matching self.experiment"""
        return metob.to_dataframe(self.lcmsruns)

    def get_lcmsruns_short_names(
        self, fields: Optional[Dict[str, Union[List[int], str]]] = None
    ) -> pd.DataFrame:
        """
        Querys DB for lcms filenames from self.experiment and returns
        a pandas DataFrame containing identifiers for each file
        inputs:
            fields: optional dict with column names as key
                    and list of lcms filename metadata fields positions or 'all' as value
        """
        if fields is None:
            fields = {
                "full_filename": "all",
                "sample_treatment": [12],
                "short_filename": [0, 2, 4, 5, 7, 9, 14],
                "short_samplename": [9, 12, 13, 14],
            }
        out = pd.DataFrame(columns=fields.keys())
        for i, lcms_file in enumerate(self.lcmsruns):
            stem = lcms_file.name.split(".")[0]
            tokens = stem.split("_")
            for name, idxs in fields.items():
                out.loc[i, name] = stem if idxs == "all" else "_".join([tokens[n] for n in idxs])
            out.loc[i, "last_modified"] = pd.to_datetime(lcms_file.last_modified, unit="s")
        if out.empty:
            return out
        out.sort_values(by="last_modified", inplace=True)
        out.drop(columns=["last_modified"], inplace=True)
        out.drop_duplicates(subset=["full_filename"], keep="last", inplace=True)
        out.set_index("full_filename", inplace=True)
        return out.sort_values(by="full_filename")

    lcmsruns_short_names: pd.DataFrame = property(get_lcmsruns_short_names)

    def write_lcmsruns_short_names(self) -> None:
        """Write short names and raise error if exists and differs from current data"""
        short_names = self.lcmsruns_short_names
        short_names["full_filename"] = short_names.index
        write_utils.export_dataframe_die_on_diff(
            short_names,
            self.additional_output_dir / "short_names.csv",
            "LCMS runs short names",
            index=False,
        )

    @property
    def groups(self) -> GroupList:
        """Return the currently selected groups"""
        return groups.get_groups_and_runs(
            self.execution,
            self.groups_controlled_vocab,
            or_default(self.include_lcmsruns, []),
            or_default(self.exclude_lcmsruns, []),
            or_default(self.include_groups, []),
            or_default(self.exclude_groups, []),
            self.all_lcmsruns,
            self._all_groups,
        )[0]

    @property
    def all_groups(self) -> GroupList:
        """Return the all groups"""
        return groups.get_groups_and_runs(
            self.execution, self.groups_controlled_vocab, [], [], [], [], self.all_lcmsruns, self._all_groups
        )[0]

    @observe("include_lcmsruns")
    def _observe_include_lcmsruns(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._exec_group_invalidation_callbacks("include_lcmsruns")

    @observe("exclude_lcmsruns")
    def _observe_exclude_lcmsruns(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._exec_group_invalidation_callbacks("exclude_lcmsruns")

    @observe("include_groups")
    def _observe_include_groups(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._exec_group_invalidation_callbacks("include_groups")

    @observe("exclude_groups")
    def _observe_exclude_groups(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._exec_group_invalidation_callbacks("exclude_groups")

    def _exec_group_invalidation_callbacks(self, source: str) -> None:
        for callback in self.groups_invalidation_callbacks:
            logger.debug(
                "Change to %s triggers invalidation callback function %s()", source, callback.__name__
            )
            callback()

    @property
    def existing_groups(self) -> List[metob.Group]:
        """Get your own groups that are prefixed by self.experiment"""
        return metob.retrieve("Groups", name=f"{self.experiment}%{self.execution}_%", username=self.username)

    @property
    def chromatography(self) -> str:
        """returns the type of chromatography used"""
        alternatives = {t.name: t.aliases for t in self.configuration.chromatography_types}
        chrom_field = self.all_lcmsruns[0].name.split("_")[7]
        chrom_type = chrom_field.split("-")[0]
        for name, alt_list in alternatives.items():
            if chrom_type == name or chrom_type in alt_list:
                return name
        logger.warning("Unknown chromatography field '%s'.", chrom_type)
        return chrom_type

    def store_all_groups(self, exist_ok: bool = False) -> None:
        """
        Save self.object_list to DB
        inputs:
            exist_ok: if False, store nothing and raise ValueError if any of the group names
                      have already been saved to the DB by you.
        """
        if not exist_ok:
            db_names = {group.name for group in self.existing_groups}
            new_names = {group.name for group in self.all_groups}
            overlap = db_names.intersection(new_names)
            try:
                if overlap:
                    raise ValueError(
                        "Not saving groups as you have already saved groups"
                        f'with these names: {", ".join(overlap)}.'
                    )
            except ValueError as err:
                logger.exception(err)
                raise err
        logger.debug("Storing %d groups in the database", len(self.all_groups))
        metob.store(self.all_groups)

    def set_output_state(self, parameters: BaseNotebookParameters, state: str):
        """set include/exclude lcmsruns/groups for the given state"""
        assert state in vars(OutputLists())
        for out_list in ["include_lcmsruns", "exclude_lcmsruns", "include_groups", "exclude_groups"]:
            current = getattr(self, out_list)
            new = getattr(getattr(parameters, out_list), state)
            if new != current:
                setattr(self, out_list, new)
                logger.info("Setting %s list to %s.", out_list, new)

    def register_groups_invalidation_callback(self, func: Callable) -> None:
        """Register a function to call when groups get invalidated"""
        self.groups_invalidation_callbacks.append(func)

    def display_groups(self):
        current_groups = self.groups
        total = len(self.all_groups)
        num = len(current_groups)
        display(HTML(f"<H4>Groups: {num} selected and displayed (total: {total})</H4>"))
        display(metob.to_dataframe(current_groups)[["name", "short_name"]])

    def display_lcmsruns(self):
        current_runs = self.lcmsruns
        total = len(self.all_lcmsruns)
        num = len(current_runs)
        display(HTML(f"<H4>LCMS runs: {num} selected and displayed (total: {total})</H4>"))
        display(metob.to_dataframe(current_runs)[["name"]])
