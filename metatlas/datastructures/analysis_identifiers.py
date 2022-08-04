""" AnalysisIdentifiers object for use with MetatlasDataset """
# pylint: disable=too-many-lines

import getpass
import logging
import os

from functools import lru_cache
from pathlib import Path
from typing import cast, Dict, List, Optional, Union

import pandas as pd
import traitlets

from traitlets import observe, validate, Bool, HasTraits, Int, Instance, TraitError, Unicode
from traitlets.traitlets import ObserveHandler

from metatlas.datastructures.id_types import (
    IterationNumber,
    Experiment,
    FileMatchList,
    GroupList,
    GroupMatchList,
    LcmsRunDict,
    LcmsRunsList,
    PathString,
    Polarity,
    POLARITIES,
    Proposal,
    ShortPolarity,
    SHORT_POLARITIES,
)
import metatlas.datastructures.metatlas_objects as metob
import metatlas.plots.dill2plots as dp

from metatlas.datastructures.utils import AtlasName, get_atlas, Username
from metatlas.io import write_utils
from metatlas.tools.config import Config
from metatlas.tools.util import or_default


logger = logging.getLogger(__name__)

MSMS_REFS_PATH = PathString("/global/cfs/cdirs/metatlas/projects/spectral_libraries/msms_refs_v3.tab")


class AnalysisIdentifiers(HasTraits):
    """Names used in generating an analysis"""

    project_directory: PathString = Unicode(read_only=True)
    experiment: Experiment = Unicode(read_only=True)
    polarity: Polarity = Unicode(default_value="positive", read_only=True)
    analysis_number: IterationNumber = Int(default_value=0, read_only=True)
    rt_alignment_number: IterationNumber = Int(default_value=0, read_only=True)
    google_folder: str = Unicode(read_only=True)
    source_atlas_unique_id: str = Unicode(read_only=True)
    copy_atlas: bool = Bool(default_value=True, read_only=True)
    username: Username = Unicode(default_value=getpass.getuser(), read_only=True)
    exclude_files: FileMatchList = traitlets.List(trait=Unicode(), default_value=[], read_only=True)
    include_groups: Optional[GroupMatchList] = traitlets.List(
        allow_none=True, default_value=None, read_only=True
    )
    exclude_groups: GroupMatchList = traitlets.List(default_value=[], read_only=True)
    groups_controlled_vocab: GroupMatchList = traitlets.List(
        trait=Unicode(), default_value=[], read_only=True
    )
    configuration: Config = Instance(klass=Config)
    workflow: str = Unicode(read_only=True)
    analysis: str = Unicode(read_only=True)
    _lcmsruns: LcmsRunsList = traitlets.List(allow_none=True, default_value=None, read_only=True)
    _all_groups: GroupList = traitlets.List(allow_none=True, default_value=None, read_only=True)
    _groups: GroupList = traitlets.List(allow_none=True, default_value=None, read_only=True)

    # pylint: disable=no-self-use,too-many-arguments,too-many-locals
    def __init__(
        self,
        project_directory,
        experiment,
        analysis_number,
        google_folder,
        polarity=None,
        source_atlas_unique_id=None,
        copy_atlas=True,
        username=None,
        exclude_files=None,
        include_groups=None,
        exclude_groups=None,
        groups_controlled_vocab=None,
        lcmsruns=None,
        all_groups=None,
        rt_alignment_number=0,
        configuration=None,
        workflow=None,
        analysis="RT_Alignment",
    ) -> None:
        super().__init__()
        self.set_trait("project_directory", project_directory)
        self.set_trait("experiment", experiment)
        self.set_trait("polarity", Polarity(or_default(polarity, "positive")))
        self.set_trait("analysis_number", analysis_number)
        self.set_trait("rt_alignment_number", rt_alignment_number)
        self.set_trait("google_folder", google_folder)
        self.set_trait("source_atlas_unique_id", source_atlas_unique_id)
        self.set_trait("copy_atlas", copy_atlas)
        self.set_trait("username", or_default(username, getpass.getuser()))
        self.set_trait("exclude_files", or_default(exclude_files, []))
        self.set_trait("include_groups", include_groups)
        self.set_trait("exclude_groups", or_default(exclude_groups, []))
        self.set_trait("groups_controlled_vocab", or_default(groups_controlled_vocab, []))
        self.set_trait("_lcmsruns", lcmsruns)
        self.set_trait("_all_groups", all_groups)
        self.set_trait("configuration", configuration)
        self.set_trait("workflow", workflow)
        self.set_trait("analysis", analysis)
        logger.info(
            "IDs: source_atlas_unique_id=%s, atlas=%s, output_dir=%s",
            self.source_atlas_unique_id,
            self.atlas,
            self.output_dir,
        )
        self.store_all_groups(exist_ok=True)

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
                get_atlas(proposal["value"])
            except ValueError as err:
                raise TraitError(str(err)) from err
            return proposal["value"]
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
    def output_dir(self) -> PathString:
        """Creates the output directory and returns the path as a string"""
        sub_dirs = [
            self.experiment_id,
            self.workflow,
            str(self.rt_alignment_number),
            str(self.analysis_number),
            "Targeted",
            self.workflow,
            self.analysis,
        ]
        out = os.path.join(self.project_directory, *sub_dirs)
        os.makedirs(out, exist_ok=True)
        return PathString(out)

    @property
    def notebook_dir(self) -> PathString:
        """Directoy where notebooks are saved"""
        return PathString(str(Path(self.output_dir).resolve().parent.parent.parent.parent))

    @property
    def cache_dir(self) -> PathString:
        """Creates directory for storing cache files and returns the path as a string"""
        out = os.path.join(self.project_directory, self.experiment_id, "cache")
        os.makedirs(out, exist_ok=True)
        return PathString(out)

    @property
    def lcmsruns(self) -> List[metob.LcmsRun]:
        """Get LCMS runs from DB matching experiment"""
        if self._lcmsruns is not None:
            return self._lcmsruns
        all_lcmsruns = dp.get_metatlas_files(experiment=self.experiment, name="%")
        if self.exclude_files is not None and len(self.exclude_files) > 0:
            self.set_trait(
                "_lcmsruns",
                [
                    r
                    for r in all_lcmsruns
                    if not any(map(r.name.__contains__, or_default(self.exclude_files, [])))
                ],
            )
            if self._lcmsruns:
                logger.info(
                    "Excluding %d LCMS runs containing any of: %s",
                    len(all_lcmsruns) - len(self._lcmsruns),
                    self.exclude_files,
                )
        else:
            self.set_trait("_lcmsruns", all_lcmsruns)
        if self._lcmsruns:
            for run in self._lcmsruns:
                logger.info("Run: %s", run.name)
        logger.info(
            "Number of LCMS output files matching '%s' is: %d.", self.experiment, len(self._lcmsruns or [])
        )
        return self._lcmsruns or []

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
            os.path.join(self.output_dir, "short_names.csv"),
            "LCMS runs short names",
            index=False,
        )

    @property
    def _files_dict(self) -> Dict[str, LcmsRunDict]:
        """
        Queries DB for all lcmsruns matching the class properties.
        Returns a dict of dicts where keys are filenames minus extensions and values are
        dicts with keys: object, group, and short_name
        """
        file_dict: Dict[str, LcmsRunDict] = {}
        for lcms_file in self.lcmsruns:
            base_name: str = lcms_file.name.split(".")[0]
            file_dict[base_name] = cast(LcmsRunDict, {"object": lcms_file, **self.group_name(base_name)})
        return file_dict

    @property
    def groups(self) -> List[metob.Group]:
        """Return the currently selected groups"""
        if self._groups is not None:
            return self._groups
        out = dp.filter_metatlas_objects_to_most_recent(self.all_groups, "name")
        if self.include_groups is not None and len(self.include_groups) > 0:
            out = dp.filter_metatlas_objects_by_list(out, "name", self.include_groups)
        if self.exclude_groups is not None and len(self.exclude_groups) > 0:
            out = dp.remove_metatlas_objects_by_list(out, "name", self.exclude_groups)
        sorted_out = sorted(dp.filter_empty_metatlas_objects(out, "items"), key=lambda x: x.name)
        self.set_trait("_groups", sorted_out)
        return self._groups

    @observe("_all_groups")
    def _observe_all_groups(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self.set_trait("_groups", None)
            logger.debug("Change to all_groups invalidates groups")

    @observe("groups_controlled_vocab")
    def _observe_groups_controlled_vocab(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self.set_trait("_lcmsruns", None)
            logger.debug("Change to groups_controlled_vocab invalidates lcmsruns")

    @observe("include_groups")
    def _observe_include_groups(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self.set_trait("_groups", None)
            logger.debug("Change to include_groups invalidates groups")

    @observe("exclude_groups")
    def _observe_exclude_groups(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self.set_trait("_groups", None)
            logger.debug("Change to exclude_groups invalidates groups")

    @observe("exclude_files")
    def _observe_exclude_files(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self.set_trait("_lcmsruns", None)
            logger.debug("Change to exclude_files invalidates lcmsruns")

    @observe("_lcmsruns")
    def _observe_lcmsruns(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self.set_trait("_all_groups", None)
            logger.debug("Change to lcmsruns invalidates all_groups")

    @property
    def existing_groups(self) -> List[metob.Group]:
        """Get your own groups that are prefixed by self.experiment"""
        return metob.retrieve("Groups", name=f"{self.experiment}%{self.execution}_%", username=self.username)

    def group_name(self, base_filename: str) -> Dict[str, str]:
        """Returns dict with keys group and short_name corresponding to base_filename"""
        tokens = base_filename.split("_")
        prefix = "_".join(tokens[:11])
        indices = [
            i
            for i, s in enumerate(or_default(self.groups_controlled_vocab, []))
            if s.lower() in base_filename.lower()
        ]
        suffix = self.groups_controlled_vocab[indices[0]].lstrip("_") if indices else tokens[12]
        group_name = f"{prefix}_{self.execution}_{suffix}"
        short_name = f"{tokens[9]}_{suffix}"  # Prepending POL to short_name
        return {"group": group_name, "short_name": short_name}

    @property
    def all_groups_dataframe(self) -> pd.DataFrame:
        """Returns pandas Dataframe with one row per file"""
        out = pd.DataFrame(self._files_dict).T
        if out.empty:
            return out
        out.drop(columns=["object"], inplace=True)
        out.index.name = "filename"
        return out.reset_index()

    @property
    def all_groups(self) -> List[metob.Group]:
        """Returns a list of Group objects"""
        if self._all_groups is not None:
            return self._all_groups
        unique_groups = self.all_groups_dataframe[["group", "short_name"]].drop_duplicates()
        self.set_trait("_all_groups", [])
        assert self._all_groups is not None  # needed for mypy
        for values in unique_groups.to_dict("index").values():
            self._all_groups.append(
                metob.Group(
                    name=values["group"],
                    short_name=values["short_name"],
                    items=[
                        file_value["object"]
                        for file_value in self._files_dict.values()
                        if file_value["group"] == values["group"]
                    ],
                )
            )
        return sorted(self._all_groups, key=lambda x: x.name)

    @property
    def chromatography(self) -> str:
        """returns the type of chromatography used"""
        alternatives = {t.name: t.aliases for t in self.configuration.chromatography_types}
        chrom_field = self.lcmsruns[0].name.split("_")[7]
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
            new_names = set(self.all_groups_dataframe["group"].to_list())
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
