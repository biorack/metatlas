""" AnalysisIdentifiers object for use with MetatlasDataset """
# pylint: disable=too-many-lines

import getpass
import logging
import os

from typing import cast, Dict, List, NewType, Optional, TypedDict

import pandas as pd
import traitlets

from traitlets import Bool, TraitError, observe, validate
from traitlets import HasTraits, Int, TraitType, Unicode
from traitlets.traitlets import ObserveHandler

from metatlas.datastructures import metatlas_objects as metob
from metatlas.datastructures.utils import AtlasName, get_atlas, Username
from metatlas.io import write_utils
from metatlas.plots import dill2plots as dp
from metatlas.tools.util import or_default

GroupList = Optional[List[metob.Group]]
LcmsRunsList = Optional[List[metob.LcmsRun]]
FileMatchList = List[str]
GroupMatchList = List[str]

Polarity = NewType("Polarity", str)
ShortPolarity = NewType("ShortPolarity", str)
Experiment = NewType("Experiment", str)
OutputType = NewType("OutputType", str)
AnalysisNumber = NewType("AnalysisNumber", int)
PathString = NewType("PathString", str)

DEFAULT_GROUPS_CONTROLLED_VOCAB = cast(GroupMatchList, ["QC", "InjBl", "ISTD"])
OUTPUT_TYPES = [
    OutputType("ISTDsEtc"),
    OutputType("FinalEMA-HILIC"),
    OutputType("FinalEMA-C18"),
    OutputType("data_QC"),
    OutputType("other"),
]
POLARITIES = [Polarity("positive"), Polarity("negative"), Polarity("fast-polarity-switching")]
SHORT_POLARITIES = {
    Polarity("positive"): ShortPolarity("POS"),
    Polarity("negative"): ShortPolarity("NEG"),
    Polarity("fast-polarity-switching"): ShortPolarity("FPS"),
}

logger = logging.getLogger(__name__)


class Proposal(TypedDict):
    """for use with traitlets.validate"""

    owner: HasTraits
    value: object
    trait: TraitType


class _LcmsRunDict(TypedDict):
    """part of return type for AnalysisIdentifiers._files_dict"""

    object: metob.LcmsRun
    group: str
    short_name: str


class AnalysisIdentifiers(HasTraits):
    """Names used in generating an analysis"""

    project_directory: PathString = Unicode(read_only=True)
    experiment: Experiment = Unicode(read_only=True)
    output_type: OutputType = Unicode(read_only=True)
    polarity: Polarity = Unicode(default_value="positive", read_only=True)
    analysis_number: AnalysisNumber = Int(default_value=0, read_only=True)
    google_folder: str = Unicode(read_only=True)

    source_atlas: Optional[AtlasName] = Unicode(allow_none=True, default_value=None, read_only=True)
    copy_atlas: bool = Bool(default_value=True, read_only=True)
    username: Username = Unicode(default_value=getpass.getuser(), read_only=True)
    exclude_files: FileMatchList = traitlets.List(trait=Unicode(), default_value=[], read_only=True)
    include_groups: GroupMatchList = traitlets.List(read_only=True)
    exclude_groups: GroupMatchList = traitlets.List(read_only=True)
    groups_controlled_vocab: GroupMatchList = traitlets.List(
        trait=Unicode(), default_value=DEFAULT_GROUPS_CONTROLLED_VOCAB, read_only=True
    )

    _lcmsruns: LcmsRunsList = traitlets.List(allow_none=True, default_value=None, read_only=True)
    _all_groups: GroupList = traitlets.List(allow_none=True, default_value=None, read_only=True)
    _groups: GroupList = traitlets.List(allow_none=True, default_value=None, read_only=True)

    # pylint: disable=no-self-use,too-many-arguments,too-many-locals
    def __init__(
        self,
        project_directory,
        experiment,
        output_type,
        polarity,
        analysis_number,
        google_folder,
        source_atlas=None,
        copy_atlas=True,
        username=None,
        exclude_files=None,
        include_groups=None,
        exclude_groups=None,
        groups_controlled_vocab=None,
        lcmsruns=None,
        all_groups=None,
    ) -> None:
        super().__init__()
        self.set_trait("project_directory", project_directory)
        self.set_trait("experiment", experiment)
        self.set_trait("output_type", output_type)
        self.set_trait("polarity", polarity)
        self.set_trait("analysis_number", analysis_number)
        self.set_trait("google_folder", google_folder)
        self.set_trait("source_atlas", source_atlas)
        self.set_trait("copy_atlas", copy_atlas)
        self.set_trait("username", or_default(username, getpass.getuser()))
        self.set_trait("exclude_files", or_default(exclude_files, []))
        self.set_trait("include_groups", or_default(include_groups, self._default_include_groups))
        self.set_trait("exclude_groups", or_default(exclude_groups, self._default_exclude_groups))
        self.set_trait(
            "groups_controlled_vocab", or_default(groups_controlled_vocab, DEFAULT_GROUPS_CONTROLLED_VOCAB)
        )
        self.set_trait("_lcmsruns", lcmsruns)
        self.set_trait("_all_groups", all_groups)
        logger.info(
            "IDs: source_atlas=%s, atlas=%s, short_experiment_analysis=%s, output_dir=%s",
            self.source_atlas,
            self.atlas,
            self.short_experiment_analysis,
            self.output_dir,
        )
        self.store_all_groups(exist_ok=True)
        self.set_trait("exclude_groups", append_inverse(self.exclude_groups, self.polarity))

    @property
    def _default_include_groups(self) -> GroupMatchList:
        if self.output_type == "data_QC":
            return ["QC"]
        return []

    def _get_default_exclude_groups(self, polarity: Polarity) -> GroupMatchList:
        out: GroupMatchList = ["InjBl", "InjBL"]
        if self.output_type in ["ISTDsEtc", "FinalEMA-HILIC"]:
            out.append("QC")
        return append_inverse(out, polarity)

    @property
    def _default_exclude_groups(self) -> GroupMatchList:
        return self._get_default_exclude_groups(self.polarity)

    @validate("polarity")
    def _valid_polarity(self, proposal: Proposal) -> Polarity:
        if proposal["value"] not in POLARITIES:
            raise TraitError(f"Parameter polarity must be one of {', '.join(POLARITIES)}")
        return cast(Polarity, proposal["value"])

    @validate("output_type")
    def _valid_output_type(self, proposal: Proposal) -> OutputType:
        if proposal["value"] not in OUTPUT_TYPES:
            raise TraitError(f"Parameter output_type must be one of {', '.join(OUTPUT_TYPES)}")
        return cast(OutputType, proposal["value"])

    @validate("source_atlas")
    def _valid_source_atlas(self, proposal: Proposal) -> Optional[AtlasName]:
        if proposal["value"] is not None:
            proposed_name = cast(AtlasName, proposal["value"])
            try:
                get_atlas(proposed_name, cast(Username, "*"))  # raises error if not found or matches multiple
            except ValueError as err:
                raise TraitError(str(err)) from err
            return proposed_name
        return None

    @validate("analysis_number")
    def _valid_analysis_number(self, proposal: Proposal) -> AnalysisNumber:
        value = cast(AnalysisNumber, proposal["value"])
        if value < 0:
            raise TraitError("Parameter analysis_number cannot be negative.")
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
    def project(self) -> int:
        """Returns project number (proposal id)"""
        return int(self._exp_tokens[3])

    @property
    def atlas(self) -> AtlasName:
        """Atlas identifier (name)"""
        if self.source_atlas is None or (self.copy_atlas and self.source_atlas is not None):
            return AtlasName(
                f"{'_'.join(self._exp_tokens[3:6])}_{self.output_type}_{self.short_polarity}_{self.analysis}"
            )
        return self.source_atlas

    @property
    def analysis(self) -> str:
        """Analysis identifier"""
        return f"{self.username}{self.analysis_number}"

    @property
    def short_experiment_analysis(self) -> str:
        """Short experiment analysis identifier"""
        return f"{self._exp_tokens[0]}_{self._exp_tokens[3]}_{self.output_type}_{self.analysis}"

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
        sub_dirs = [self.experiment, self.analysis, self.output_type]
        if self.output_type in ["ISTDsEtc", "FinalEMA-HILIC"]:
            sub_dirs.append(self.short_polarity)
        out = os.path.join(self.project_directory, *sub_dirs)
        os.makedirs(out, exist_ok=True)
        return PathString(out)

    @property
    def cache_dir(self) -> PathString:
        """Creates directory for storing cache files and returns the path as a string"""
        out = os.path.join(self.project_directory, self.experiment, "cache")
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

    def get_lcmsruns_short_names(self, fields: Optional[Dict[str, List[int]]] = None) -> pd.DataFrame:
        """
        Querys DB for lcms filenames from self.experiment and returns
        a pandas DataFrame containing identifiers for each file
        inputs:
            fields: optional dict with column names as key
                    and list of lcms filename metadata fields positions as value
        """
        if fields is None:
            fields = {
                "full_filename": list(range(16)),
                "sample_treatment": [12],
                "short_filename": [0, 2, 4, 5, 7, 9, 14],
                "short_samplename": [9, 12, 13, 14],
            }
        out = pd.DataFrame(columns=fields.keys())
        for i, lcms_file in enumerate(self.lcmsruns):
            tokens = lcms_file.name.split(".")[0].split("_")
            for name, idxs in fields.items():
                out.loc[i, name] = "_".join([tokens[n] for n in idxs])
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
    def _files_dict(self) -> Dict[str, _LcmsRunDict]:
        """
        Queries DB for all lcmsruns matching the class properties.
        Returns a dict of dicts where keys are filenames minus extensions and values are
        dicts with keys: object, group, and short_name
        """
        file_dict: Dict[str, _LcmsRunDict] = {}
        for lcms_file in self.lcmsruns:
            base_name: str = lcms_file.name.split(".")[0]
            file_dict[base_name] = cast(_LcmsRunDict, {"object": lcms_file, **self.group_name(base_name)})
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
        self.set_trait("_groups", dp.filter_empty_metatlas_objects(out, "items"))
        return self._groups or []

    @observe("polarity")
    def _observe_polarity(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self.set_trait("exclude_groups", self._get_default_exclude_groups(signal.new))
            logger.debug("Change to polarity invalidates exclude_groups")

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
        return metob.retrieve("Groups", name=f"{self.experiment}%{self.analysis}_%", username=self.username)

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
        group_name = f"{prefix}_{self.analysis}_{suffix}"
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
        for values in unique_groups.to_dict("index").values():
            if self._all_groups is not None:  # needed for mypy
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
        return self._all_groups or []

    @property
    def chromatography(self) -> str:
        """returns the type of chromatography used"""
        return self.lcmsruns[0].name.split("_")[7]

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

    def remove_from_exclude_groups(self, remove_groups: List[str]) -> None:
        """Remove items in remove_groups from exclude_groups"""
        self.set_trait("exclude_groups", remove_items(self.exclude_groups, remove_groups))


def append_inverse(in_list: List[str], polarity: Polarity) -> List[str]:
    """appends short version of inverse of polarity to and retuns the list"""
    inverse = {"positive": "NEG", "negative": "POS"}
    return in_list + [inverse[polarity]] if polarity in inverse else in_list


def remove_items(edit_list: List[str], remove_list: List[str], ignore_case: bool = True) -> List[str]:
    """Returns list of items in edit_list but not in remove_list"""
    if ignore_case:
        lower_remove_list = [x.lower() for x in remove_list]
        return [x for x in edit_list if x.lower() not in lower_remove_list]
    return [x for x in edit_list if x not in remove_list]
