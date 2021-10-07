""" object oriented interface to metatlas_dataset """
# pylint: disable=too-many-lines

import datetime
import getpass
import glob
import logging
import os
import pickle
import shutil
import uuid

from pathlib import Path
from typing import cast, Any, Dict, List, NewType, Optional, Tuple, TypedDict, Union

import humanize
import pandas as pd
import traitlets

from traitlets import TraitError, default, observe, validate
from traitlets import Bool, Float, HasTraits, Instance, Int, TraitType, Unicode
from traitlets.traitlets import ObserveHandler

from metatlas.datastructures import metatlas_objects as metob
from metatlas.datastructures import object_helpers as metoh
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.io import targeted_output
from metatlas.io import write_utils
from metatlas.plots import dill2plots as dp
from metatlas.tools import parallel
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
AtlasName = NewType("AtlasName", str)
PathString = NewType("PathString", str)
Username = NewType("Username", str)

MSMS_REFS_PATH = PathString(
    "/global/project/projectdirs/metatlas/projects/spectral_libraries/msms_refs_v3.tab"
)
DEFAULT_GROUPS_CONTROLLED_VOCAB = cast(GroupMatchList, ["QC", "InjBl", "ISTD"])
OUTPUT_TYPES = [OutputType("ISTDsEtc"), OutputType("FinalEMA-HILIC"), OutputType("data_QC")]
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
    """part of return type for AnalysisIds._files_dict"""

    object: metob.LcmsRun
    group: str
    short_name: str


class MsSummary(TypedDict):
    """part of MetatlasDataset._data"""

    num_ms1_datapoints: int
    mz_peak: float
    rt_peak: float
    mz_centroid: float
    rt_centroid: float
    peak_height: float
    peak_area: float


class Eic(TypedDict):
    """part of MetatlasDataset._data"""

    mz: List[float]
    rt: List[float]
    intensity: List[float]


class MsmsDataDict(TypedDict):
    """part of MetatlasDataset._data"""

    mz: List[float]
    i: List[float]
    rt: List[float]
    polarity: List[float]
    precursor_MZ: List[float]
    precursor_intensity: List[float]
    collision_energy: List[float]


class MsmsDict(TypedDict):
    """part of MetatlasDataset._data"""

    data: MsmsDataDict


class MsDataDict(TypedDict):
    """part of MetatlasDataset._data"""

    msms: MsmsDict
    eic: Eic
    ms1_summary: MsSummary


class CompoundDict(TypedDict):
    """part of MetatlasDataset._data"""

    atlas_name: AtlasName
    atlas_unique_id: str
    lcmsrun: metob.LcmsRun
    group: metob.Group
    identification: metob.CompoundIdentification
    data: MsDataDict


class AnalysisIdentifiers(HasTraits):
    """Names used in generating an analysis"""

    source_atlas: Optional[AtlasName] = Unicode(allow_none=True, default_value=None)
    experiment: Experiment = Unicode()
    output_type: OutputType = Unicode()
    polarity: Polarity = Unicode(default_value="positive")
    analysis_number: AnalysisNumber = Int(default_value=0)
    username: Username = Unicode(default_value=getpass.getuser())
    project_directory: PathString = Unicode()
    google_folder: str = Unicode()
    exclude_files: FileMatchList = traitlets.List(trait=Unicode(), default_value=[])
    include_groups: GroupMatchList = traitlets.List()
    exclude_groups: GroupMatchList = traitlets.List()
    groups_controlled_vocab: GroupMatchList = traitlets.List(
        trait=Unicode(), default_value=DEFAULT_GROUPS_CONTROLLED_VOCAB
    )
    _lcmsruns: LcmsRunsList = traitlets.List(allow_none=True, default_value=None)
    _all_groups: GroupList = traitlets.List(allow_none=True, default_value=None)
    _groups: GroupList = traitlets.List(allow_none=True, default_value=None)

    # pylint: disable=no-self-use
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.exclude_groups = append_inverse(self.exclude_groups, self.polarity)
        logger.info(
            "IDs: source_atlas=%s, atlas=%s, short_experiment_analysis=%s, output_dir=%s",
            self.source_atlas,
            self.atlas,
            self.short_experiment_analysis,
            self.output_dir,
        )
        self.store_all_groups(exist_ok=True)

    @default("include_groups")
    def _default_include_groups(self) -> List[OutputType]:
        if self.output_type == "data_QC":
            return [OutputType("QC")]
        return []

    @default("exclude_groups")
    def _default_exclude_groups(self) -> GroupMatchList:
        out: GroupMatchList = ["InjBl", "InjBL"]
        if self.output_type != "data_QC":
            out.append("QC")
        return append_inverse(out, self.polarity)

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
        if len(value.split("_")) != 9:
            raise TraitError('Parameter experiment does contain 9 fields when split on "_".')
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
        return AtlasName(
            f"{'_'.join(self._exp_tokens[3:6])}_{self.output_type}_{self.short_polarity}_{self.analysis}"
        )

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
        if self.output_type != "data_QC":
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
            self._lcmsruns = [
                r
                for r in all_lcmsruns
                if not any(map(r.name.__contains__, or_default(self.exclude_files, [])))
            ]
            logger.info(
                "Excluding %d LCMS runs containing any of: %s",
                len(all_lcmsruns) - len(self._lcmsruns),
                self.exclude_files,
            )
        else:
            self._lcmsruns = all_lcmsruns
        for run in self._lcmsruns:
            logger.info("Run: %s", run.name)
        logger.info("Number of LCMS output files matching '%s' is: %d.", self.experiment, len(self._lcmsruns))
        return self._lcmsruns

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
        self._groups = dp.filter_empty_metatlas_objects(out, "items")
        return self._groups

    @observe("polarity")
    def _observe_polarity(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self.exclude_groups = append_inverse(self.exclude_groups, signal.new)
            logger.debug("Change to polarity invalidates exclude_groups")

    @observe("_all_groups")
    def _observe_all_groups(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._groups = None
            logger.debug("Change to all_groups invalidates groups")

    @observe("groups_controlled_vocab")
    def _observe_groups_controlled_vocab(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._lcmsruns = None
            logger.debug("Change to groups_controlled_vocab invalidates lcmsruns")

    @observe("include_groups")
    def _observe_include_groups(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._groups = None
            logger.debug("Change to include_groups invalidates groups")

    @observe("exclude_groups")
    def _observe_exclude_groups(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._groups = None
            logger.debug("Change to exclude_groups invalidates groups")

    @observe("exclude_files")
    def _observe_exclude_files(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._lcmsruns = None
            logger.debug("Change to exclude_files invalidates lcmsruns")

    @observe("_lcmsruns")
    def _observe_lcmsruns(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._all_groups = None
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
        self._all_groups = []
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
        return self._all_groups

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


class MetatlasSample:
    """
    Object oriented interface to second level of metatlas_dataset. Each instance is one sample (LCMS run).
    """

    def __init__(self, compounds: Tuple[CompoundDict, ...]) -> None:
        self.compounds: Tuple[CompoundDict, ...] = compounds

    def __getitem__(self, idx: int) -> CompoundDict:
        """get sample at idx"""
        return self.compounds[idx]

    def __len__(self) -> int:
        """len is from data"""
        return len(self.compounds)


class MetatlasDataset(HasTraits):
    """
    Like the non-object oriented metatlas_dataset, you can index into this class by file_idx and compound_idx:
    metatlas_dataset = MetatlasDataset(analysis_ids)
    metatlas_dataset[0][0]['identification'].compound[0].inchi_key

    But MetatlasDataset adds additional functionality, such as:
    metatlas_dataset.hits returns the msms hits dataframe
    metatlas_dataset.atlas returns the atlas
    metatlas_dataset.atlas_df returns the atlas dataframe

    If you change a class property that another property is dependent on, then the second property
    automatically is re-calculated the next time you access the second property. For example:
    metatlas_dataset.extra_time = 0.5  # this invalidates the current hits property
    metatlas_dataset.hits  # this re-generates the hits before returning them

    MetatlasDataset also has methods for updating RT values and identification notes while keeping
    the atlas, atlas_df, metatlas_dataset, and database in sync. This removes the need to do kernel
    restarts between steps in the workflow.


    ids: AnalysisIdentifiers instance defining the analysis
    save_metadata: if True, write metadata files containing data sources and LCMS runs short name
    """

    extra_time: float = Float(default_value=0.75)
    extra_mz: float = Float(default_value=0)
    frag_mz_tolerance: float = Float(default_value=0.01)
    max_cpus: int = Int(default_value=1)
    save_metadata: bool = Bool(default_value=True)
    keep_nonmatches: bool = Bool(default_value=True)
    msms_refs_loc: PathString = Unicode(default_value=MSMS_REFS_PATH)
    ids: AnalysisIdentifiers = Instance(klass=AnalysisIdentifiers)
    atlas: metob.Atlas = Instance(klass=metob.Atlas)
    _atlas_df: Optional[pd.DataFrame] = Instance(klass=pd.DataFrame, allow_none=True, default_value=None)
    _data: Optional[Tuple[MetatlasSample, ...]] = traitlets.Tuple(allow_none=True, default_value=None)
    _hits: Optional[pd.DataFrame] = Instance(klass=pd.DataFrame, allow_none=True, default_value=None)

    # pylint: disable=too-many-instance-attributes, too-many-arguments, too-many-public-methods, no-self-use
    def __init__(self, **kwargs) -> None:
        """Constructor"""
        super().__init__(**kwargs)
        logger.debug("Creating new MetatlasDataset instance...")
        self._hits_valid_for_rt_bounds = False  # based only on RT min/max changes
        self._data_valid_for_rt_bounds = False  # based only on RT min/max changes
        if self.ids.source_atlas is not None:
            self._get_atlas()
        if self.save_metadata:
            logger.debug("Writing MetatlasDataset metadata files")
            self.write_data_source_files()
            self.ids.write_lcmsruns_short_names()

    def _save_to_cache(self, data: Any, metadata: dict) -> None:
        assert "_variable_name" in metadata.keys()
        metadata = metadata.copy()
        name = metadata["_variable_name"]
        base_name = f"{name}_{uuid.uuid4()}"
        metadata["_pickle_file_name"] = f"{base_name}.pkl"  # relative to metadata file
        metadata_file_name = os.path.join(self.ids.cache_dir, f"{base_name}.metadata")
        pickle_path = os.path.join(self.ids.cache_dir, metadata["_pickle_file_name"])
        with open(pickle_path, "wb") as pickle_fh:
            pickle.dump(data, pickle_fh)
        with open(metadata_file_name, "wb") as metadata_fh:
            pickle.dump(metadata, metadata_fh)
        logger.info("Caching %s in %s.", name, pickle_path)

    def _query_cache(self, required_metadata: dict) -> Optional[Any]:
        assert "_variable_name" in required_metadata.keys()
        name = required_metadata["_variable_name"]
        files_matching_metadata = []
        for metadata_file in Path(self.ids.cache_dir).glob(f"{name}_*.metadata"):
            with open(metadata_file, "rb") as metadata_fh:
                potential_metadata = pickle.load(metadata_fh)
                pickle_file_name = Path(self.ids.cache_dir) / potential_metadata["_pickle_file_name"]
                # require_metadata does not have a '_pickle_file_name' key, so remove before equality test
                del potential_metadata["_pickle_file_name"]
                if required_metadata == potential_metadata:
                    files_matching_metadata.append(pickle_file_name)
        if len(files_matching_metadata) == 0:
            return None
        newest_file = sorted(files_matching_metadata, key=lambda x: x.stat().st_mtime)[-1]
        with open(newest_file, "rb") as pickle_fh:
            logger.info("Loading cached %s from %s.", name, pickle_file_name)
            return pickle.load(pickle_fh)

    def write_data_source_files(self) -> None:
        """Write the data source files if they don't already exist"""
        data_sources_dir = os.path.join(self.ids.output_dir, f"{self.ids.short_polarity}_data_sources")
        if len(glob.glob(os.path.join(data_sources_dir, "*"))) >= 4:
            logger.warning(
                (
                    "Data sources directory already populated from previous work on this analysis. "
                    "Not overwriting."
                )
            )
        else:
            shutil.rmtree(data_sources_dir, ignore_errors=True)
            logger.info("Writing data source files to %s.", data_sources_dir)
            ma_data.make_data_sources_tables(
                self.ids.groups, self.atlas, self.ids.output_dir, self.ids.short_polarity
            )

    def _get_atlas(self) -> None:
        """
        Copy source atlas from database into current analysis atlas
        If the atlas does not yet exist, it will be copied from source_atlas and there will be an
        an additional side effect that all mz_tolerances in the resulting atlas
        get their value from source_atlas' atlas.compound_identifications[0].mz_references[0].mz_tolerance
        """
        atlases = metob.retrieve("Atlas", name=self.ids.atlas, username=self.ids.username)
        if len(atlases) == 1:
            logger.warning(
                (
                    "Destination atlas, %s, already exists, so not copying source atlas, "
                    "%s, to destination. Not overwriting."
                ),
                self.ids.atlas,
                self.ids.source_atlas,
            )
            self.atlas = atlases[0]
        elif len(atlases) > 1:
            try:
                raise ValueError(
                    (
                        f"{len(atlases)} atlases with name {self.ids.atlas} "
                        f"and owned by {self.ids.username} already exist."
                    )
                )
            except ValueError as err:
                logger.exception(err)
                raise err
        elif self.ids.source_atlas is not None:
            self.atlas = self._clone_source_atlas()
        else:
            try:
                raise ValueError("Could not load atlas as source_atlas is None.")
            except ValueError as err:
                logger.exception(err)
                raise err

    def _clone_source_atlas(self) -> metob.Atlas:
        logger.info("Retriving source atlas: %s", self.ids.source_atlas)
        source_atlas = get_atlas(cast(AtlasName, self.ids.source_atlas), cast(Username, "*"))
        source_atlas_df = ma_data.make_atlas_df(source_atlas)
        logger.info("Cloning atlas %s", self.ids.source_atlas)
        return dp.make_atlas_from_spreadsheet(
            source_atlas_df,
            self.ids.atlas,
            filetype="dataframe",
            sheetname="",
            polarity=self.ids.polarity,
            store=True,
            mz_tolerance=source_atlas.compound_identifications[0].mz_references[0].mz_tolerance,
        )

    def _build(self) -> None:
        """Populate self._data from database and h5 files."""
        start_time = datetime.datetime.now()
        files = [
            (h5_file, group, self.atlas_df, self.atlas, self.extra_time, self.extra_mz)
            for group in self.ids.groups
            for h5_file in group.items
        ]
        try:
            if len(files) == 0:
                raise ValueError("No matching h5 files were found")
        except ValueError as err:
            logger.exception(err)
            raise err
        logger.info("Generating MetatlasDataset by reading MSMS data from h5 files")
        samples = parallel.parallel_process(
            ma_data.get_data_for_atlas_df_and_file, files, self.max_cpus, unit="sample", spread_args=False
        )
        self._data = tuple(MetatlasSample(x) for x in samples)
        logger.info(
            "MetatlasDataset with %d files built in %s.",
            len(files),
            _duration_since(start_time),
        )

    def _remove_compound_id(self, idx: int) -> None:
        """
        Remove compound identification at index idx from both in db and self.atlas
        Does not invalidate _data or _hits or _atlas_df
        This bypasses several ORM layers and therefore is a hack, but I couldn't get it to work with the ORM.
        """
        cid_id = self.atlas.compound_identifications[idx].unique_id
        del self.atlas.compound_identifications[idx]
        atlas_id = self.atlas.unique_id
        link_table = "atlases_compound_identifications"
        target = f"target_id='{cid_id}'"
        workspace = metob.Workspace.get_instance()
        workspace.get_connection()
        workspace.db.begin()
        try:
            workspace.db.query(f"delete from {link_table} where ({target} and source_id='{atlas_id}')")
            links = workspace.db.query(f"select source_id from {link_table} where {target}")
            if len(list(links)) == 0:  # other atlases are not linked to this CompoundIdentification
                workspace.db.query(f"delete from compoundidentifications where unique_id='{cid_id}'")
            workspace.db.commit()
        except Exception as err:  # pylint: disable=broad-except
            metoh.rollback_and_log(workspace.db, err)
            raise Exception from err
        workspace.close_connection()

    def filter_compounds(
        self, keep_idxs: Optional[List[int]] = None, remove_idxs: Optional[List[int]] = None
    ) -> None:
        """
        inputs:
            keep_idxs: the indexes of compounds to keep
            remove_idxs: the indexes of compounds to remove
                Exactly one of keep_idxs or remove_idxs must be None
        output:
            If keep_idxs is not None then update self.atlas to contain only the compound_identifications at
            keep_idxs. If remove_idxs is not None then update self.atlas to contain only the compound
            identifications not at remove_idxs. Raises ValueError if both keep_idxs and remove_idxs are None.

            Does not invalidate _data or _hits
        """
        if (keep_idxs is None) == (remove_idxs is None):
            raise ValueError("Exactly one of keep_idxs and remove_idxs should be None")
        _error_if_bad_idxs(self.atlas_df, keep_idxs or remove_idxs)
        start_len = len(self.atlas_df)
        in_idxs: List[int] = keep_idxs or self.atlas_df.index.difference(remove_idxs)
        if len(in_idxs) == start_len:
            return
        out_idxs: List[int] = remove_idxs or self.atlas_df.index.difference(keep_idxs)
        self._atlas_df = self.atlas_df.iloc[in_idxs].copy().reset_index(drop=True)
        if self._data is not None:
            self._data = tuple(
                MetatlasSample(
                    tuple(compound for idx, compound in enumerate(sample.compounds) if idx in in_idxs)
                )
                for sample in self._data
            )
        for i in sorted(out_idxs, reverse=True):
            self._remove_compound_id(i)
        logger.info(
            "Filtering reduced atlas from %d to %d compounds (%d removed).",
            start_len,
            len(self.atlas_df),
            start_len - len(self.atlas_df),
        )
        if self._hits is not None:
            self.filter_hits_by_atlas()

    def filter_hits_by_atlas(self) -> None:
        """Remove any hits that do not have a corresponding inchi_key-adduct pair in atlas_df"""
        start_len = len(self.hits)
        keep_adducts = self.atlas_df.loc[:, ["inchi_key", "adduct"]].drop_duplicates()
        logger.info("Number of inchi_key-adduct pairs is %d.", len(keep_adducts))
        hits_plus = self.hits.copy()
        hits_plus["copy_index"] = hits_plus.index
        new_hits = hits_plus.merge(keep_adducts, on=["inchi_key", "adduct"], how="inner")
        logger.info("Number rows in new_hits is %d.", len(new_hits))
        new_hits.index = pd.MultiIndex.from_tuples(new_hits["copy_index"], names=self.hits.index.names)
        new_hits.drop(["copy_index"], axis=1)
        self._hits = new_hits
        logger.info(
            "Filtering reduced number of MSMS hits from %d to %d (%d removed).",
            start_len,
            len(self.hits),
            start_len - len(self.hits),
        )

    def filter_compounds_ms1_notes_remove(self) -> None:
        """
        output:
            updates self.atlas to contain only the compound_identifications that do not have ms1_notes
            starting with 'remove' (case insensitive)
            There is an additional side effect that all mz_tolerances in the returned atlas
            get their value from self.atlas.compound_identifications[0].mz_references[0].mz_tolerance
        """
        logger.debug("Filtering atlas to exclude ms1_notes=='remove'.")
        self.filter_compounds(remove_idxs=self.compound_indices_marked_remove())

    def filter_compounds_by_signal(self, num_points: int, peak_height: float) -> None:
        """
        inputs:
            num_points: number of points in EIC that must be exceeded in one or more samples
                        in order for the compound to remain in the atlas
            peak_height: max intensity in the EIC that must be exceeded in one or more samples
                         in order for the compound to remain in the atlas
        """
        logger.debug("Filtering atlas on num_points=%d, peak_height=%d.", num_points, peak_height)
        keep_idxs = dp.strong_signal_compound_idxs(self, num_points, peak_height)
        self.filter_compounds(keep_idxs=keep_idxs)

    def store_atlas(self, even_if_exists: bool = False) -> None:
        """
        inputs:
            even_if_exists: if True, will save the atlas even if the atlas name already is in the database
                            with your username
        side effects:
            Saves the altas to the database.
            Raises ValueError if even_if_exists==False and name is already in the database with your username
        """
        start_time = datetime.datetime.now()
        name = self.atlas.name
        username = self.ids.username
        try:
            if not even_if_exists and len(metob.retrieve("Atlas", name=name, username=username)) > 0:
                raise ValueError(f"An atlas with name {name} and owned by {username} already exists.")
        except ValueError as err:
            logger.exception(err)
            raise err
        metob.store(self.atlas)
        logger.info(
            "Atlas %s stored in database with owner %s in %s.",
            self.ids.atlas,
            self.ids.username,
            _duration_since(start_time),
        )

    def export_atlas_to_csv(self, filename: Optional[str] = None) -> None:
        """
        save atlas, including ms1_notes, ms2_notes, identification_notes, rt_min, rt_max to filename
        if filename is not provided, then the export is saved to the working directory with filename
        atlas.name + '.csv'
        """
        filename = f"{self.atlas.name}.csv" if filename is None else filename
        dp.export_atlas_to_spreadsheet(self, filename)

    def __getitem__(self, idx: int) -> MetatlasSample:
        """get sample at idx"""
        return self.data[idx]

    @property
    def data(self) -> Tuple[MetatlasSample, ...]:
        """data getter, update ._data if necessary"""
        if self._data is None:
            self._build()
            self._data_valid_for_rt_bounds = True
        return cast(Tuple[MetatlasSample, ...], self._data)

    @property
    def atlas_df(self) -> pd.DataFrame:
        """atlas_df getter, update ._atlas_df if necessary"""
        if self._atlas_df is None:
            start_time = datetime.datetime.now()
            logger.info("Generating atlas_df")
            self._atlas_df = ma_data.make_atlas_df(self.atlas)
            logger.info(
                "Generated atlas_df with %d rows in %s.",
                len(self.atlas_df),
                _duration_since(start_time),
            )
        return self._atlas_df

    @observe("atlas")
    def _observe_atlas(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._atlas_df = None
            self._data = None
            logger.debug("Change to atlas invalidates atlas_df, data")

    @observe("_atlas_df")
    def _observe_atlas_df(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._data = None
            logger.debug("Change to atlas_df invalidates data")

    @property
    def polarity(self) -> Polarity:
        """
        polarity getter assumes all polarities within class are the same
        returns 'positive' if there are no samples or no compound identifications
        """
        try:
            cid = self.data[0][0]["identification"]
        except IndexError:
            return Polarity("positive")
        return Polarity(cid.mz_references[0].detected_polarity)

    @observe("extra_time")
    def _observe_extra_time(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._hits = None
            self._data = None
            logger.debug("Change to extra_time invalidates hits, data")

    @observe("extra_mz")
    def _observe_extra_mz(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._hits = None
            self._data = None
            logger.debug("Change to extra_mz invalidates hits, data")

    @observe("keep_nonmatches")
    def _observe_keep_nonmatches(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._hits = None
            logger.debug("Change to keep_nonmatches invalidates hits")

    @observe("frag_mz_tolerance")
    def _observe_frag_mz_tolerance(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._hits = None
            logger.debug("Change to frag_mz_tolerance invalidates hits")

    @observe("msms_refs_loc")
    def _observe_msms_refs_loc(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._hits = None
            logger.debug("Change to msms_refs_loc invalidates hits")

    @property
    def hits(self) -> pd.DataFrame:
        """get msms hits DataFrame"""
        if self._hits:
            return self._hits
        metadata = self._get_hits_metadata()
        self._hits = self._query_cache(metadata)
        if self._hits:
            self._hits_valid_for_rt_bounds = False  # unsure, so assume False
            return self._hits
        _ = self.atlas_df  # regenerate if needed before logging hits generation
        _ = self.data  # regenerate if needed before logging hits generation
        logger.info(
            "Generating hits with extra_time=%.3f, frag_mz_tolerance=%.4f, msms_refs_loc=%s.",
            self.extra_time,
            self.frag_mz_tolerance,
            self.msms_refs_loc,
        )
        start_time = datetime.datetime.now()
        self._hits = dp.get_msms_hits(
            self.data,
            extra_time=self.extra_time > 0,
            keep_nonmatches=self.keep_nonmatches,
            frag_mz_tolerance=self.frag_mz_tolerance,
            ref_loc=self.msms_refs_loc,
        )
        logger.info("Generated %d hits in %s.", len(self._hits), _duration_since(start_time))
        self._hits_valid_for_rt_bounds = True
        self._save_to_cache(self._hits, metadata)
        return self._hits

    def _get_hits_metadata(self) -> Dict[str, Any]:
        return {
            "_variable_name": "hits",
            "polarity": self.ids.polarity,
            "extra_time": self.extra_time,
            "keep_nonmatches": self.keep_nonmatches,
            "frag_mz_tolerance": self.frag_mz_tolerance,
            "ref_loc": self.msms_refs_loc,
            "extra_mz": self.extra_mz,
            "output_type": self.ids.output_type,
        }

    def __len__(self) -> int:
        """len is from data"""
        return len(self.data)

    @property
    def rts(self) -> Tuple[metob.RtReference, ...]:
        """
        Allow Rt_Reference objects to be accessed
        use set_rt() if you want to modify the RT values held by this class.
        """
        if self.atlas is None:
            return tuple()  # noqa: C408
        return tuple(cid.rt_references[0] for cid in self.atlas.compound_identifications)

    def set_rt(self, compound_idx: int, which: str, time: float) -> None:
        """
        inputs:
            compound_idx: index of of compound to update
            which: 'rt_min', 'rt_max', or 'rt_peak'
            time: a floating point value for the number of minutes
        updates the RT value in database, self.atlas, self.atlas_df, self.data
        so that no datastructures need to be invalidated
        """
        try:
            if self.atlas is None:
                raise ValueError("Cannot set RTs when atlas is None.")
        except ValueError as err:
            logger.exception(err)
            raise err
        assert which in ["rt_min", "rt_peak", "rt_max"]
        atlas_rt_ref = self.atlas.compound_identifications[compound_idx].rt_references[0]
        setattr(atlas_rt_ref, which, time)
        for sample in self.data:
            setattr(sample[compound_idx]["identification"].rt_references[0], which, time)
        self.atlas_df.loc[compound_idx, which] = time
        metob.store(atlas_rt_ref)
        if which in ["rt_min", "rt_max"]:
            self._hits_valid_for_rt_bounds = False
            self._data_valid_for_rt_bounds = False

    def set_note(self, compound_idx: int, which: str, value: str) -> None:
        """
        inputs:
            compound_idx: index of of compound to update
            which: 'ms1_notes', 'ms2_notes' or 'identification_notes'
            value: a string with the note content
        updates the notes value in database, self.atlas, self.atlas_df, self.data
        so that no datastructures need to be invalidated
        """
        try:
            if self.atlas is None:
                raise ValueError("Cannot set notes when atlas is None.")
        except ValueError as err:
            logger.exception(err)
            raise err
        assert which in ["ms1_notes", "ms2_notes", "identification_notes"]
        atlas_cid = self.atlas.compound_identifications[compound_idx]
        setattr(atlas_cid, which, value)
        data_cid = self.data[0][compound_idx]["identification"]
        setattr(data_cid, which, value)
        self.atlas_df.loc[compound_idx, which] = value
        metob.store(atlas_cid)

    def compound_indices_marked_remove(self) -> List[int]:
        """
        outputs:
            list of compound_idx of the compound identifications with ms1_notes to remove
        """
        ids = ["identification", "ms1_notes"]
        return [i for i, j in enumerate(self.data[0].compounds) if _is_remove(ma_data.extract(j, ids))]

    def compound_idxs_not_evaluated(self) -> List[int]:
        """
        Returns list of compound indices where ms1 note is not 'remove' and
        ms2 note is None or 'no selection'
        """
        out = []
        for i, compound in enumerate(self.data[0].compounds):
            ms1_note = ma_data.extract(compound, ["identification", "ms1_notes"])
            ms2_note = ma_data.extract(compound, ["identification", "ms2_notes"])
            if (not _is_remove(ms1_note)) and (not _has_selection(ms2_note)):
                out.append(i)
        return out

    def error_if_not_all_evaluated(self) -> None:
        """Raises ValueError if there are compounds that have not been evaluated"""
        not_evaluated = self.compound_idxs_not_evaluated()
        try:
            if len(not_evaluated) != 0:
                raise ValueError(
                    (
                        "Compounds with the following indices need notes selected via radio "
                        f"buttons before continuing: {','.join([str(i) for i in not_evaluated])}"
                    )
                )
        except ValueError as err:
            logger.exception(err)
            raise err

    def annotation_gui(
        self, compound_idx: int = 0, width: float = 15, height: float = 3, alpha: float = 0.5, colors=""
    ) -> dp.adjust_rt_for_selected_compound:
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
        return dp.adjust_rt_for_selected_compound(
            self,
            msms_hits=self.hits,
            color_me=colors,
            compound_idx=compound_idx,
            alpha=alpha,
            width=width,
            height=height,
        )

    def generate_all_outputs(self, msms_fragment_ions: bool = False, overwrite: bool = False) -> None:
        """
        Generates the default set of outputs for a targeted experiment
        inputs:
            msms_fragment_ions: if True, generate msms fragment ions report
            overwrite: if False, throw error if any output files already exist
        """
        if not self._hits_valid_for_rt_bounds:
            self._hits = None  # force hits to be regenerated
        if not self._data_valid_for_rt_bounds:
            self._data = None  # force data to be regenerated
        self.extra_time = 0.5
        logger.info("extra_time set to 0.5 minutes for output generation.")
        logger.info("Removing InjBl from exclude_groups.")
        self.ids.exclude_groups = remove_items(self.ids.exclude_groups, ["InjBl"])
        targeted_output.write_atlas_to_spreadsheet(self, overwrite=overwrite)
        targeted_output.write_stats_table(self, overwrite=overwrite)
        targeted_output.write_chromatograms(self, overwrite=overwrite, max_cpus=self.max_cpus)
        targeted_output.write_identification_figure(self, overwrite=overwrite)
        targeted_output.write_metrics_and_boxplots(self, overwrite=overwrite, max_cpus=self.max_cpus)
        if msms_fragment_ions:
            targeted_output.write_msms_fragment_ions(self, overwrite=overwrite)
        logger.info("Generation of output files completed sucessfully.")
        targeted_output.archive_outputs(self.ids)
        targeted_output.copy_outputs_to_google_drive(self.ids)


def _duration_since(start: datetime.datetime) -> str:
    """
    inputs:
        start: a datetime object of when the duration started
    returns:
        string with humanized duration of start to now
    """
    return humanize.precisedelta(datetime.datetime.now() - start)


def _is_remove(obj: object) -> bool:
    """is obj a string that starts with 'remove' (case insensitive)?"""
    return isinstance(obj, str) and obj.lower().startswith("remove")


def _has_selection(obj: object) -> bool:
    """is obj a string that is not None, '', or 'no selection' (case insensitive)?"""
    if obj is None or not isinstance(obj, str):
        return False
    return obj.lower() not in ["", "no selection"]


def _set_nested(data: Any, ids: List[Union[int, str, Tuple[str]]], value: Any):
    """
    inputs:
        data: hierarchical data structure consisting of lists, dicts, and objects with attributes.
        ids: a list of idices, key names, and attribute names
        value: object
    output:
        modifies data in place so that the value is stored at the location indicated by the ids list

    Strings in ids are first tried as key name and if no such key name exists, then they are
    tried as attribute names. To designate that a member of ids should be used as an attribute
    and not a key name, make it a tuple with the attribute name string as the first member, such
    as: ('attribute_name',). If you want to make it more explict to the reader, you can add a
    second member to the tuple, which will not be used, such as ('attribute_name', 'as attribute')
    """
    try:
        if len(ids) == 0:
            raise ValueError("ids cannot be empty")
    except ValueError as err:
        logger.exception(err)
        raise err
    if len(ids) == 1:
        if isinstance(ids[0], tuple):
            setattr(data, ids[0][0], value)
        elif isinstance(ids[0], str) and hasattr(data, ids[0]):
            setattr(data, ids[0], value)
        else:
            data[ids[0]] = value  # works for list or dict
    else:
        if isinstance(ids[0], tuple):
            _set_nested(getattr(data, ids[0][0]), ids[1:], value)
        elif isinstance(ids[0], str) and hasattr(data, ids[0]):
            _set_nested(getattr(data, ids[0]), ids[1:], value)
        else:
            _set_nested(data[ids[0]], ids[1:], value)


def _error_if_bad_idxs(dataframe: pd.DataFrame, test_idx_list: List[int]) -> None:
    """Raise IndexError if any members of of test_idx_list are not in dataframe's index"""
    bad = set(test_idx_list) - set(dataframe.index)
    try:
        if len(bad) > 0:
            raise IndexError(f"Invalid index values: {bad}.")
    except IndexError as err:
        logger.exception(err)
        raise err


def get_atlas(name: AtlasName, username: Username) -> metob.Atlas:
    """Load atlas from database"""
    atlases = metob.retrieve("Atlas", name=name, username=username)
    try:
        if len(atlases) == 0:
            raise ValueError(f"Database does not contain an atlas {name} owned by {username}.")
    except ValueError as err:
        logger.exception(err)
        raise err
    try:
        if len(atlases) > 1:
            raise ValueError(f"Database contains more than one atlas {name} owned by {username}.")
    except ValueError as err:
        logger.exception(err)
        raise err
    return atlases[0]


def quoted_string_list(strings: List[str]) -> str:
    """Adds double quotes around each string and seperates with ', '."""
    return ", ".join([f'"{x}"' for x in strings])


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


# pylint: disable=too-many-arguments
def pre_annotation(
    source_atlas: AtlasName,
    experiment: Experiment,
    output_type: OutputType,
    polarity: Polarity,
    analysis_number: AnalysisNumber,
    project_directory: PathString,
    google_folder: str,
    groups_controlled_vocab: GroupMatchList,
    exclude_files: FileMatchList,
    num_points: int,
    peak_height: float,
    max_cpus: int,
    username: Username = None,
) -> MetatlasDataset:
    """All data processing that needs to occur before the annotation GUI in Targeted notebook"""
    ids = AnalysisIdentifiers(
        source_atlas=source_atlas,
        experiment=experiment,
        output_type=output_type,
        polarity=polarity,
        analysis_number=analysis_number,
        project_directory=project_directory,
        google_folder=google_folder,
        groups_controlled_vocab=groups_controlled_vocab,
        exclude_files=exclude_files,
        username=getpass.getuser() if username is None else username,
    )
    metatlas_dataset = MetatlasDataset(ids=ids, max_cpus=max_cpus)
    if metatlas_dataset.ids.output_type in ["FinalEMA-HILIC"]:
        metatlas_dataset.filter_compounds_by_signal(num_points=num_points, peak_height=peak_height)
    return metatlas_dataset


def post_annotation(metatlas_dataset: MetatlasDataset) -> None:
    """All data processing that needs to occur after the annotation GUI in Targeted notebook"""
    if metatlas_dataset.ids.output_type in ["FinalEMA-HILIC"]:
        metatlas_dataset.error_if_not_all_evaluated()
        metatlas_dataset.filter_compounds_ms1_notes_remove()
    metatlas_dataset.generate_all_outputs()
    logger.info("DONE - execution of notebook is complete.")
