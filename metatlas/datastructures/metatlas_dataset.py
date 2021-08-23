""" object oriented interface to metatlas_dataset """
import datetime
import getpass
import glob
import logging
import os
import shutil

import humanize
import pandas as pd

from traitlets import HasTraits, TraitError, default, observe, validate
from traitlets import Bool, Float, Instance, Int, List, Tuple, Unicode

from metatlas.datastructures import metatlas_objects as metob
from metatlas.datastructures import object_helpers as metoh
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.io import targeted_output
from metatlas.io import write_utils
from metatlas.plots import dill2plots as dp
from metatlas.tools import parallel

MSMS_REFS_PATH = "/global/project/projectdirs/metatlas/projects/spectral_libraries/msms_refs_v3.tab"
DEFAULT_GROUPS_CONTROLLED_VOCAB = ["QC", "InjBl", "ISTD"]
DEFAULT_EXCLUDE_GROUPS = ["InjBl", "InjBL", "QC"]
OUTPUT_TYPES = ["ISTDsEtc", "FinalEMA-HILIC", "data_QC"]
POLARITIES = ["positive", "negative", "fast-polarity-switching"]
SHORT_POLARITIES = {"positive": "POS", "negative": "NEG", "fast-polarity-switching": "FPS"}

logger = logging.getLogger(__name__)


class AnalysisIdentifiers(HasTraits):
    """Names used in generating an analysis"""

    source_atlas = Unicode(allow_none=True, default_value=None)
    experiment = Unicode()
    output_type = Unicode()
    polarity = Unicode(default_value="positive")
    analysis_number = Int(default_value=0)
    username = Unicode(default_value=getpass.getuser())
    project_directory = Unicode()
    google_folder = Unicode()
    exclude_files = List(trait=Unicode(), allow_none=True, default_value=[])
    include_groups = List(allow_none=True, default_value=None)
    exclude_groups = List(allow_none=True, default_value=DEFAULT_EXCLUDE_GROUPS)
    groups_controlled_vocab = List(
        trait=Unicode(), allow_none=True, default_value=DEFAULT_GROUPS_CONTROLLED_VOCAB
    )
    _lcmsruns = List(allow_none=True, default_value=None)
    _all_groups = List(allow_none=True, default_value=None)
    _groups = List(allow_none=True, default_value=None)

    # pylint: disable=no-self-use
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.polarity == "positive":
            self.exclude_groups.append("NEG")
        elif self.polarity == "negative":
            self.exclude_groups.append("POS")
        logger.info(
            "IDs: source_atlas=%s, atlas=%s, short_experiment_analysis=%s, output_dir=%s",
            self.source_atlas,
            self.atlas,
            self.short_experiment_analysis,
            self.output_dir,
        )
        self.store_all_groups(exist_ok=True)

    @default("include_groups")
    def _default_include_groups(self):
        if self.output_type == "data_QC":
            return ["QC"]
        return []

    @validate("polarity")
    def _valid_polarity(self, proposal):
        if proposal["value"] not in POLARITIES:
            raise TraitError(f"Parameter polarity must be one of {', '.join(POLARITIES)}")
        return proposal["value"]

    @validate("output_type")
    def _valid_output_type(self, proposal):
        if proposal["value"] not in OUTPUT_TYPES:
            raise TraitError(f"Parameter output_type must be one of {', '.join(OUTPUT_TYPES)}")
        return proposal["value"]

    @validate("source_atlas")
    def _valid_source_atlas(self, proposal):
        if proposal["value"] is not None:
            try:
                get_atlas(proposal["value"], self.username)  # raises error if not found or matches multiple
            except ValueError as err:
                raise TraitError(str(err)) from err
        return proposal["value"]

    @validate("analysis_number")
    def _valid_analysis_number(self, proposal):
        if proposal["value"] < 0:
            raise TraitError("Parameter analysis_number cannot be negative.")
        return proposal["value"]

    @validate("experiment")
    def _valid_experiment(self, proposal):
        if len(proposal["value"].split("_")) != 9:
            raise TraitError('Parameter experiment does contain 9 fields when split on "_".')
        return proposal["value"]

    @property
    def _exp_tokens(self):
        """Returns list of strings from the experiment name"""
        return self.experiment.split("_")

    @property
    def project(self):
        """Returns project number (proposal id)"""
        return self._exp_tokens[3]

    @property
    def atlas(self):
        """Atlas identifier (name)"""
        return f"{'_'.join(self._exp_tokens[3:6])}_{self.output_type}_{self.short_polarity}_{self.analysis}"

    @property
    def analysis(self):
        """Analysis identifier"""
        return f"{self.username}{self.analysis_number}"

    @property
    def short_experiment_analysis(self):
        """Short experiment analysis identifier"""
        return f"{self._exp_tokens[0]}_{self._exp_tokens[3]}_{self.output_type}_{self.analysis}"

    @property
    def short_polarity(self):
        """Short polarity identifier: 3 letters, upper case"""
        return SHORT_POLARITIES[self.polarity]

    @property
    def short_polarity_inverse(self):
        """Returns the short_polarity values not used in this analysis"""
        return list(set(SHORT_POLARITIES.values()) - {self.short_polarity})

    @property
    def output_dir(self):
        """Creates the output directory and returns the path as a string"""
        out = os.path.join(self.project_directory, self.experiment, self.analysis, self.output_type)
        os.makedirs(out, exist_ok=True)
        return out

    @property
    def lcmsruns(self):
        """Get LCMS runs from DB matching experiment"""
        if self._lcmsruns is not None:
            return self._lcmsruns
        all_lcmsruns = dp.get_metatlas_files(experiment=self.experiment, name="%")
        if len(self.exclude_files) > 0:
            self._lcmsruns = [
                r for r in all_lcmsruns if not any(map(r.name.__contains__, self.exclude_files))
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
    def lcmsruns_dataframe(self):
        """Returns a pandas DataFrame with lcmsrun matching self.experiment"""
        return metob.to_dataframe(self.lcmsruns)

    def get_lcmsruns_short_names(self, fields=None):
        """
        Querys DB for lcms filenames from self.experiment and returns
        a pandas DataFrame containing identifiers for each file
        inputs:
            fields: optional dict with column names as key
                    and list of lcms filename metadata fields positions as value
        """
        if fields is None:
            fields = {
                "full_filename": range(16),
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

    lcmsruns_short_names = property(get_lcmsruns_short_names)

    def write_lcmsruns_short_names(self):
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
    def _files_dict(self):
        """
        Queries DB for all lcmsruns matching the class properties.
        Returns a dict of dicts where keys are filenames minus extensions and values are
        dicts with keys: object, group, and short_name
        """
        file_dict = {}
        for lcms_file in self.lcmsruns:
            base_name = lcms_file.name.split(".")[0]
            file_dict[base_name] = {"object": lcms_file, **self.group_name(base_name)}
        return file_dict

    @property
    def groups(self):
        """Return the currently selected groups"""
        if self._groups is not None:
            return self._groups
        out = dp.filter_metatlas_objects_to_most_recent(self.all_groups, "name")
        if len(self.include_groups) > 0:
            out = dp.filter_metatlas_objects_by_list(out, "name", self.include_groups)
        if len(self.exclude_groups) > 0:
            out = dp.remove_metatlas_objects_by_list(out, "name", self.exclude_groups)
        self._groups = dp.filter_empty_metatlas_objects(out, "items")
        return self._groups

    @observe("polarity")
    def _observe_polarity(self, signal):
        if signal.type == "change":
            if signal.new == "positive":
                self.exclude_groups.append("NEG")
            elif signal.new == "negative":
                self.exclude_groups.append("POS")
            logger.debug("Change to polarity invalidates exclude_groups")

    @observe("_all_groups")
    def _observe_all_groups(self, signal):
        if signal.type == "change":
            self._groups = None
            logger.debug("Change to all_groups invalidates groups")

    @observe("groups_controlled_vocab")
    def _observe_groups_controlled_vocab(self, signal):
        if signal.type == "change":
            self._lcmsruns = None
            logger.debug("Change to groups_controlled_vocab invalidates lcmsruns")

    @observe("include_groups")
    def _observe_include_groups(self, signal):
        if signal.type == "change":
            self._groups = None
            logger.debug("Change to include_groups invalidates groups")

    @observe("exclude_groups")
    def _observe_exclude_groups(self, signal):
        if signal.type == "change":
            self._groups = None
            logger.debug("Change to exclude_groups invalidates groups")

    @observe("exclude_files")
    def _observe_exclude_files(self, signal):
        if signal.type == "change":
            self._lcmsruns = None
            logger.debug("Change to exclude_files invalidates lcmsruns")

    @observe("_lcmsruns")
    def _observe_lcmsruns(self, signal):
        if signal.type == "change":
            self._all_groups = None
            logger.debug("Change to lcmsruns invalidates all_groups")

    @property
    def existing_groups(self):
        """Get your own groups that are prefixed by self.experiment"""
        return metob.retrieve("Groups", name=f"{self.experiment}%{self.analysis}_%", username=self.username)

    def group_name(self, base_filename):
        """Returns dict with keys group and short_name corresponding to base_filename"""
        tokens = base_filename.split("_")
        prefix = "_".join(tokens[:11])
        indices = [
            i for i, s in enumerate(self.groups_controlled_vocab) if s.lower() in base_filename.lower()
        ]
        suffix = self.groups_controlled_vocab[indices[0]].lstrip("_") if indices else tokens[12]
        group_name = f"{prefix}_{self.analysis}_{suffix}"
        short_name = f"{tokens[9]}_{suffix}"  # Prepending POL to short_name
        return {"group": group_name, "short_name": short_name}

    @property
    def all_groups_dataframe(self):
        """Returns pandas Dataframe with one row per file"""
        out = pd.DataFrame(self._files_dict).T
        if out.empty:
            return out
        out.drop(columns=["object"], inplace=True)
        out.index.name = "filename"
        return out.reset_index()

    @property
    def all_groups(self):
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

    def store_all_groups(self, exist_ok=False):
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
                        "Not saving groups as you have already saved groups with these names: %s."
                        % ", ".join(overlap),
                    )
            except ValueError as err:
                logger.exception(err)
                raise err
        logger.debug("Storing %d groups in the database", len(self.all_groups))
        metob.store(self.all_groups)


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

    extra_time = Float(default_value=0.75)
    extra_mz = Float(default_value=0)
    frag_mz_tolerance = Float(default_value=0.01)
    max_cpus = Int(default_value=1)
    save_metadata = Bool(default_value=True)
    keep_nonmatches = Bool(default_value=True)
    msms_refs_loc = Unicode(default_value=MSMS_REFS_PATH)
    ids = Instance(klass=AnalysisIdentifiers)
    atlas = Instance(klass=metob.Atlas, allow_none=True, default_value=None)
    _atlas_df = Instance(klass=pd.DataFrame, allow_none=True, default_value=None)
    _data = Tuple(allow_none=True, default_value=None)
    _hits = Instance(klass=pd.DataFrame, allow_none=True, default_value=None)

    # pylint: disable=too-many-instance-attributes, too-many-arguments, too-many-public-methods, no-self-use
    def __init__(self, **kwargs):
        """Constructor"""
        super().__init__(**kwargs)
        logger.debug("Creating new MetatlasDataset instance...")
        self._hits_valid_for_rt_bounds = False  # based only on RT min/max changes
        if self.ids.source_atlas is not None:
            self._get_atlas()
        if self.save_metadata:
            logger.debug("Writing MetatlasDataset metadata files")
            self.write_data_source_files()
            self.ids.write_lcmsruns_short_names()

    def write_data_source_files(self):
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

    def _get_atlas(self):
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
        else:
            logger.info("Retriving source atlas: %s", self.ids.source_atlas)
            source_atlas = get_atlas(self.ids.source_atlas, self.ids.username)
            source_atlas_df = ma_data.make_atlas_df(source_atlas)
            logger.info("Cloning source atlas")
            self.atlas = dp.make_atlas_from_spreadsheet(
                source_atlas_df,
                self.ids.atlas,
                filetype="dataframe",
                sheetname="",
                polarity=self.ids.polarity,
                store=True,
                mz_tolerance=source_atlas.compound_identifications[0].mz_references[0].mz_tolerance,
            )

    def _build(self):
        """Populate self._data from database and h5 files."""
        start_time = datetime.datetime.now()
        files = []
        for group in self.ids.groups:
            for h5_file in group.items:
                files.append(
                    (
                        h5_file,
                        group,
                        self.atlas_df,
                        self.atlas,
                        self.extra_time,
                        self.extra_mz,
                    )
                )
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

    def _remove_compound_id(self, idx):
        """
        Remove compound identification at index idx from both in db and self.atlas
        Does not invalidate _data or _hits or _atlas_df
        This bypasses several ORM layers and therefore is a hack, but I couldn't get it to work with the ORM.
        """
        cid_id = self.atlas.compound_identifications[idx].unique_id
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
            del self.atlas.compound_identifications[idx]
        except Exception as err:  # pylint: disable=broad-except
            metoh.rollback_and_log(workspace.db, err)
        workspace.close_connection()

    def filter_compounds(self, keep_idxs=None, remove_idxs=None):
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
        start_len = len(self.atlas_df)
        if remove_idxs is not None:
            _error_if_bad_idxs(self.atlas_df, remove_idxs)
            keep_idxs = self.atlas_df.index.difference(remove_idxs)
        self._atlas_df = self.atlas_df.iloc[keep_idxs].copy().reset_index(drop=True)
        if self._data is not None:
            self._data = [
                [compound for idx, compound in enumerate(sample) if idx in keep_idxs] for sample in self._data
            ]
        if remove_idxs is None:
            remove_idxs = [
                idx for idx, _ in enumerate(self.atlas.compound_identifications) if idx not in keep_idxs
            ]
        _ = [self._remove_compound_id(idx) for idx in sorted(remove_idxs, reverse=True)]
        logger.info(
            "Filtering reduced atlas from %d to %d compounds (%d removed).",
            start_len,
            len(self.atlas_df),
            start_len - len(self.atlas_df),
        )
        if self._hits is not None:
            self.filter_hits_by_atlas()

    def filter_hits_by_atlas(self):
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

    def filter_compounds_ms1_notes_remove(self):
        """
        output:
            updates self.atlas to contain only the compound_identifications that do not have ms1_notes
            starting with 'remove' (case insensitive)
            There is an additional side effect that all mz_tolerances in the returned atlas
            get their value from self.atlas.compound_identifications[0].mz_references[0].mz_tolerance
        """
        logger.debug("Filtering atlas to exclude ms1_notes=='remove'.")
        self.filter_compounds(remove_idxs=self.compound_indices_marked_remove())

    def filter_compounds_by_signal(self, num_points, peak_height):
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

    def store_atlas(self, even_if_exists=False):
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

    def export_atlas_to_csv(self, filename=None):
        """
        save atlas, including ms1_notes, ms2_notes, identification_notes, rt_min, rt_max to filename
        if filename is not provided, then the export is saved to the working directory with filename
        atlas.name + '.csv'
        """
        filename = f"{self.atlas.name}.csv" if filename is None else filename
        dp.export_atlas_to_spreadsheet(self, filename)

    def __getitem__(self, idx):
        """get sample at idx"""
        return self.data[idx]

    @property
    def data(self):
        """data getter, update ._data if necessary"""
        if self._data is None:
            self._build()
        return self._data

    @property
    def atlas_df(self):
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
    def _observe_atlas(self, signal):
        if signal.type == "change":
            self._atlas_df = None
            self._data = None
            logger.debug("Change to atlas invalidates atlas_df, data")

    @observe("_atlas_df")
    def _observe_atlas_df(self, signal):
        if signal.type == "change":
            self._data = None
            logger.debug("Change to atlas_df invalidates data")

    @property
    def polarity(self):
        """
        polarity getter assumes all polarities within class are the same
        returns 'positive' if there are no samples or no compound identifications
        """
        try:
            cid = self.data[0][0]["identification"]
        except IndexError:
            return "positive"
        return cid.mz_references[0].detected_polarity

    @observe("extra_time")
    def _observe_extra_time(self, signal):
        if signal.type == "change":
            self._hits = None
            self._data = None
            logger.debug("Change to extra_time invalidates hits, data")

    @observe("extra_mz")
    def _observe_extra_mz(self, signal):
        if signal.type == "change":
            self._hits = None
            self._data = None
            logger.debug("Change to extra_mz invalidates hits, data")

    @observe("keep_nonmatches")
    def _observe_keep_nonmatches(self, signal):
        if signal.type == "change":
            self._hits = None
            logger.debug("Change to keep_nonmatches invalidates hits")

    @observe("frag_mz_tolerance")
    def _observe_frag_mz_tolerance(self, signal):
        if signal.type == "change":
            self._hits = None
            logger.debug("Change to frag_mz_tolerance invalidates hits")

    @observe("msms_refs_loc")
    def _observe_msms_refs_loc(self, signal):
        if signal.type == "change":
            self._hits = None
            logger.debug("Change to msms_refs_loc invalidates hits")

    @property
    def hits(self):
        """get msms hits DataFrame"""
        _ = self.atlas_df  # regenerate if needed before logging hits generation
        _ = self.data  # regenerate if needed before logging hits generation
        if self._hits is None:
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
        return self._hits

    def __len__(self):
        """len is from data"""
        return len(self.data)

    @property
    def rts(self):
        """
        Allow Rt_Reference objects to be accessed
        use set_rt() if you want to modify the RT values held by this class.
        """
        return tuple(cid.rt_references[0] for cid in self.atlas.compound_identifications)

    def set_rt(self, compound_idx, which, time):
        """
        inputs:
            compound_idx: index of of compound to update
            which: 'rt_min', 'rt_max', or 'rt_peak'
            time: a floating point value for the number of minutes
        updates the RT value in database, self.atlas, self.atlas_df, self.data
        so that no datastructures need to be invalidated
        """
        assert which in ["rt_min", "rt_peak", "rt_max"]
        atlas_rt_ref = self.atlas.compound_identifications[compound_idx].rt_references[0]
        setattr(atlas_rt_ref, which, time)
        for sample in self.data:
            setattr(sample[compound_idx]["identification"].rt_references[0], which, time)
        self.atlas_df.loc[compound_idx, which] = time
        metob.store(atlas_rt_ref)
        if which in ["rt_min", "rt_max"]:
            self._hits_valid_for_rt_bounds = False

    def set_note(self, compound_idx, which, value):
        """
        inputs:
            compound_idx: index of of compound to update
            which: 'ms1_notes', 'ms2_notes' or 'identification_notes'
            value: a string with the note content
        updates the notes value in database, self.atlas, self.atlas_df, self.data
        so that no datastructures need to be invalidated
        """
        assert which in ["ms1_notes", "ms2_notes", "identification_notes"]
        atlas_cid = self.atlas.compound_identifications[compound_idx]
        setattr(atlas_cid, which, value)
        data_cid = self.data[0][compound_idx]["identification"]
        setattr(data_cid, which, value)
        self.atlas_df.loc[compound_idx, which] = value
        metob.store(atlas_cid)

    def compound_indices_marked_remove(self):
        """
        outputs:
            list of compound_idx of the compound identifications with ms1_notes to remove
        """
        ids = ["identification", "ms1_notes"]
        return [i for i, j in enumerate(self.data[0]) if _is_remove(ma_data.extract(j, ids))]

    def compound_idxs_not_evaluated(self):
        """NOT YET IMPLEMENTED"""
        for compound_idx, _ in enumerate(self.data[0]):
            print(compound_idx)
        return []

    def annotation_gui(self, compound_idx=0, width=15, height=3, alpha=0.5, colors=""):
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

    def generate_all_outputs(self, msms_fragment_ions=False, overwrite=False):
        """
        Generates the default set of outputs for a targeted experiment
        inputs:
            msms_fragment_ions: if True, generate msms fragment ions report
            overwrite: if False, throw error if any output files already exist
        """
        if not self._hits_valid_for_rt_bounds:
            self._hits = None  # force hits to be regenerated
        self.extra_time = 0.5
        logger.info("extra_time set to 0.5 minutes for output generation.")
        targeted_output.write_atlas_to_spreadsheet(self, overwrite)
        targeted_output.write_stats_table(self, overwrite)
        targeted_output.write_chromatograms(self, overwrite, max_cpus=self.max_cpus)
        targeted_output.write_identification_figure(self, overwrite)
        targeted_output.write_metrics_and_boxplots(self, overwrite, max_cpus=self.max_cpus)
        if msms_fragment_ions:
            targeted_output.write_msms_fragment_ions(self, overwrite)
        logger.info("Generation of output files completed sucessfully.")
        targeted_output.archive_outputs(self.ids)
        targeted_output.copy_outputs_to_google_drive(self.ids)


class MetatlasSample:
    """
    Object oriented interface to second level of metatlas_dataset. Each instance is one sample (LCMS run).
    """

    def __init__(self, data):
        self._data = data

    def __getitem__(self, idx):
        """get sample at idx"""
        return self._data[idx]

    def __len__(self):
        """len is from data"""
        return len(self._data)


def _duration_since(start):
    """
    inputs:
        start: a datetime object of when the duration started
    returns:
        string with humanized duration of start to now
    """
    return humanize.precisedelta(datetime.datetime.now() - start)


def _is_remove(obj):
    """is obj a string that starts with 'remove' (case insensitive)?"""
    return isinstance(obj, str) and obj.lower().startswith("remove")


def _set_nested(data, ids, value):
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


def _error_if_bad_idxs(dataframe, test_idx_list):
    """Raise IndexError if any members of of test_idx_list are not in dataframe's index"""
    bad = set(test_idx_list) - set(dataframe.index)
    try:
        if len(bad) > 0:
            raise IndexError(f"Invalid index values: {bad}.")
    except IndexError as err:
        logger.exception(err)
        raise err


def get_atlas(name, username):
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


def quoted_string_list(strings):
    """Adds double quotes around each string and seperates with ', '."""
    return ", ".join([f'"{x}"' for x in strings])
