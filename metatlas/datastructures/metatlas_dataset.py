""" object oriented interface to metatlas_dataset """

import datetime
import logging
import pickle
import shutil
import uuid

from collections.abc import Iterable
from pathlib import Path
from typing import cast, Any, Dict, List, Optional, Sequence, Tuple, TypedDict, Union

import humanize
import pandas as pd
import traitlets

from traitlets import default, observe
from traitlets import Bool, Float, HasTraits, Instance, Int
from traitlets.traitlets import ObserveHandler

import metatlas.plots.dill2plots as dp
import metatlas.datastructures.analysis_identifiers as analysis_ids

from metatlas.datastructures.id_types import Polarity
from metatlas.datastructures import metatlas_objects as metob
from metatlas.datastructures import object_helpers as metoh
from metatlas.datastructures.utils import AtlasName, get_atlas, set_atlas_mz_tolerance
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.io.metatlas_get_data_helper_fun import extract
from metatlas.tools import parallel
from metatlas.tools.config import NotebookParameters_co
from metatlas.tools.util import is_atlas_df_subset

logger = logging.getLogger(__name__)


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


class MetatlasSample:
    """
    Object oriented interface to second level of metatlas_dataset. Each instance is one sample (LCMS run).
    """

    def __init__(self, compounds: Tuple[CompoundDict, ...]) -> None:
        self.compounds: Tuple[CompoundDict, ...] = compounds

    def __getitem__(self, idx: int) -> CompoundDict:
        """get sample at idx"""
        return self.compounds[idx]

    def __iter__(self) -> Iterable:
        """allow use as an iterator"""
        yield from self.compounds

    def __len__(self) -> int:
        """len is from data"""
        return len(self.compounds)


SampleSet = Tuple[MetatlasSample, ...]


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


    ids: analysis_ids.AnalysisIdentifiers instance defining the analysis
    save_metadata: if True, write metadata files containing data sources and LCMS runs short name
    """

    extra_time: float = Float()
    extra_mz: float = Float(default_value=0)
    max_cpus: int = Int(default_value=1)
    save_metadata: bool = Bool(default_value=True)
    keep_nonmatches: bool = Bool(default_value=True)
    ids: analysis_ids.AnalysisIdentifiers = Instance(klass=analysis_ids.AnalysisIdentifiers)
    atlas: metob.Atlas = Instance(klass=metob.Atlas)
    rt_min_delta = Int(allow_none=True, default_value=None)
    rt_max_delta = Int(allow_none=True, default_value=None)
    _atlas_df: Optional[pd.DataFrame] = Instance(klass=pd.DataFrame, allow_none=True, default_value=None)
    # _all_data contanis all data in experiement before any filtering
    _all_data: Optional[SampleSet] = traitlets.Tuple(allow_none=True, default_value=None)
    # _data contains the post-filtering data that is current being utilized
    _data: Optional[SampleSet] = traitlets.Tuple(allow_none=True, default_value=None)
    _hits: Optional[pd.DataFrame] = Instance(klass=pd.DataFrame, allow_none=True, default_value=None)

    # pylint: disable=too-many-instance-attributes, too-many-arguments, too-many-public-methods
    def __init__(self, **kwargs) -> None:
        """Constructor"""
        super().__init__(**kwargs)
        logger.debug("Creating new MetatlasDataset instance...")
        self._hits_valid_for_rt_bounds = False  # based only on RT min/max changes
        self._data_valid_for_rt_bounds = False  # based only on RT min/max changes
        self._get_atlas()
        if self.save_metadata:
            logger.debug("Writing MetatlasDataset metadata files")
            self.write_data_source_files()
            self.ids.write_lcmsruns_short_names()
        self.ids.register_groups_invalidation_callback(self.invalidate_data)

    def _save_to_cache(self, data: Any, metadata: dict) -> None:
        assert "_variable_name" in metadata.keys()
        metadata = metadata.copy()
        name = metadata["_variable_name"]
        base_name = f"{name}_{uuid.uuid4()}"
        metadata["_pickle_file_name"] = f"{base_name}.pkl"  # relative to metadata file
        metadata_file_name = self.ids.cache_dir / f"{base_name}.metadata"
        pickle_path = self.ids.cache_dir / metadata["_pickle_file_name"]
        with open(pickle_path, "wb") as pickle_fh:
            pickle.dump(data, pickle_fh)
        with open(metadata_file_name, "wb") as metadata_fh:
            pickle.dump(metadata, metadata_fh)
        logger.info("Caching %s in %s.", name, pickle_path)

    def _query_cache(self, required_metadata: dict) -> Optional[Any]:
        assert "_variable_name" in required_metadata.keys()
        required_metadata = required_metadata.copy()
        name = required_metadata["_variable_name"]
        if name == "hits" and "atlas_df" in required_metadata:
            del required_metadata["atlas_df"]
        files_matching_metadata = []
        for metadata_file in Path(self.ids.cache_dir).glob(f"{name}_*.metadata"):
            with open(metadata_file, "rb") as metadata_fh:
                potential_metadata = pickle.load(metadata_fh)
                pickle_file_name = Path(self.ids.cache_dir) / potential_metadata["_pickle_file_name"]
                if name == "hits":
                    if "atlas_df" not in potential_metadata:
                        continue  # do not use cache files pre-dating atlas_df in metadata
                    cached_atlas_df = potential_metadata["atlas_df"]
                    if not is_atlas_df_subset(cached_atlas_df, self.atlas_df):
                        continue
                    del potential_metadata["atlas_df"]
                # require_metadata does not have a '_pickle_file_name' key
                # so remove before equality test
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
        data_sources_dir = self.ids.additional_output_dir / f"{self.ids.short_polarity}_data_sources"
        if len(list(data_sources_dir.glob("*"))) >= 4:
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
                self.ids.groups, self.atlas, self.ids.additional_output_dir, self.ids.short_polarity
            )

    def _get_atlas(self) -> None:
        """
        Copy source atlas from database into current analysis atlas.
        If the atlas does not yet exist, it will be copied from source_atlas.
        Adjusts rt_min and rt_max if rt_min_delta or rt_max_delta are not None.
        Sets mz_tolerance based on parameters mz_tolerance_override and mz_tolerance_default
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
            temp_atlas = atlases[0]
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
            temp_atlas = self._clone_source_atlas()
        self.atlas = metob.adjust_atlas_rt_range(temp_atlas, self.rt_min_delta, self.rt_max_delta)
        if self.parameters.mz_tolerance_override is not None:
            set_atlas_mz_tolerance(self.atlas, self.parameters.mz_tolerance_override, True)
        else:
            set_atlas_mz_tolerance(self.atlas, self.parameters.mz_tolerance_default)

    def _clone_source_atlas(self) -> metob.Atlas:
        logger.info("Retriving source atlas with unique_id: %s", self.ids.source_atlas_unique_id)
        source_atlas = get_atlas(self.ids.source_atlas_unique_id)
        logger.info("Cloning atlas %s", source_atlas.name)
        logger.info("Traits of source atlas: %s", source_atlas.traits())
        new_atlas = source_atlas.clone(recursive=False)
        logger.info("Traits of new atlas: %s", new_atlas.traits())
        new_atlas.name = self.ids.atlas
        logger.info("Storing new atlas %s", new_atlas.name)
        metob.store(new_atlas)
        logger.info("Returning new atlas %s", new_atlas.name)
        return new_atlas

    def _get_all_data(self) -> SampleSet:
        """Populate self._all_data from database and h5 files."""
        start_time = datetime.datetime.now()
        files = [
            (h5_file, group, self.atlas_df, self.atlas, self.extra_time, self.extra_mz)
            for group in self.ids.all_groups
            for h5_file in group.items
        ]
        if len(files) == 0:
            logger.warning("No matching h5 files were found!")
        logger.info("Reading MSMS data from h5 files...")
        samples = parallel.parallel_process(
            ma_data.get_data_for_atlas_df_and_file, files, self.max_cpus, unit="sample"
        )
        logger.info(
            "Loaded MS data from %d h5 files in %s.",
            len(files),
            _duration_since(start_time),
        )
        return tuple(MetatlasSample(x) for x in samples)

    def _get_data(self) -> SampleSet:
        """Generate _data from _all_data based on current attributes"""

        def contains(inchi_key: str, adduct: str) -> bool:
            adf = self.atlas_df
            return len(adf[(adf["inchi_key"] == inchi_key) & (adf["adduct"] == adduct)]) > 0

        adduct_indexes = ["compounds", 0, "identification", "mz_references", 0, "adduct"]
        inchi_key_indexes = ["compounds", 0, "identification", "compound", 0, "inchi_key"]
        h5_files = {run.hdf5_file for group in self.ids.groups for run in group.items}
        return tuple(
            x
            for x in self.all_data
            if (
                extract(x, [0, "lcmsrun", "hdf5_file"]) in h5_files
                and contains(extract(x, inchi_key_indexes), extract(x, adduct_indexes))
            )
        )

    def _remove_compound_id(self, idx: int) -> None:
        """
        Remove compound identification at index idx from db, self.atlas, and self._data
        Does not invalidate _data or _hits or _atlas_df
        This bypasses several ORM layers and therefore is a hack, but I couldn't get it to work with the ORM.
        """
        cid_id = self.atlas.compound_identifications[idx].unique_id
        del self.atlas.compound_identifications[idx]
        self._filter_data(list(set(range(len(extract(self._data, [0], [])))) - {idx}))
        if self._atlas_df is not None:
            self._atlas_df.drop(index=idx, inplace=True)
            self._atlas_df.reset_index(drop=True, inplace=True)
        atlas_id = self.atlas.unique_id
        link_table = "atlases_compound_identifications"
        target = f"target_id='{cid_id}'"
        workspace = metob.Workspace.get_instance()
        db_conn = workspace.get_connection()
        db_conn.begin()
        try:
            db_conn.query(f"delete from {link_table} where ({target} and source_id='{atlas_id}')")
            links = db_conn.query(f"select source_id from {link_table} where {target}")
            if len(list(links)) == 0:  # other atlases are not linked to this CompoundIdentification
                db_conn.query(f"delete from compoundidentifications where unique_id='{cid_id}'")
            db_conn.commit()
        except Exception as err:  # pylint: disable=broad-except
            metoh.rollback_and_log(db_conn, err)
            raise Exception from err
        finally:
            metoh.close_db_connection(db_conn)

    def _filter_data(self, keep_idxs: Sequence[int]) -> None:
        """Removes any compounds from _all_data and _data that do not have their index in keep_idxs"""
        if self._data is None:
            return
        for sample in self.data:
            sample.compounds = tuple(
                compound for idx, compound in enumerate(sample.compounds) if idx in keep_idxs
            )

    def _assert_consistent_num_compounds(self):
        """Test that internal data structures are in sync"""
        atlas_len = len(self.atlas.compound_identifications)
        assert len(self.atlas_df) == atlas_len
        assert len(self.rts) == atlas_len
        if self._data is not None and len(self._data) > 0:
            assert len(self[0]) == atlas_len

    def filter_compounds(
        self, keep_idxs: Optional[Sequence[int]] = None, remove_idxs: Optional[Sequence[int]] = None
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
        _error_if_bad_idxs(self.atlas_df, cast(List[int], remove_idxs if keep_idxs is None else keep_idxs))
        self._assert_consistent_num_compounds()
        start_len = len(self.atlas_df)
        in_idxs = cast(
            List[int], keep_idxs if remove_idxs is None else self.atlas_df.index.difference(remove_idxs)
        )
        if len(in_idxs) != start_len:
            out_idxs = cast(
                List[int], remove_idxs if keep_idxs is None else self.atlas_df.index.difference(keep_idxs)
            )
            for i in sorted(out_idxs, reverse=True):
                self._remove_compound_id(i)
        logger.info(
            "Filtering reduced atlas from %d to %d compounds (%d removed).",
            start_len,
            len(self.atlas_df),
            start_len - len(self.atlas_df),
        )
        if len(in_idxs) != start_len and self._hits is not None:
            self.filter_hits_by_atlas()
            self._save_to_cache(self._hits, self._get_hits_metadata())
        self._assert_consistent_num_compounds()

    def filter_hits_by_atlas(self) -> None:
        """Remove any hits that do not have a corresponding inchi_key-adduct pair in atlas_df"""
        start_len = len(self.hits)
        keep_adducts = self.atlas_df.loc[:, ["inchi_key", "adduct"]].drop_duplicates()
        logger.info("Number of inchi_key-adduct pairs is %d.", len(keep_adducts))
        hits_plus = self.hits.copy()
        hits_plus["copy_index"] = hits_plus.index.to_numpy()
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
        """
        logger.debug("Filtering atlas to exclude ms1_notes=='remove'.")
        self.filter_compounds(remove_idxs=self.compound_indices_marked_remove)

    def filter_compounds_by_signal(
        self,
        num_points: Optional[int] = None,
        peak_height: Optional[float] = None,
        msms_score: Optional[float] = None,
    ) -> None:
        """
        inputs:
            num_points: number of points in EIC that must be exceeded in one or more samples
                        in order for the compound to remain in the atlas
            peak_height: max intensity in the EIC that must be exceeded in one or more samples
                         in order for the compound to remain in the atlas
            msms_score: spectra similarity score that must be exceeded in one or more samples
                         in order for the compound to remain in the atlas
        """
        if num_points is not None:
            logger.info("Filtering atlas on num_points=%d.", num_points)
        if peak_height is not None:
            logger.info("Filtering atlas on peak_height=%.2f.", peak_height)
        if msms_score is not None:
            logger.info("Filtering atlas on msms_score=%.2f.", msms_score)
        keep_idxs = dp.strong_signal_compound_idxs(self, num_points, peak_height, msms_score)
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

    def __iter__(self) -> Iterable:
        """allow use as an iterator"""
        yield from self.data

    @default("extra_time")
    def get_extra_time_default(self) -> float:
        """set extra_time based on chromatography type"""
        # TODO this should get moved into the config file
        return 0.2 if self.ids.chromatography == "C18" else 0.75

    @property
    def all_data(self) -> SampleSet:
        """all_dat getter, update ._all_data if necessary"""
        if self._all_data is None:
            self._all_data = self._get_all_data()
        return self._all_data

    @property
    def data(self) -> SampleSet:
        """data getter, update ._data if necessary"""
        if self._data is None:
            self._data = self._get_data()
            self._data_valid_for_rt_bounds = True
        return self._data

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
            self._all_data = None
            logger.debug("Change to atlas invalidates atlas_df, all_data")

    @observe("_atlas_df")
    def _observe_atlas_df(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._data = None
            logger.debug("Change to atlas_df invalidates data")

    @observe("extra_time")
    def _observe_extra_time(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._hits = None
            self._all_data = None
            logger.debug("Change to extra_time invalidates hits, all_data")

    @observe("_data")
    def _observe_data(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._hits = None
            logger.debug("Change to data invalidates hits")

    @observe("_all_data")
    def _observe_all_data(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._data = None
            logger.debug("Change to all_data invalidates data")

    @observe("extra_mz")
    def _observe_extra_mz(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._all_data = None
            logger.debug("Change to extra_mz invalidates all_data")

    @observe("keep_nonmatches")
    def _observe_keep_nonmatches(self, signal: ObserveHandler) -> None:
        if signal.type == "change":
            self._hits = None
            logger.debug("Change to keep_nonmatches invalidates hits")

    @property
    def hits(self) -> pd.DataFrame:
        """get msms hits DataFrame"""
        if self._hits is not None:
            return self._hits
        metadata = self._get_hits_metadata()
        self._hits = self._query_cache(metadata)
        if self._hits is not None:
            self._hits_valid_for_rt_bounds = False  # unsure, so assume False
            return self._hits
        _ = self.data  # regenerate if needed before logging hits generation
        logger.info(
            "Generating hits with extra_time=%.3f, frag_mz_tolerance=%.4f, msms_refs=%s.",
            self.extra_time,
            self.parameters.frag_mz_tolerance,
            self.parameters.msms_refs,
        )
        start_time = datetime.datetime.now()
        self._hits = dp.get_msms_hits(
            self.data,
            extra_time=self.extra_time > 0,
            keep_nonmatches=self.keep_nonmatches,
            frag_mz_tolerance=self.parameters.frag_mz_tolerance,
            ref_loc=self.parameters.msms_refs,
            resolve_by=self.parameters.resolve_msms_matches_by,
        )
        logger.info("Generated %d hits in %s.", len(self._hits), _duration_since(start_time))
        self._hits_valid_for_rt_bounds = True
        self._save_to_cache(self._hits, metadata)
        return self._hits

    def _get_hits_metadata(self) -> Dict[str, Any]:
        return {
            "_variable_name": "hits",
            "extra_time": self.extra_time,
            "keep_nonmatches": self.keep_nonmatches,
            "frag_mz_tolerance": self.parameters.frag_mz_tolerance,
            "ref_loc": self.parameters.msms_refs,
            "extra_mz": self.extra_mz,
            "source_atlas": self.ids.source_atlas,
            "include_lcmsruns": self.ids.include_lcmsruns,
            "exclude_lcmsruns": self.ids.exclude_lcmsruns,
            "include_groups": self.ids.include_groups,
            "exclude_groups": self.ids.exclude_groups,
            "atlas_df": self.atlas_df,
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

    @property
    def compound_indices_marked_remove(self) -> List[int]:
        """
        outputs:
            list of compound_idx of the compound identifications with ms1_notes to remove
        """
        ids = ["identification", "ms1_notes"]
        return [i for i, j in enumerate(self.data[0].compounds) if _is_remove(ma_data.extract(j, ids))]

    @property
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

    @property
    def parameters(self) -> NotebookParameters_co:
        workflow = self.ids.configuration.get_workflow(self.ids.workflow)
        return workflow.get_analysis(self.ids.analysis).parameters

    def error_if_not_all_evaluated(self) -> None:
        """Raises ValueError if there are compounds that have not been evaluated"""
        not_evaluated = self.compound_idxs_not_evaluated
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

    def update(self) -> None:
        """update hits and data if they no longer are based on current rt bounds"""
        if not self._hits_valid_for_rt_bounds:
            self._hits = None  # force hits to be regenerated
        if not self._data_valid_for_rt_bounds:
            self._all_data = None
            self._data = None

    def invalidate_data(self) -> None:
        """Force data to be reloaded from disk on next usage"""
        self._data = None


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


def _set_nested(data: Any, ids: Sequence[Union[int, str, Tuple[str]]], value: Any):
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


def _error_if_bad_idxs(dataframe: pd.DataFrame, test_idx_list: Sequence[int]) -> None:
    """Raise IndexError if any members of of test_idx_list are not in dataframe's index"""
    bad = set(test_idx_list) - set(dataframe.index)
    try:
        if len(bad) > 0:
            raise IndexError(f"Invalid index values: {bad}.")
    except IndexError as err:
        logger.exception(err)
        raise err


def quoted_string_list(strings: Sequence[str]) -> str:
    """Adds double quotes around each string and seperates with ', '."""
    return ", ".join([f'"{x}"' for x in strings])


def _filter_data_by_atlas_df(data: SampleSet, atlas_df: pd.DataFrame) -> SampleSet:
    """Remove compounds-adduct pairs from data that are not in atlas_df"""

    def contains(adf: pd.DataFrame, inchi_key: str, adduct: str) -> bool:
        return len(adf[(adf["inchi_key"] == inchi_key) & (adf["adduct"] == adduct)]) > 0

    adduct_indexes = [0, "compounds", 0, "identification", "mz_references", 0, "adduct"]
    inchi_key_indexes = [0, "compounds", 0, "identification", "compound", 0, "inchi_key"]
    return tuple(
        x for x in data if contains(atlas_df, extract(x, inchi_key_indexes), extract(x, adduct_indexes))
    )
