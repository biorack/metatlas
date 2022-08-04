""" object oriented interface to metatlas_dataset """

import datetime
import glob
import logging
import os
import pickle
import shutil
import uuid

from pathlib import Path
from typing import cast, Any, Dict, List, Optional, Tuple, TypedDict, Union

import humanize
import pandas as pd
import traitlets

from traitlets import default, observe
from traitlets import Bool, Float, HasTraits, Instance, Int, Unicode
from traitlets.traitlets import ObserveHandler

import metatlas.plots.dill2plots as dp
import metatlas.datastructures.analysis_identifiers as analysis_ids

from metatlas.datastructures.id_types import PathString, Polarity
from metatlas.datastructures import metatlas_objects as metob
from metatlas.datastructures import object_helpers as metoh
from metatlas.datastructures.utils import AtlasName, get_atlas
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.tools import parallel

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


    ids: analysis_ids.AnalysisIdentifiers instance defining the analysis
    save_metadata: if True, write metadata files containing data sources and LCMS runs short name
    """

    extra_time: float = Float()
    extra_mz: float = Float(default_value=0)
    frag_mz_tolerance: float = Float(default_value=0.01)
    max_cpus: int = Int(default_value=1)
    save_metadata: bool = Bool(default_value=True)
    keep_nonmatches: bool = Bool(default_value=True)
    msms_refs_loc: PathString = Unicode(default_value=analysis_ids.MSMS_REFS_PATH)
    ids: analysis_ids.AnalysisIdentifiers = Instance(klass=analysis_ids.AnalysisIdentifiers)
    atlas: metob.Atlas = Instance(klass=metob.Atlas)
    rt_min_delta = Int(allow_none=True, default_value=None)
    rt_max_delta = Int(allow_none=True, default_value=None)
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
        Adjusts rt_min and rt_max if rt_min_delta or rt_max_delta are not None
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

    def _clone_source_atlas(self) -> metob.Atlas:
        logger.info("Retriving source atlas with unique_id: %s", self.ids.source_atlas_unique_id)
        source_atlas = get_atlas(self.ids.source_atlas_unique_id)
        source_atlas_df = ma_data.make_atlas_df(source_atlas)
        logger.info("Cloning atlas %s", source_atlas.name)
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
        if len(files) == 0:
            logger.warning("No matching h5 files were found!")
        logger.info("Generating MetatlasDataset by reading MSMS data from h5 files")
        samples = parallel.parallel_process(
            ma_data.get_data_for_atlas_df_and_file, files, self.max_cpus, unit="sample"
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

    def _filter_data(self, keep_idxs: List[int]) -> None:
        """Removes any compounds from _data that do not have their index in keep_idxs"""
        if self._data is None:
            return
        self._data = tuple(
            MetatlasSample(
                tuple(compound for idx, compound in enumerate(sample.compounds) if idx in keep_idxs)
            )
            for sample in self._data
        )

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
        _error_if_bad_idxs(self.atlas_df, cast(List[int], remove_idxs if keep_idxs is None else keep_idxs))
        start_len = len(self.atlas_df)
        in_idxs = cast(
            List[int], keep_idxs if remove_idxs is None else self.atlas_df.index.difference(remove_idxs)
        )
        if len(in_idxs) != start_len:
            out_idxs = cast(
                List[int], remove_idxs if keep_idxs is None else self.atlas_df.index.difference(keep_idxs)
            )
            self._atlas_df = self.atlas_df.iloc[in_idxs].copy().reset_index(drop=True)
            self._filter_data(in_idxs)
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
            There is an additional side effect that all mz_tolerances in the returned atlas
            get their value from self.atlas.compound_identifications[0].mz_references[0].mz_tolerance
        """
        logger.debug("Filtering atlas to exclude ms1_notes=='remove'.")
        self.filter_compounds(remove_idxs=self.compound_indices_marked_remove())

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

    @default("extra_time")
    def get_extra_time_default(self) -> float:
        return 0.2 if self.ids.chromatography == "C18" else 0.75

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
        if self._hits is not None:
            return self._hits
        metadata = self._get_hits_metadata()
        self._hits = self._query_cache(metadata)
        if self._hits is not None:
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
            "extra_time": self.extra_time,
            "keep_nonmatches": self.keep_nonmatches,
            "frag_mz_tolerance": self.frag_mz_tolerance,
            "ref_loc": self.msms_refs_loc,
            "extra_mz": self.extra_mz,
            "source_atlas": self.ids.source_atlas,
            "exclude_files": self.ids.exclude_files,
            "exclude_groups": self.ids.exclude_groups,
            "include_groups": self.ids.include_groups,
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

    def update(self) -> None:
        """update hits and data if they no longer are based on current rt bounds"""
        if not self._hits_valid_for_rt_bounds:
            self._hits = None  # force hits to be regenerated
        if not self._data_valid_for_rt_bounds:
            self._data = None  # force data to be regenerated


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


def quoted_string_list(strings: List[str]) -> str:
    """Adds double quotes around each string and seperates with ', '."""
    return ", ".join([f'"{x}"' for x in strings])
