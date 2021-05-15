""" object oriented interface to metatlas_dataset """
import datetime
import getpass
import logging
import multiprocessing

import humanize
import pandas as pd

from metatlas.datastructures import metatlas_objects as metob
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.plots import dill2plots as dp

MSMS_REFS_PATH = "/global/project/projectdirs/metatlas/projects/spectral_libraries/msms_refs_v3.tab"
logger = logging.getLogger(__name__)


class MetatlasDataset:
    """
    Like the non-object oriented metatlas_dataset, you can index into this class by file_idx and compound_idx:
    metatlas_dataset = MetatlasDataset(atlas, groups)
    metatlas_dataset[0][0]['identification'].compound[0].inchi_key

    But MetatlasDataset adds additional functionality, such as:
    metatlas_dataset.hits returns the msms hits dataframe
    metatlas_dataset.atlas returns the atlas
    metatlas_dataset.atlas_df returns the atlas dataframe

    If you change a class property that another property is dependent on, then the second property
    automatically is re-calculated the next time you access the second property. For example:
    metatlas_dataset.extra_time = 0.5  # this invalidates the current hits property
    metatlas_dataset.hits  # this re-generates the hits before returning them
    """

    # pylint: disable=too-many-instance-attributes, too-many-arguments, too-many-public-methods
    def __init__(
        self,
        atlas,
        groups,
        extra_time=0.75,
        extra_mz=0,
        keep_nonmatches=True,
        frag_mz_tolerance=0.01,
        msms_refs_loc=MSMS_REFS_PATH,
        max_cpus=1,
    ):
        self._atlas = atlas
        self._atlas_df = None
        self._atlas_df_valid = False
        self._data = None
        self._data_valid = False
        self._hits = None
        self._hits_valid = False
        self._groups = groups
        self._extra_time = extra_time
        self._extra_mz = extra_mz
        self._keep_nonmatches = keep_nonmatches
        self._frag_mz_tolerance = frag_mz_tolerance
        self._msms_refs_loc = msms_refs_loc
        self.max_cpus = max_cpus

    def _build(self):
        """Populate self._data from database and h5 files."""
        start_time = datetime.datetime.now()
        files = []
        for group in self._groups:
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
        if self.max_cpus > 1:
            with multiprocessing.Pool(processes=min(self.max_cpus, len(files))) as pool:
                samples = pool.map(ma_data.get_data_for_atlas_df_and_file, files)
        else:  # skip multiprocessing as this makes for easier debugging
            samples = [ma_data.get_data_for_atlas_df_and_file(i) for i in files]
        self._data = [MetatlasSample(x) for x in samples]
        logger.info(
            "MetatlasDataset with %d files built in %s.",
            len(files),
            _duration_since(start_time),
        )

    def filter_compounds(self, keep_idxs=None, remove_idxs=None, name=None):
        """
        inputs:
            keep_idxs: the indexes of compounds to keep
            remove_idxs: the indexes of compounds to remove
                Exactly one of keep_idxs or remove_idxs must be None
            name: the name for the new atlas, defaults to current name + '_compound_filtered'
        output:
            If keep_idxs is not None then update self.atlas to contain only the compound_identifications at
            keep_idxs. If remove_idxs is not None then update self.atlas to contain only the compound
            identifications not at remove_idxs. Raises ValueError if both keep_idxs and remove_idxs are None.

            There is an additional side effect that all mz_tolerances in the returned atlas
            get their value from self.atlas.compound_identifications[0].mz_references[0].mz_tolerance

            Does not invalidate _data or _hits
        """
        if (keep_idxs is None) == (remove_idxs is None):
            raise ValueError("Exactly one of keep_idxs and remove_idxs should be None")
        start_len = len(self.atlas_df)
        keep_idxs = keep_idxs if remove_idxs is None else self.atlas_df.index.difference(remove_idxs)
        self._atlas_df = self.atlas_df.iloc[keep_idxs].copy().reset_index(drop=True)
        self._atlas_df_valid = True
        name = f"{self.atlas.name}_compound_filtered" if name is None else name
        mz_tolerance = self.atlas.compound_identifications[0].mz_references[0].mz_tolerance
        if self._data_valid:
            self._data = [
                [compound for idx, compound in enumerate(sample) if idx in keep_idxs] for sample in self._data
            ]
        self._atlas = dp.make_atlas_from_spreadsheet(
            self.atlas_df,
            name,
            filetype="dataframe",
            polarity=self.polarity,
            store=False,
            mz_tolerance=mz_tolerance,
        )
        logger.info(
            "Filtering reduced atlas from %d to %d compounds (%d removed).",
            start_len,
            len(self.atlas_df),
            start_len - len(self.atlas_df),
        )
        if self._hits_valid:
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

    def filter_compounds_ms1_notes_remove(self, name=None):
        """
        inputs:
            name: the name for the new atlas, defaults to current name + '_kept'
        output:
            updates self.atlas to contain only the compound_identifications that do not have ms1_notes
            starting with 'remove' (case insensitive)
            There is an additional side effect that all mz_tolerances in the returned atlas
            get their value from self.atlas.compound_identifications[0].mz_references[0].mz_tolerance
        """
        logger.debug("Filtering atlas to exclude ms1_notes=='remove'.")
        name = f"{self.atlas.name}_kept" if name is None else name
        self.filter_compounds(remove_idxs=self.compound_indices_marked_remove(), name=name)

    def filter_compounds_by_signal(self, num_points, peak_height, name=None):
        """
        inputs:
            num_points: number of points in EIC that must be exceeded in one or more samples
                        in order for the compound to remain in the atlas
            peak_height: max intensity in the EIC that must be exceeded in one or more samples
                         in order for the compound to remain in the atlas
        """
        logger.debug("Filtering atlas on num_points=%d, peak_height=%d.")
        name = f"{self.atlas.name}_strong" if name is None else name
        keep_idxs = dp.strong_signal_compound_idxs(self, num_points, peak_height)
        self.filter_compounds(keep_idxs=keep_idxs, name=name)

    def store_atlas(self, name=None, even_if_exists=False):
        """
        inputs:
            name: name to save to database, if None then use self.atlas.name
            even_if_exists: if True, will save the atlas even if the atlas name already is in the database
                            with your username
        side effects:
            Saves the altas to the database.
            Raises ValueError if even_if_exists==False and name is already in the database with your username
        """
        name = self.atlas.name if name is None else name
        username = getpass.getuser()
        if not even_if_exists and len(metob.retrieve("atlases", name=name, username=username)) > 0:
            raise ValueError(f"An atlas with name {name} and owned by {username} already exists.")
        metob.store(self.atlas)

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

    def __setitem__(self, idx, value):
        """assign value for sample at idx"""
        self._data[idx] = value

    def _set_and_invalidate_properties(self, attribute_name, new_value, property_names):
        """
        inputs:
            attribute_name: name of the class attribute being modified
            new_value: value to assign to attribute
            propert_names: list of names of the class propertys that are dependent on the attribute's value
        side effects:
            If the property is valid and new_value is different from previous value, then invalidate.
            And set attribute to new_value
        """
        for prop in property_names:
            valid_attr_name = f"_{prop}_valid"
            setattr(
                self,
                valid_attr_name,
                getattr(self, valid_attr_name) and new_value == getattr(self, attribute_name),
            )
        setattr(self, f"_{attribute_name}", new_value)

    @property
    def data(self):
        """data getter, update ._data if necessary"""
        if not self._data_valid:
            self._build()
            self._data_valid = True
        return self._data

    @property
    def atlas_df(self):
        """atlas_df getter, update ._atlas_df if necessary"""
        if not self._atlas_df_valid:
            start_time = datetime.datetime.now()
            self._atlas_df = ma_data.make_atlas_df(self.atlas)
            self._atlas_df_valid = True
            logger.info(
                "Generated atlas_df with %d rows in %s.",
                len(self.atlas_df),
                _duration_since(start_time),
            )
        return self._atlas_df

    @property
    def atlas(self):
        """atlas getter"""
        return self._atlas

    @atlas.setter
    def atlas(self, atlas):
        """atlas setter, invalidate atlas_df and data"""
        self._set_and_invalidate_properties("atlas", atlas, ["atlas_df", "data"])

    @property
    def groups(self):
        """groups getter"""
        return self._groups

    @groups.setter
    def groups(self, groups):
        """groups setter, invalidate data"""
        self._set_and_invalidate_properties("groups", groups, ["data"])

    @property
    def polarity(self):
        """polarity getter assumes all polarities within class are the same"""
        return self.data[0][0]["identification"].mz_references[0].detected_polarity

    @property
    def extra_time(self):
        """extra_time getter"""
        return self._extra_time

    @extra_time.setter
    def extra_time(self, extra_time):
        """extra_time setter, invalidates data and hits"""
        self._set_and_invalidate_properties("extra_time", extra_time, ["data", "hits"])

    @property
    def extra_mz(self):
        """extra_mz getter"""
        return self._extra_mz

    @extra_mz.setter
    def extra_mz(self, extra_mz):
        """extra_mz setter, invalidates data and hits"""
        self._set_and_invalidate_properties("extra_mz", extra_mz, ["data", "hits"])

    @property
    def keep_nonmatches(self):
        """keep_nonmatches getter"""
        return self._keep_nonmatches

    @keep_nonmatches.setter
    def keep_nonmatches(self, keep_nonmatches):
        """keep_nonmatches setter, invalidates hits"""
        self._set_and_invalidate_properties("keep_nonmatches", keep_nonmatches, ["hits"])

    @property
    def frag_mz_tolerance(self):
        """frag_mz_tolerance getter"""
        return self._frag_mz_tolerance

    @frag_mz_tolerance.setter
    def frag_mz_tolerance(self, frag_mz_tolerance):
        """frag_mz_tolerance setter, invlidates hits"""
        self._set_and_invalidate_properties("frag_mz_tolerance", frag_mz_tolerance, ["hits"])

    @property
    def msms_refs_loc(self):
        """msms_refs_loc getter"""
        return self._msms_refs_loc

    @msms_refs_loc.setter
    def msms_refs_loc(self, msms_refs_loc):
        """msms_refs_loc setter, invalidates hits"""
        self._set_and_invalidate_properties("msms_refs_loc", msms_refs_loc, ["hits"])

    @property
    def hits(self):
        """get msms hits DataFrame"""
        if not self._hits_valid:
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
            self._hits_valid = True
        return self._hits

    def __len__(self):
        """len is from data"""
        return len(self.data)

    def set_data(self, ids, value):
        """update a value within self._data"""
        self._atlas_df_valid = False
        self._data_valid = False
        _set_nested(self._data, ids, value)

    @property
    def rts(self):
        """
        Allow Rt_Reference objects to be accessed
        Because this returns a list, the return value is mutable add will modify
        data internal to this class, but atlas_df and data would get out of sync with atlas.
        So don't modify the values returned by this property!
        """
        return [cid.rt_references[0] for cid in self.atlas.compound_identifications]

    def set_rt(self, compound_idx, which, time):
        """
        inputs:
            compound_idx: index of of compound to update
            which: 'rt_min', 'rt_max', or 'rt_peak'
            time: a floating point value for the number of minutes
        updates the RT value in 3 places so that no datastructures need to be invalidated
        """
        assert which in ["rt_min", "rt_peak", "rt_max"]
        atlas_rt_ref = self.atlas.compound_identifications[compound_idx].rt_references[0]
        setattr(atlas_rt_ref, which, time)
        data_rt_ref = self._data[0][compound_idx]["identification"].rt_references[0]
        setattr(data_rt_ref, which, time)
        self._atlas_df.loc[compound_idx, which] = time

    def set_note(self, compound_idx, which, value):
        """
        inputs:
            compound_idx: index of of compound to update
            which: 'ms1_notes', 'ms2_notes' or 'identification_notes'
            value: a string with the note content
        updates the RT value in 3 places so that no datastructures need to be invalidated
        """
        assert which in ["ms1_notes", "ms2_notes", "identification_notes"]
        atlas_cid = self.atlas.compound_identifications[compound_idx]
        setattr(atlas_cid, which, value)
        data_cid = self._data[0][compound_idx]["identification"]
        setattr(data_cid, which, value)
        self._atlas_df.loc[compound_idx, which] = value

    def compound_indices_marked_remove(self):
        """
        outputs:
            list of compound_idx of the compound identifications with ms1_notes to remove
        """
        ids = ["identification", "ms1_notes"]
        return [i for i, j in enumerate(self.data[0]) if _is_remove(ma_data.extract(j, ids))]


class MetatlasSample:
    """
    Object oriented interface to second level of metatlas_dataset. Each instance is one sample (LCMS run).
    """

    def __init__(self, data):
        self._data = data

    def __getitem__(self, idx):
        """get sample at idx"""
        return self._data[idx]

    def __setitem__(self, idx, value):
        """assign value for sample at idx"""
        self._data[idx] = value

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
    if len(ids) == 0:
        raise ValueError("ids cannot be empty")
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
