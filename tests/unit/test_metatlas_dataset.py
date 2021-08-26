""" tests for MetatasDataset """
# pylint: disable=missing-function-docstring,protected-access,unused-argument,too-many-arguments

import datetime
import glob
import logging
import os

import pandas as pd
import pytest
import traitlets

from metatlas.datastructures import metatlas_dataset as mads
from metatlas.datastructures import metatlas_objects as metob
from metatlas.datastructures import object_helpers as metoh
from metatlas.io import metatlas_get_data_helper_fun as ma_data


def test_metatlas_dataset_build01(metatlas_dataset):
    assert len(metatlas_dataset) == 1
    assert len(metatlas_dataset[0]) == 1
    assert metatlas_dataset[0][0]["identification"].compound[0].inchi_key == "OLXZPDWKRNYJJZ-RRKCRQDMSA-N"
    assert metatlas_dataset[0][0]["data"]["ms1_summary"]["rt_peak"] == 2.2922415733
    assert (
        metatlas_dataset[0][0]["lcmsrun"].experiment
        == "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583"
    )


@pytest.mark.xfail
def test_metatlas_dataset_build02(mocker, atlas, group_with_2_lcmsruns, df_container, analysis_ids):
    # need to mock multiprocessing for this to work
    mocker.patch(
        "metatlas.io.metatlas_get_data_helper_fun.df_container_from_metatlas_file", return_value=df_container
    )
    metatlas_dataset = mads.MetatlasDataset(ids=analysis_ids, max_cpus=2)
    assert len(metatlas_dataset) == 2
    assert len(metatlas_dataset[0]) == 1


def test_filter_compounds_ms1_notes_remove01(mocker, metatlas_dataset_with_2_cids, compound):
    mocker.patch("metatlas.datastructures.metatlas_objects.retrieve", return_value=[compound])
    metatlas_dataset = metatlas_dataset_with_2_cids
    metatlas_dataset.filter_compounds_ms1_notes_remove()
    assert len(metatlas_dataset[0]) == 2
    metatlas_dataset.set_note(1, "ms1_notes", "Remove")
    metatlas_dataset.filter_compounds_ms1_notes_remove()
    assert len(metatlas_dataset[0]) == 1


def test_filter_compounds01(metatlas_dataset_with_2_cids):
    # mocker.patch("metatlas.datastructures.metatlas_objects.retrieve", return_value=[compound])
    metatlas_dataset = metatlas_dataset_with_2_cids
    metatlas_dataset.filter_compounds(remove_idxs=[])
    assert len(metatlas_dataset[0]) == 2
    metatlas_dataset.filter_compounds(keep_idxs=[0, 1])
    assert len(metatlas_dataset[0]) == 2
    metatlas_dataset.filter_compounds(keep_idxs=[])
    assert len(metatlas_dataset[0]) == 0
    with pytest.raises(ValueError):
        metatlas_dataset.filter_compounds()


def test_filter_compounds02(mocker, metatlas_dataset_with_2_cids, compound):
    mocker.patch("metatlas.datastructures.metatlas_objects.retrieve", return_value=[compound])
    metatlas_dataset = metatlas_dataset_with_2_cids
    with pytest.raises(ValueError):
        metatlas_dataset.filter_compounds(keep_idxs=[0], remove_idxs=[1])


def test_filter_compounds03(mocker, metatlas_dataset, compound):
    mocker.patch("metatlas.datastructures.metatlas_objects.retrieve", return_value=[compound])
    with pytest.raises(IndexError):
        metatlas_dataset.filter_compounds(keep_idxs=[999])


def test_filter_compounds04(mocker, metatlas_dataset, compound):
    mocker.patch("metatlas.datastructures.metatlas_objects.retrieve", return_value=[compound])
    with pytest.raises(IndexError):
        metatlas_dataset.filter_compounds(remove_idxs=[999])


def test_filter_compounds05(mocker, metatlas_dataset_with_2_cids, username):
    original_rt_min = metatlas_dataset_with_2_cids.rts[1].rt_min
    print([r.rt_min for r in metatlas_dataset_with_2_cids.rts])
    updated_rt_min = 9.99
    metatlas_dataset_with_2_cids.set_rt(1, "rt_min", updated_rt_min)
    metatlas_dataset_with_2_cids.filter_compounds(remove_idxs=[0])
    atlas = metob.retrieve("Atlas", name=metatlas_dataset_with_2_cids.atlas.name, username=username)[0]
    assert atlas.compound_identifications[0].rt_references[0].rt_min != original_rt_min
    assert atlas.compound_identifications[0].rt_references[0].rt_min == updated_rt_min


def test_filter_hits_by_atlas01(mocker, metatlas_dataset_with_2_cids, hits, compound):
    mocker.patch("metatlas.plots.dill2plots.get_msms_hits", return_value=hits)
    mocker.patch("metatlas.datastructures.metatlas_objects.retrieve", return_value=[compound])
    hits = metatlas_dataset_with_2_cids.hits
    start_num = len(hits)
    metatlas_dataset_with_2_cids.filter_compounds(keep_idxs=[0])
    assert start_num > len(metatlas_dataset_with_2_cids.hits)
    metatlas_dataset_with_2_cids.filter_compounds(remove_idxs=[0])
    assert len(metatlas_dataset_with_2_cids.hits) == 0


def test_polarity(metatlas_dataset):
    assert metatlas_dataset.polarity == "positive"
    metatlas_dataset.filter_compounds(remove_idxs=[0])
    assert len(metatlas_dataset[0]) == 0
    assert metatlas_dataset.polarity == "positive"


def test_extra_time_setter(metatlas_dataset, hits, mocker):
    mocker.patch("metatlas.plots.dill2plots.get_msms_hits", return_value=hits)
    metatlas_dataset.hits  # pylint: disable=pointless-statement
    assert metatlas_dataset._hits is not None
    metatlas_dataset.extra_time = 0.3
    assert metatlas_dataset._hits is None
    metatlas_dataset.hits  # pylint: disable=pointless-statement
    assert metatlas_dataset._hits is not None


def test_rts01(metatlas_dataset):
    metatlas_dataset.set_rt(0, "rt_min", 9.99)
    assert metatlas_dataset.rts[0].rt_min == 9.99
    assert len(metatlas_dataset.rts) == 1


def test_rts02(metatlas_dataset):
    metatlas_dataset._atlas_df = None
    metatlas_dataset.set_rt(0, "rt_max", 9.99)
    assert metatlas_dataset.rts[0].rt_max == 9.99
    assert len(metatlas_dataset.rts) == 1


def test_rts03(metatlas_dataset, analysis_ids):
    assert metatlas_dataset.rts[0].rt_max != 9.99
    metatlas_dataset.set_rt(0, "rt_max", 9.99)
    second_metatlas_dataset = mads.MetatlasDataset(ids=analysis_ids)
    assert second_metatlas_dataset.rts[0].rt_max == 9.99
    assert len(second_metatlas_dataset.rts) == 1


def test_rts04(analysis_ids, sqlite_with_atlas, mocker, lcmsrun, df_container):
    mocker.patch(
        "metatlas.io.metatlas_get_data_helper_fun.df_container_from_metatlas_file", return_value=df_container
    )
    mocker.patch("metatlas.plots.dill2plots.get_metatlas_files", return_value=[lcmsrun])
    first = mads.MetatlasDataset(ids=analysis_ids)
    first.set_rt(0, "rt_max", 1.11)
    second = mads.MetatlasDataset(ids=analysis_ids)
    assert second.rts[0].rt_max == 1.11
    second.set_rt(0, "rt_max", 2.22)
    third = mads.MetatlasDataset(ids=analysis_ids)
    assert third.rts[0].rt_max == 2.22


def test_set_note01(metatlas_dataset, sqlite):
    metatlas_dataset.set_note(0, "ms2_notes", "Foobar")
    assert metatlas_dataset[0][0]["identification"].ms2_notes == "Foobar"


def test_set_note02(metatlas_dataset):
    metatlas_dataset._atlas_df = None
    metatlas_dataset.set_note(0, "ms1_notes", "keeper")
    assert metatlas_dataset[0][0]["identification"].ms1_notes == "keeper"


def test_compound_indices_marked_remove01(sqlite, metatlas_dataset):
    assert len(metatlas_dataset.compound_indices_marked_remove()) == 0
    metatlas_dataset.set_note(0, "ms1_notes", "REMOVE")
    assert len(metatlas_dataset.compound_indices_marked_remove()) == 1


def test_set_nested01():
    with pytest.raises(ValueError):
        mads._set_nested([], [], 0)


def test_set_nested02(atlas):
    mads._set_nested(atlas, ["compound_identifications", 0, ("compound",), 0, ("inchi_key",)], "FOOBAR")
    assert atlas.compound_identifications[0].compound[0].inchi_key == "FOOBAR"


def test_set_nested03(atlas):
    mads._set_nested(atlas, ["name"], "My Atlas")
    assert atlas.name == "My Atlas"


def test_set_nested04(atlas):
    with pytest.raises(TypeError):
        mads._set_nested(atlas, ["zoop"], None)


def test_set_nested05():
    my_dict = {}
    mads._set_nested(my_dict, ["zoop"], None)
    assert my_dict["zoop"] is None


def test_error_if_bad_idxs():
    data = pd.DataFrame(data={"a": [1, 2], "b": [3, 4]})
    mads._error_if_bad_idxs(data, [0])
    with pytest.raises(IndexError):
        mads._error_if_bad_idxs(data, [2])


def test_is_remove():
    assert not mads._is_remove([])
    assert not mads._is_remove("foobar")
    assert mads._is_remove("Remove")
    assert mads._is_remove("REMOVE AND MORE")


def test_duration_since():
    assert mads._duration_since(datetime.datetime.now()) == "0.00 seconds"


def test_filter_compounds_by_signal01(mocker, metatlas_dataset_with_2_cids, df_container, compound):
    mocker.patch(
        "metatlas.io.metatlas_get_data_helper_fun.df_container_from_metatlas_file", return_value=df_container
    )
    mocker.patch("metatlas.datastructures.metatlas_objects.retrieve", return_value=[compound])
    assert len(metatlas_dataset_with_2_cids[0]) == 2
    metatlas_dataset_with_2_cids.filter_compounds_by_signal(73, 2.30e6)
    assert len(metatlas_dataset_with_2_cids[0]) == 1
    metatlas_dataset_with_2_cids.filter_compounds_by_signal(73, 2.36e6)
    assert len(metatlas_dataset_with_2_cids[0]) == 0


def test_filter_compounds_by_signal02(mocker, metatlas_dataset_with_2_cids, df_container):
    mocker.patch(
        "metatlas.io.metatlas_get_data_helper_fun.df_container_from_metatlas_file", return_value=df_container
    )
    assert len(metatlas_dataset_with_2_cids[0]) == 2
    metatlas_dataset_with_2_cids.filter_compounds_by_signal(74, 1e5)
    assert len(metatlas_dataset_with_2_cids[0]) == 0


def test_export_atlas_to_csv01(metatlas_dataset, tmp_path):
    out_file = tmp_path / "export.csv"
    metatlas_dataset.export_atlas_to_csv(out_file)
    in_df = pd.read_csv(out_file)
    assert list(in_df.columns) == [
        "Unnamed: 0",
        "chebi_id",
        "chebi_url",
        "creation_time",
        "description",
        "formula",
        "head_id",
        "hmdb_id",
        "hmdb_url",
        "img_abc_id",
        "inchi",
        "inchi_key",
        "iupac_name",
        "kegg_id",
        "kegg_url",
        "last_modified",
        "lipidmaps_id",
        "lipidmaps_url",
        "metacyc_id",
        "mono_isotopic_molecular_weight",
        "name",
        "neutralized_2d_inchi",
        "neutralized_2d_inchi_key",
        "neutralized_inchi",
        "neutralized_inchi_key",
        "num_free_radicals",
        "number_components",
        "permanent_charge",
        "prev_uid",
        "pubchem_compound_id",
        "pubchem_url",
        "source",
        "synonyms",
        "unique_id",
        "username",
        "wikipedia_url",
        "label",
        "id_notes",
        "ms1_notes",
        "ms2_notes",
        "identification_notes",
        "rt_min",
        "rt_max",
        "rt_peak",
        "mz",
        "mz_tolerance",
        "adduct",
        "polarity",
    ]
    assert len(in_df) == 1
    assert in_df.loc[0, "inchi_key"] == "OLXZPDWKRNYJJZ-RRKCRQDMSA-N"


def test_atlas_setter01(metatlas_dataset, atlas_with_2_cids):
    metatlas_dataset.data  # pylint: disable=pointless-statement
    metatlas_dataset.atlas = atlas_with_2_cids
    assert metatlas_dataset._data is None
    assert len(metatlas_dataset[0]) == 2


def test_atlas_setter02(metatlas_dataset):
    with pytest.raises(traitlets.traitlets.TraitError):
        metatlas_dataset.atlas = [1, 2]


def test_groups01(metatlas_dataset):
    assert metatlas_dataset.ids.groups[0].short_name == "POS_Cone-S1"


def test_set_extra_mz_setter(metatlas_dataset, mocker, hits):
    mocker.patch("metatlas.plots.dill2plots.get_msms_hits", return_value=hits)
    metatlas_dataset.data  # pylint: disable=pointless-statement
    metatlas_dataset.hits  # pylint: disable=pointless-statement
    metatlas_dataset.extra_mz = 0.43
    assert metatlas_dataset._data is None
    assert metatlas_dataset._hits is None
    assert metatlas_dataset.extra_mz == 0.43


def test_set_keep_nonmatches_setter(metatlas_dataset, mocker, hits):
    mocker.patch("metatlas.plots.dill2plots.get_msms_hits", return_value=hits)
    metatlas_dataset.hits  # pylint: disable=pointless-statement
    metatlas_dataset.keep_nonmatches = False
    assert metatlas_dataset._hits is None
    assert not metatlas_dataset.keep_nonmatches


def test_set_frag_mz_tolerance_setter(metatlas_dataset, mocker, hits):
    mocker.patch("metatlas.plots.dill2plots.get_msms_hits", return_value=hits)
    metatlas_dataset.hits  # pylint: disable=pointless-statement
    metatlas_dataset.frag_mz_tolerance = 1e-4
    assert metatlas_dataset._hits is None
    assert metatlas_dataset.frag_mz_tolerance == 1e-4


def test_set_msms_refs_loc_setter(metatlas_dataset, mocker, hits):
    mocker.patch("metatlas.plots.dill2plots.get_msms_hits", return_value=hits)
    metatlas_dataset.hits  # pylint: disable=pointless-statement
    metatlas_dataset.msms_refs_loc = "/tmp/some_file.tab"
    assert metatlas_dataset._hits is None
    assert metatlas_dataset.msms_refs_loc == "/tmp/some_file.tab"


def test_store_atlas01(atlas, sqlite, username):
    atlas.name = "test_store_atlas01"
    atlas_list = metob.retrieve("Atlas", name=atlas.name, username=username)
    assert len(atlas_list) == 0
    metob.store(atlas)
    second = metob.retrieve("Atlas", name=atlas.name, username=username)
    assert len(second) == 1


def test_store_atlas02(metatlas_dataset, username):
    atlas_list = metob.retrieve("Atlas", name=metatlas_dataset.ids.source_atlas, username=username)
    assert len(atlas_list) == 1
    second = metob.retrieve("Atlas", name=metatlas_dataset.atlas.name, username=username)
    assert len(second) == 1
    metatlas_dataset.store_atlas(even_if_exists=True)
    second = metob.retrieve("Atlas", name=metatlas_dataset.atlas.name, username=username)
    assert len(second) == 1


def test_store_atlas03(metatlas_dataset, atlas, sqlite, username):
    metatlas_dataset.atlas.name = "test_store_atlas01"
    atlas_list = metob.retrieve("Atlas", name=metatlas_dataset.atlas.name, username=username)
    assert len(atlas_list) == 0
    metatlas_dataset.store_atlas()
    second = metob.retrieve("Atlas", name=metatlas_dataset.atlas.name, username=username)
    assert len(second) == 1


def test_store_atlas04(metatlas_dataset, sqlite, username):
    metatlas_dataset.atlas.name = "test_store_atlas01"
    atlas_list = metob.retrieve("Atlas", name=metatlas_dataset.atlas.name, username=username)
    assert len(atlas_list) == 0
    metatlas_dataset.store_atlas()
    second = metob.retrieve("Atlas", name=metatlas_dataset.atlas.name, username=username)
    assert len(second) == 1
    metatlas_dataset.store_atlas(even_if_exists=True)
    with pytest.raises(ValueError):
        metatlas_dataset.store_atlas()


def test_store_atlas05(atlas, sqlite, username):
    atlas.name = "test atlas"
    metob.store(atlas)
    second = metob.retrieve("Atlas", name=atlas.name, username=username)
    assert len(second) == 1


def test_store_atlas06(atlas, sqlite_with_atlas, username):
    atlas.name = "test atlas 06"
    metob.store(atlas)
    second = metob.retrieve("Atlas", name=atlas.name, username=username)
    assert len(second) == 1


def test_store_atlas07(atlas, sqlite, username):
    atlas.name = "test_store_atlas07"
    metob.store(atlas)
    metoh.Workspace.get_instance().close_connection()
    metoh.Workspace.instance = None
    atlases = metob.retrieve("Atlas", name=atlas.name, username=username)
    assert len(atlases) == 1


def test_analysis_identifiers01(sqlite):
    with pytest.raises(traitlets.traitlets.TraitError, match=r"Database does not contain an atlas.*"):
        mads.AnalysisIdentifiers(
            source_atlas="Not_A_Real_Atlas_Name",
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            output_type="FinalEMA-HILIC",
            polarity="positive",
            analysis_number=1,
            project_directory="/foo/bar",
        )


def test_analysis_identifiers02(sqlite_with_atlas, username):
    with pytest.raises(
        traitlets.traitlets.TraitError,
        match="Parameter output_type must be one of ISTDsEtc, FinalEMA-HILIC, data_QC",
    ):
        mads.AnalysisIdentifiers(
            source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            output_type="output_type_not_valid",
            polarity="positive",
            analysis_number=1,
            project_directory="/foo/bar",
        )


def test_analysis_identifiers03(username, sqlite_with_atlas):
    with pytest.raises(
        traitlets.traitlets.TraitError,
        match="Parameter polarity must be one of positive, negative, fast-polarity-switching",
    ):
        mads.AnalysisIdentifiers(
            source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            output_type="FinalEMA-HILIC",
            polarity="not a polarity value",
            analysis_number=1,
            project_directory="/foo/bar",
        )


def test_analysis_identifiers04(username, sqlite_with_atlas):
    with pytest.raises(
        traitlets.traitlets.TraitError,
        match="The 'analysis_number' trait of an AnalysisIdentifiers instance expected an int, not",
    ):
        mads.AnalysisIdentifiers(
            source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            output_type="FinalEMA-HILIC",
            polarity="positive",
            analysis_number="this is a string",
            project_directory="/foo/bar",
        )


def test_analysis_identifiers05(username, sqlite_with_atlas):
    with pytest.raises(
        traitlets.traitlets.TraitError,
        match="The 'analysis_number' trait of an AnalysisIdentifiers instance expected an int, not",
    ):
        mads.AnalysisIdentifiers(
            source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            output_type="FinalEMA-HILIC",
            polarity="positive",
            analysis_number="1",
            project_directory="/foo/bar",
        )


def test_analysis_identifiers06(username, sqlite_with_atlas):
    with pytest.raises(traitlets.traitlets.TraitError, match="Parameter analysis_number cannot be negative."):
        mads.AnalysisIdentifiers(
            source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
            experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
            output_type="FinalEMA-HILIC",
            polarity="positive",
            analysis_number=-9,
            project_directory="/foo/bar",
        )


def test_analysis_identifiers07(username, sqlite_with_atlas):
    with pytest.raises(
        traitlets.traitlets.TraitError, match='Parameter experiment does contain 9 fields when split on "_".'
    ):
        mads.AnalysisIdentifiers(
            source_atlas=f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
            experiment="experiment_name_not_valid",
            output_type="FinalEMA-HILIC",
            polarity="positive",
            analysis_number=0,
            project_directory="/foo/bar",
        )


def test_analysis_identifiers_atlas01(analysis_ids, username):
    assert analysis_ids.atlas == f"505892_OakGall_final_FinalEMA-HILIC_POS_{username}0"


def test_analysis_identifiers_atlas02(analysis_ids, username):
    # call .atlas twice to get cached value
    analysis_ids.atlas  # pylint: disable=pointless-statement
    assert analysis_ids.atlas == f"505892_OakGall_final_FinalEMA-HILIC_POS_{username}0"


def test_write_data_source_files01(metatlas_dataset, mocker, caplog):
    mocker.patch("glob.glob", return_value=range(10))
    metatlas_dataset.write_data_source_files()
    assert "Data sources directory already populated" in caplog.text


def test_write_data_source_files02(metatlas_dataset, mocker, caplog):
    mocker.patch("glob.glob", return_value=range(3))
    mocker.patch("shutil.rmtree")
    mocker.patch("metatlas.io.metatlas_get_data_helper_fun.make_data_sources_tables")
    caplog.set_level(logging.INFO)
    metatlas_dataset.write_data_source_files()
    assert "Writing data source files to" in caplog.text
    assert ma_data.make_data_sources_tables.called  # pylint: disable=no-member


def test_get_atlas01(mocker, analysis_ids, df_container, lcmsrun, atlas, username):
    mocker.patch(
        "metatlas.io.metatlas_get_data_helper_fun.df_container_from_metatlas_file", return_value=df_container
    )
    mocker.patch("metatlas.plots.dill2plots.get_metatlas_files", return_value=[lcmsrun])
    mocker.patch("glob.glob", return_value=range(10))
    metatlas_dataset = mads.MetatlasDataset(ids=analysis_ids)
    assert metatlas_dataset.atlas.name == f"505892_OakGall_final_FinalEMA-HILIC_POS_{username}0"


def test_get_atlas02(mocker, analysis_ids, caplog):
    mocker.patch("metatlas.datastructures.metatlas_objects.retrieve", return_value=[])
    caplog.set_level(logging.INFO)
    with pytest.raises(ValueError):
        mads.MetatlasDataset(ids=analysis_ids)
    assert "Database does not contain an atlas" in caplog.text


def test_get_atlas03(mocker, analysis_ids, caplog, username):
    mocker.patch("metatlas.datastructures.metatlas_objects.retrieve", return_value=[0, 0])
    with pytest.raises(ValueError):
        mads.MetatlasDataset(ids=analysis_ids)
    atlas = f"505892_OakGall_final_FinalEMA-HILIC_POS_{username}0"
    assert f"2 atlases with name {atlas} and owned by {username} already exist." in caplog.text


def test_get_atlas04(metatlas_dataset, username):
    atlases = metob.retrieve("Atlas", name="This_atlas_does_not_exists", username=username)
    assert len(atlases) == 0


def test_existing_groups(mocker, metatlas_dataset):
    """This test has little value, but is needed for coverage"""
    mocker.patch("metatlas.datastructures.metatlas_objects.retrieve", return_value=[])
    assert metatlas_dataset.ids.existing_groups == []


def test_lcmsruns_dataframe(metatlas_dataset):
    assert metatlas_dataset.ids.lcmsruns_dataframe.shape == (1, 15)


def test_store_groups01(metatlas_dataset, mocker):
    mocker.patch("metatlas.datastructures.metatlas_objects.retrieve", return_value=[])
    mocker.patch("metatlas.datastructures.metatlas_objects.store")
    metatlas_dataset.ids.store_all_groups()
    assert metob.store.called  # pylint: disable=no-member


def test_store_groups02(metatlas_dataset, mocker, username):
    def group():
        pass

    group.name = (
        f"20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_{username}0_Cone-S1"
    )
    mocker.patch("metatlas.datastructures.metatlas_objects.retrieve", return_value=[group])
    with pytest.raises(ValueError):
        metatlas_dataset.ids.store_all_groups()


def test_annotation_gui01(metatlas_dataset, hits, mocker):
    mocker.patch("metatlas.plots.dill2plots.get_msms_hits", return_value=hits)
    agui = metatlas_dataset.annotation_gui()
    agui.compound_idx = 0
    agui.set_msms_flag("1, co-isolated precursor but all reference ions are in sample spectrum")
    agui.set_peak_flag("remove")
    agui.data.set_rt(0, "rt_min", 2.1245)
    agui.data.set_rt(0, "rt_max", 2.4439)
    assert metatlas_dataset.rts[0].rt_min == 2.1245
    assert metatlas_dataset.rts[0].rt_max == 2.4439
    assert metatlas_dataset.data[0][0]["identification"].ms1_notes == "remove"
    assert (
        metatlas_dataset.data[0][0]["identification"].ms2_notes
        == "1, co-isolated precursor but all reference ions are in sample spectrum"
    )


def test_generate_all_outputs01(metatlas_dataset, hits, mocker):
    mocker.patch("metatlas.plots.dill2plots.get_msms_hits", return_value=hits)
    metatlas_dataset.generate_all_outputs()
    assert len(glob.glob(metatlas_dataset.ids.output_dir + "/*")) == 10
    assert len(glob.glob(metatlas_dataset.ids.output_dir + "/*/*")) == 19


def test_short_polarity_inverse01(analysis_ids):
    assert set(analysis_ids.short_polarity_inverse) == {"NEG", "FPS"}


def test_access_data_compound_name(metatlas_dataset):
    assert metatlas_dataset.data[0][0]["identification"].name == "2'-deoxyadenosine"


def test_cid_type01(atlas):
    assert isinstance(atlas.compound_identifications[0], metob.CompoundIdentification)


def test_load_atlas01(atlas, sqlite_with_atlas, username):
    atlases = metob.retrieve("Atlas", name=atlas.name, username=username)
    assert isinstance(atlases[0].compound_identifications[0], metob.CompoundIdentification)


def test_load_atlas02(atlas, sqlite_with_atlas, username):
    results = metob.retrieve("Atlas", name=atlas.name, username=username)
    assert isinstance(results[0].compound_identifications[0], metob.CompoundIdentification)


def test_load_atlas03(sqlite_with_atlas, atlas, username):
    results = metob.retrieve("Atlas", name=atlas.name, username=username)
    assert results[0].compound_identifications[0].rt_references[0].rt_peak == 2.1964640053707174


def test_invalidation01(analysis_ids):
    _ = analysis_ids.groups
    assert analysis_ids._groups is not None
    analysis_ids.exclude_files = ["Cone-S1"]
    assert analysis_ids._groups is None


def test_negative_polarity01(sqlite_with_atlas, username, lcmsrun, mocker, groups_controlled_vocab):
    mocker.patch("metatlas.plots.dill2plots.get_metatlas_files", return_value=[lcmsrun])
    ids = mads.AnalysisIdentifiers(
        experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
        output_type="FinalEMA-HILIC",
        polarity="negative",
        analysis_number=0,
        project_directory=str(os.getcwd()),
        groups_controlled_vocab=groups_controlled_vocab,
    )
    assert "POS" in ids.exclude_groups


def test_include_groups01(sqlite_with_atlas, username, lcmsrun, mocker, groups_controlled_vocab):
    mocker.patch("metatlas.plots.dill2plots.get_metatlas_files", return_value=[lcmsrun])
    ids = mads.AnalysisIdentifiers(
        experiment="20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
        output_type="data_QC",
        polarity="negative",
        analysis_number=0,
        project_directory=str(os.getcwd()),
        groups_controlled_vocab=groups_controlled_vocab,
    )
    assert "QC" in ids.include_groups


def test_project01(analysis_ids):
    assert analysis_ids.project == 505892


def test_exclude_files01(analysis_ids):
    analysis_ids.exclude_files = ["POS"]
    assert len(analysis_ids.lcmsruns) == 0
    assert analysis_ids.lcmsruns_short_names.empty


def test_invlidate_groups_controlled_vocab01(analysis_ids):
    _ = analysis_ids.lcmsruns
    assert analysis_ids._lcmsruns is not None
    analysis_ids.groups_controlled_vocab = ["FOOBAR"]
    assert analysis_ids._lcmsruns is None
