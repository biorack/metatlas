"""
per-directory pytest configuration
fixtures used across multiple files should go in here
"""
# pylint: disable=missing-function-docstring,unused-argument,line-too-long,too-many-lines,too-many-arguments

import getpass
import logging
import os
import sqlite3

from importlib import reload

import pytest
import numpy as np
import pandas as pd

from sqlalchemy.orm import close_all_sessions

from metatlas.datastructures import metatlas_dataset as mads
from metatlas.datastructures import metatlas_objects as metob
from metatlas.datastructures import object_helpers as metoh


logger = logging.getLogger(__name__)


@pytest.fixture(name="username", scope="session")
def fixture_username():
    return getpass.getuser()


@pytest.fixture(name="analysis_ids")
def fixture_analysis_ids(sqlite_with_atlas, username):
    return mads.AnalysisIdentifiers(
        f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0",
        "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
        "FinalEMA-HILIC",
        "positive",
        0,
        str(os.getcwd()),
    )


@pytest.fixture(name="analysis_ids_with_2_cids")
def fixture_analysis_ids_with_2_cids(sqlite_with_atlas_with_2_cids, username):
    return mads.AnalysisIdentifiers(
        f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}1",
        "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
        "FinalEMA-HILIC",
        "positive",
        0,
        str(os.getcwd()),
    )


@pytest.fixture(name="sqlite")
def fixture_sqlite(username, change_test_dir, atlas):
    logging.debug("creating database file in %s", os.getcwd())
    assert not os.path.exists(f"{username}_workspace.db")
    sqlite3.connect(f"{username}_workspace.db").close()
    logger.debug("reloading metoh")
    reload(metoh)
    logger.debug("Storing empty objects to create tables")
    metob.store(metob.Atlas())
    metob.store(metob.CompoundIdentification())
    metob.store(metob.Compound())
    metob.store(metob.MzReference())
    metob.store(metob.RtReference())
    metob.store(metob.LcmsRun())
    logger.debug("Done storing empty objects to create tables")
    yield
    metob.workspace.close_connection()
    close_all_sessions()


@pytest.fixture(name="sqlite_with_atlas")
def fixture_sqlite_with_atlas(sqlite, atlas, username):
    atlas.name = f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}0"
    logger.debug("Saving atlas %s", atlas.name)
    metob.store(atlas)


@pytest.fixture(name="sqlite_with_atlas_with_2_cids")
def fixture_sqlite_with_atlas_with_2_cids(sqlite, atlas_with_2_cids, username):
    atlas_with_2_cids.name = f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}1"
    logger.debug("Saving atlas %s", atlas_with_2_cids.name)
    metob.store(atlas_with_2_cids)


@pytest.fixture(name="change_test_dir", scope="function", autouse=True)
def fixture_change_test_dir(request, tmp_path):
    os.chdir(tmp_path)
    logger.debug("changing dir to %s", tmp_path)
    yield
    os.chdir(request.config.invocation_dir)


@pytest.fixture(name="ms1_pos")
def fixture_ms1_pos():
    return pd.DataFrame(
        data={
            "mz": {
                "0": 252.1089324951,
                "1": 252.1090087891,
                "2": 252.1088104248,
                "3": 252.1090087891,
                "4": 252.10887146,
                "5": 252.1089324951,
                "6": 252.1089324951,
                "7": 252.1088256836,
                "8": 252.1088867188,
                "9": 252.1090393066,
                "10": 252.1089782715,
                "11": 252.1089630127,
                "12": 252.1089630127,
                "13": 252.1089782715,
                "14": 252.1090240479,
                "15": 252.1089782715,
                "16": 252.1090240479,
                "17": 252.1089324951,
                "18": 252.1090393066,
                "19": 252.1088867188,
                "20": 252.10887146,
                "21": 252.1089324951,
                "22": 252.1089630127,
                "23": 252.1089935303,
                "24": 252.1089172363,
                "25": 252.1089477539,
                "26": 252.1090545654,
                "27": 252.1089630127,
                "28": 252.1090240479,
                "29": 252.1090087891,
                "30": 252.1090393066,
                "31": 252.1090240479,
                "32": 252.1089935303,
                "33": 252.1090240479,
                "34": 252.1089630127,
                "35": 252.1090087891,
                "36": 252.1090240479,
                "37": 252.1089172363,
                "38": 252.1089019775,
                "39": 252.1089477539,
                "40": 252.1089324951,
                "41": 252.1089477539,
                "42": 252.1089477539,
                "43": 252.1089477539,
                "44": 252.1089782715,
                "45": 252.1088867188,
                "46": 252.1089172363,
                "47": 252.1089324951,
                "48": 252.1089782715,
                "49": 252.1089477539,
                "50": 252.1089172363,
                "51": 252.1089324951,
                "52": 252.1089630127,
                "53": 252.1088867188,
                "54": 252.1089630127,
                "55": 252.1085205078,
                "56": 252.1090545654,
                "57": 252.1089935303,
                "58": 252.1088104248,
                "59": 252.1086578369,
                "60": 252.1089935303,
                "61": 252.1085510254,
                "62": 252.1082763672,
                "63": 252.1082458496,
                "64": 252.1084136963,
                "65": 252.1092224121,
                "66": 252.1091766357,
                "67": 252.1092834473,
                "68": 252.1087493896,
                "69": 252.1112518311,
                "70": 252.1088409424,
                "71": 252.1086425781,
                "72": 252.1091766357,
                "73": 252.1094055176,
            },
            "i": {
                "0": 312203.5,
                "1": 387914.59375,
                "2": 308308.5,
                "3": 334653.59375,
                "4": 339521.625,
                "5": 345527.21875,
                "6": 292437.34375,
                "7": 413614.53125,
                "8": 300285.28125,
                "9": 383848.71875,
                "10": 404313.21875,
                "11": 377231.34375,
                "12": 453965.5625,
                "13": 431327.0,
                "14": 523180.0625,
                "15": 510239.8125,
                "16": 631459.1875,
                "17": 807419.5,
                "18": 842647.5625,
                "19": 1053031.625,
                "20": 1082361.625,
                "21": 1198966.625,
                "22": 1109162.375,
                "23": 1126347.125,
                "24": 1373071.5,
                "25": 1589018.375,
                "26": 1281309.875,
                "27": 1660166.75,
                "28": 1492912.25,
                "29": 2029801.5,
                "30": 2029874.125,
                "31": 2035966.625,
                "32": 2010867.875,
                "33": 2036981.375,
                "34": 2148879.25,
                "35": 2359861.25,
                "36": 2054066.125,
                "37": 1691976.0,
                "38": 1778159.125,
                "39": 1776166.125,
                "40": 1752154.125,
                "41": 1575676.875,
                "42": 1199910.625,
                "43": 1259708.25,
                "44": 1087384.375,
                "45": 826077.125,
                "46": 802296.875,
                "47": 547785.125,
                "48": 545340.0625,
                "49": 584624.4375,
                "50": 468524.8125,
                "51": 305931.1875,
                "52": 330310.34375,
                "53": 309740.625,
                "54": 289212.71875,
                "55": 230440.9375,
                "56": 210549.390625,
                "57": 169972.390625,
                "58": 140521.234375,
                "59": 116637.953125,
                "60": 117197.625,
                "61": 84652.1171875,
                "62": 117615.578125,
                "63": 103500.921875,
                "64": 89320.9453125,
                "65": 76313.9296875,
                "66": 55575.00390625,
                "67": 76784.6796875,
                "68": 28829.162109375,
                "69": 26051.6171875,
                "70": 42957.18359375,
                "71": 50342.6953125,
                "72": 37611.33984375,
                "73": 38202.83203125,
            },
            "rt": {
                "0": 2.1030805111,
                "1": 2.1084616184,
                "2": 2.1139531136,
                "3": 2.1193552017,
                "4": 2.1248509884,
                "5": 2.1302509308,
                "6": 2.135682106,
                "7": 2.1411821842,
                "8": 2.1459801197,
                "9": 2.1513926983,
                "10": 2.1568279266,
                "11": 2.1622362137,
                "12": 2.1676549911,
                "13": 2.1730883121,
                "14": 2.179015398,
                "15": 2.1845297813,
                "16": 2.1900422573,
                "17": 2.1949694157,
                "18": 2.20002985,
                "19": 2.2055358887,
                "20": 2.2110378742,
                "21": 2.2165191174,
                "22": 2.2219588757,
                "23": 2.2273921967,
                "24": 2.2328462601,
                "25": 2.2382712364,
                "26": 2.2437169552,
                "27": 2.2492566109,
                "28": 2.2547125816,
                "29": 2.2601687908,
                "30": 2.2656960487,
                "31": 2.2704958916,
                "32": 2.2758042812,
                "33": 2.2813498974,
                "34": 2.2868082523,
                "35": 2.2922415733,
                "36": 2.2976748943,
                "37": 2.3031060696,
                "38": 2.308131218,
                "39": 2.313628912,
                "40": 2.3185498714,
                "41": 2.3239560127,
                "42": 2.3293914795,
                "43": 2.3349123001,
                "44": 2.3403663635,
                "45": 2.346799612,
                "46": 2.3522267342,
                "47": 2.3576600552,
                "48": 2.3631224632,
                "49": 2.3685662746,
                "50": 2.3740911484,
                "51": 2.3794057369,
                "52": 2.3848536015,
                "53": 2.3903660774,
                "54": 2.3953785896,
                "55": 2.4006638527,
                "56": 2.4062638283,
                "57": 2.411709547,
                "58": 2.4171659946,
                "59": 2.4226117134,
                "60": 2.4302260876,
                "61": 2.4357616901,
                "62": 2.4407405853,
                "63": 2.4461927414,
                "64": 2.451615572,
                "65": 2.4571509361,
                "66": 2.4627010822,
                "67": 2.4681572914,
                "68": 2.4735822678,
                "69": 2.4735822678,
                "70": 2.4787945747,
                "71": 2.4842174053,
                "72": 2.4896612167,
                "73": 2.495146513,
            },
            "polarity": {
                "0": 1,
                "1": 1,
                "2": 1,
                "3": 1,
                "4": 1,
                "5": 1,
                "6": 1,
                "7": 1,
                "8": 1,
                "9": 1,
                "10": 1,
                "11": 1,
                "12": 1,
                "13": 1,
                "14": 1,
                "15": 1,
                "16": 1,
                "17": 1,
                "18": 1,
                "19": 1,
                "20": 1,
                "21": 1,
                "22": 1,
                "23": 1,
                "24": 1,
                "25": 1,
                "26": 1,
                "27": 1,
                "28": 1,
                "29": 1,
                "30": 1,
                "31": 1,
                "32": 1,
                "33": 1,
                "34": 1,
                "35": 1,
                "36": 1,
                "37": 1,
                "38": 1,
                "39": 1,
                "40": 1,
                "41": 1,
                "42": 1,
                "43": 1,
                "44": 1,
                "45": 1,
                "46": 1,
                "47": 1,
                "48": 1,
                "49": 1,
                "50": 1,
                "51": 1,
                "52": 1,
                "53": 1,
                "54": 1,
                "55": 1,
                "56": 1,
                "57": 1,
                "58": 1,
                "59": 1,
                "60": 1,
                "61": 1,
                "62": 1,
                "63": 1,
                "64": 1,
                "65": 1,
                "66": 1,
                "67": 1,
                "68": 1,
                "69": 1,
                "70": 1,
                "71": 1,
                "72": 1,
                "73": 1,
            },
        }
    )


@pytest.fixture(name="ms2_pos")
def fixture_ms2_pos():
    return pd.DataFrame(
        data={
            "mz": {
                "0": 252.1081695557,
                "1": 252.1564941406,
                "2": 252.1087036133,
                "3": 252.1572875977,
                "4": 252.1089019775,
                "5": 252.1550292969,
                "6": 252.1090698242,
                "7": 252.1557617188,
            },
            "i": {
                "0": 32103.3515625,
                "1": 6470.0009765625,
                "2": 93112.0859375,
                "3": 7624.11328125,
                "4": 131062.0,
                "5": 6535.4560546875,
                "6": 76976.7265625,
                "7": 6090.6440429688,
            },
            "rt": {
                "0": 2.0097544193,
                "1": 2.0097544193,
                "2": 2.2203779221,
                "3": 2.2203779221,
                "4": 2.327804327,
                "5": 2.327804327,
                "6": 2.3452186584,
                "7": 2.3452186584,
            },
            "polarity": {"0": 1, "1": 1, "2": 1, "3": 1, "4": 1, "5": 1, "6": 1, "7": 1},
            "precursor_MZ": {
                "0": 252.0195159912,
                "1": 252.0195159912,
                "2": 252.10887146,
                "3": 252.10887146,
                "4": 252.0194854736,
                "5": 252.0194854736,
                "6": 252.1089477539,
                "7": 252.1089477539,
            },
            "precursor_intensity": {
                "0": 2748235.5,
                "1": 2748235.5,
                "2": 2872807.5,
                "3": 2872807.5,
                "4": 3536752.25,
                "5": 3536752.25,
                "6": 3046732.75,
                "7": 3046732.75,
            },
            "collision_energy": {
                "0": 23.3333339691,
                "1": 23.3333339691,
                "2": 23.3333339691,
                "3": 23.3333339691,
                "4": 23.3333339691,
                "5": 23.3333339691,
                "6": 23.3333339691,
                "7": 23.3333339691,
            },
        }
    )


@pytest.fixture(name="ms1_neg_empty")
def fixture_ms1_neg_empty():
    return pd.DataFrame(data={"mz": {}, "i": {}, "rt": {}, "polarity": {}})


@pytest.fixture(name="ms2_neg_empty")
def fixture_ms2_neg_empty():
    return pd.DataFrame(
        data={
            "mz": {},
            "i": {},
            "rt": {},
            "polarity": {},
            "precursor_MZ": {},
            "precursor_intensity": {},
            "collision_energy": {},
        }
    )


@pytest.fixture(name="df_container")
def fixture_df_container(ms1_pos, ms2_pos, ms1_neg_empty, ms2_neg_empty):
    return {"ms1_neg": ms1_neg_empty, "ms1_pos": ms1_pos, "ms2_neg": ms2_neg_empty, "ms2_pos": ms2_pos}


@pytest.fixture(name="ms1_summary")
def fixture_ms1_summary():
    return {
        "num_ms1_datapoints": 85.0,
        "mz_peak": 252.1092987060547,
        "rt_peak": 2.2775044441223145,
        "mz_centroid": 252.10915042669814,
        "rt_centroid": 2.218492414487913,
        "peak_height": 304761.90625,
        "peak_area": 7696977.46875,
    }


@pytest.fixture(name="msms")
def fixture_msms():
    return {
        "data": {
            "mz": np.array([], dtype=np.float64),
            "i": np.array([], dtype=np.float64),
            "rt": np.array([], dtype=np.float64),
            "polarity": np.array([], dtype=np.float64),
            "precursor_MZ": np.array([], dtype=np.float64),
            "precursor_intensity": np.array([], dtype=np.float64),
            "collision_energy": np.array([], dtype=np.float64),
        }
    }


@pytest.fixture(name="groups_controlled_vocab")
def fixture_groups_controlled_vocab():
    return ["QC", "InjBl", "ISTD"]


@pytest.fixture(name="metatlas_dataset")
def fixture_metatlas_dataset(
    mocker, df_container, analysis_ids, groups_controlled_vocab, lcmsrun, sqlite_with_atlas
):
    mocker.patch(
        "metatlas.io.metatlas_get_data_helper_fun.df_container_from_metatlas_file", return_value=df_container
    )
    mocker.patch("metatlas.plots.dill2plots.get_metatlas_files", return_value=[lcmsrun])
    return mads.MetatlasDataset(analysis_ids, groups_controlled_vocab, save_metadata=False)


@pytest.fixture(name="metatlas_dataset_with_2_cids")
def fixture_metatlas_dataset_with_2_cids(
    mocker,
    df_container,
    analysis_ids_with_2_cids,
    groups_controlled_vocab,
    lcmsrun,
    sqlite_with_atlas_with_2_cids,
):
    mocker.patch(
        "metatlas.io.metatlas_get_data_helper_fun.df_container_from_metatlas_file", return_value=df_container
    )
    mocker.patch("metatlas.plots.dill2plots.get_metatlas_files", return_value=[lcmsrun])
    return mads.MetatlasDataset(analysis_ids_with_2_cids, groups_controlled_vocab, save_metadata=False)


@pytest.fixture(name="eic")
def fixture_eic():
    return {
        "mz": [
            252.1089324951172,
            252.10943603515625,
            252.10926818847656,
            252.109375,
            252.10923767089844,
            252.10910034179688,
            252.10914611816406,
            252.1089630126953,
            252.10971069335938,
            252.1093292236328,
            252.10934448242188,
            252.109130859375,
            252.10935974121094,
            252.10939025878906,
            252.1090545654297,
            252.10916137695312,
            252.10946655273438,
            252.10923767089844,
            252.1093292236328,
            252.10919189453125,
            252.10914611816406,
            252.10897827148438,
            252.10934448242188,
            252.10928344726562,
            252.10888671875,
            252.10926818847656,
            252.109130859375,
            252.1090087890625,
            252.10934448242188,
            252.10939025878906,
            252.1093292236328,
            252.1091766357422,
            252.109130859375,
            252.1095428466797,
            252.10890197753906,
            252.1095428466797,
            252.109130859375,
            252.10911560058594,
            252.1091766357422,
            252.1088409423828,
            252.10916137695312,
            252.10935974121094,
            252.10928344726562,
            252.10922241210938,
            252.10914611816406,
            252.10922241210938,
            252.10894775390625,
            252.10906982421875,
            252.10914611816406,
            252.10916137695312,
            252.10910034179688,
            252.10916137695312,
            252.10934448242188,
            252.10899353027344,
            252.10928344726562,
            252.10897827148438,
            252.10916137695312,
            252.10928344726562,
            252.1092987060547,
            252.1089324951172,
            252.10914611816406,
            252.1090545654297,
            252.10914611816406,
            252.1090850830078,
            252.10894775390625,
            252.10914611816406,
            252.10911560058594,
            252.1090850830078,
            252.109130859375,
            252.10903930664062,
            252.10890197753906,
            252.109130859375,
            252.10885620117188,
            252.10914611816406,
            252.10926818847656,
            252.10888671875,
            252.109619140625,
            252.10922241210938,
            252.1092529296875,
            252.1099853515625,
            252.10972595214844,
            252.10910034179688,
            252.10935974121094,
            252.1088409423828,
            252.10838317871094,
            252.11212158203125,
        ],
        "rt": [
            1.7180122137069702,
            1.8222843408584595,
            1.838305115699768,
            1.8444031476974487,
            1.8705799579620361,
            1.875998616218567,
            1.8913277387619019,
            1.9020838737487793,
            1.9127358198165894,
            1.9397128820419312,
            1.9451169967651367,
            1.9505127668380737,
            1.955920934677124,
            1.966427206993103,
            1.9718105792999268,
            1.9769750833511353,
            1.9823375940322876,
            1.987752079963684,
            1.9932082891464233,
            1.9986457824707031,
            2.0094456672668457,
            2.019866466522217,
            2.030582904815674,
            2.036003589630127,
            2.0568389892578125,
            2.062201499938965,
            2.0675911903381348,
            2.0834577083587646,
            2.088857650756836,
            2.0939910411834717,
            2.099109649658203,
            2.104536771774292,
            2.1208388805389404,
            2.1262447834014893,
            2.1420176029205322,
            2.152921676635742,
            2.15836763381958,
            2.163788318634033,
            2.169198751449585,
            2.1755259037017822,
            2.180954933166504,
            2.18635892868042,
            2.191038131713867,
            2.1964569091796875,
            2.2018840312957764,
            2.2069132328033447,
            2.21236515045166,
            2.2177650928497314,
            2.2228589057922363,
            2.2283151149749756,
            2.2338151931762695,
            2.239321231842041,
            2.244842052459717,
            2.250317096710205,
            2.255610704421997,
            2.261033535003662,
            2.2665293216705322,
            2.2720251083374023,
            2.2775044441223145,
            2.28295636177063,
            2.288454294204712,
            2.29386043548584,
            2.299298048019409,
            2.304720878601074,
            2.310127019882202,
            2.3155603408813477,
            2.320981025695801,
            2.326420545578003,
            2.33160400390625,
            2.3370935916900635,
            2.3428516387939453,
            2.3483099937438965,
            2.3535475730895996,
            2.3589975833892822,
            2.364443302154541,
            2.3699119091033936,
            2.375347375869751,
            2.3808369636535645,
            2.3862972259521484,
            2.3917577266693115,
            2.397282600402832,
            2.402780294418335,
            2.4081971645355225,
            2.419055461883545,
            2.457223892211914,
            3.3080079555511475,
        ],
        "intensity": [
            34249.71484375,
            28511.658203125,
            41718.13671875,
            33448.546875,
            40190.94140625,
            32525.16015625,
            37058.60546875,
            51132.91015625,
            36473.0546875,
            42659.0859375,
            45187.6171875,
            51186.30078125,
            58456.5859375,
            43299.24609375,
            52062.02734375,
            42501.8671875,
            39734.91015625,
            41848.02734375,
            48979.640625,
            42957.48046875,
            54214.27734375,
            63583.64453125,
            38661.046875,
            47146.54296875,
            36974.3046875,
            37674.35546875,
            37412.4609375,
            47036.44921875,
            32295.888671875,
            39751.12109375,
            47359.0,
            57496.41796875,
            33690.4765625,
            36853.53515625,
            33045.0703125,
            33235.64453125,
            52481.1015625,
            48210.37109375,
            62178.734375,
            73049.2109375,
            52741.03125,
            88225.1953125,
            101593.296875,
            127965.625,
            124079.859375,
            134410.46875,
            148749.0,
            134068.8125,
            141625.515625,
            202721.015625,
            204341.703125,
            172160.484375,
            185859.765625,
            195729.234375,
            216657.453125,
            239248.65625,
            172232.296875,
            195105.046875,
            304761.90625,
            181052.265625,
            222467.5625,
            251571.53125,
            205874.765625,
            224279.0625,
            173697.359375,
            236325.078125,
            153999.28125,
            156835.59375,
            118963.8046875,
            105766.234375,
            103081.484375,
            97180.5625,
            95681.4140625,
            74239.0703125,
            69208.8984375,
            60604.1484375,
            37020.84765625,
            32874.484375,
            24641.875,
            23305.75,
            23413.94140625,
            42582.77734375,
            35980.16796875,
            25743.97265625,
            21777.99609375,
            59454.40234375,
        ],
    }


@pytest.fixture(name="atlas_df")
def fixture_atlas_df(metatlas_dataset):
    return metatlas_dataset.atlas_df


@pytest.fixture(name="compound")
def fixture_compound(username):
    compound = metob.Compound()
    compound.unique_id = "60cd6743e56545c6a6cb066ec3553450"
    compound.mono_isotopic_molecular_weight = 251.101839276
    compound.creation_time = 1466212395
    compound.synonyms = "2'-deoxyadenosine"  # value was pruned down
    compound.inchi_key = "OLXZPDWKRNYJJZ-RRKCRQDMSA-N"
    compound.chebi_url = "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:17256"
    compound.permanent_charge = 0
    compound.img_abc_id = ""
    compound.neutralized_2d_inchi = "InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)"  # noqa: E501
    compound.lipidmaps_url = ""
    compound.source = "gnps///chebi///metacyc///hmdb"
    compound.kegg_url = "http://www.genome.jp/dbget-bin/www_bget?C00559"
    compound.hmdb_url = "http://www.hmdb.ca/metabolites/HMDB00101"
    compound.wikipedia_url = ""
    compound.head_id = "60cd6743e56545c6a6cb066ec3553450"
    compound.formula = "C10H13N5O3"
    compound.number_components = 1
    compound.iupac_name = ""
    compound.username = username
    compound.pubchem_compound_id = "13730"
    compound.description = "A purine 2'-deoxyribonucleoside having adenine as the nucleobase."
    compound.metacyc_id = "DEOXYADENOSINE"
    compound.kegg_id = "C00559"
    compound.hmdb_id = "HMDB00101"
    compound.chebi_id = "CHEBI:17256"
    compound.inchi = "InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1"  # noqa: E501
    compound.neutralized_inchi_key = "OLXZPDWKRNYJJZ-RRKCRQDMSA-N"
    compound.prev_uid = "origin"
    compound.neutralized_inchi = "InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1"  # noqa: E501
    compound.name = "2'-deoxyadenosine"
    compound.neutralized_2d_inchi_key = "OLXZPDWKRNYJJZ-UHFFFAOYSA-N"
    compound.num_free_radicals = 0
    compound.lipidmaps_id = ""
    compound.last_modified = 1612996604
    compound.pubchem_url = "http://pubchem.ncbi.nlm.nih.gov/compound/13730"
    return compound


@pytest.fixture(name="rt_reference")
def fixture_rt_reference(username):
    rt_ref = metob.RtReference()
    rt_ref.unique_id = "a845ddfdf8ef4713bcef3bdb84999030"
    rt_ref.username = username
    rt_ref.rt_units = "min"
    rt_ref.description = "No description"
    rt_ref.rt_peak = "2.1964640053707174"
    rt_ref.enabled = True
    rt_ref.creation_time = 1613002850
    rt_ref.lcms_run = None
    rt_ref.rt_min = 1.6964640053707174
    rt_ref.last_modified = 1613002979
    rt_ref.ref_type = ""
    rt_ref.prev_uid = "origin"
    rt_ref.rt_max = 2.6964640053707174
    rt_ref.name = "Untitled"
    rt_ref.head_id = "a845ddfdf8ef4713bcef3bdb84999030"
    return rt_ref


@pytest.fixture(name="mz_reference")
def fixture_mz_reference(username):
    mz_ref = metob.MzReference()
    mz_ref.unique_id = "eb6d03c9ef574051b92dad7b2fc259a2"
    mz_ref.username = username
    mz_ref.adduct = "[M+H]+"
    mz_ref.description = "No description"
    mz_ref.mz_tolerance_units = "ppm"
    mz_ref.enabled = True
    mz_ref.mz = 252.1091393
    mz_ref.creation_time = 1613002850
    mz_ref.lcms_run = None
    mz_ref.mz_tolerance = 20.0
    mz_ref.last_modified = 1613002979
    mz_ref.detected_polarity = "positive"
    mz_ref.modification = ""
    mz_ref.ref_type = ""
    mz_ref.observed_formula = ""
    mz_ref.prev_uid = "origin"
    mz_ref.name = "Untitled"
    mz_ref.head_id = "eb6d03c9ef574051b92dad7b2fc259a2"
    return mz_ref


@pytest.fixture(name="compound_identification")
def fixture_compound_identification(compound, rt_reference, mz_reference, username):
    ident = metob.CompoundIdentification()
    ident.unique_id = "18737c7141cc4efaa4545bead13ac751"
    ident.username = username
    ident.description = "No description"
    ident.creation_time = 1613002849
    ident.last_modified = 1613002979
    ident.identification_grade = None
    ident.prev_uid = "origin"
    ident.name = "2'-deoxyadenosine"
    ident.head_id = "18737c7141cc4efaa4545bead13ac751"
    ident.internal_standard_to_use = ""
    ident.internal_standard_id = ""
    ident.do_normalization = False
    ident.identification_notes = "my id note"
    ident.ms2_notes = "bad match to ref"
    ident.ms1_notes = "keep"
    ident.frag_references = []
    ident.intensity_references = []
    ident.compound = [compound]
    ident.mz_references = [mz_reference]
    ident.rt_references = [rt_reference]
    return ident


@pytest.fixture(name="atlas")
def fixture_atlas(compound_identification):
    small_atlas = metob.Atlas()
    small_atlas.compound_identifications = [compound_identification]
    return small_atlas


@pytest.fixture(name="compound_2")
def fixture_compound_2(username):
    compound = metob.Compound()
    compound.chebi_id = "CHEBI:16335"
    compound.chebi_url = "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:16335"
    compound.creation_time = 1466212384
    compound.description = "A ribonucleoside composed of a molecule of adenine attached to a ribofuranose moiety via a beta1N9-glycosidic bond."
    compound.formula = "C10H13N5O4"
    compound.head_id = "1ad02275f47b4033a451e99874f4764f"
    compound.hmdb_id = "HMDB00050"
    compound.hmdb_url = "http://www.hmdb.ca/metabolites/HMDB00050"
    compound.img_abc_id = ""
    compound.inchi = "InChI=1S/C10H13N5O4/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(18)6(17)4(1-16)19-10/h2-4,6-7,10,16-18H,1H2,(H2,11,12,13)/t4-,6-,7-,10-/m1/s1"
    compound.inchi_key = "OIRDTQYFTABQOQ-KQYNXXCUSA-N"
    compound.iupac_name = ""
    compound.kegg_id = "C00212"
    compound.kegg_url = "http://www.genome.jp/dbget-bin/www_bget?C00212"
    compound.last_modified = 1612996604
    compound.lipidmaps_id = ""
    compound.lipidmaps_url = ""
    compound.metacyc_id = "ADENOSINE"
    compound.mono_isotopic_molecular_weight = 267.096753896
    compound.name = "adenosine"
    compound.neutralized_2d_inchi = "InChI=1S/C10H13N5O4/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(18)6(17)4(1-16)19-10/h2-4,6-7,10,16-18H,1H2,(H2,11,12,13)"
    compound.neutralized_2d_inchi_key = "OIRDTQYFTABQOQ-UHFFFAOYSA-N"
    compound.neutralized_inchi = "InChI=1S/C10H13N5O4/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(18)6(17)4(1-16)19-10/h2-4,6-7,10,16-18H,1H2,(H2,11,12,13)/t4-,6-,7-,10-/m1/s1"
    compound.neutralized_inchi_key = "OIRDTQYFTABQOQ-KQYNXXCUSA-N"
    compound.num_free_radicals = 0
    compound.number_components = 1
    compound.permanent_charge = 0
    compound.prev_uid = "origin"
    compound.pubchem_compound_id = "60961"
    compound.pubchem_url = "http://pubchem.ncbi.nlm.nih.gov/compound/60961"
    compound.source = "chebi///wikidata///metacyc///gnps///hmdb"
    compound.synonyms = "adenosine///58-61-7///Adenocard///Adenoscan"  # this value was pruned down
    compound.unique_id = "1ad02275f47b4033a451e99874f4764f"
    compound.username = username
    compound.wikipedia_url = ""
    return compound


@pytest.fixture(name="rt_reference_2")
def fixture_rt_reference_2(username):
    rt_ref = metob.RtReference()
    rt_ref.creation_time = 1613002857
    rt_ref.description = "No description"
    rt_ref.enabled = True
    rt_ref.head_id = "f74622bcef924f5390ba6e127633e731"
    rt_ref.last_modified = 1613002980
    rt_ref.lcms_run = None
    rt_ref.name = "Untitled"
    rt_ref.prev_uid = "origin"
    rt_ref.ref_type = ""
    rt_ref.rt_max = 3.5233184079926665
    rt_ref.rt_min = 2.5233184079926665
    rt_ref.rt_peak = 3.0233184079926665
    rt_ref.rt_units = "min"
    rt_ref.unique_id = "f74622bcef924f5390ba6e127633e731"
    rt_ref.username = username
    return rt_ref


@pytest.fixture(name="mz_reference_2")
def fixture_mz_reference_2(username):
    mz_ref = metob.MzReference()
    mz_ref.adduct = "[M+H]+"
    mz_ref.creation_time = 1613002857
    mz_ref.description = "No description"
    mz_ref.detected_polarity = "positive"
    mz_ref.enabled = True
    mz_ref.head_id = "b0e3cf0df44a4079be7908c6b525d3ac"
    mz_ref.last_modified = 1613002980
    mz_ref.lcms_run = None
    mz_ref.modification = ""
    mz_ref.mz = 268.1040539
    mz_ref.mz_tolerance = 20.0
    mz_ref.mz_tolerance_units = "ppm"
    mz_ref.name = "Untitled"
    mz_ref.observed_formula = ""
    mz_ref.prev_uid = "origin"
    mz_ref.ref_type = ""
    mz_ref.unique_id = "b0e3cf0df44a4079be7908c6b525d3ac"
    mz_ref.username = username
    return mz_ref


@pytest.fixture(name="compound_identification_2")
def fixture_compound_identification_2(compound_2, rt_reference_2, mz_reference_2, username):
    ident = metob.CompoundIdentification()
    ident.creation_time = 1613002856
    ident.description = "No description"
    ident.do_normalization = False
    ident.frag_references = []
    ident.head_id = "6cca7aa44c0e4a109f695ba980d69472"
    ident.identification_grade = None
    ident.identification_notes = ""
    ident.intensity_references = []
    ident.internal_standard_id = ""
    ident.internal_standard_to_use = ""
    ident.last_modified = 1613002980
    ident.ms1_notes = ""
    ident.ms2_notes = ""
    ident.name = "adenosine"
    ident.prev_uid = "origin"
    ident.unique_id = "6cca7aa44c0e4a109f695ba980d69472"
    ident.username = username
    ident.frag_references = []
    ident.intensity_references = []
    ident.compound = [compound_2]
    ident.mz_references = [mz_reference_2]
    ident.rt_references = [rt_reference_2]
    return ident


@pytest.fixture(name="atlas_with_2_cids")
def fixture_atlas_with_2_cids(compound_identification, compound_identification_2):
    small_atlas = metob.Atlas()
    small_atlas.compound_identifications = [
        compound_identification,
        compound_identification_2,
    ]
    return small_atlas


@pytest.fixture(name="lcmsrun")
def fixture_lcmsrun(username):
    run = metob.LcmsRun()
    run.unique_id = "7ce51039cfca4426b4e51999ac45d018"
    run.username = username
    run.hdf5_file = "/project/projectdirs/metatlas/raw_data/akuftin/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5"  # noqa: E501
    run.description = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.mzML"  # noqa: E501
    run.creation_time = 1605311923
    run.sample = None
    run.last_modified = 1620101765
    run.mzml_file = "/project/projectdirs/metatlas/raw_data/akuftin/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.mzML"  # noqa: E501
    run.prev_uid = "28323058b6e84a9db0f9e802544764e3"
    run.method = None
    run.name = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.mzML"  # noqa: E501
    run.head_id = "7ce51039cfca4426b4e51999ac45d018"
    run.experiment = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583"
    run.injection_volume = 0.0
    run.injection_volume_units = "uL"
    run.acquisition_time = 1604770080
    run.pass_qc = False
    return run


@pytest.fixture(name="group")
def fixture_group(lcmsrun, username):
    grp = metob.Group()
    grp.items = [lcmsrun]
    grp.unique_id = "61041d07b5a24ca5b88efbda8f319654"
    grp.username = username
    grp.description = "No description"
    grp.creation_time = 1620146477
    grp.last_modified = 1620146477
    grp.prev_uid = "origin"
    grp.name = (
        f"20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_{username}0_Cone-S1"
    )
    grp.head_id = "61041d07b5a24ca5b88efbda8f319654"
    grp.short_name = "POS_Cone-S1"
    return grp


@pytest.fixture(name="group_with_2_lcmsruns")
def fixture_group_with_2_lcmsruns(lcmsrun, username):
    grp = metob.Group()
    grp.items = [lcmsrun, lcmsrun]
    grp.unique_id = "61041d07b5a24ca5b88efbda8f319654"
    grp.username = username
    grp.description = "No description"
    grp.creation_time = 1620146477
    grp.last_modified = 1620146477
    grp.prev_uid = "origin"
    grp.name = (
        f"20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_{username}0_Cone-S1"
    )
    grp.head_id = "61041d07b5a24ca5b88efbda8f319654"
    grp.short_name = "POS_Cone-S1"
    return grp


@pytest.fixture(name="hits")
def fixture_hits():
    """
    the 'data' parameter to pd.DataFrame is generated by:
    1. Running the docker testing image docker/local_jupyter.sh
    2. open /src/notebooks/reference/Targeted.ipynba
    3. Put the following in the second code block:
        source_atlas = 'HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_root0'
        metatlas_repo_path = '/src'
        project_directory = '/out'
        max_cpus = 2
        experiment = '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583'
    4. After metatlas_dataset has been created, add a code block:
        import json
        import pandas as pd
        temp_df = metatlas_dataset.hits
        temp_df['copy_index'] = temp_df.index
        slice_df = temp_df.loc[:,:,"20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5"]
        slice_df.index = pd.MultiIndex.from_tuples(
            slice_df["copy_index"], names=["database", "id", "file_name", "msms_scan"]
        )
        parsed = json.loads(slice_df.iloc[:4].to_json())
        print(json.dumps(parsed, indent=4, sort_keys=True).replace('null', 'np.nan'))
    5. copy the output from the code block into 'data' parameter of the DataFrame definition below
    """
    hits_plus = pd.DataFrame(
        data={
            "adduct": {
                "('metatlas', '29247268c3cf4acfb649ebce7b0c9e0c', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": "[M+H]+",
                "('metatlas', '50334867a31f4cab973459a59d5731c4', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": "[M+H]+",
                "('metatlas', '8ba70c0f245247eeb6ba90011026763a', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": "[M+H]+",
                "('metatlas', '9d53a44c42004e16a468e92e2b0a7009', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": "[M+H]+",
            },
            "copy_index": {
                "('metatlas', '29247268c3cf4acfb649ebce7b0c9e0c', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": [
                    "metatlas",
                    "29247268c3cf4acfb649ebce7b0c9e0c",
                    "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5",
                    2.6239302158,
                ],
                "('metatlas', '50334867a31f4cab973459a59d5731c4', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": [
                    "metatlas",
                    "50334867a31f4cab973459a59d5731c4",
                    "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5",
                    2.6239302158,
                ],
                "('metatlas', '8ba70c0f245247eeb6ba90011026763a', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": [
                    "metatlas",
                    "8ba70c0f245247eeb6ba90011026763a",
                    "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5",
                    2.6239302158,
                ],
                "('metatlas', '9d53a44c42004e16a468e92e2b0a7009', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": [
                    "metatlas",
                    "9d53a44c42004e16a468e92e2b0a7009",
                    "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5",
                    2.6239302158,
                ],
            },
            "inchi_key": {
                "('metatlas', '29247268c3cf4acfb649ebce7b0c9e0c', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": "GFFGJBXGBJISGV-UHFFFAOYSA-N",
                "('metatlas', '50334867a31f4cab973459a59d5731c4', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": "GFFGJBXGBJISGV-UHFFFAOYSA-N",
                "('metatlas', '8ba70c0f245247eeb6ba90011026763a', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": "GFFGJBXGBJISGV-UHFFFAOYSA-N",
                "('metatlas', '9d53a44c42004e16a468e92e2b0a7009', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": "GFFGJBXGBJISGV-UHFFFAOYSA-N",
            },
            "measured_precursor_intensity": {
                "('metatlas', '29247268c3cf4acfb649ebce7b0c9e0c', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 1779719.0,
                "('metatlas', '50334867a31f4cab973459a59d5731c4', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 1779719.0,
                "('metatlas', '8ba70c0f245247eeb6ba90011026763a', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 1779719.0,
                "('metatlas', '9d53a44c42004e16a468e92e2b0a7009', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 1779719.0,
            },
            "measured_precursor_mz": {
                "('metatlas', '29247268c3cf4acfb649ebce7b0c9e0c', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 136.06199646,
                "('metatlas', '50334867a31f4cab973459a59d5731c4', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 136.06199646,
                "('metatlas', '8ba70c0f245247eeb6ba90011026763a', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 136.06199646,
                "('metatlas', '9d53a44c42004e16a468e92e2b0a7009', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 136.06199646,
            },
            "msv_query_aligned": {
                "('metatlas', '29247268c3cf4acfb649ebce7b0c9e0c', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": [
                    [
                        np.nan,
                        53.2601699829,
                        59.4822044373,
                        65.2955932617,
                        66.7956771851,
                        75.0065155029,
                        75.0689544678,
                        75.4281921387,
                        84.2779464722,
                        91.0504608154,
                        94.0367355347,
                        102.1198806763,
                        108.4924850464,
                        119.0352630615,
                        121.0889511108,
                        123.1165771484,
                        135.7551269531,
                        136.0224761963,
                        136.0620117188,
                        136.1121368408,
                        136.3276824951,
                        137.046295166,
                    ],
                    [
                        np.nan,
                        2901.2893066406,
                        3058.2041015625,
                        2817.9626464844,
                        3278.6765136719,
                        3068.3347167969,
                        8541.603515625,
                        2778.4802246094,
                        2839.1333007812,
                        4060.1638183594,
                        5292.673828125,
                        3443.1560058594,
                        3947.8520507812,
                        8919.974609375,
                        5798.638671875,
                        3330.2827148438,
                        2859.4689941406,
                        18918.111328125,
                        625742.3125,
                        91467.8984375,
                        4438.6645507812,
                        11957.54296875,
                    ],
                ],
                "('metatlas', '50334867a31f4cab973459a59d5731c4', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": [
                    [
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        53.2601699829,
                        59.4822044373,
                        65.2955932617,
                        66.7956771851,
                        75.0065155029,
                        75.0689544678,
                        75.4281921387,
                        84.2779464722,
                        91.0504608154,
                        94.0367355347,
                        102.1198806763,
                        108.4924850464,
                        119.0352630615,
                        121.0889511108,
                        123.1165771484,
                        135.7551269531,
                        136.0224761963,
                        136.0620117188,
                        136.1121368408,
                        136.3276824951,
                        137.046295166,
                    ],
                    [
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        2901.2893066406,
                        3058.2041015625,
                        2817.9626464844,
                        3278.6765136719,
                        3068.3347167969,
                        8541.603515625,
                        2778.4802246094,
                        2839.1333007812,
                        4060.1638183594,
                        5292.673828125,
                        3443.1560058594,
                        3947.8520507812,
                        8919.974609375,
                        5798.638671875,
                        3330.2827148438,
                        2859.4689941406,
                        18918.111328125,
                        625742.3125,
                        91467.8984375,
                        4438.6645507812,
                        11957.54296875,
                    ],
                ],
                "('metatlas', '8ba70c0f245247eeb6ba90011026763a', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": [
                    [
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        53.2601699829,
                        59.4822044373,
                        65.2955932617,
                        66.7956771851,
                        75.0065155029,
                        75.0689544678,
                        75.4281921387,
                        84.2779464722,
                        91.0504608154,
                        94.0367355347,
                        102.1198806763,
                        108.4924850464,
                        119.0352630615,
                        121.0889511108,
                        123.1165771484,
                        135.7551269531,
                        136.0224761963,
                        136.0620117188,
                        136.1121368408,
                        136.3276824951,
                        137.046295166,
                    ],
                    [
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        2901.2893066406,
                        3058.2041015625,
                        2817.9626464844,
                        3278.6765136719,
                        3068.3347167969,
                        8541.603515625,
                        2778.4802246094,
                        2839.1333007812,
                        4060.1638183594,
                        5292.673828125,
                        3443.1560058594,
                        3947.8520507812,
                        8919.974609375,
                        5798.638671875,
                        3330.2827148438,
                        2859.4689941406,
                        18918.111328125,
                        625742.3125,
                        91467.8984375,
                        4438.6645507812,
                        11957.54296875,
                    ],
                ],
                "('metatlas', '9d53a44c42004e16a468e92e2b0a7009', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": [
                    [
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        53.2601699829,
                        59.4822044373,
                        65.2955932617,
                        66.7956771851,
                        75.0065155029,
                        75.0689544678,
                        75.4281921387,
                        84.2779464722,
                        91.0504608154,
                        94.0367355347,
                        102.1198806763,
                        108.4924850464,
                        119.0352630615,
                        121.0889511108,
                        123.1165771484,
                        135.7551269531,
                        136.0224761963,
                        136.0620117188,
                        136.1121368408,
                        136.3276824951,
                        137.046295166,
                    ],
                    [
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        2901.2893066406,
                        3058.2041015625,
                        2817.9626464844,
                        3278.6765136719,
                        3068.3347167969,
                        8541.603515625,
                        2778.4802246094,
                        2839.1333007812,
                        4060.1638183594,
                        5292.673828125,
                        3443.1560058594,
                        3947.8520507812,
                        8919.974609375,
                        5798.638671875,
                        3330.2827148438,
                        2859.4689941406,
                        18918.111328125,
                        625742.3125,
                        91467.8984375,
                        4438.6645507812,
                        11957.54296875,
                    ],
                ],
            },
            "msv_ref_aligned": {
                "('metatlas', '29247268c3cf4acfb649ebce7b0c9e0c', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": [
                    [
                        51.3947,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        91.0548,
                        94.0404,
                        np.nan,
                        np.nan,
                        119.035,
                        np.nan,
                        np.nan,
                        np.nan,
                        136.022,
                        136.062,
                        136.112,
                        np.nan,
                        137.046,
                    ],
                    [
                        1870.1,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        3051.11,
                        13543.2,
                        np.nan,
                        np.nan,
                        28284.0,
                        np.nan,
                        np.nan,
                        np.nan,
                        55585.3,
                        1607820.0,
                        17469.6,
                        np.nan,
                        43758.8,
                    ],
                ],
                "('metatlas', '50334867a31f4cab973459a59d5731c4', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": [
                    [
                        52.1001,
                        53.5537,
                        54.6096,
                        57.8238,
                        63.3067,
                        64.108,
                        82.7587,
                        93.0862,
                        94.6115,
                        111.471,
                        113.584,
                        115.21,
                        137.067,
                        137.476,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        94.0407,
                        np.nan,
                        np.nan,
                        119.036,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        136.062,
                        np.nan,
                        np.nan,
                        137.046,
                    ],
                    [
                        491091.0,
                        614205.0,
                        486992.0,
                        569335.0,
                        2513570.0,
                        554436.0,
                        577010.0,
                        580100.0,
                        930338.0,
                        567270.0,
                        515519.0,
                        616418.0,
                        17234000.0,
                        693366.0,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        2437690.0,
                        np.nan,
                        np.nan,
                        7680000.0,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        514804000.0,
                        np.nan,
                        np.nan,
                        4940020.0,
                    ],
                ],
                "('metatlas', '8ba70c0f245247eeb6ba90011026763a', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": [
                    [
                        59.3596,
                        62.4513,
                        63.2027,
                        76.4601,
                        86.8208,
                        115.912,
                        115.975,
                        123.375,
                        137.067,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        94.0407,
                        np.nan,
                        np.nan,
                        119.036,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        136.062,
                        np.nan,
                        np.nan,
                        137.046,
                    ],
                    [
                        55769.1,
                        43616.3,
                        118692.0,
                        54358.0,
                        48393.1,
                        45996.2,
                        55157.9,
                        61623.1,
                        1357390.0,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        121260.0,
                        np.nan,
                        np.nan,
                        306316.0,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        41864400.0,
                        np.nan,
                        np.nan,
                        370525.0,
                    ],
                ],
                "('metatlas', '9d53a44c42004e16a468e92e2b0a7009', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": [
                    [
                        55.0301,
                        56.3854,
                        66.7513,
                        67.0298,
                        81.1529,
                        82.4076,
                        92.0251,
                        92.3892,
                        104.302,
                        109.051,
                        112.051,
                        135.054,
                        135.653,
                        136.227,
                        136.474,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        94.0405,
                        np.nan,
                        np.nan,
                        119.035,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        136.062,
                        np.nan,
                        np.nan,
                        137.046,
                    ],
                    [
                        246689.0,
                        186484.0,
                        198526.0,
                        974057.0,
                        232546.0,
                        306008.0,
                        388476.0,
                        265393.0,
                        246201.0,
                        1625240.0,
                        1318880.0,
                        345780.0,
                        925801.0,
                        254046.0,
                        715569.0,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        7436560.0,
                        np.nan,
                        np.nan,
                        23732500.0,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        884493000.0,
                        np.nan,
                        np.nan,
                        23845700.0,
                    ],
                ],
            },
            "name": {
                "('metatlas', '29247268c3cf4acfb649ebce7b0c9e0c', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": "adenine",
                "('metatlas', '50334867a31f4cab973459a59d5731c4', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": "adenine",
                "('metatlas', '8ba70c0f245247eeb6ba90011026763a', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": "adenine",
                "('metatlas', '9d53a44c42004e16a468e92e2b0a7009', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": "adenine",
            },
            "num_matches": {
                "('metatlas', '29247268c3cf4acfb649ebce7b0c9e0c', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 7,
                "('metatlas', '50334867a31f4cab973459a59d5731c4', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 4,
                "('metatlas', '8ba70c0f245247eeb6ba90011026763a', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 4,
                "('metatlas', '9d53a44c42004e16a468e92e2b0a7009', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 4,
            },
            "precursor_mz": {
                "('metatlas', '29247268c3cf4acfb649ebce7b0c9e0c', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 136.0617952,
                "('metatlas', '50334867a31f4cab973459a59d5731c4', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 136.0617952,
                "('metatlas', '8ba70c0f245247eeb6ba90011026763a', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 136.0617952,
                "('metatlas', '9d53a44c42004e16a468e92e2b0a7009', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 136.0617952,
            },
            "score": {
                "('metatlas', '29247268c3cf4acfb649ebce7b0c9e0c', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 0.7861480398,
                "('metatlas', '50334867a31f4cab973459a59d5731c4', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 0.8248297009,
                "('metatlas', '8ba70c0f245247eeb6ba90011026763a', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 0.8078499983,
                "('metatlas', '9d53a44c42004e16a468e92e2b0a7009', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5', 2.6239302158355713)": 0.8274397807,
            },
        }
    )
    hits_plus.index = pd.MultiIndex.from_tuples(
        hits_plus["copy_index"], names=["database", "id", "file_name", "msms_scan"]
    )
    hits_plus.drop(columns=["copy_index"], inplace=True)
    return hits_plus
