"""unit tests of feature tools functions"""

from metatlas.io import feature_tools as feature_tools
import numpy as np
import pandas as pd


def test_group_consecutive(metatlas_dataset_with_2_cids):

    data = metatlas_dataset_with_2_cids.atlas_df['mz'].values[:]
    stepsize=10.0
    do_ppm=True

    group_indices = feature_tools.group_consecutive(data, stepsize, do_ppm)

    #assert isinstance(group_indices, np.ndarray)  ## Assert by type
    assert (group_indices == np.array([0, 1], dtype=np.int64)).all() ## Assert by value

def test_setup_file_slicing_parameters001(mocker, lcmsrun, metatlas_dataset_with_2_cids):  ## Test positive polarity return

    mocker.patch("metatlas.io.feature_tools.group_consecutive", return_value=np.array([0, 1], dtype=np.int64))
    filenames = [lcmsrun.hdf5_file]
    atlas = metatlas_dataset_with_2_cids.atlas_df

    slicing_dicts = feature_tools.setup_file_slicing_parameters(atlas, filenames, extra_time=0.1,
                                                                ppm_tolerance=20,
                                                                polarity='positive',
                                                                project_dir=False,
                                                                base_dir = '/project/projectdirs/metatlas/projects/',
                                                                overwrite=True)

    assert all(isinstance(ele, dict) for ele in slicing_dicts)  ## Assert by type
    #assert slicing_dicts[0].get('lcmsrun') == filenames[0]  ## Assert by a value that's not explicitly passed


def test_setup_file_slicing_parameters002(mocker, lcmsrun, metatlas_dataset_with_2_cids):  ## Test negative polarity return

    mocker.patch("metatlas.io.feature_tools.group_consecutive", return_value=np.array([0, 1], dtype=np.int64))
    filenames = [lcmsrun.hdf5_file]
    atlas = metatlas_dataset_with_2_cids.atlas_df

    slicing_dicts = feature_tools.setup_file_slicing_parameters(atlas, filenames, extra_time=0.1,
                                                                ppm_tolerance=20,
                                                                polarity='negative',
                                                                project_dir=False,
                                                                base_dir = '/project/projectdirs/metatlas/projects/',
                                                                overwrite=True)
    
    assert all(isinstance(ele, dict) for ele in slicing_dicts)  ## Assert by type
    #assert slicing_dicts[0].get('lcmsrun') == filenames[0]  ## Assert by a value that's not explicitly passed

def test_calculate_ms1_summary001(df_container, feature_filter=True):  ## Test when feature_filter is true (removes Compound2)

    desired_key = 'ms1_pos'    

    unfiltered_data = df_container[desired_key]
    unfiltered_data['label'] = np.concatenate([np.array(np.repeat("Compound1", 37)), np.array(np.repeat("Compound2", 37))])  ## Fake two compounds
    unfiltered_data['in_feature'] = np.concatenate([np.array(np.repeat(True, 37)), np.array(np.repeat(False, 37))])  ## Fake rt_window (aka extra_time)

    summary_df = feature_tools.calculate_ms1_summary(unfiltered_data, feature_filter = feature_filter)

    assert summary_df.shape == (1,6)  ## MS1 data is split by compound only, resulting in two groups (rows in summary)


def test_calculate_ms1_summary002(df_container, feature_filter=False):  ## Test when feature_filter is false (nothing removed)

    desired_key = 'ms1_pos'    

    unfiltered_data = df_container[desired_key]
    unfiltered_data['label'] = np.concatenate([np.array(np.repeat("Compound1", 37)), np.array(np.repeat("Compound2", 37))])  ## Fake two compounds
    unfiltered_data['in_feature'] = np.concatenate([np.array(np.repeat(True, 37)), np.array(np.repeat(False, 37))])  ## Fake rt_window (aka extra_time)

    summary_df = feature_tools.calculate_ms1_summary(unfiltered_data, feature_filter = feature_filter)

    assert summary_df.shape == (2,6)  ## MS1 data is split by compound only, resulting in two groups (rows in summary)

def test_calculate_ms2_summary001(df_container, feature_filter=True):  ## Test when feature_filter is true (removes Compound2)

    desired_key = 'ms2_pos'    

    unfiltered_data = df_container[desired_key]
    unfiltered_data['label'] = np.concatenate([np.array(np.repeat("Compound1", 4)), np.array(np.repeat("Compound2", 4))])  ## Fake two compounds
    unfiltered_data['in_feature'] = np.concatenate([np.array(np.repeat(True, 4)), np.array(np.repeat(False, 4))])  ## Fake rt_window (aka extra_time)

    summary_df = feature_tools.calculate_ms2_summary(unfiltered_data, feature_filter = feature_filter)

    assert summary_df.shape == (2,5)  ## MS2 df is split by compound and by rt, resulting in two groups of two (rows in summary)


def test_calculate_ms2_summary002(df_container, feature_filter=False):  ## Test when feature_filter is false (nothing removed)

    desired_key = 'ms2_pos'    

    unfiltered_data = df_container[desired_key]
    unfiltered_data['label'] = np.concatenate([np.array(np.repeat("Compound1", 4)), np.array(np.repeat("Compound2", 4))])  ## Fake two compounds
    unfiltered_data['in_feature'] = np.concatenate([np.array(np.repeat(True, 4)), np.array(np.repeat(False, 4))])  ## Fake rt_window (aka extra_time)

    summary_df = feature_tools.calculate_ms2_summary(unfiltered_data, feature_filter = feature_filter)

    assert summary_df.shape == (4,5)  ## MS2 df is split by compound and by rt, resulting in two groups of two (rows in summary)


def test_map_mzgroups_to_data(metatlas_dataset_with_2_cids, eic):
    
    mz_atlas = metatlas_dataset_with_2_cids.atlas_df['mz'].values[:]
    mz_group_indices = np.array([0, 1], dtype=np.int64)
    mz_data = np.array(eic.get('mz'), dtype=float)

    map_mzgroups_to_data_results = feature_tools.map_mzgroups_to_data(mz_atlas, mz_group_indices, mz_data)

    assert map_mzgroups_to_data_results.shape == mz_data.shape ## Assert by shape


def test_df_container_from_metatlas_file001(mocker, df_container, lcmsrun): ## Test MS1 return

    desired_key = 'ms1_pos'
    mocker.patch("pandas.read_hdf", return_value = df_container[desired_key])
    hdf5_filename = [lcmsrun.hdf5_file][0]

    df_container_out = feature_tools.df_container_from_metatlas_file(hdf5_filename, desired_key)

    assert {'mz','i','rt'}.issubset(df_container_out.columns) == True


def test_df_container_from_metatlas_file002(mocker, df_container, lcmsrun): ## Test MS1 return

    desired_key = 'ms2_pos'
    mocker.patch("pandas.read_hdf", return_value = df_container[desired_key])
    hdf5_filename = [lcmsrun.hdf5_file][0]

    df_container_out = feature_tools.df_container_from_metatlas_file(hdf5_filename, desired_key)

    assert {'mz','i','rt','precursor_MZ','precursor_intensity','collision_energy'}.issubset(df_container_out.columns) == True


def test_filter_raw_data_using_atlas(df_container, metatlas_dataset_with_2_cids):

    desired_key = 'ms1_pos'
    msdata = df_container[desired_key]
    atlas = metatlas_dataset_with_2_cids.atlas_df
    
    atlas['group_index'] = np.array([0, 1]) ## This mocks group_consecutive() because it isn't actually called in the function
    atlas['extra_time'] = 0.1 ## This mocks setup_file_slicing_parameters() because it isn't actually called in the function
    atlas['ppm_tolerance'] = 20 ## This mocks setup_file_slicing_parameters() because it isn't actually called in the function
    msdata['group_index'] = np.array(np.repeat(0, 74))  ## This mocks map_mzgroups_to_data() because it isn't actually called in the function

    filtered_msdata = feature_tools.filter_raw_data_using_atlas(atlas, msdata)

    assert {'mz'}.issubset(filtered_msdata.columns) == True


def test_get_atlas_data_from_file001(mocker, lcmsrun, df_container, metatlas_dataset_with_2_cids):  ## Test MS1 return

    desired_key = 'ms1_pos'    

    mocker.patch("metatlas.io.feature_tools.df_container_from_metatlas_file", return_value=df_container[desired_key])
    mocker.patch("metatlas.io.feature_tools.map_mzgroups_to_data", return_value=np.array(np.repeat(0, 74)))
    
    atlas = metatlas_dataset_with_2_cids.atlas_df
    atlas['group_index'] = np.array([0, 1]) ## This mocks group_consecutive() because it isn't actually called in the function
    atlas['extra_time'] = 0.1 ## This mocks setup_file_slicing_parameters() because it isn't actually called in the function
    atlas['ppm_tolerance'] = 20 ## This mocks setup_file_slicing_parameters() because it isn't actually called in the function
    filename = [lcmsrun.hdf5_file]

    filtered_msdata = feature_tools.get_atlas_data_from_file(filename,atlas,desired_key=desired_key)

    assert {'label','rt','mz','i','in_feature'}.issubset(filtered_msdata.columns) == True


def test_get_atlas_data_from_file002(mocker, lcmsrun, df_container, metatlas_dataset_with_2_cids):  ## Test MS2 return

    desired_key = 'ms2_pos'

    mocker.patch("metatlas.io.feature_tools.df_container_from_metatlas_file", return_value=df_container[desired_key])
    mocker.patch("metatlas.io.feature_tools.map_mzgroups_to_data", return_value=np.array(np.repeat(0, 4)))  ## Four identical rts from the 8 frags
    
    atlas = metatlas_dataset_with_2_cids.atlas_df
    atlas['group_index'] = np.array([0, 1]) ## This mocks group_consecutive() because it isn't actually called in the function
    atlas['extra_time'] = 0.1 ## This mocks setup_file_slicing_parameters() because it isn't actually called in the function
    atlas['ppm_tolerance'] = 20 ## This mocks setup_file_slicing_parameters() because it isn't actually called in the function   
    filename = [lcmsrun.hdf5_file]

    filtered_msdata = feature_tools.get_atlas_data_from_file(filename,atlas,desired_key=desired_key)

    assert {'rt','i','mz','precursor_MZ','precursor_intensity','collision_energy'}.issubset(filtered_msdata.columns) == True


def test_get_data001(mocker, metatlas_dataset, df_container):  ## Test negative polarity return

    fake_input_dict = {'atlas': metatlas_dataset.atlas_df, 
                       'lcmsrun': "/path/to/lcmsrun",
                       'polarity': 'negative'}
        
    fake_summary_df = pd.DataFrame.from_dict({})

    mocker.patch("metatlas.io.feature_tools.get_atlas_data_from_file", side_effects=[df_container['ms1_%s'%fake_input_dict['polarity'][:3]], 
                                                                                     df_container['ms2_%s'%fake_input_dict['polarity'][:3]]])
    
    mocker.patch("metatlas.io.feature_tools.calculate_ms1_summary", return_value=fake_summary_df)

    out_data = feature_tools.get_data(fake_input_dict, return_data=True,save_file=False)

    assert all(keys in out_data for keys in ('atlas','ms1_data','ms1_summary','ms2_data')) == True


def test_get_data002(mocker, metatlas_dataset, df_container):  ## Test positive polarity return

    fake_input_dict = {'atlas': metatlas_dataset.atlas_df, 
                       'lcmsrun': "/path/to/lcmsrun",
                       'polarity': 'positive'}
        
    fake_summary_df = pd.DataFrame.from_dict({})

    mocker.patch("metatlas.io.feature_tools.get_atlas_data_from_file", side_effects=[df_container['ms1_%s'%fake_input_dict['polarity'][:3]], 
                                                                                     df_container['ms2_%s'%fake_input_dict['polarity'][:3]]])
    
    mocker.patch("metatlas.io.feature_tools.calculate_ms1_summary", return_value=fake_summary_df)

    out_data = feature_tools.get_data(fake_input_dict, return_data=True,save_file=False)

    assert all(keys in out_data for keys in ('atlas','ms1_data','ms1_summary','ms2_data')) == True