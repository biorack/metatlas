import sys
import os
import multiprocessing as mp
import pprint

from metatlas.helpers import metatlas_get_data_helper_fun as ma_data
from metatlas import metatlas_objects as metob
from metatlas.helpers import dill2plots as dp
from metatlas.helpers import chromplotplus as cpp
from metatlas.helpers import spectralprocessing as sp

import numpy as np
import pandas as pd

loose_param = {'min_intensity': 1e3,
               'rt_tolerance': .25,
               'mz_tolerance': 25,
               'min_msms_score': 0.3, 'allow_no_msms': True,
               'min_num_frag_matches': 1,  'min_relative_frag_intensity': .01}

strict_param = {'min_intensity': 1e5,
                'rt_tolerance': .25,
                'mz_tolerance': 5,
                'min_msms_score': .6, 'allow_no_msms': False,
                'min_num_frag_matches': 3, 'min_relative_frag_intensity': .1}

def make_stats_table(input_fname = '', input_dataset = [],
                     include_lcmsruns = [], exclude_lcmsruns = [], include_groups = [], exclude_groups = [],
                     output_loc = None,
                     min_peak_height=0, rt_tolerance=np.inf, ppm_tolerance=np.inf,
                     min_msms_score=0, min_num_frag_matches=0,
                     allow_no_msms=False, min_relative_frag_intensity=None,
                     use_labels=False, return_all=False,
                     msms_refs_loc='/project/projectdirs/metatlas/projects/spectral_libraries/msms_refs_v2.tab',
                     dependencies = {'peak_height': [],
                                     'peak_area': ['peak_height'],
                                     'rt_peak': ['peak_height', 'rt_delta'],
                                     'rt_delta': ['peak_height'],
                                     'mz_centroid': ['peak_height', 'mz_ppm'],
                                     'mz_ppm': ['peak_height'],
                                     'msms_score': ['peak_height', 'num_frag_matches'],
                                     'num_frag_matches': ['peak_height', 'msms_score']}):

    assert output_loc is not None or return_all

    if not input_dataset:
        metatlas_dataset = ma_data.get_dill_data(os.path.expandvars(input_fname))
    else:
        metatlas_dataset = input_dataset

    if output_loc is not None and not os.path.exists(output_loc):
        os.mkdir(output_loc)

    # filter runs from the metatlas dataset
    if include_lcmsruns:
        metatlas_dataset = dp.filter_lcmsruns_in_dataset_by_include_list(metatlas_dataset,'lcmsrun',include_lcmsruns)
    if include_groups:
        metatlas_dataset = dp.filter_lcmsruns_in_dataset_by_include_list(metatlas_dataset,'group',include_groups)

    if exclude_lcmsruns:
        metatlas_dataset = dp.filter_lcmsruns_in_dataset_by_exclude_list(metatlas_dataset,'lcmsrun',exclude_lcmsruns)
    if exclude_groups:
        metatlas_dataset = dp.filter_lcmsruns_in_dataset_by_exclude_list(metatlas_dataset,'group',exclude_groups)

    file_names = ma_data.get_file_names(metatlas_dataset)
    compound_names = ma_data.get_compound_names(metatlas_dataset,use_labels=use_labels)[0]

    metrics = ['msms_score', 'num_frag_matches', 'mz_centroid', 'mz_ppm', 'rt_peak', 'rt_delta', 'peak_height', 'peak_area']

    dfs = {m:None for m in metrics}
    passing = {m:np.ones((len(compound_names), len(file_names))).astype(float) for m in metrics}

    for metric in ['peak_height', 'peak_area', 'rt_peak', 'mz_centroid']:
        dfs[metric] = dp.make_output_dataframe(input_dataset=metatlas_dataset, fieldname=metric, use_labels=use_labels)

    dfs['mz_ppm'] = dfs['peak_height'].copy()
    dfs['mz_ppm'] *= np.nan

    dfs['msms_score'] = dfs['mz_ppm'].copy()
    dfs['num_frag_matches'] = dfs['mz_ppm'].copy()
    dfs['rt_delta'] = dfs['mz_ppm'].copy()

    passing['peak_height'] = (np.nan_to_num(dfs['peak_height'].values) >= min_peak_height).astype(float)

    msms_hits_df = dp.get_msms_hits(metatlas_dataset, use_labels, ref_index=['database', 'id', 'inchi_key', 'precursor_mz'])
    msms_hits_df.reset_index(inplace=True)

    for compound_idx, compound_name in enumerate(compound_names):

        ref_rt_peak = metatlas_dataset[0][compound_idx]['identification'].rt_references[0].rt_peak
        ref_mz = metatlas_dataset[0][compound_idx]['identification'].mz_references[0].mz

        dfs['rt_delta'].iloc[compound_idx] = abs(ref_rt_peak - dfs['rt_peak'].iloc[compound_idx])
        passing['rt_delta'][compound_idx] = (abs(ref_rt_peak - np.nan_to_num(dfs['rt_peak'].iloc[compound_idx].values)) <= rt_tolerance).astype(float)

        dfs['mz_ppm'].iloc[compound_idx] = 1e6*(abs(ref_mz - dfs['mz_centroid'].iloc[compound_idx]) / ref_mz)
        passing['mz_ppm'][compound_idx] = (dfs['mz_ppm'].iloc[compound_idx].values <= ppm_tolerance).astype(float)

        inchi_key = metatlas_dataset[0][compound_idx]['identification'].compound[0].inchi_key

        for file_idx, file_name in enumerate(file_names):
            rows = msms_hits_df[(msms_hits_df['inchi_key'] == inchi_key) &
                                (msms_hits_df['file_name'] == file_name) &
                                (abs(msms_hits_df['precursor_mz'].values.astype(float) - data[0][compound_idx]['identification'].mz_references[0].mz) \
                                   <= data[0][compound_idx]['identification'].mz_references[0].mz_tolerance*1e-6)]

            if len(rows) == 0:
                dfs['msms_score'].iat[compound_idx, file_idx] = np.nan
                dfs['num_frag_matches'].iat[compound_idx, file_idx] = np.nan
            else:
                dfs['msms_score'].iat[compound_idx, file_idx] = rows.loc[rows['score'].idxmax()]['score']
                dfs['num_frag_matches'].iat[compound_idx, file_idx] = rows.loc[rows['score'].idxmax()]['num_matches']

    passing['msms_score'] = (np.nan_to_num(dfs['msms_score'].values) >= min_msms_score).astype(float)
    passing['num_frag_matches'] = (np.nan_to_num(dfs['num_frag_matches'].values) >= min_num_frag_matches).astype(float)

    for metric in metrics:
        passing[metric][passing[metric] == 0] = np.nan

    stats_table = []

    for metric in metrics:
        test = np.product(np.array([passing[dep] for dep in dependencies[metric]]), axis=0)
        # group_df = (dfs[metric] * test).T.groupby('group').describe()
        if output_loc is not None:
            (dfs[metric] * test).to_csv(os.path.join(output_loc, 'filtered_%s.tab'%metric), sep='\t')
        stats_df = (dfs[metric] * test * passing[metric]).T.describe().T
        stats_df['range'] = stats_df['max'] - stats_df['min']
        stats_df.columns = pd.MultiIndex.from_product([['filtered'], [metric], stats_df.columns])
        stats_table.append(stats_df)

    for metric in metrics:
        if output_loc is not None:
            dfs[metric].to_csv(os.path.join(output_loc, 'unfiltered_%s.tab'%metric), sep='\t')
        stats_df = dfs[metric].T.describe().T
        stats_df['range'] = stats_df['max'] - stats_df['min']
        stats_df.columns = pd.MultiIndex.from_product([['unfiltered'], [metric], stats_df.columns])
        stats_table.append(stats_df)

    stats_table = pd.concat(stats_table, axis=1)

    if output_loc is not None:
        stats_table.to_csv(os.path.join(output_loc, 'stats_table.tab'), sep='\t')

        with open(os.path.join(output_loc, 'stats_table.readme'), 'w') as readme:
            for var in ['dependencies', 'min_peak_height', 'rt_tolerance', 'ppm_tolerance', 'min_msms_score', 'min_num_frag_matches']:
                readme.write('%s\n'%var)
                try:
                    if np.isinf(eval(var)):
                        pprint.pprint('default', readme)
                    else:
                        pprint.pprint(eval(var), readme)
                except TypeError:
                    pprint.pprint(eval(var), readme)
                readme.write('\n')

    if return_all:
        return stats_table, dfs, passing


def make_scores_df(metatlas_dataset):
    """
    Returns pandas dataframe with columns 'max_intensity', 'median_rt_shift','median_mz_ppm', 'max_msms_score',
    'num_frag_matches', and 'max_relative_frag_intensity', rows of compounds in metatlas_dataset, and values
    of the best "score" for a given compound across all files.

    'max_intensity': highest intensity across all files for given compound
    'median_rt_shift': median shift of RT across all files for given compound to reference
    'median_mz_ppm': median ppm of mz across all files for given compound relative to reference
    'max_msms_score': highest compound dot-product score across all files for given compound relative to reference
    'num_frag_matches': number of matching mzs when calculating max_msms_score
    'max_relative_frag_intensity': ratio of second highest to first highest intensity of matching sample mzs

    :param metatlas_dataset:

    :return scores_df: pandas dataframe
    """

    file_names = ma_data.get_file_names(metatlas_dataset)
    compound_names = ma_data.get_compound_names(metatlas_dataset)[0]

    scores = []

    msms_hits_df = dp.get_msms_hits(metatlas_dataset, ref_index=['database', 'id', 'inchi_key', 'precursor_mz'])
    msms_hits_df.reset_index(inplace=True)

    for compound_idx in range(len(compound_names)):
        intensities = []
        rt_shifts = []
        mz_ppms = []
        max_msms_score = np.nan
        num_frag_matches = np.nan
        max_relative_frag_intensity = np.nan

        compound_ref_rt_peak = metatlas_dataset[0][compound_idx]['identification'].rt_references[0].rt_peak
        compound_ref_mz = metatlas_dataset[0][compound_idx]['identification'].mz_references[0].mz
        inchi_key = metatlas_dataset[0][compound_idx]['identification'].compound[0].inchi_key

        if len(msms_hits_df) == 0:
            comp_msms_hits = msms_hits_df
        else:
            comp_msms_hits = msms_hits_df[(msms_hits_df['inchi_key'] == data[0][compound_idx]['identification'].compound[0].inchi_key) \
                                          & (abs(msms_hits_df['precursor_mz'].values.astype(float) - data[0][compound_idx]['identification'].mz_references[0].mz) \
                                             <= data[0][compound_idx]['identification'].mz_references[0].mz_tolerance*1e-6)]

        for file_idx in range(len(file_names)):
            try:
                assert(metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['peak_height'] > 0)
                intensities.append(metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['peak_height'])
            except:# AssertionError:
                pass

            try:
                assert(metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['num_ms1_datapoints'] > 0)
                rt_shifts.append(abs(compound_ref_rt_peak - metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['rt_centroid']))
                mz_ppms.append(1e6*(abs(compound_ref_mz - metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['mz_centroid']) / compound_ref_mz))
            except:# AssertionError:
                pass

        if len(comp_msms_hits['score']) > 0:
            row = comp_msms_hits.loc[comp_msms_hits['score'].idxmax()]
            max_msms_score = row['score']
            num_frag_matches = row['num_matches']

            if num_frag_matches > 1:
                msv_sample_matches = sp.partition_aligned_ms_vectors(row['msv_query_aligned'],
                                                                     row['msv_ref_aligned'])[0]
                msv_sample_matches = msv_sample_matches[:, msv_sample_matches[1].argsort()[::-1]]
                msv_sample_matches_by_intensity = msv_sample_matches[:, msv_sample_matches[1].argsort()]

                max_relative_frag_intensity = msv_sample_matches_by_intensity[1,-2] / msv_sample_matches_by_intensity[1,-1]

        try:
            max_intensity = np.nanmax(intensities)
        except ValueError:
            max_intensity = np.nan
        try:
            median_rt_shift = np.nanmedian(rt_shifts)
        except ValueError:
            median_rt_shift = np.nan
        try:
            median_mz_ppm = np.nanmedian(mz_ppms)
        except ValueError:
            median_mz_ppm = np.nan

        # assign scores
        scores.append([metatlas_dataset[0][compound_idx]['identification'].compound[0].name,
                       metatlas_dataset[0][compound_idx]['identification'].compound[0].inchi_key,
                       max_intensity,
                       median_rt_shift,
                       median_mz_ppm,
                       max_msms_score,
                       num_frag_matches,
                       max_relative_frag_intensity])

    scores_df = pd.DataFrame(scores,
                             columns=['name',
                                      'inchi_key',
                                      'max_intensity',
                                      'median_rt_shift',
                                      'median_mz_ppm',
                                      'max_msms_score',
                                      'num_frag_matches',
                                      'max_relative_frag_intensity'])

    return scores_df


def test_scores_df(scores_df,
                   min_intensity, rt_tolerance, mz_tolerance,
                   min_msms_score, allow_no_msms, min_num_frag_matches, min_relative_frag_intensity):
    """
    Returns pandas series containing boolean values for each
    compound in scores_df describing if it passes minimum requirements set by:
    'min_intensity', 'rt_tolerance','mz_tolerance', 'min_msms_score',
    'min_num_frag_matches', and 'min_relative_frag_intensity'.

    'min_intensity' <= highest intensity across all files for given compound
    'rt_tolerance' >= shift of median RT across all files for given compound to reference
    'mz_tolerance' >= ppm of median mz across all files for given compound relative to reference
    'min_msms_score' <= highest compound dot-product score across all files for given compound relative to reference
    'min_num_frag_matches' <= number of matching mzs when calculating max_msms_score
    'min_relative_frag_intensity' <= ratio of second highest to first highest intensity of matching sample mzs

    :param scores_df:

    :return test_series:
    """

    def test_row(row):
        if (row.max_intensity >= min_intensity):
            if (row.median_rt_shift <= rt_tolerance):
                if (row.median_mz_ppm <= mz_tolerance):
                    if allow_no_msms:
                        if np.isnan(row.max_msms_score):
                            return True
                    if (row.max_msms_score >= min_msms_score):
                        if (row.num_frag_matches >= min_num_frag_matches):
                            if (row.num_frag_matches <= 1):
                                return True
                            if (row.max_relative_frag_intensity >= min_relative_frag_intensity):
                                return True
        return False

    return scores_df.apply(test_row, axis=1)



def filter_atlas_and_dataset(scores_df, atlas_df, metatlas_dataset,
                             column='passing'):
    """
    Splits atlas and metatlas_dataset by compound according to if it
    passes/fails minimum requirements set by:
    'min_intensity', 'rt_tolerance','mz_tolerance', 'min_msms_score',
    'min_num_frag_matches', and 'min_relative_frag_intensity'.

    'min_intensity' <= highest intensity across all files for given compound
    'rt_tolerance' >= shift of median RT across all files for given compound to reference
    'mz_tolerance' >= ppm of median mz across all files for given compound relative to reference
    'min_msms_score' <= highest compound dot-product score across all files for given compound relative to reference
    'min_num_frag_matches' <= number of matching mzs when calculating max_msms_score
    'min_relative_frag_intensity' <= ratio of second highest to first highest intensity of matching sample mzs

    :param scores_df:
    :param atlas:
    :param metatlas_dataset:

    :return pass_atlas_df, fail_atlas_df, pass_dataset, fail_dataset:
    """

    try:
        assert(column in scores_df)
    except AssertionError:
        print 'Error: ' + column + ' not in scores_df. Either set column where pass/fail boolean values are or run test_scores_df().'

    pass_atlas_df = atlas_df[atlas_df.inchi_key.isin(scores_df[scores_df['passing']].inchi_key.values)]
    fail_atlas_df = atlas_df[atlas_df.inchi_key.isin(scores_df[~scores_df['passing']].inchi_key.values)]

    pass_dataset = []
    fail_dataset = []

    for file_idx in range(len(metatlas_dataset)):
        pass_temp = []
        fail_temp = []

        for compound_idx in range(len(metatlas_dataset[0])):
            if metatlas_dataset[0][compound_idx]['identification'].compound[0].inchi_key in scores_df[scores_df['passing']].inchi_key.values:
                pass_temp.append(metatlas_dataset[file_idx][compound_idx])
            else:
                fail_temp.append(metatlas_dataset[file_idx][compound_idx])

        pass_dataset.append(pass_temp)
        fail_dataset.append(fail_temp)

    return pass_atlas_df, fail_atlas_df, pass_dataset, fail_dataset


def filter_and_output(atlas_df, metatlas_dataset, output_dir,
                      min_intensity,
                      rt_tolerance,
                      mz_tolerance,
                      min_msms_score, allow_no_msms,
                      min_num_frag_matches,  min_relative_frag_intensity,
                      num_threads=4,
                      output_pass=True, output_fail=False,
                      compress=False):

    """
    Splits atlas and metatlas_dataset by compound according to if it
    passes/fails minimum requirements set by:
    'min_intensity', 'rt_tolerance','mz_tolerance', 'min_msms_score',
    'min_num_frag_matches', and 'min_relative_frag_intensity' and
    creates error bars, chromatograms, and identification figures in output_dir.

    'min_intensity' <= highest intensity across all files for given compound
    'rt_tolerance' >= shift of median RT across all files for given compound to reference
    'mz_tolerance' >= ppm of median mz across all files for given compound relative to reference
    'min_msms_score' <= highest compound dot-product score across all files for given compound relative to reference
    'min_num_frag_matches' <= number of matching mzs when calculating max_msms_score
    'min_relative_frag_intensity' <= ratio of second highest to first highest intensity of matching sample mzs
    'num_threads' = number of threads to use in multiprocessing

    Returns the unfiltered metatlas dataset and filtered dataset that can be used for downstream processing steps.

    :param atlas:
    :param groups:
    :param output_dir:
    """

    with open(os.path.join(output_dir, 'test_parameters.txt'), 'w') as f:
        f.write('min_intensity=' + str(min_intensity) + '\n' +
                'rt_tolerance=' + str(rt_tolerance) + '\n' +
                'mz_tolerance=' + str(mz_tolerance) + '\n' +
                'min_msms_score=' + str(min_msms_score) + '\n' +
                'allow_no_msms=' + str(allow_no_msms) + '\n' +
                'min_num_frag_matches=' + str(min_num_frag_matches) + '\n' +
                'min_relative_frag_intensity=' + str(min_relative_frag_intensity))

    print 'making scores_df'
    # scores compounds in metatlas dataset
    scores_df = make_scores_df(metatlas_dataset)

    print 'testing and making compound_scores.csv'
    # scores dataframe
    scores_df['passing'] = test_scores_df(scores_df,
                                          min_intensity, rt_tolerance, mz_tolerance,
                                          min_msms_score, allow_no_msms, min_num_frag_matches, min_relative_frag_intensity)
    scores_df.to_csv(os.path.join(output_dir, 'compound_scores.csv'))

    print 'filtering atlas and dataset'
    # filter dataset by scores
    pass_atlas_df, fail_atlas_df, pass_dataset, fail_dataset = filter_atlas_and_dataset(scores_df, atlas_df, metatlas_dataset)

    outputs = []

    if output_pass:
        try:
            pass_dataset[0][0]['data']
            outputs.append((pass_atlas_df, pass_dataset, os.path.join(output_dir, 'pass')))
        except:
            pass

    if output_fail:
        try:
            fail_dataset[0][0]['data']
            outputs.append((fail_atlas_df, fail_dataset, os.path.join(output_dir, 'fail')))
        except:
            pass

    for atlas_df, filtered_dataset, output_dir in outputs:

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        print 'saving atlas'
        atlas_df.to_csv(os.path.join(output_dir, 'filtered_atlas_export.csv'))

        print 'making info tables'
        peak_height = dp.make_output_dataframe(input_fname = '',input_dataset = filtered_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='peak_height' , output_loc=os.path.join(output_dir,'sheets'))
        peak_area = dp.make_output_dataframe(input_fname = '',input_dataset = filtered_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='peak_area' , output_loc=os.path.join(output_dir,'sheets'))
        mz_peak = dp.make_output_dataframe(input_fname = '',input_dataset = filtered_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='mz_peak' , output_loc=os.path.join(output_dir,'sheets'))
        rt_peak = dp.make_output_dataframe(input_fname = '', input_dataset = filtered_dataset,include_lcmsruns = [],exclude_lcmsruns = [],fieldname='rt_peak' , output_loc=os.path.join(output_dir,'sheets'))
        mz_centroid = dp.make_output_dataframe(input_fname = '',input_dataset = filtered_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='mz_centroid' , output_loc=os.path.join(output_dir,'sheets'))
        rt_centroid = dp.make_output_dataframe(input_fname = '',input_dataset = filtered_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='rt_centroid' , output_loc=os.path.join(output_dir,'sheets'))

        print 'making error bars'
        #Error bars
        peak_height = dp.make_output_dataframe(input_fname='', input_dataset=filtered_dataset, include_lcmsruns=[], exclude_lcmsruns=[], fieldname='peak_height')
        dp.plot_errorbar_plots(peak_height, output_loc=os.path.join(output_dir, 'error_bar_peak_height'))

        print 'making identification figures'
        #Identification figures
        dp.make_identification_figure_v2(input_dataset=filtered_dataset, include_lcmsruns = [],exclude_lcmsruns=[], output_loc=os.path.join(output_dir, 'identification'))

        print 'making chromatograms'
        # Chromatograms
        group = 'sort'  # 'page' or 'index' or 'sort' or None
        save = True
        share_y = True

        file_names = ma_data.get_file_names(filtered_dataset)
        compound_names = ma_data.get_compound_names(filtered_dataset)[0]
        args_list = []

        chromatogram_str = 'compound_chromatograms'

        if not os.path.exists(os.path.join(output_dir, chromatogram_str)):
            os.makedirs(os.path.join(output_dir, chromatogram_str))

        for compound_idx, my_compound in enumerate(compound_names):
            my_data = list()
            for file_idx, my_file in enumerate(file_names):
                my_data.append(filtered_dataset[file_idx][compound_idx])
            kwargs = {'data': my_data,
                      'file_name': os.path.join(output_dir, chromatogram_str, my_compound+'.pdf'),
                      'group': group,
                      'save': save,
                      'share_y': share_y,
                      'names': file_names}
            args_list.append(kwargs)

        pool = mp.Pool(processes=min(num_threads, len(filtered_dataset[0])))
        pool.map(cpp.chromplotplus, args_list)
        pool.close()
        pool.terminate()

    print 'done'
    return pass_atlas_df, fail_atlas_df, pass_dataset, fail_dataset
