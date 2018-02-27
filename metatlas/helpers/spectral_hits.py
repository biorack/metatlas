#!/global/common/software/m2650/python-cori/bin/python

import sys
import os
import argparse
import ConfigParser

import mzml_loader as ml
import spectralprocessing as sp

import pandas as pd
import numpy as np


def find_spectral_hits(mzml_loc, tab_loc=None, **kwargs):

    if os.path.isfile(mzml_loc):
        mzml_loc = os.path.realpath(mzml_loc)
    else:
        sys.stderr.write(mzml_loc + ' not a file')
        sys.exit(1)

    if tab_loc is None or tab_loc == '':
        tab_loc = os.path.splitext(mzml_loc)[0] + '_spectral-hits.tab.gz'

    ms_types = kwargs.pop('ms_types', ['ms2_pos', 'ms2_neg'])

    ref_loc = kwargs.pop('ref_loc', '/global/project/projectdirs/metatlas/projects/spectral_libraries/msms_refs.tab')
    ref_dtypes = kwargs.pop('ref_dtypes', {'database':str, 'id':str, 'name':str,
                                           'spectrum':object,'decimal':float, 'precursor_mz':float,
                                           'polarity':str, 'adduct':str, 'fragmentation_method':str,
                                           'collision_energy':str, 'instrument':str, 'instrument_type':str,
                                           'formula':str, 'exact_mass':float,
                                           'inchi_key':str, 'inchi':str, 'smiles':str})
    ref_index = kwargs.pop('ref_index', ['database', 'id'])

    pre_query = kwargs.pop('pre_query', 'index == index or index == @pd.NaT')
    post_query = kwargs.pop('post_query', 'index == index or index == @pd.NaT')

    if 'ref_df' in kwargs:
        ref_df = kwargs.pop('ref_df')
    else:
        ref_df = pd.read_csv(ref_loc,
                             sep='\t',
                             dtype=ref_dtypes
                            ).set_index(ref_index)

    ref_df = ref_df.query(pre_query,
                          local_dict=dict(locals(), **kwargs))

    if ref_df['spectrum'].apply(type).eq(str).all():
        ref_df['spectrum'] = ref_df['spectrum'].apply(lambda s: eval(s))

    mzml_dfs = ml.mzml_to_df(mzml_loc)
    try:
        mzml_df = pd.concat([df
                             for ms_type, df in mzml_dfs.items()
                             if ms_type in ms_types
                             and isinstance(df, pd.DataFrame)], axis=1)
    except ValueError:
        return

    mzml_rt_group = mzml_df.groupby('rt')

    data_df = mzml_df[['rt', 'polarity', 'precursor_mz', 'precursor_intensity', 'collision_energy']].drop_duplicates()

    data_df['spectrum'] = map(lambda s: np.array(list(s)), zip(mzml_rt_group['mz'].apply(list),
                                                               mzml_rt_group['i'].apply(list)))
    data_df.set_index(['rt'], inplace=True)

    spectral_hits_dfs = []

    for rt, scan in data_df.iterrows():

        hit_df = sp.search_ms_refs(scan['spectrum'],
                                   ref_df=ref_df,
                                   **dict(scan.to_dict(),
                                          **kwargs))

        if len(hit_df.index) > 0:
            hit_df['rt'] = rt
            hit_df['polarity'] = scan['polarity']
            hit_df['precursor_mz'] = scan['precursor_mz']
            hit_df['precursor_intensity'] = scan['precursor_intensity']

            spectral_hits_dfs.append(hit_df)

    if len(spectral_hits_dfs) == 0:
        with open(tab_loc, 'w') as blank:
            blank.write('\t'.join(['database', 'id', 'rt',
                                   'precursor_mz', 'precursor_intensity',
                                   'score', 'num_matches',
                                   'msv_query_aligned', 'msv_ref_aligned']))
        return

    spectral_hits_df = pd.concat(spectral_hits_dfs
                                 ).query(post_query,
                                         local_dict=dict(locals(), **kwargs))

    spectral_hits_df.set_index('rt', append=True, inplace=True)
    spectral_hits_df.reorder_levels(['rt'] + ref_index)

    spectral_hits_df['msv_query_aligned'] = spectral_hits_df['msv_query_aligned'].apply(lambda a: a.tolist())
    spectral_hits_df['msv_ref_aligned'] = spectral_hits_df['msv_ref_aligned'].apply(lambda a: a.tolist())

    spectral_hits_df.to_csv(tab_loc, columns=['precursor_mz', 'precursor_intensity',
                                              'score', 'num_matches',
                                              'msv_query_aligned', 'msv_ref_aligned'],
                            sep='\t', compression='gzip')


def generate_worklist(worklist_loc, mzml_dir, tab_dir=None):

    worklist = []

    assert os.path.isdir(mzml_dir)
    mzml_dir = os.path.realpath(mzml_dir)

    if tab_dir is None or tab_dir == '':
        tab_dir = mzml_dir
    else:
        tab_dir = os.path.realpath(tab_dir)

    if not os.path.exists(tab_dir):
        os.makedirs(tab_dir)

    for root, sub_dirs, filenames in os.walk(mzml_dir):
        for filename in filenames:
            name, ext = os.path.splitext(filename)
            if ext.lower() == '.mzml':
                out_path = os.path.join(tab_dir,
                                        os.path.relpath(root, mzml_dir)
                                        if root != mzml_dir
                                        else '')
                if not os.path.exists(out_path):
                    os.makedirs(out_path)

                out_file = os.path.join(out_path,
                                        name + '_spectral-hits.tab.gz')

                if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
                    worklist.append(os.path.realpath(__file__) +
                                    ' -t \"' + os.path.join(root, filename) +
                                    '\" \"' + out_file + '\"')

    with open(worklist_loc, 'w') as worklist_file:
        worklist_file.write('\n'.join(worklist))


def arg_parser():
    parser = argparse.ArgumentParser()
    task_type = parser.add_mutually_exclusive_group(required=True)
    task_type.add_argument('-t', '--tab', nargs='+', type=str,
                           help='mzml_loc [tab_loc]')
    task_type.add_argument('-w', '--worklist', nargs='+', type=str,
                           help='worklist_loc mzml_dir [tab_dir]')

    return parser

def main():
    config = ConfigParser.SafeConfigParser()
    config.read(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             'spectral_hits.config'))
    kwargs = dict(config.items('DEFAULT'))

    kwargs = dict([(k, eval(v))
                   for k, v in kwargs.items()])

    parser = arg_parser()
    args, cmd_kwargs = parser.parse_known_args()

    if args.tab:
        find_spectral_hits(*args.tab, **kwargs)

    if args.worklist:
        generate_worklist(*args.worklist)


if __name__ == '__main__':
    main()
