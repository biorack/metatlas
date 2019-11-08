#!/usr/bin/env python
# python2

import sys
import os
import argparse
import configparser

import mzml_loader as ml
import spectralprocessing as sp

import pandas as pd
import numpy as np
from scipy.stats import rankdata

def make_nistified(msv):
    msv[0,~np.isnan(msv[0])] = -1
    msv[1] = rankdata(msv[1], method='dense')*(msv[1]-msv[1]+1)

def find_spectral_hits(mzml_loc, tab_loc=None, force=False, nistify=False, **kwargs):

    if os.path.isfile(mzml_loc):
        mzml_loc = os.path.realpath(mzml_loc)
    else:
        sys.stderr.write('error: ' + mzml_loc + ' not found\n')
        sys.exit(1)

    if tab_loc is None or tab_loc == '':
        tab_loc = os.path.splitext(mzml_loc)[0] + '_spectral-hits.tab.gz'

    if not force and os.path.exists(tab_loc) and os.path.getsize(tab_loc) != 0:
        sys.stderr.write('error: ' + tab_loc + ' already exists\n')
        sys.exit(1)

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
        mzml_df = mzml_df[mzml_df['mz'] < mzml_df['precursor_mz'] + 2.5]

    except ValueError:
        sys.stderr.write('error: ' + mzml_loc + ' not convertable to dataframe\n')
        sys.exit(1)

    mzml_rt_group = mzml_df.groupby('rt')

    data_df = mzml_df[['rt', 'polarity', 'precursor_mz', 'precursor_intensity', 'collision_energy']].drop_duplicates()

    data_df['spectrum'] = map(lambda s: np.array(list(s)), zip(mzml_rt_group['mz'].apply(list),
                                                               mzml_rt_group['i'].apply(list)))


    data_df.set_index(['rt'], inplace=True)

    spectral_hits_dfs = []

    for rt, row in data_df.iterrows():

        try:
            hit_df = sp.search_ms_refs(row['spectrum'],
                                       ref_df=ref_df,
                                       **dict(row.to_dict(),
                                              **kwargs))
        except:
            sys.stderr.write('error on: \n' + repr(row))
            sys.exit(1)

        if len(hit_df.index) > 0:
            hit_df['rt'] = rt
            hit_df['polarity'] = row['polarity']
            hit_df['precursor_mz'] = row['precursor_mz']
            hit_df['precursor_intensity'] = row['precursor_intensity']

            spectral_hits_dfs.append(hit_df)

    if len(spectral_hits_dfs) == 0:
        with open(tab_loc, 'w') as blank:
            blank.write('\t'.join(['database', 'id', 'rt',
                                   'precursor_mz', 'precursor_intensity', 'polarity',
                                   'score', 'num_matches',
                                   'msv_query_aligned', 'msv_ref_aligned']))
        return

    spectral_hits_df = pd.concat(spectral_hits_dfs
                                 ).query(post_query,
                                         local_dict=dict(locals(), **kwargs))

    spectral_hits_df.set_index('rt', append=True, inplace=True)
    spectral_hits_df.reorder_levels(['rt'] + ref_index)

    spectral_hits_df['msv_query_aligned'] = spectral_hits_df['msv_query_aligned'].apply(lambda a: a.tolist())
    if nistify:
        spectral_hits_df.xs('nist', level='database')['msv_ref_aligned'].apply(make_nistified)
    spectral_hits_df['msv_ref_aligned'] = spectral_hits_df['msv_ref_aligned'].apply(lambda a: a.tolist())

    spectral_hits_df.to_csv(tab_loc, columns=['precursor_mz', 'precursor_intensity', 'polarity',
                                              'score', 'num_matches',
                                              'msv_query_aligned', 'msv_ref_aligned'],
                            sep='\t', compression='gzip')

    sys.stdout.write('created: ' + tab_loc + '\n')



def generate_worklist(worklist_loc, mzml_dir, tab_dir=None, force=False, nistify=False):

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
                out_file = os.path.join(name + '_spectral-hits.tab.gz')

                if force or not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
                    command = os.path.realpath(__file__)
                    if force:
                        command += ' -f'
                    if nistify:
                        command += ' -n'
                    command += ' -m \"' + os.path.join(root, filename) + '\" -l \"' + out_file + '\"'
                    worklist.append(command)

    with open(worklist_loc, 'w') as worklist_file:
        worklist_file.write('\n'.join(worklist))


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', required=True,
                         help='the config file')
    parser.add_argument('-f', '--force', action='store_true', required=False,
                         help='forces file(s) to be remade if they exist')
    parser.add_argument('-n', '--nistify', action='store_true', required=False,
                         help='obfuscates nist reference msms')
    parser.add_argument('-l', '--location', type=str, required=False, metavar='out_dir/loc',
                         help='changes output location for files')

    task_type = parser.add_mutually_exclusive_group(required=True)
    task_type.add_argument('-m', '--make', type=str, metavar='mzml_loc',
                           help='makes spectral_hits file')
    task_type.add_argument('-w', '--worklist', nargs=2, type=str, metavar=('worklist_loc', 'mzml_dir'),
                           help='creates worklist to make many spectral_hits files')


    return parser

def main():
    parser = arg_parser()
    args, cmd_kwargs = parser.parse_known_args()

    config = configparser.SafeConfigParser()
    config.read(os.path.join(os.path.dirname(os.path.realpath(__file__)), args.config))
    kwargs = dict(config.items('DEFAULT'))

    kwargs = dict([(k, eval(v))
                   for k, v in kwargs.items()])

    print (kwargs)

    if args.make:
        find_spectral_hits(args.make, args.location, args.force, args.nistify, **kwargs)

    if args.worklist:
        generate_worklist(args.worklist[0], args.worklist[1], args.location, args.force, args.nistify)


if __name__ == '__main__':
    main()

