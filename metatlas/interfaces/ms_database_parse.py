#!/global/common/software/m2650/python-cori/bin/python

import sys
import os
import io
import numpy as np
import multiprocessing as mp
import pandas as pd
import json
import ijson
import re
import argparse
from decimal import Decimal


class DecimalEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, Decimal):
            return float(o)
        return super(DecimalEncoder, self).default(o)


COLUMNS = ['database', 'id', 'name',
           'spectrum', 'decimal',
           'precursor_mz', 'polarity', 'adduct', 'fragmentation_method', 'collision_energy',
           'instrument', 'instrument_type',
           'formula', 'exact_mass', 'inchikey', 'inchi', 'smiles']


def tuplify(value):
    if isinstance(value, list):
        return tuple(tuplify(x) for x in value)
    else:
        return value


def msp_parse(msp_file):
    """
    Yields data delimitted by blank new lines from MSP format.
    """
    entry = []
    for line in msp_file.readlines():
        if line.strip() == '':
            yield ''.join(entry)
            entry = []
        else:
            entry.append(line)


def split_file(file_type, in_file_loc, out_dir_loc,
               files_per_folder=5000):
    """
    """

    assert(file_type in ['metatlas', 'mona', 'nist'])

    file_to_split = io.open(in_file_loc, 'r', encoding='utf-8')

    if not os.path.exists(out_dir_loc):
        os.makedirs(out_dir_loc)

    file_number = 1
    folder_number = 1

    subdir_loc = os.path.join(out_dir_loc,
                              str(file_number) +
                              '-' +
                              str(folder_number*files_per_folder))

    if not os.path.exists(subdir_loc):
        os.makedirs(subdir_loc)

    if file_type == 'metatlas' or file_type == 'mona':
        file_iter = ijson.items(file_to_split, 'item')
    if file_type == 'nist':
        file_iter = msp_parse(file_to_split)
        nist_id_parse = re.compile(r'(?m)NIST#: (\d+)')

    for item in file_iter:

        file_number += 1

        if file_number % files_per_folder == 0:
            folder_number += 1

            subdir_loc = os.path.join(out_dir_loc,
                                      str(file_number) +
                                      '-' +
                                      str(folder_number*files_per_folder))

            if not os.path.exists(subdir_loc):
                os.makedirs(subdir_loc)

        if file_type == 'metatlas':
            with open(os.path.join(subdir_loc,
                                   item['_trait_values']['unique_id']
                                   + '.json'), 'w') as file:
                json.dump(item, file, cls=DecimalEncoder)

        if file_type == 'mona':
            with open(os.path.join(subdir_loc,
                                   item['id']+'.json'), 'w') as file:
                json.dump(item, file, cls=DecimalEncoder)

        if file_type == 'nist':
            name = nist_id_parse.search(item).group(1)
            with open(os.path.join(subdir_loc,
                                   name+'.msp'), 'w') as file:
                file.write(item)

    file_to_split.close()


def convert_mona_to_csv(in_file_loc, out_file_loc, column_names=False):
    list_dict = []

    json_file = io.open(in_file_loc, 'r', encoding='utf-8')

    item = json.load(json_file)

    # identifier
    database = 'mona'
    id = item['id']
    name = item['compound'][0]['names'][0]['name']

    # spectrum
    spectrum = [[], []]
    mona_spectrum_parse = re.compile(r'([0-9\.]+):([0-9\.]+)')
    for mz, i in mona_spectrum_parse.findall(item['spectrum']):
        spectrum[0].append(Decimal(mz))
        spectrum[1].append(Decimal(i))

    # decimal
    decimal = max([-mz.as_tuple().exponent for mz in spectrum[0]])
    spectrum = str(spectrum)

    # experimental
    precursor_mz = None
    polarity = None
    adduct = None
    fragmentation_method = None
    collision_energy = None
    instrument = None
    instrument_type = None
    for subitem in item['metaData']:
        if subitem['name'] == 'precursor m/z':
            precursor_mz = str(subitem['value'])
        if subitem['name'] == 'ionization mode':
            polarity = subitem['value']
        if subitem['name'] == 'precursor type':
            adduct = subitem['value']
        if subitem['name'] == 'fragmentation mode':
            fragmentation_method = subitem['value']
        if subitem['name'] == 'fragmentation method':
            fragmentation_method = subitem['value']
        if subitem['name'] == 'collision energy':
            collision_energy = subitem['value']
        if subitem['name'] == 'instrument':
            json_instrument = subitem['value']
        if subitem['name'] == 'instrument type':
            instrument_type = subitem['value']

    # compound
    formula = None
    exact_mass = None
    inchikey = None
    inchi = None
    smiles = None
    for subitem in item['compound'][0]['metaData']:
        if subitem['name'] == 'molecular formula':
            formula = subitem['value']
        if subitem['name'] == 'total exact mass':
            exact_mass = str(subitem['value'])
        if subitem['name'] == 'InChIKey':
            inchikey = subitem['value']
        if subitem['name'] == 'InChI':
            inchi = subitem['value']
        if subitem['name'] == 'SMILES':
            smiles = subitem['value']

    list_dict.append({})
    for column in COLUMNS:
        list_dict[-1][column] = eval(column)

    df = pd.DataFrame(list_dict)
    df.set_index(['database', 'id'], inplace=True)
    df.to_csv(out_file_loc, sep='\t',
              columns=COLUMNS[2:], header=column_names,
              encoding='utf-8')

    json_file.close()


def convert_metatlas_to_csv(in_file_loc, out_file_loc, column_names=False):
    list_dict = []

    json_file = io.open(in_file_loc, 'r', encoding='utf-8')

    item = json.load(json_file)

    # identifier
    database = 'metatlas'
    id = item['_trait_values']['unique_id']
    name = item['_trait_values']['name']

    # spectrum
    spectrum = None

    try:
        spectrum = [(Decimal(str(subitem['_trait_values']['mz'])),
                    Decimal(str(subitem['_trait_values']['intensity'])))
                    for subitem in
                    item['_trait_values']['frag_references'][0]['_trait_values']['mz_intensities']]
        spectrum = list(map(list, zip(*spectrum)))
    except Exception:
        pass

    # decimal
    decimal = None

    try:
        decimal = max([-mz.as_tuple().exponent for mz in spectrum[0]])
    except Exception:
        pass
    spectrum = str(spectrum)

    # experimental
    precursor_mz = None
    polarity = None
    adduct = None
    fragmentation_method = None
    collision_energy = None
    instrument = None
    instrument_type = None

    try:
        precursor_mz = str(item['_trait_values']['frag_references'][0]['_trait_values']['precursor_mz'])
    except Exception:
        pass
    try:
        polarity = item['_trait_values']['frag_references'][0]['_trait_values']['polarity']
    except Exception:
        pass
    try:
        fragmentation_method = item['_trait_values']['frag_references'][0]['_trait_values']['technique']
    except Exception:
        pass
    try:
        collision_energy = item['_trait_values']['frag_references'][0]['_trait_values']['collision_energy']
    except Exception:
        pass

    # compound
    formula = None
    exact_mass = None
    inchikey = None
    inchi = None
    smiles = None

    try:
        formula = item['_trait_values']['compound'][0]['_trait_values']['formula']
    except Exception:
        pass
    try:
        exact_mass = str(item['_trait_values']['compound'][0]['_trait_values']['mono_isotopic_molecular_weight'])
    except Exception:
        pass
    try:
        inchikey = item['_trait_values']['compound'][0]['_trait_values']['inchi_key']
    except Exception:
        pass
    try:
        inchi = item['_trait_values']['compound'][0]['_trait_values']['inchi']
    except Exception:
        pass

    list_dict.append({})
    for column in COLUMNS:
        list_dict[-1][column] = eval(column)

    df = pd.DataFrame(list_dict)
    df.set_index(['database', 'id'], inplace=True)
    df.to_csv(out_file_loc, sep='\t',
              columns=COLUMNS[2:], header=column_names,
              encoding='utf-8')

    json_file.close()


def convert_nist_to_csv(in_file_loc, out_file_loc, column_names=False):
    list_dict = []

    msp_file = io.open(in_file_loc, 'r', encoding='utf-8')

    item = msp_file.read()

    nist_id_parse = re.compile(r'NIST#: (\d+)')
    nist_name_parse = re.compile(r'Name: (.+)')
    nist_spectrum_parse = re.compile(r'(?m)(?:(^[0-9\.]+) ([0-9\.]+)\s)')
    nist_precursor_mz_parse = re.compile(r'PrecursorMZ: ([0-9\.]+)')
    nist_polarity_parse = re.compile(r'Ion_mode: ([P|N])')
    nist_adduct_parse = re.compile(r'Precursor_type: (\[.+\]\d*[+-])')
    nist_fragmentation_method_parse = re.compile(r'Instrument_type: (.+)')
    nist_collision_energy_parse = re.compile(r'Collision_energy: (.+)')
    nist_instrument_parse = re.compile(r'Instrument: (.+)')
    nist_formula_parse = re.compile(r'Formula: ([A-Za-z0-9]+)')
    nist_exact_mass_parse = re.compile(r'ExactMass: ([0-9\.]+)')
    nist_inchikey_parse = re.compile(r'InChIKey: (.+)')

    # identifier
    database = 'nist'
    id = nist_id_parse.search(item).group(1)
    name = None
    try:
        name = nist_name_parse.search(item).group(1)
    except Exception:
        pass

    # spectrum
    spectrum = [[], []]
    for mz, i in nist_spectrum_parse.findall(item.split('Num Peaks: ', 1)[1]):
        spectrum[0].append(Decimal(mz))
        spectrum[1].append(Decimal(i))

    # decimal
    decimal = max([-mz.as_tuple().exponent for mz in spectrum[0]])
    spectrum = str(spectrum)

    # experimental
    precursor_mz = None
    polarity = None
    adduct = None
    fragmentation_method = None
    collision_energy = None
    instrument = None
    instrument_type = None

    try:
        precursor_mz = nist_precursor_mz_parse.search(item).group(1)
    except Exception:
        pass
    try:
        polarity = nist_polarity_parse.search(item).group(1)
    except Exception:
        pass
    try:
        adduct = nist_adduct_parse.search(item).group(1)
    except Exception:
        pass
    try:
        fragmentation_method = nist_fragmentation_method_parse.search(item).group(1)
    except Exception:
        pass
    try:
        collision_energy = nist_collision_energy_parse.search(item).group(1)
    except Exception:
        pass
    try:
        instrument = nist_instrument_parse.search(item).group(1)
    except Exception:
        pass

    # compound
    formula = None
    exact_mass = None
    inchikey = None
    inchi = None
    smiles = None

    try:
        formula = nist_formula_parse.search(item).group(1)
    except Exception:
        pass
    try:
        exact_mass = nist_exact_mass_parse.search(item).group(1)
    except Exception:
        pass
    try:
        inchikey = nist_inchikey_parse.search(item).group(1)
    except Exception:
        pass

    list_dict.append({})
    for column in COLUMNS:
        list_dict[-1][column] = eval(column)

    df = pd.DataFrame(list_dict)
    df.set_index(['database', 'id'], inplace=True)
    df.to_csv(out_file_loc, sep='\t',
              columns=COLUMNS[2:], header=column_names,
              encoding='utf-8')

    msp_file.close()


def generate_worklist(in_dir_loc, out_dir_loc, out_worklist_loc, database):

    valid_file_type = {'metatlas': '.json', 'mona': '.json', 'nist': '.msp'}
    worklist = []

    if not os.path.exists(out_dir_loc):
        os.makedirs(out_dir_loc)

    in_dir_loc = os.path.abspath(in_dir_loc)
    out_dir_loc = os.path.abspath(out_dir_loc)

    for root, sub_dirs, filenames in os.walk(in_dir_loc):

        for sub_dir in sub_dirs:
            out_sub_dir = os.path.join(out_dir_loc, sub_dir)
            if not os.path.exists(out_sub_dir):
                os.makedirs(out_sub_dir)

        for filename in filenames:
            name, ext = os.path.splitext(filename)
            if valid_file_type[database] == ext:
                out_file = os.path.join(out_dir_loc,
                                        os.path.relpath(root, in_dir_loc),
                                        name + '.csv')
                if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
                    worklist.append(os.path.abspath(__file__) +
                                    ' -f ' + database +
                                    ' -c \"' + os.path.join(root, filename) +
                                    '\" \"' + out_file + '\"')

    with open(out_worklist_loc, 'w') as worklist_file:
        worklist_file.write('\n'.join(worklist))


def arg_parser():
    parser = argparse.ArgumentParser()
    task_type = parser.add_mutually_exclusive_group(required=True)
    task_type.add_argument('-s', '--split', nargs='+', type=str,
                           help='[in_file_loc] [out_dir_loc]')
    task_type.add_argument('-c', '--convert', nargs='+', type=str,
                           help='[in_file_loc] [out_file_loc]')
    task_type.add_argument('-w', '--worklist', nargs='+', type=str,
                           help='[in_dir_loc] [out_dir_loc] [out_worklist_loc]')
    file_type = parser.add_mutually_exclusive_group()

    parser.add_argument('-f', '--filetype',
                        choices=['metatlas', 'mona', 'nist'],
                        required=True,
                        help='Select which filetype to work on.')

    return parser


def main():
    parser = arg_parser()
    args = parser.parse_args()

    if args.split:
        split_file(args.filetype, args.split[0], args.split[1],
                   files_per_folder=5000)

    if args.convert:
        if args.filetype == 'metatlas':
            convert_metatlas_to_csv(args.convert[0], args.convert[1],
                                    column_names=False)
        if args.filetype == 'nist':
            convert_nist_to_csv(args.convert[0], args.convert[1],
                                column_names=False)
        if args.filetype == 'mona':
            convert_mona_to_csv(args.convert[0], args.convert[1],
                                column_names=False)

    if args.worklist:
        generate_worklist(args.worklist[0], args.worklist[1],
                          args.worklist[2], args.filetype)


if __name__ == '__main__':
    main()
