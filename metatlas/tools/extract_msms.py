import logging
import uuid
import sys
import os

import pandas as pd
import numpy as np

import numpy.typing as npt
from typing import Optional

import pymzml

from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from IPython.display import display
from tqdm.auto import tqdm
tqdm.pandas()

from metatlas.plots import dill2plots as dp

EASYIC_MZS = {'positive': 202.07770, 'negative': 202.07880}
EASYIC_MZ_TOLERANCE = 0.001

MS2Scan = npt.NDArray[npt.NDArray[float]]

logger = logging.getLogger(__name__)


def get_file_polarity(file_name: str) -> str:
    """
    Get polarity of file from file name conforming to Northen lab conventions.
    """
    pol = file_name.split('_')[9]
    
    if pol == 'NEG':
        return 'negative'
    elif pol == 'POS':
        return 'positive'
    else:
        return 'fps'


def get_scan_data(scan_num: int, file: str, exp_dir: str, raw_data_dir: str = '/global/cfs/cdirs/metatlas/raw_data/jgi') -> MS2Scan:
    """
    Extracts MS2 scan data for a given file and scan number.
    """
    
    file_path = os.path.join(raw_data_dir, exp_dir, "{}.mzML".format(file))
    
    reader = pymzml.run.Reader(file_path, build_index_from_scratch=True)
    mz = None
    i = None
    
    for spec in reader:
        if spec.ID == scan_num:
            if spec.ms_level <= 1:
                raise Exception("MS scan level must be greater than 1 for file: {} - scan num: {}".format(file, scan_num))
            mz = spec.mz
            i = spec.i
            
    reader.close()
    scan_data = np.array([mz, i])
    
    return scan_data


def sanitize_ms2_scan(ms2_scan: MS2Scan , precursor_mz: float, polarity: str) -> MS2Scan:
    """
    Removes fluoranthene signal from MS2 scan.
    """
    assert polarity == 'positive' or polarity == 'negative'
    
    sanitized_mzs = []
    sanitized_is = []
    
    for i, mz in enumerate(ms2_scan[0]):
    
        if abs(mz - EASYIC_MZS[polarity]) < EASYIC_MZ_TOLERANCE:
            continue
        if (mz - 2.0) > precursor_mz:
            continue

        sanitized_mzs.append(mz)
        sanitized_is.append(ms2_scan[1][i])
    
    sanitzed_ms2_scan = np.array([sanitized_mzs, sanitized_is])
    
    return sanitzed_ms2_scan


def make_text_spectrum(ms2_scan: MS2Scan) -> list[list[str]]:
    """
    Convert numpy array MS2 scan into text for storage in CSV format.
    """
    mzs = ms2_scan[0]
    intensities = ms2_scan[1]
    
    idx = np.argsort(mzs).flatten()
    mzs = [mzs[i] for i in idx]
    intensities = [intensities[i] for i in idx]
    
    spectrum = [['%.5f'%m for m in mzs], ['%.5f'%x for x in intensities]]
    spectrum = str(spectrum)
    spectrum = spectrum.replace('\'','')
    
    return spectrum


def is_ordered(ms2_scan: MS2Scan) -> bool:
    """
    Checks that MS2 scan has mz values in acending order.
    
    If not, these scans will break MatchMS during downstream operations.
    """
    
    for i in range(len(ms2_scan[0])):
        if i == 0:
            continue
        if ms2_scan[0][i] < ms2_scan[0][i-1]:
            return False
        
    return True


def enrich_metadata(msms_refs: pd.DataFrame, frag_method: str, instrument_type: str, 
                    decimal: float, id_prefix: Optional[str]) -> None:
    """
    Add metadata to MSMS refs file.
    """
    msms_refs['database'] = 'metatlas'
    msms_refs['fragmentation_method'] = frag_method
    msms_refs['instrument_type'] = instrument_type
    msms_refs['decimal'] = decimal
    msms_refs['name'] = msms_refs['label']
    
    msms_refs['instrument'] = msms_refs['file'].str.split('_', expand=True)[6]
    msms_refs['id'] = [id_prefix + str(uuid.uuid4().hex) for i in range(msms_refs.shape[0])]
    msms_refs['polarity'] = msms_refs['file'].apply(get_file_polarity)
    msms_refs['collision_energy'] = msms_refs.apply(lambda x: "{}-{}".format(x.ce_type.lower(), str(x.ce)), axis=1)
    msms_refs['precursor_mz'] = np.round(msms_refs['mz'].tolist(), 4)
    msms_refs['smiles'] = msms_refs['inchi'].apply(lambda x: Chem.MolToSmiles(Chem.MolFromInchi(x)))
    msms_refs['exact_mass'] = msms_refs['inchi'].apply(lambda x: ExactMolWt(Chem.MolFromInchi(x)))
    msms_refs['formula'] = msms_refs['inchi'].apply(lambda x: CalcMolFormula(Chem.MolFromInchi(x)))
    
    
def generate_msms_refs(msms_data_file: str, msms_refs_file: str, id_prefix: Optional[str], 
                       raw_data_dir: str, frag_method: str, instrument_type: str, decimal: float) -> None:
    """
    Generate tab delimited MSMS refs file.
    """
    display(dp.LOGGING_WIDGET)
    
    if id_prefix is None:
        id_prefix = ''
    
    msms_data = pd.read_csv(msms_data_file)
    msms_data = msms_data.dropna(subset='scan_num')
    msms_data = msms_data.astype({'scan_num': int})
    
    logger.info('Collecting scan data.')
    msms_data['spectrum'] = msms_data.progress_apply(lambda x: get_scan_data(x.scan_num, x.file, x.exp_dir, raw_data_dir=raw_data_dir), axis=1)
    enrich_metadata(msms_data, frag_method, instrument_type, decimal, id_prefix)
    
    logger.info('Sanitizing raw scan data.')
    msms_data['spectrum'] = msms_data.progress_apply(lambda x: sanitize_ms2_scan(x.spectrum, x.precursor_mz, x.polarity), axis=1)
    
    test_order = msms_data['spectrum'].apply(is_ordered)
    if not test_order.all():
        raise Exception("1 or more spectra are disordered")
        
    logger.info('Converting spectra to text.')
    msms_data['spectrum'] = msms_data['spectrum'].progress_apply(make_text_spectrum)
    msms_data.drop(columns=['label', 'ce_type', 'ce', 'mz', 'scan_num', 'file', 'exp_dir'], inplace=True)
    
    msms_data.to_csv(msms_refs_file, sep='\t')
    logger.info('DONE - execution of notebook%s is complete.')
