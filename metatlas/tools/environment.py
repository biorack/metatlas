"""Environment setup functions"""

import getpass
import json
import logging
import os
import re
import shutil
import sys

from pathlib import Path

logger = logging.getLogger(__name__)


def install_kernel():
    """
    Copies kernel.json from repo to active location under home directory.
    Only for use on NERC!
    """
    logger.info('Installing kernel.json for "Metatlas Targeted".')
    repo_path = Path(__file__).resolve().parent.parent.parent
    source = repo_path / "notebooks" / "kernels" / "metatlas-targeted.kernel.json"
    dest_dir = Path.home() / ".local" / "share" / "jupyter" / "kernels" / "metatlas-targeted"
    os.makedirs(dest_dir, exist_ok=True)
    shutil.copyfile(source, dest_dir / "kernel.json")
    logger.info(
            'Kernel installation complete. Reload Jupyter notebook page to see new kernel". '
    )


def repo_dir():
    """Returns a string with the path to the root of the Metatlas git repo"""
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def create_all_notebooks(output_type, base_output_dir, experiment_id, analysis_number):
    """
    Creates Jupyter notebooks with appropriate filename and pre-populated parameters
    inputs:
        output_type: one of 'RT-Predict', 'ISTDsEtc', 'FinalEMA-HILIC'
        base_output_dir: project directory containing the experiment directories
        experiment_id: '_' delimited experiment identifier
        analysis_number: increment to not overwrite existing analysis outputs
    """
    possible_outputs = ['RT-Predict', 'ISTDsEtc', 'FinalEMA-HILIC']
    outputs = possible_outputs[:(1+possible_outputs.index(output_type))]
    source_notebook = {'RT-Predict': 'Workflow_Notebook_VS_Auto_RT_Predict_V2.ipynb',
                       'ISTDsEtc': 'Targeted.ipynb',
                       'FinalEMA-HILIC': 'Targeted.ipynb'}
    source_atlas_prefix = {
            'RT-Predict': None,
            'ISTDsEtc': 'HILICz150_ANT20190824_PRD_IS_LabUnlab2_',
            'FinalEMA-HILIC': 'HILICz150_ANT20190824_TPL_EMA_Unlab_'}
    parameters = {
            'experiment': experiment_id,
            'metatlas_repo_path': repo_dir(),
            'output_directory': base_output_dir,
            'analysis_number': analysis_number}
    analysis_id = f"{getpass.getuser()}{parameters['analysis_number']}"
    tokens = parameters['experiment'].split('_')
    output_dir = os.path.join(base_output_dir, experiment_id)
    os.makedirs(output_dir, exist_ok=True)
    for output in outputs:
        parameters['output_type'] = output
        for polarity in (['positive', 'negative'] if output != 'RT-Predict' else [None]):
            source = os.path.join(repo_dir(), 'notebooks', 'reference', source_notebook[output])
            if polarity is not None:
                parameters['polarity'] = polarity
                pol = polarity[:3].upper()
                parameters['source_atlas'] = f"{source_atlas_prefix[output]}_{pol}_{tokens[3]}_{analysis_id}"
            generate_notebook(source, output_dir, parameters)


def generate_notebook(source, output_dir, parameters):
    """
    Creates a notebook from source in output_dir that has updated parameters.
    inputs:
        source: path of input Jupyter notebook
        output_dir: directory to write output Jupyter notebook
        parameters: dict of parameters to update in the notebook
    parameters must have atleast the following keys: analysis_number, experiment, output_type
    """
    if 'polarity' in parameters:
        pol = parameters['polarity'][:3].upper()
        suffix = f"{parameters['output_type']}_{pol}"
    else:
        suffix = 'RT-Predict'
    tokens = parameters['experiment'].split('_')
    dest = os.path.join(output_dir, '_'.join(tokens[3:5]+[suffix])+'.ipynb')
    create_notebook_with_parameters(source, dest, parameters)


def create_notebook_with_parameters(source, dest, parameters):
    with open(source) as source_fh:
        data = json.load(source_fh)
    eq_pat = re.compile('^[^#]')
    for line in data['cells'][1]['source']:
        if '=' in line:
            print(line)


def validate_data_dir(base_data_dir, experiment_id):
    """Raise FileNotFoundError if base_data_dir / experiment_id is not an existing directory"""
    experiment_dir = os.path.join(base_data_dir, experiment_id)
    try:
        if not os.path.isdir(experiment_dir):
            raise FileNotFoundError(f"Data directory does not exist at {experiment_dir}.")
    except FileNotFoundError as err:
        logger.exception(err)
        raise err

source = '/global/homes/w/wjholtz/metatlas-dev/notebooks/reference/Targeted.ipynb'
create_notebook_with_parameters(source, 'foo', {})
