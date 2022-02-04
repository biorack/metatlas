"""Functions for multi-CPU execution"""

import logging
import multiprocessing
from tqdm.notebook import tqdm

logger = logging.getLogger(__name__)


def parallel_process(function, data, max_cpus, unit='it'):
    """
    performs imap_unordered(function, data) using multiprocessing module but
    adds a progress bar and bypasses multiprocessing in the 1 cpu case as this makes debugging easier
    inputs:
        function: the function to apply
        data: iterater containing the inputs to function
        max_cpus: number of cpus to use
        unit: string label for what is processed in one iteration, default 'it'
    """
    if max_cpus > 1 and len(data) > 1:
        logger.debug("Starting parallel processing of %s with %d cpus.", function.__name__, max_cpus)
        with multiprocessing.Pool(processes=min(max_cpus, len(data))) as pool:
            return list(tqdm(pool.imap_unordered(function, data), total=len(data), unit=unit))
    logger.debug("Processing of %s with 1 cpu.", function.__name__)
    return list(map(function, data))
