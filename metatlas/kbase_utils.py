"""Selected functions from https://github.com/kbase/transform/blob/master/lib/biokbase/Transform/script_utils.py
For use in Travis testing
"""

import logging
import re
import sys
import time
import requests
import os
import io
from requests_toolbelt import MultipartEncoder


def stderrlogger(name, level=logging.INFO):
    """
    Return a standard python logger with a stderr handler attached and using a prefix
    format that will make logging consistent between scripts.
    """

    logger = logging.getLogger(name)
    logger.setLevel(level)

    # send messages to sys.stderr
    streamHandler = logging.StreamHandler(sys.stderr)

    formatter = logging.Formatter("%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
    formatter.converter = time.gmtime
    streamHandler.setFormatter(formatter)

    logger.addHandler(streamHandler)

    return logger


def parse_docs(docstring=None):
    """
    Parses the docstring of a function and returns a dictionary of the elements.
    """

    # TODO, revisit this, probably can use other ways of doing this
    script_details = dict()

    keys = ["Authors","Returns","Args"]

    remainder = docstring[:]
    for k in keys:
        remainder, script_details[k] = remainder.split(k+":",1)
        script_details[k] = script_details[k].strip()


    script_details["Description"] = remainder

    # special treatment for Args since we want a dict, split on :, then cleanup whitespace
    # keep the : in the keys to do clean splits when getting the values
    argument_keys = [x.strip() for x in re.findall(".*:",script_details["Args"])]

    # split on the string in reverse by the keys, then wash out the extra whitespace
    remainder = script_details["Args"]
    argument_values = list()
    for k in reversed(argument_keys):
        remainder, _, value = remainder.rpartition(k)
        argument_values.append(" ".join([x.strip() for x in value.split("\n")]))

    # create the dict using they keys without :, then get the values in the correct order
    script_details["Args"] = dict(zip([x.replace(":","") for x in argument_keys], reversed(argument_values)))

    return script_details


def upload_file_to_shock(logger = stderrlogger(__file__),
                         shock_service_url = None,
                         filePath = None,
                         ssl_verify = True,
                         token = None):
    """
    Use HTTP multi-part POST to save a file to a SHOCK instance.
    """

    if token is None:
        raise Exception("Authentication token required!")
    
    #build the header
    header = dict()
    header["Authorization"] = "Oauth {0}".format(token)

    if filePath is None:
        raise Exception("No file given for upload to SHOCK!")

    dataFile = open(os.path.abspath(filePath), 'rb')
    m = MultipartEncoder(fields={'upload': (os.path.split(filePath)[-1], dataFile)})
    header['Content-Type'] = m.content_type

    logger.info("Sending {0} to {1}".format(filePath,shock_service_url))

    try:
        response = requests.post(shock_service_url + "/node", headers=header, data=m, allow_redirects=True, verify=ssl_verify)
        dataFile.close()
    except:
        dataFile.close()
        raise    

    if not response.ok:
        response.raise_for_status()

    result = response.json()

    if result['error']:     
        raise Exception(result['error'][0])
    else:
        return result["data"]   


def download_file_from_shock(logger = stderrlogger(__file__),
                             shock_service_url = None,
                             shock_id = None,
                             filename = None,
                             directory = None,
                             token = None):
    """
    Given a SHOCK instance URL and a SHOCK node id, download the contents of that node
    to a file on disk.
    """

    header = dict()
    header["Authorization"] = "Oauth {0}".format(token)

    logger.info("Downloading shock node {0}/node/{1}".format(shock_service_url,shock_id))

    metadata_response = requests.get("{0}/node/{1}?verbosity=metadata".format(shock_service_url, shock_id), headers=header, stream=True, verify=True)
    shock_metadata = metadata_response.json()['data']
    shockFileName = shock_metadata['file']['name']
    shockFileSize = shock_metadata['file']['size']
    metadata_response.close()
        
    download_url = "{0}/node/{1}?download_raw".format(shock_service_url, shock_id)
        
    data = requests.get(download_url, headers=header, stream=True, verify=True)

    if filename is not None:
        shockFileName = filename

    if directory is not None:
        filePath = os.path.join(directory, shockFileName)
    else:
        filePath = shockFileName

    chunkSize = shockFileSize/4
    
    maxChunkSize = 2**30
    
    if chunkSize > maxChunkSize:
        chunkSize = maxChunkSize
    
    f = io.open(filePath, 'wb')
    try:
        for chunk in data.iter_content(chunkSize):
            if chunk:                
                f.write(chunk)
                f.flush()            
    finally:
        data.close()
        f.close()


def upload_from_nersc(user, relative_path):
    """Load a data file from NERSC to an H5 file

    Parameters
    ----------
    user : str
        NERSC user account
    relative_path : str
        Path to file from "/project/projectdirs/metatlas/original_data/<user>/"
    """
    import pexpect
    from IPython.display import clear_output

    cmd = 'scp -o StrictHostKeyChecking=no %s@edisongrid.nersc.gov:/project/projectdirs/metatlas/original_data/%s/%s . && echo "Download Complete"'
    cmd = cmd % (user, user, relative_path)
    print(cmd)
    return
    proc = pexpect.spawn(cmd)
    proc.expect("assword:*")
    if sys.version.startswith('3'):
        passwd = input()
    else:
        passwd = raw_input()
    clear_output()
    proc.send(passwd)
    proc.send('\r')
    proc.expect('Download Complete')
    proc.close()
    return os.path.abspath(os.path.basename(relative_path))
