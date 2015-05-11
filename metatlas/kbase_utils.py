"""Selected functions from https://github.com/kbase/transform/blob/master/lib/biokbase/Transform/script_utils.py
For use in Travis testing
"""

import logging
import re
import sys
import time
import os
import io

import requests
from requests_toolbelt import MultipartEncoder
import tables


SHOCK_URL = 'https://ci.kbase.us/services/shock-api'
WS_URL = 'https://ci.kbase.us/services/ws/'
HANDLE_URL = 'https://ci.kbase.us/services/handle_service'


def stderrlogger(name, level=logging.INFO):
    """
    Return a standard python logger with a stderr handler attached and using a
    prefix format that will make logging consistent between scripts.
    """

    logger = logging.getLogger(name)
    logger.setLevel(level)

    # send messages to sys.stderr
    streamHandler = logging.StreamHandler(sys.stderr)

    formatter = logging.Formatter(
        "%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
    formatter.converter = time.gmtime
    streamHandler.setFormatter(formatter)

    logger.addHandler(streamHandler)

    return logger


def parse_docs(docstring=None):
    """
    Parses the docstring of a function and returns a dictionary of the
    elements.
    """

    # TODO, revisit this, probably can use other ways of doing this
    script_details = dict()

    keys = ["Authors", "Returns", "Args"]

    remainder = docstring[:]
    for k in keys:
        remainder, script_details[k] = remainder.split(k + ":", 1)
        script_details[k] = script_details[k].strip()

    script_details["Description"] = remainder

    # special treatment for Args since we want a dict, split on :, then
    # cleanup whitespace
    # keep the : in the keys to do clean splits when getting the values
    argument_keys = [x.strip() for x in re.findall(".*:",
                                                   script_details["Args"])]

    # split on the string in reverse by the keys, then wash out the extra
    # whitespace
    remainder = script_details["Args"]
    argument_values = list()
    for k in reversed(argument_keys):
        remainder, _, value = remainder.rpartition(k)
        argument_values.append(" ".join([x.strip()
                                         for x in value.split("\n")]))

    # create the dict using they keys without :, then get the values in the
    # correct order
    script_details["Args"] = dict(zip([x.replace(":", "")
                                       for x in argument_keys],
                                      reversed(argument_values)))

    return script_details


def upload_file_to_shock(logger=stderrlogger(__file__),
                         shock_service_url=SHOCK_URL,
                         filePath=None,
                         ssl_verify=True,
                         token=None):
    """
    Use HTTP multi-part POST to save a file to a SHOCK instance.
    """

    token = token or os.environ.get('KB_AUTH_TOKEN', None)

    if token is None:
        raise Exception("Authentication token required!")

    # build the header
    header = dict()
    header["Authorization"] = "Oauth {0}".format(token)

    if filePath is None:
        raise Exception("No file given for upload to SHOCK!")

    dataFile = open(os.path.abspath(filePath), 'rb')
    m = MultipartEncoder(fields={'upload': (os.path.split(filePath)[-1],
                                            dataFile)})
    header['Content-Type'] = m.content_type

    logger.info("Sending {0} to {1}".format(filePath, shock_service_url))

    try:
        response = requests.post(shock_service_url + "/node", headers=header,
                                 data=m, allow_redirects=True,
                                 verify=ssl_verify)
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


def delete_file_from_shock(logger=stderrlogger(__file__),
                           shock_service_url=SHOCK_URL,
                           shock_id=None,
                           token=None):

    token = token or os.environ.get('KB_AUTH_TOKEN', None)

    if token is None:
        raise Exception("Authentication token required!")

    # build the header
    header = dict()
    header["Authorization"] = "Oauth {0}".format(token)

    logger.info("Deleting {0} from {1}".format(shock_id, shock_service_url))

    response = requests.delete(shock_service_url + "/node/" + str(shock_id),
                               headers=header, allow_redirects=True,
                               verify=True)

    if not response.ok:
        response.raise_for_status()

    result = response.json()

    if result['error']:
        raise Exception(result['error'][0])
    else:
        return result


def download_file_from_shock(logger=stderrlogger(__file__),
                             shock_service_url=SHOCK_URL,
                             shock_id=None,
                             filename=None,
                             directory=None,
                             token=None):
    """
    Given a SHOCK instance URL and a SHOCK node id, download the contents of
    that node to a file on disk.
    """
    token = token or os.environ.get('KB_AUTH_TOKEN', None)

    if token is None:
        raise Exception("Authentication token required!")

    header = dict()
    header["Authorization"] = "Oauth {0}".format(token)

    logger.info(
        "Downloading shock node {0}/node/{1}".format(shock_service_url,
                                                     shock_id))

    metadata_response = requests.get("{0}/node/{1}?verbosity=metadata".format(
        shock_service_url, shock_id), headers=header, stream=True, verify=True)
    shock_metadata = metadata_response.json()['data']
    shockFileName = shock_metadata['file']['name']
    shockFileSize = shock_metadata['file']['size']
    metadata_response.close()

    download_url = "{0}/node/{1}?download_raw".format(
        shock_service_url, shock_id)

    data = requests.get(download_url, headers=header, stream=True, verify=True)

    if filename is not None:
        shockFileName = filename

    if directory is not None:
        filePath = os.path.join(directory, shockFileName)
    else:
        filePath = shockFileName

    chunkSize = shockFileSize / 4

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


def getHandles(logger=stderrlogger(__file__),
               shock_service_url=SHOCK_URL,
               handle_service_url=HANDLE_URL,
               shock_ids=None,
               handle_ids=None,
               token=None):
    """
    Retrieve KBase handles for a list of shock ids or a list of handle ids.
    """
    from .handle_service import AbstractHandle as HandleService

    token = token or os.environ.get('KB_AUTH_TOKEN', None)
    if token is None:
        raise Exception("Authentication token required!")

    hs = HandleService(url=handle_service_url, token=token)

    handles = list()
    if shock_ids is not None:
        header = dict()
        header["Authorization"] = "Oauth {0}".format(token)

        for sid in shock_ids:
            info = None

            try:
                logger.info(
                    "Found shock id {0}, retrieving information about the data.".format(sid))

                response = requests.get(
                    "{0}/node/{1}".format(shock_service_url, sid), headers=header, verify=True)
                info = response.json()["data"]
            except:
                logger.error("There was an error retrieving information about the shock node id {0} from url {1}".format(
                    sid, shock_service_url))

            try:
                logger.info("Retrieving a handle id for the data.")
                handle = hs.persist_handle({"id": sid,
                                            "type": "shock",
                                            "url": shock_service_url,
                                            "file_name": info["file"]["name"],
                                            "remote_md5": info["file"]["checksum"]["md5"]})
                handles.append(handle)
            except:
                try:
                    handle_id = hs.ids_to_handles([sid])[0]["hid"]
                    single_handle = hs.hids_to_handles([handle_id])

                    assert len(single_handle) != 0

                    if info is not None:
                        single_handle[0]["file_name"] = info["file"]["name"]
                        single_handle[0]["remote_md5"] = info[
                            "file"]["checksum"]["md5"]
                        logger.debug(single_handle)

                    handles.append(single_handle[0])
                except:
                    logger.error(
                        "The input shock node id {} is already registered or could not be registered".format(sid))

                    hs = HandleService(url=handle_service_url, token=token)
                    all_handles = hs.list_handles()

                    for x in all_handles:
                        if x[0] == sid:
                            logger.info("FOUND shock id as existing handle")
                            logger.info(x)
                            break
                    else:
                        logger.info(
                            "Unable to find a handle containing shock id")

                        logger.info(
                            "Trying again to get a handle id for the data.")
                        handle_id = hs.persist_handle({"id": sid,
                                                       "type": "shock",
                                                       "url": shock_service_url,
                                                       "file_name": info["file"]["name"],
                                                       "remote_md5": info["file"]["checksum"]["md5"]})
                        handles.append(handle_id)

                    raise
    elif handle_ids is not None:
        for hid in handle_ids:
            try:
                single_handle = hs.hids_to_handles([hid])

                assert len(single_handle) != 0

                handles.append(single_handle[0])
            except:
                logger.error("Invalid handle id {0}".format(hid))
                raise

    return handles


def save_ws_object(obj):
    """Save an object to the workspace

    Parameters
    ----------
    obj : dict
        Object with the fields: type, data, and name.
        The type must be the full typespec
        (e.g. 'MetaboliteAtlas.Compound-0.3')

    Returns
    -------
    id : str
        Object workspace id
    """
    from biokbase.workspace.client import Workspace
    ws = Workspace(WS_URL)
    obj.setdefault('hidden', 0)
    wks = ws.list_workspaces({'excludeGlobal': 1})
    ws_id = [wk[-1] for wk in wks if wk[0] == os.environ['KB_WORKSPACE_ID']][0]
    save_objects_params = {'id': ws_id, 'objects': [obj]}
    return ws.save_objects(save_objects_params)[0][-2]


def get_object_uid(name):
    WS_URL = 'https://ci.kbase.us/services/ws/'
    from biokbase.workspace.client import Workspace
    ws = Workspace(WS_URL)
    info = ws.get_objects([dict(workspace=os.environ['KB_WORKSPACE_ID'],
                                name=name)])[0]['info']
    return '%s/%s/%s' % (info[6], info[0], info[4])


def get_object(name):
    WS_URL = 'https://ci.kbase.us/services/ws/'
    from biokbase.workspace.client import Workspace
    ws = Workspace(WS_URL)
    return ws.get_objects([dict(workspace=os.environ['KB_WORKSPACE_ID'],
                                name=name)])[0]['data']


def get_object_from_ref(ref):
    objid = int(ref.split('/')[1])
    WS_URL = 'https://ci.kbase.us/services/ws/'
    from biokbase.workspace.client import Workspace
    ws = Workspace(WS_URL)
    return ws.get_objects([dict(workspace=os.environ['KB_WORKSPACE_ID'],
                                objid=objid)])[0]['data']
