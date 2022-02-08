
import os
from app import app
import ftputil
import credentials
import json
import requests
from werkzeug.utils import secure_filename

ALLOWED_EXTENSIONS = set(['mgf', 'mzxml', 'mzml', 'csv', 'txt', 'raw', 'msp'])


def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


def thermoconvert_localfile(input_filename, save_dir):
    extension = input_filename.rsplit('.', 1)[-1].lower()
    print(extension)

    """Do Nothing"""
    if extension != "raw":
        return input_filename

    """Perform Conversion"""
    cmd = "mono /src/bin/x64/Debug/ThermoRawFileParser.exe -i=%s -o=%s -f=1" % (input_filename, save_dir)
    os.system(cmd)
    os.remove(input_filename)
    output_filename = os.path.join(save_dir, os.path.basename(input_filename).replace(".raw", ".mzML"))
    return output_filename

def upload_single_file(request, group):
    sessionid = request.cookies.get('sessionid')

    filename = ""

    if 'file' not in request.files:
        return "{}"
    request_file = request.files['file']

    return upload_single_file_push(request_file, sessionid, group)

def upload_single_file_push(request_file, uuid_folder, collection_name):
    if request_file.filename == '':
        return "{}"
    if request_file and allowed_file(request_file.filename):
        filename = secure_filename(request_file.filename)
        save_dir = os.path.join(app.config['UPLOAD_FOLDER'], uuid_folder, collection_name)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        local_filename = os.path.join(save_dir, filename)
        request_file.save(local_filename)

        """If we need to convert raw file, we do it here"""
        local_filename = thermoconvert_localfile(local_filename, save_dir)

        #Uploading to FTP
        upload_to_gnps(local_filename, uuid_folder, collection_name)

        #Remove local file
        os.remove(local_filename)
    else:
        print("not allowed")
        return json.dumps({"status": "Invalid File Type"})

    return json.dumps({"filename": filename})


def check_ftp_folders(username):
    url = "ccms-ftp01.ucsd.edu"
    present_folders = []

    with ftputil.FTPHost(url, credentials.USERNAME, credentials.PASSWORD) as ftp_host:
        names = ftp_host.listdir(ftp_host.curdir)
        if not username in names:
            return present_folders

        ftp_host.chdir(username)

        return ftp_host.listdir(ftp_host.curdir)

    return present_folders


def upload_to_gnps(input_filename, folder_for_spectra, group_name, username=credentials.USERNAME, password=credentials.PASSWORD):
    url = "ccms-ftp01.ucsd.edu"

    with ftputil.FTPHost(url, username, password) as ftp_host:
        names = ftp_host.listdir(ftp_host.curdir)
        try:
            if not folder_for_spectra in names:
                print("MAKING DIR")
                ftp_host.mkdir(folder_for_spectra)
        except:
            print("Cannot Make Folder", folder_for_spectra)

        ftp_host.chdir(folder_for_spectra)
        try:
            if not group_name in ftp_host.listdir(ftp_host.curdir):
                print("MAKING Group DIR")
                ftp_host.mkdir(group_name)
        except:
            print("Cannot Make Folder", group_name)
        ftp_host.chdir(group_name)

        ftp_host.upload(input_filename, os.path.basename(input_filename))

def get_classic_networking_lowres_parameters():
    invokeParameters = {}
    invokeParameters["workflow"] = "METABOLOMICS-SNETS-V2"
    invokeParameters["protocol"] = "None"
    invokeParameters["workflow_version"] = "release_22"
    invokeParameters["library_on_server"] = "d.speclibs;"
    invokeParameters["tolerance.PM_tolerance"] = "2.0"
    invokeParameters["tolerance.Ion_tolerance"] = "0.5"
    invokeParameters["PAIRS_MIN_COSINE"] = "0.70"
    invokeParameters["MIN_MATCHED_PEAKS"] = "6"
    invokeParameters["TOPK"] = "10"
    invokeParameters["CLUSTER_MIN_SIZE"] = "2"
    invokeParameters["RUN_MSCLUSTER"] = "on"
    invokeParameters["MAXIMUM_COMPONENT_SIZE"] = "100"
    invokeParameters["MIN_MATCHED_PEAKS_SEARCH"] = "6"
    invokeParameters["SCORE_THRESHOLD"] = "0.7"
    invokeParameters["ANALOG_SEARCH"] = "0"
    invokeParameters["MAX_SHIFT_MASS"] = "100.0"
    invokeParameters["FILTER_STDDEV_PEAK_datasetsINT"] = "0.0"
    invokeParameters["MIN_PEAK_INT"] = "0.0"
    invokeParameters["FILTER_PRECURSOR_WINDOW"] = "1"
    invokeParameters["FILTER_LIBRARY"] = "1"
    invokeParameters["WINDOW_FILTER"] = "1"
    invokeParameters["CREATE_CLUSTER_BUCKETS"] = "1"
    invokeParameters["CREATE_ILI_OUTPUT"] = "0"
    invokeParameters["FILTER_G6_BLANKS"] = "0"
    

    return invokeParameters

def get_classic_networking_highres_parameters():
    invokeParameters = {}
    invokeParameters["workflow"] = "METABOLOMICS-SNETS-V2"
    invokeParameters["protocol"] = "None"
    invokeParameters["workflow_version"] = "release_22"
    invokeParameters["library_on_server"] = "d.speclibs;"
    invokeParameters["tolerance.PM_tolerance"] = "0.05"
    invokeParameters["tolerance.Ion_tolerance"] = "0.05"
    invokeParameters["PAIRS_MIN_COSINE"] = "0.70"
    invokeParameters["MIN_MATCHED_PEAKS"] = "6"
    invokeParameters["TOPK"] = "10"
    invokeParameters["CLUSTER_MIN_SIZE"] = "2"
    invokeParameters["RUN_MSCLUSTER"] = "on"
    invokeParameters["MAXIMUM_COMPONENT_SIZE"] = "100"
    invokeParameters["MIN_MATCHED_PEAKS_SEARCH"] = "6"
    invokeParameters["SCORE_THRESHOLD"] = "0.7"
    invokeParameters["ANALOG_SEARCH"] = "0"
    invokeParameters["MAX_SHIFT_MASS"] = "100.0"
    invokeParameters["FILTER_STDDEV_PEAK_datasetsINT"] = "0.0"
    invokeParameters["MIN_PEAK_INT"] = "0.0"
    invokeParameters["FILTER_PRECURSOR_WINDOW"] = "1"
    invokeParameters["FILTER_LIBRARY"] = "1"
    invokeParameters["WINDOW_FILTER"] = "1"
    invokeParameters["CREATE_CLUSTER_BUCKETS"] = "1"
    invokeParameters["CREATE_ILI_OUTPUT"] = "0"
    invokeParameters["FILTER_G6_BLANKS"] = "0"

    return invokeParameters

def launch_GNPS_workflow(ftp_path, job_description, username, password, groups_present, email, preset):
    invokeParameters = {}

    if preset == "LOWRES":
        invokeParameters = get_classic_networking_lowres_parameters()
    elif preset == "HIGHRES":
        invokeParameters = get_classic_networking_highres_parameters()
    else:
        return "Error No Preset"

    invokeParameters["desc"] = job_description
    invokeParameters["spec_on_server"] = "d." + ftp_path + "/G1;"
    if "G2" in groups_present:
        invokeParameters["spec_on_server_group2"] = "d." + ftp_path + "/G2;"
    if "G3" in groups_present:
        invokeParameters["spec_on_server_group3"] = "d." + ftp_path + "/G3;"

    invokeParameters["email"] = email


    task_id = invoke_workflow("gnps.ucsd.edu", invokeParameters, username, password)

    return task_id

def launch_GNPS_featurenetworking_workflow(ftp_path, job_description, username, password, email, featuretool, present_folders, preset,metadata_filename,mgf_filename,quant_filename):
    invokeParameters = {}

    if preset == "LOWRES":
        invokeParameters = get_featurenetworking_lowres_parameters()
    elif preset == "HIGHRES":
        invokeParameters = get_featurenetworking_highres_parameters()
    else:
        return "Error No Preset"

    #Specific Parameters Update
    invokeParameters["desc"] = job_description

#     invokeParameters["quantification_table"] = "d." + ftp_path + "/featurequantification;"
#     invokeParameters["spec_on_server"] = "d." + ftp_path + "/featurems2;"
#     if "samplemetadata" in present_folders:
#         invokeParameters["metadata_table"] = "d." + ftp_path + "/samplemetadata;"
    invokeParameters["quantification_table"] = "d.%s/"%username + ftp_path + "/featurequantification/%s;"%quant_filename
    invokeParameters["spec_on_server"] = "d.%s/"%username + ftp_path + "/featurems2/%s;"%mgf_filename
    invokeParameters["metadata_table"] = "d.%s/"%username + ftp_path + "/samplemetadata/%s;"%metadata_filename

    # This is a folder of mzml files 
    # make a subfolder called "mzml_pos"
    # make a subfolder 
    # It will look just like this, but a differnt folder 
        # invokeParameters["metadata_table"] = "d.%s/"%username + ftp_path + "/samplemetadata/%s;"%metadata_filename
    invokeParameters["raw_spectra"] = "d.%s/"%username + ftp_path + "/rawdata/;"
    
    #Quant
    invokeParameters["QUANT_TABLE_SOURCE"] = featuretool

    #Additional Pairs
    if "additionalpairs" in present_folders:
        invokeParameters["additional_pairs"] = "d." + ftp_path + "/additionalpairs;"

    invokeParameters["email"] = email

    task_id = invoke_workflow("gnps.ucsd.edu", invokeParameters, username, password)

    return task_id

def get_featurenetworking_lowres_parameters():
    invokeParameters = {}
    invokeParameters["workflow"] = "FEATURE-BASED-MOLECULAR-NETWORKING"
    invokeParameters["protocol"] = "None"
    invokeParameters["workflow_version"] = "release_27"
    invokeParameters["desc"] = "Job Description"
    invokeParameters["library_on_server"] = "d.speclibs;"

    #Networking
    invokeParameters["tolerance.PM_tolerance"] = "2.0"
    invokeParameters["tolerance.Ion_tolerance"] = "0.5"
    invokeParameters["PAIRS_MIN_COSINE"] = "0.70"
    invokeParameters["MIN_MATCHED_PEAKS"] = "6"
    invokeParameters["TOPK"] = "10"
    invokeParameters["MAX_SHIFT"] = "500"

    #Network Pruning
    invokeParameters["MAXIMUM_COMPONENT_SIZE"] = "100"

    #Library Search
    invokeParameters["MIN_MATCHED_PEAKS_SEARCH"] = "6"
    invokeParameters["SCORE_THRESHOLD"] = "0.7"
    invokeParameters["TOP_K_RESULTS"] = "1"
    invokeParameters["ANALOG_SEARCH"] = "0"
    invokeParameters["MAX_SHIFT_MASS"] = "100.0"
    invokeParameters["FILTER_STDDEV_PEAK_datasetsINT"] = "0.0"
    invokeParameters["MIN_PEAK_INT"] = "0.0"
    invokeParameters["FILTER_PRECURSOR_WINDOW"] = "1"
    invokeParameters["FILTER_LIBRARY"] = "1"
    invokeParameters["WINDOW_FILTER"] = "1"

    #Quant
    invokeParameters["QUANT_TABLE_SOURCE"] = ""
    invokeParameters["GROUP_COUNT_AGGREGATE_METHOD"] = "Mean"
    invokeParameters["QUANT_FILE_NORM"] = "RowSum"

    # Stats
    invokeParameters["RUN_STATS"] = "No"
    invokeParameters["METADATA_COLUMN"] = "None"
    invokeParameters["METADATA_COLUMN_FACET"] = "None"
    invokeParameters["METADATA_CONDITION_ONE"] = "None"
    invokeParameters["METADATA_CONDITION_TWO"] = "None"

    #External tools
    invokeParameters["RUN_DEREPLICATOR"] = "1"

    # Qiime2
    invokeParameters["QIIME2_PCOA_DISTANCE"] = "cosine"

    # Metadata
    invokeParameters["googlesheetsmetadata"] = "None"

    invokeParameters["email"] = "ccms.web@gmail.com"
    invokeParameters["uuid"] = "1DCE40F7-1211-0001-979D-15DAB2D0B500"

    return invokeParameters

def get_featurenetworking_highres_parameters():
    invokeParameters = {}
    invokeParameters["workflow"] = "FEATURE-BASED-MOLECULAR-NETWORKING"
    invokeParameters["protocol"] = "None"
    invokeParameters["workflow_version"] = "release_28.2"
    invokeParameters["desc"] = "Job Description"
    invokeParameters["library_on_server"] = "d.speclibs;d.staticlibs/NIST_17.1;"

    #Networking
    invokeParameters["tolerance.PM_tolerance"] = "0.01"
    invokeParameters["tolerance.Ion_tolerance"] = "0.02"
    invokeParameters["PAIRS_MIN_COSINE"] = "0.70"
    invokeParameters["MIN_MATCHED_PEAKS"] = "3"
    invokeParameters["TOPK"] = "10"
    invokeParameters["MAX_SHIFT"] = "500"

    #Network Pruning
    invokeParameters["MAXIMUM_COMPONENT_SIZE"] = "0"

    #Library Search
    invokeParameters["MIN_MATCHED_PEAKS_SEARCH"] = "3"
    invokeParameters["SCORE_THRESHOLD"] = "0.4"
    invokeParameters["TOP_K_RESULTS"] = "20"
    invokeParameters["ANALOG_SEARCH"] = "0"
    invokeParameters["MAX_SHIFT_MASS"] = "200.0"
    invokeParameters["FILTER_STDDEV_PEAK_datasetsINT"] = "0.0"
    invokeParameters["MIN_PEAK_INT"] = "0.0"
    invokeParameters["FILTER_PRECURSOR_WINDOW"] = "1"
    invokeParameters["FILTER_LIBRARY"] = "1"
    invokeParameters["WINDOW_FILTER"] = "1"

    #Quant
    invokeParameters["QUANT_TABLE_SOURCE"] = ""
    invokeParameters["GROUP_COUNT_AGGREGATE_METHOD"] = "Mean"
    invokeParameters["QUANT_FILE_NORM"] = "None"

    # Stats
    invokeParameters["RUN_STATS"] = "Yes"
    invokeParameters["METADATA_COLUMN"] = "ATTRIBUTE_sampletype"
    invokeParameters["METADATA_COLUMN_FACET"] = "None"
    invokeParameters["METADATA_CONDITION_ONE"] = "None"
    invokeParameters["METADATA_CONDITION_TWO"] = "None"

    #External tools
    invokeParameters["RUN_DEREPLICATOR"] = "1"

    # Qiime2
    invokeParameters["QIIME2_PCOA_DISTANCE"] = "cosine"

    # Metadata
    invokeParameters["googlesheetsmetadata"] = "None"

    invokeParameters["email"] = "ben.bowen@gmail.com"
    invokeParameters["uuid"] = "1DCE40F7-1211-0001-979D-15DAB2D0B500"

    return invokeParameters




def invoke_workflow(base_url, parameters, login, password):
    username = login
    password = password

    s = requests.Session()

    payload = {
        'user' : username,
        'password' : password,
        'login' : 'Sign in'
    }

    r = s.post('https://' + base_url + '/ProteoSAFe/user/login.jsp', data=payload, verify=False)
    r = s.post('https://' + base_url + '/ProteoSAFe/InvokeTools', data=parameters, verify=False)
    task_id = r.text

    import sys
    print(r.text, file=sys.stderr, flush=True)

    if len(task_id) > 4 and len(task_id) < 60:
        print("Launched Task: : " + r.text)
        return task_id
    else:
        print(task_id)
        return None
