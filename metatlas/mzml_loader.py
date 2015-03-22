import pymzml
import os, pwd
import csv
import sys
import struct
from scidbpy import interface, SciDBQueryError, SciDBArray
import random, string
from pymongo import MongoClient
from pyteomics import mzml
from tempfile import NamedTemporaryFile
import datetime

DEBUG=False
# MONGO_URL = "mongodb://metatlas_admin:b3igeISthenewbl%40ck@localhost:27017/metatlas"
MONGO_URL = "mongodb://metatlas_admin:b3igeISthenewbl%40ck@mongodb01.nersc.gov:27017/metatlas"
SCIDB_ARRAYNAME="lcms_test_1"
SCIDB_IP="http://128.55.57.21:22800"

try:
    from local_settings import *
except ImportError:
    pass

def load_bin_to_scidb(file_name, scidb_arrayname):
    """Loads binary data into SciDB; Returns the name of the SciDB array

    Keyword Arguments:
    file_name -- name of binary file to be loaded
    scidb_arrayname -- array in scidb to load to
    """

    if DEBUG: print "STATUS: Loading %s to SciDB (%s)..." % (file_name, scidb_arrayname)
    # Connect to SciDB
    sdb = interface.SciDBShimInterface(SCIDB_IP)

    tmparrayname = 'tmp_load_'+''.join(random.choice(string.lowercase+string.digits) for i in range(10))

    # Create flat array for initial data load
    if DEBUG: print "Creating flat loading array: %s"  % tmparrayname
    querystr = "CREATE ARRAY %(arr_name)s <fileid:int32,MZ:double, RT:double, I:double,P:int32,L:int32,precursorMZ:double, precursorIntensity:double,collisionEnergy:double> [i=0:*, 500000,0]" % {
        "arr_name": tmparrayname,
    }
    sdb._execute_query(querystr)

    # Load initial data into flat array
    if DEBUG: print "Loading %s to flat array: %s" % (file_name, tmparrayname)
    querystr = "load(%(arr_name)s,'%(bin_file)s',-2, '(int64,double,double,double,int32,int32,double,double,double)')" % {
        "arr_name": tmparrayname,
        "bin_file": file_name,
    }

    sdb._execute_query(querystr)

    arrname = tmparrayname + "_loaded"
    # Create new, redimensioned array template
    if DEBUG: print "Creating correctly dimensioned array: %s" % arrname
    querystr = "CREATE ARRAY %(arr_name)s <MZ:double,RT:double,I:double,precursorMZ:double,precursorIntensity:double,collisionEnergy:double>[D_L=1:2,1,0, D_P=0:1,1,0, D_file=0:100000,10,0, D_RT=%(D_RT_lower)d:%(D_RT_upper)d,%(D_RT_chunk_size)d,0, D_LOGI=%(D_LOGI_lower)d:%(D_LOGI_upper)d,%(D_LOGI_chunk_size)d,0, D_LOGMZ=%(D_LOGMZ_lower)d:%(D_LOGMZ_upper)d,%(D_LOGMZ_chunk_size)d,0]" % {
        "arr_name": arrname,
        "D_LOGMZ_lower": 1000000,
        "D_LOGMZ_upper": 3500000,
        "D_LOGMZ_chunk_size": 200000,
        "D_LOGI_lower": 0,
        "D_LOGI_upper": 800000000,
        "D_LOGI_chunk_size": 2000000,
        "D_RT_lower": 0,
        "D_RT_upper": 80000000,
        "D_RT_chunk_size": 2000000,
    }
    sdb._execute_query(querystr)

    # Redimension flat array into new array
    if DEBUG: print "Redimensioning flat array"
    querystr = "redimension_store(apply(%(source_arr)s, D_LOGMZ, int64(log10(MZ)*1000000), D_LOGI, int64(log10(I)*1000000), D_RT, int64(RT*1000000), D_L, L,D_P,P,D_file,fileid), %(dest_arr)s)" % {
        "source_arr": tmparrayname,
        "dest_arr": arrname,
    }
    sdb._execute_query(querystr)

    # Check if large array exist (otherwise create it)
    # This query will throw an exception if the array doesn't exist
    querystr = "show(%s)" % scidb_arrayname
    array_exists=True
    try:
        sdb._execute_query(querystr)
    except SciDBQueryError:
        array_exists=False

    if array_exists:
        # Merge into large array    
        if DEBUG: print "Inserting %s into %s" % (arrname, scidb_arrayname)
        querystr = "insert(%s, %s)" % (arrname, scidb_arrayname)
        sdb._execute_query(querystr)
        # Remove loaded array from scidb
        querystr = "remove(%s)" % arrname
        sdb._execute_query(querystr)
    else:
        if DEBUG: print "Inserting %s into %s" % (arrname, scidb_arrayname)
        querystr = "rename(%(source_arr)s,%(dest_arr)s)" % {
            "source_arr": arrname,
            "dest_arr": scidb_arrayname,
        }
        sdb._execute_query(querystr)
    
    if DEBUG: print "Removing temporary arrays"
    # Remove temporary flat array from SciDB
    querystr = "remove(%s)" % tmparrayname
    sdb._execute_query(querystr)

    if DEBUG: print "STATUS: Done"

    return arrname

def mzml_to_bin(in_file_name, out_file_name, fileid):
    """Converts in_file (mzml) to binary and stores it in out_file
    """
    out_file = open(out_file_name, "w+b")
    if DEBUG: print "STATUS: Converting %s to %s (mzML to binary)" % (in_file_name, out_file_name)

    # Extra accessions for pymzml to read
    extraAccessions = [
        ('MS:1000016',['value','unitName']), # scan start time
        ('MS:1000129',['name']), # negative scan
        ('MS:1000130',['name']), # positive scan
        ('MS:1000744',['name','value']), # selected ion m/z
        ('MS:1000042',['name','value']), # peak intensity
        ('MS:1000045',['name','value']), # collision energy
    ]

    mzml_reader = pymzml.run.Reader(in_file_name, extraAccessions=extraAccessions)

    for spectrum in mzml_reader:
        try:
            polarity = 1 if 'positive scan' in spectrum.keys() else 0
            ms_level = spectrum['ms level']
            #check if scan start time exists. some thermo spectra 
        #are missing this value
        if 'scan start time' in spectrum.keys():
            scan_time = spectrum['scan start time'][0]
        else:
            scan_time = spectrum['scan start time'][0]
        except KeyError:
            continue

        for mz, i in spectrum.peaks:
            precursor_MZ = 0.0
            precursor_intensity = 0.0
            collision_energy = 0.0

            if ms_level == 2:
                collision_energy=spectrum['collision energy'][1]
                if 'peak intensity' in spectrum.keys():
                    precursor_intensity=spectrum['peak intensity'][1]
                else:
                    precursor_intensity=0.0
                precursor_MZ= spectrum['selected ion m/z'][0]
            
            # Bundle information into binary struct
            package = struct.pack('qdddiiddd',
                                  int(fileid),
                                  float(mz), 
                                  float(scan_time),
                                  float(i),
                                  int(polarity),
                                  int(ms_level),
                                  float(precursor_MZ),
                                  float(precursor_intensity),
                                  float(collision_energy))
            out_file.write(package)
    out_file.close()
    os.chown(out_file_name, os.getuid(), 60734)
    os.chmod(out_file_name, 0770)
    if DEBUG: print "STATUS: Finished mzML to binary conversion"


def from_mongo(fileid, array_name, output_name=None, input_name=None):
    conn = MongoClient(MONGO_URL)
    db = conn['metatlas']
    col = db['metadata']
    if output_name is None:
        output = NamedTemporaryFile(dir="/global/project/projectdirs/metatlas/loading/", delete=False)
        output.close()
        output_name = output.name
    if input_name is None:
        input_name = col.find_one({"_id.file_id": fileid, "_id.array_name":array_name})['in_file']

    # Add some metadata to MongoDB
    col.update({"_id.file_id": fileid, "_id.array_name":array_name}, {
        "$set": {
            "out_file": output_name,
        }
    })

    # Load mzml file to the temporary output file
    mzml_to_bin(input_name, output_name, fileid)
    # Load output file to scidb
    load_bin_to_scidb(output_name, array_name)
    # Remove temporary file
    if DEBUG: print "STATUS: Removing temporary files"
    os.remove(output.name)

    # Update MongoDB 'pending' status
    col.update({"_id.file_id": fileid, "_id.array_name":array_name}, {
        "$set": {
            "pending": 0,
            "upload_date": datetime.datetime.utcnow(),
            "uploaded_by": pwd.getpwuid(os.getuid())[0],
        }
    })

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Load mzml files to SciDB")
    parser.add_argument("fileid", type=int, nargs=1, help="File ID for SciDB")
    parser.add_argument("array_name", nargs=1, type=str, help="SciDB array name")
    parser.add_argument("-o", "--output", nargs=1, type=str, help="Output file name", required=False)
    parser.add_argument("-i", "--input", type=str, nargs=1, help="Input mzML file", required=False)
    parser.add_argument("-d", "--debug", help="Sets debug mode", action="store_true")

    args = parser.parse_args()
    kwargs = {}
    if args.output:
        kwargs['output_name'] = args.output[0]
    if args.input:
        kwargs['input_name'] = args.array[0]
    
    # Toggles debug mode base on --debug flag
    DEBUG = args.debug
    
    from_mongo(args.fileid[0], args.array_name[0], **kwargs)
