
from __future__ import print_function

import pandas as pd
import re
import glob as glob
import os
from collections import defaultdict
from xml.etree import cElementTree as ET

try:
    basestring
except NameError:  # python3
    basestring = str

BATCH_FILE_PATH = '/project/projectdirs/metatlas/projects/mzmine_parameters/batch_files/'
BINARY_PATH = '/project/projectdirs/metatlas/projects/mzmine_parameters/MZmine'

def configure_mass_detection(new_d,noise_floor):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'MassDetectionModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Mass detector' in d['@name']][0]
    idx3 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module']) if 'Centroid' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']['#text'] = '%.2f'%(noise_floor)
    return new_d

def configure_chromatogram_builder(new_d,min_peak_duration,min_peak_height,mz_tolerance):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'ChromatogramBuilderModule' in d['@method']][0]
    new_d['batch']['batchstep'][idx]['parameter']

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Min time span' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(min_peak_duration)

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Min height' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(min_peak_height)

    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    
    return new_d


def configure_peak_deconvolution(new_d,min_peak_height,min_sn_ratio,min_peak_duration,max_peak_duration):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'DeconvolutionModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Algorithm' in d['@name']][0]
    idx3 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module']) if 'Local minimum search' in d['@name']][0]
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Chromatographic threshold' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%0.1
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Search minimum in RT range (min)' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%0.05
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Minimum relative height' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%0.001
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Minimum absolute height' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%min_peak_height
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Min ratio of peak top/edge' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['#text'] = '%.3f'%min_sn_ratio
    idx4 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter']) if 'Peak duration range (min)' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['min'] = '%.3f'%min_peak_duration
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['module'][idx3]['parameter'][idx4]['max'] = '%.3f'%max_peak_duration
    
    #following deconvolution, many small peaks are created.  Filter them out here
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'PeakFilterModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Height' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['min'] = '%.3f'%(min_peak_height)
    
    return new_d

def configure_isotope_adduct_fragment_search(new_d,mz_tolerance,rt_tol_perfile,polarity,min_peak_height):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'Isotope' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_perfile)

    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'Adduct' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'RT tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_perfile)
    if polarity == 'negative':
        idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Adducts' in d['@name']][0]
        #the default is setup for positive mode adducts.
        #only change them if you are in negative mode
        for i,a in enumerate(new_d['batch']['batchstep'][idx]['parameter'][idx2]['adduct']):
            if a['@selected'] == 'true':
                new_d['batch']['batchstep'][idx]['parameter'][idx2]['adduct'][i]['@selected'] = 'false'
            else:
                new_d['batch']['batchstep'][idx]['parameter'][idx2]['adduct'][i]['@selected'] = 'true'

    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'ComplexSearchModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_perfile)
    if polarity == 'negative':
        #the default is setup for positive mode adducts.
        #only change them if you are in negative mode
        idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Ionization method' in d['@name']][0]
        new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '[M-H]-'

    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'FragmentSearchModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_perfile)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Min MS2 peak height' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(min_peak_height)
    return new_d

def configure_join_aligner(new_d,mz_tolerance,rt_tol_multifile):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'JoinAlignerModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Retention time tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_multifile)
    return new_d

def configure_duplicate_filter(new_d,mz_tolerance,rt_tol_perfile):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'DuplicateFilterModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'RT tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = '%.3f'%(rt_tol_perfile)
    return new_d

def configure_gap_filling(new_d,mz_tolerance):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'SameRangeGapFillerModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'm/z tolerance' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['ppmtolerance'] = '%.3f'%(mz_tolerance)
    return new_d

def configure_csv_output(new_d,output_csv):
    """
    
    """
    idx = [i for i,d in enumerate(new_d['batch']['batchstep']) if 'CSVExportModule' in d['@method']][0]
    idx2 = [i for i,d in enumerate(new_d['batch']['batchstep'][idx]['parameter']) if 'Filename' in d['@name']][0]
    new_d['batch']['batchstep'][idx]['parameter'][idx2]['#text'] = output_csv
    return new_d

def indent_tree(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent_tree(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def get_batch_file_template(loc='do_not_change_batch_file.xml'):
    """
    return string of text from the template batch file
    """
    with open(os.path.join(BATCH_FILE_PATH,loc),'r') as fid:
        file_text = fid.read()
    return file_text

def get_latest_mzmine_binary(system='Cori',version='most_recent'):
    """
    Returns the path to the mzmine launch script.
    Default is most recent.  Alternatively specify the folder containng version you want

    for example:
        version='MZmine-2.23'
    will use the launch script in that folder

    #
    cd /project/projectdirs/metatlas/projects/mzmine_parameters/MZmine
    wget $(curl -s https://api.github.com/repos/mzmine/mzmine2/releases/latest | grep 'browser_' | cut -d\" -f4) -O mzmine_latest.zip
    unzip mzmine_latest.zip
    cp ../MZmine-2.23/startMZmine_NERSC_* .
    """
    mzmine_versions = glob.glob(os.path.join(BINARY_PATH,'*'))
    if version == 'most_recent':
        most_recent = sorted([os.path.basename(m) for m in mzmine_versions if 'MZmine-' in m])[-1]
    else:
        most_recent = 'version'
    launch_script = os.path.join(os.path.join(BINARY_PATH,most_recent),'startMZmine_NERSC_Headless_%s.sh'%system)
    if os.path.isfile(launch_script):
        return launch_script
    else:
        print('See the docstring, the launch script seems to be missing.')

def replace_files(d,file_list):
    """
    Replace files for mzmine task
    
    Inputs:
    d: an xml derived dictionary of batch commands
    file_list: a list of full paths to mzML files
    
    Outputs:
    d: an xml derived dict with new files in it
    """
    for i,step in enumerate(d['batch']['batchstep']):
        if 'RawDataImportModule' in step['@method']:
            d['batch']['batchstep'][i]['parameter']['file'] = file_list
    return d





def tree_to_xml(t,filename=None):
    """

    """
    xml_str = ET.tostring(t)
    if filename:
        with open(filename,'w') as fid:
            fid.write(xml_str)
    return xml_str

def dict_to_etree(d):
    """
    Convert a python dictionary to an xml str
    http://stackoverflow.com/questions/7684333/converting-xml-to-dictionary-using-elementtree

    Example:
    from collections import defaultdict
    from xml.etree import cElementTree as ET

    try:
        basestring
    except NameError:  # python3
        basestring = str

    #d is a python dictionary
    ET.tostring(dict_to_etree(d))

    """
    def _to_etree(d, root):
        if not d:
            pass
        elif isinstance(d, basestring):
            root.text = d
        elif isinstance(d, dict):
            for k,v in d.items():
                assert isinstance(k, basestring)
                if k.startswith('#'):
                    assert k == '#text' and isinstance(v, basestring)
                    root.text = v
                elif k.startswith('@'):
                    assert isinstance(v, basestring)
                    root.set(k[1:], v)
                elif isinstance(v, list):
                    for e in v:
                        _to_etree(e, ET.SubElement(root, k))
                else:
                    _to_etree(v, ET.SubElement(root, k))
        else: assert d == 'invalid type', (type(d), d)
    assert isinstance(d, dict) and len(d) == 1
    tag, body = next(iter(d.items()))
    node = ET.Element(tag)
    _to_etree(body, node)
    return node

def xml_to_dict(xml_str):
    """
    Convert an xml file into a python dictionary.
    http://stackoverflow.com/questions/7684333/converting-xml-to-dictionary-using-elementtree

    Example:
    from xml.etree import cElementTree as ET
    filename = '/global/homes/b/bpb/batch_params/xmlfile.xml'
    with open(filename,'r') as fid:
        xml_str = fid.read()

    d = xml_to_dict(xml_str)
    """
    t = ET.XML(xml_str)
    d = etree_to_dict(t)
    return d

def etree_to_dict(t):
    """
    Convert an xml tree into a python dictionary.
    http://stackoverflow.com/questions/7684333/converting-xml-to-dictionary-using-elementtree

    """

    d = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.iteritems():
                dd[k].append(v)
        d = {t.tag: {k:v[0] if len(v) == 1 else v for k, v in dd.iteritems()}}
    if t.attrib:
        d[t.tag].update(('@' + k, v) for k, v in t.attrib.iteritems())
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
                d[t.tag]['#text'] = text
        else:
            d[t.tag] = text
    return d


def metatlas_formatted_atlas_from_mzmine_output(filename,polarity,make_atlas=True,atlas_name=None,
    do_store=False,min_rt=None,max_rt=None,min_mz=None,mz_tolerance=8,
    max_mz=None,remove_adducts=False,remove_fragments=False,remove_clusters=False):
    # 
    '''
    Turn mzmine output into conforming metatlas_atlas input
    
    Input:
    filename: csv file from mzmine output
    polarity: (positive,negative)
    atlas_name: string describing the atlas. useful incase you want to save it later.
    
    Output:
    atlas_df: dataframe of atlas content
    myAtlas: metatlas atlas object (if make_atlas=True)
    mzmine_df: dataframe of all mzmine info (if make_atlas=False)

    '''
    
    mzmine_df = pd.read_csv(filename)
    if min_rt:
        mzmine_df = mzmine_df[mzmine_df['row retention time']>min_rt]
    if max_rt:
        mzmine_df = mzmine_df[mzmine_df['row retention time']<max_rt]
    if min_mz:
        mzmine_df = mzmine_df[mzmine_df['row m/z']>min_mz]
    if max_mz:
        mzmine_df = mzmine_df[mzmine_Df['row m/z']<max_mz]
    if remove_adducts:
        mzmine_df = mzmine_df[~mzmine_df['row identity'].str.contains('Adduct',na=False)]
    if remove_fragments:
        mzmine_df = mzmine_df[~mzmine_df['row identity'].str.contains('Fragment',na=False)]
    if remove_clusters:
        mzmine_df = mzmine_df[~mzmine_df['row identity'].str.contains('Complex',na=False)]


    def clean_adducts(x):
        x = re.sub(r'\d+\.\d{2,}', lambda m: format(float(m.group(0)), '.2f'), x)
        new_x = ';'.join(
            [s.strip() for s in pd.unique(x.split(';'))]
            )
        return x 

    metatlas_atlas = pd.DataFrame()
    metatlas_atlas['label'] = mzmine_df.apply(lambda x: '%.4f@%.2f'%(x['row m/z'],x['row retention time']),axis=1)
    mzmine_df['row identity'].fillna('',inplace=True)
    metatlas_atlas['adduct_assignments'] = mzmine_df['row identity']#.apply(clean_adducts)
    metatlas_atlas['mz'] = mzmine_df.apply(lambda x: x['row m/z'],axis=1)
    metatlas_atlas['mz_tolerance'] = mz_tolerance#metatlas_atlas.mz.apply(lambda x: mz_tolerance*float(x)/1e6)
    metatlas_atlas['rt_peak'] = mzmine_df.apply(lambda x: x['row retention time'],axis=1)
    rt_min_cols = [col for col in mzmine_df.columns if 'Peak RT start' in col]
    metatlas_atlas['rt_min'] = mzmine_df[rt_min_cols].apply(lambda x: x.min(),axis=1)
    rt_min_cols = [col for col in mzmine_df.columns if 'Peak RT end' in col]
    metatlas_atlas['rt_max'] = mzmine_df[rt_min_cols].apply(lambda x: x.max(),axis=1)
    metatlas_atlas['inchi_key'] = None
    metatlas_atlas['detected_polarity'] = polarity
    
    #tuplize the 'Identification method' and 'Name' from adducts and fragments 

    # stick on the peak height columns
    pk_height = [col for col in list(mzmine_df) if 'Peak height' in col]
    metatlas_atlas = pd.concat([metatlas_atlas,mzmine_df[pk_height]],axis=1)
    metatlas_atlas['max_intensity'] = metatlas_atlas[pk_height].max(axis=1)
    metatlas_atlas.reset_index(inplace=True)
    metatlas_atlas.drop('index',axis=1,inplace=True)
    if make_atlas:
        if not atlas_name:
            atlas_name = filename
        myAtlas = dp.make_atlas_from_spreadsheet(metatlas_atlas,
                                           atlas_name,
                                           filetype='dataframe',
                                           sheetname='',
                                           polarity = polarity,
                                           store=do_store)
        return metatlas_atlas, myAtlas
    else:
        return metatlas_atlas, mzmine_df


def make_mzmine_scripts(mzml_files, 
                        outfile='/project/projectdirs/metatlas/projects/mzmine_parameters/pathway_1_mzmine_output.csv',
                        new_batch_file = '/project/projectdirs/metatlas/projects/mzmine_parameters/Erwinia_pyrifoliae_Pathway1_job_script_parameters.xml',
                        new_sbatch_file = '/project/projectdirs/metatlas/projects/mzmine_parameters/mod_mzmine_job.sbatch',
                        new_qsub_file = '/project/projectdirs/metatlas/projects/mzmine_parameters/mod_mzmine_job.qsub',
                        base_batch_file='/project/projectdirs/metatlas/projects/mzmine_parameters/C18_MSMS_Secondary_Metabolite_Parameters.xml',
                        base_sbatch_file = '/project/projectdirs/metatlas/projects/mzmine_parameters/mzmine_job.sbatch',
                        base_qsub_file = '/project/projectdirs/metatlas/projects/mzmine_parameters/mzmine_job.qsub',
                       mzmine_runner = '/project/projectdirs/metatlas/projects/mzmine_parameters/MZmine-2.21/startMZmine_Headless.sh'):
    '''
    This function makes three things:
        1. an xml file for using the mzmine batch queue system
        2. an sbatch file for running on Cori and Edison
        3. a qsub file for running on Genepool
        
    Inputs:
        mzml_files: list of full path to files that mzmine can parse
        outfile: name of the csv file you want mzmine to generate
        new_batch_file: name of the xml file this function will generate that tells mzmine what to do
        new_sbatch_file: sbatch filename for you to run on Cori/Edison
        new_qsub_file: qsub filename for you to run on Genepool
        base_batch_file: xml file template that has the commands for mzmine to use
        base_sbatch_file: sbatch file to modify
        base_qsub_file: qsub file to modify
        
    '''
    # 1. make a batch xml file for mzmine. modify a pre-existing one:

    with open(base_batch_file,'r') as fid:
        base_batch_text = fid.read()

    file_str = []
    for f in mzml_files:
        file_str.append('<file>%s</file>'%f)
    file_str = '\n'.join(file_str)

    #find all the files
    t = re.findall('<file>[\s\S]*?<\/file>',base_batch_text)
    #remove all but the last element
    for e in t[:-1]:
        base_batch_text = base_batch_text.replace(e+'\n','')

    base_batch_text = base_batch_text.replace(t[-1],file_str)

    old_outfile = '<parameter name="Filename">output.csv</parameter>'
    new_outfile = '<parameter name="Filename">'+outfile+'</parameter>'
    base_batch_text = base_batch_text.replace(old_outfile,new_outfile)

    with open(new_batch_file,'w') as fid:
        fid.write(base_batch_text)

    # 2. make an sbatch file.  modify a preexisting one:
    with open(base_sbatch_file,'r') as fid:
        base_sbatch_text = fid.read()
    # base_sbatch_text = base_sbatch_text.replace(base_batch_file,new_batch_file)
    base_sbatch_text = base_sbatch_text + '\n\n%s %s'%(mzmine_runner,new_batch_file)
    with open(new_sbatch_file,'w') as fid:
        fid.write(base_sbatch_text)

    # 3. make a qsub file.  modify a preexisting one:
    with open(base_qsub_file,'r') as fid:
        base_qsub_text = fid.read()
    # base_qsub_text = base_qsub_text.replace(base_batch_file,new_batch_file)
    base_qsub_text = base_qsub_text + '\n\n%s %s'%(mzmine_runner,new_batch_file)
    with open(new_qsub_file,'w') as fid:
        fid.write(base_qsub_text)
    print('qsub',new_qsub_file)

    # 3. submit job with NEWT
    # from requests import Session
    # import getpass
    # computer = 'edison'

    # if not usr:
    #     usr = getpass.getuser()
    # if not pwd:
    #     pwd = getpass.getpass("enter password for user %s: " % usr)
    # s = Session()
    # r = s.post("https://newt.nersc.gov/newt/auth", {"username": usr, "password": pwd})
    # r = s.post("https://newt.nersc.gov/newt/queue/%s/"%computer, {"jobfile": new_sbatch_file})
    # r.text

