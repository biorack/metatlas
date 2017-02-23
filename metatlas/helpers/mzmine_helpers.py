
import pandas as pd

def metatlas_formatted_atlas_from_mzmine_output(filename,polarity,make_atlas=True,atlas_name='',do_store=False):
    '''
    Turn mzmine output into conforming metatlas_atlas input
    
    Input:
    filename: csv file from mzmine output
    polarity: (positive,negative)
    atlas_name: string describing the atlas. useful incase you want to save it later.
    
    Output:
    atlas_df: dataframe of atlas content
    myAtlas: metatlas atlas object

    '''
    
    mzmine_df = pd.read_csv(filename)
    metatlas_atlas = pd.DataFrame()
    metatlas_atlas['label'] = mzmine_df.apply(lambda x: '%.4f@%.2f'%(x['row m/z'],x['row retention time']),axis=1)
    metatlas_atlas['mz'] = mzmine_df.apply(lambda x: x['row m/z'],axis=1)
    metatlas_atlas['mz_tolerance'] = metatlas_atlas.mz.apply(lambda x: 0.015/float(x)*1e6)
    metatlas_atlas['rt_peak'] = mzmine_df.apply(lambda x: x['row retention time'],axis=1)
    rt_min_cols = [col for col in mzmine_df.columns if 'Peak RT start' in col]
    metatlas_atlas['rt_min'] = mzmine_df[rt_min_cols].apply(lambda x: x.min(),axis=1)
    rt_min_cols = [col for col in mzmine_df.columns if 'Peak RT end' in col]
    metatlas_atlas['rt_max'] = mzmine_df[rt_min_cols].apply(lambda x: x.max(),axis=1)
    metatlas_atlas['inchi_key'] = None
    metatlas_atlas['detected_polarity'] = polarity
    
    # stick on the peak height columns
    pk_height = [col for col in list(mzmine_df) if 'Peak height' in col]
    metatlas_atlas = pd.concat([metatlas_atlas,mzmine_df[pk_height]],axis=1)
    metatlas_atlas['max_intensity'] = metatlas_atlas[pk_height].max(axis=1)
    metatlas_atlas.reset_index(inplace=True)
    metatlas_atlas.drop('index',axis=1,inplace=True)
    if make_atlas:
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
                        base_qsub_file = '/project/projectdirs/metatlas/projects/mzmine_parameters/mzmine_job.qsub'):
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
    import re

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
    base_sbatch_text = base_sbatch_text.replace(base_batch_file,new_batch_file)
    with open(new_sbatch_file,'w') as fid:
        fid.write(base_sbatch_text)

    # 3. make a qsub file.  modify a preexisting one:
    with open(base_qsub_file,'r') as fid:
        base_qsub_text = fid.read()
    base_qsub_text = base_qsub_text.replace(base_batch_file,new_batch_file)
    with open(new_qsub_file,'w') as fid:
        fid.write(base_qsub_text)
    print 'qsub',new_qsub_file

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

