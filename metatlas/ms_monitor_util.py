#import sys
#sys.path.insert(0,'/global/project/projectdirs/metatlas/anaconda/lib/python2.7/site-packages' )
#sys.path.append('/project/projectdirs/metatlas/projects/ms_monitor_tools')
from metatlas.helpers import metatlas_get_data_helper_fun as ma_data

from metatlas import metatlas_objects as metob
from metatlas import h5_query as h5q
from metatlas import gui as mgui
import numpy as np
import time
import os

from IPython.display import display

try:
    import ipywidgets as widgets
except ImportError:
    from IPython.html import widgets

try:
    import qgrid
    qgrid.nbinstall(overwrite=True)
    qgrid.set_grid_option('defaultColumnWidth', 200)
except Exception:
    print('Could not import QGrid')

from datetime import datetime
import pandas as pd
import json
import gspread

# this behaviour is why we have pinned the version of oauth2client
from oauth2client.client import SignedJwtAssertionCredentials
#new version uses this
#from oauth2client.service_account import ServiceAccountCredentials

from matplotlib import pyplot as plt
import re

def clean_string(oldstr):
    newstr = re.sub('[\[\]]','',oldstr)
    newstr = re.sub('[^A-Za-z0-9+-]+', '_', newstr)
    newstr = re.sub('i_[A-Za-z]+_i_', '', newstr)
    return newstr

def get_rt_mz_tolerance_from_user():
    mz_tolerance = float(raw_input('Enter mz tolerance in ppm (ex "20"): ').replace('ppm',''))
    rt_tolerance = float(raw_input('Enter the retention time uncertainty in minutes (ex "0.3"): '))
    return mz_tolerance, rt_tolerance

def get_blank_qc_pos_neg_string():
    qc_widget = widgets.Text(description="QC ID: ",value='QC')
    blank_widget = widgets.Text(description="Blank ID:",value = 'BLANK')
    pos_widget = widgets.Text(description="Neg ID: ",value='NEG')
    neg_widget = widgets.Text(description="Pos ID:",value = 'POS')
    container = widgets.VBox([widgets.HBox([qc_widget, blank_widget]),widgets.HBox([pos_widget, neg_widget])])

    display(container)
    return qc_widget, blank_widget, pos_widget, neg_widget

def get_files_for_experiment(experiment_name):
    files = metob.retrieve('LcmsRun',username='*',experiment=experiment_name)
    flist = []
    for f in files:
        flist.append(f.hdf5_file)    
    flist = np.unique(flist)
    df = pd.DataFrame()
    for counter,f in enumerate(flist):
        df.loc[counter,'file'] = os.path.basename(f)
#    del df['index']   
    df.set_index('file', drop=True, append=False, inplace=True)
    #df.reset_index(drop=True,inplace=True)
    
    options = qgrid.grid.defaults.grid_options
    options['defaultColumnWidth'] = 600
    #mygrid = qgrid.show_grid(df, remote_js=True,grid_options = options)
    grid = qgrid.grid.QGridWidget(df=df,
                                  precision=6,
                                  grid_options=options,
                                  remote_js=True)

    def handle_msg(widget, content, buffers=None):
        if content['type'] == 'cell_change':
            obj = objects[content['row']]
            try:
                setattr(obj, content['column'], content['value'])
            except Exception:
                pass

    grid.on_msg(handle_msg)
    gui = widgets.Box([grid])
    display(gui)  
    return files

def get_recent_experiments(num_days):
    if not num_days:
        num_days = 5
    query = 'SELECT DISTINCT experiment,creation_time FROM lcmsruns where creation_time >= UNIX_TIMESTAMP(DATE_SUB(CURDATE(), INTERVAL %d DAY))'%num_days
    entries = [e for e in metob.database.query(query)]
    df = pd.DataFrame() 
    counter = 0
    experiments = []
    for entry in entries:
        if entry['experiment']:
            experiments.append( entry['experiment'] )
    experiments = np.unique(experiments)
    experiment_widget = widgets.Dropdown(
        options=list(experiments),
        value=experiments[0],
        description='Experiments: '
    )
    display(experiment_widget)
    #files = get_files_for_experiment(experiment_widget.value)
    #def experiment_change(trait,value):
    #    files = get_files_for_experiment(value)
    #    return files
    #experiment_widget.on_trait_change(experiment_change,'value')

    return experiment_widget

def get_files_from_recent_experiment(num_days):
    if not num_days:
        num_days = 5
    query = 'SELECT DISTINCT experiment,creation_time,username FROM lcmsruns where creation_time >= UNIX_TIMESTAMP(DATE_SUB(CURDATE(), INTERVAL %d DAY))'%num_days
    entries = [e for e in metob.database.query(query)]
    df = pd.DataFrame() 
    counter = 0
    for entry in entries:
        if entry['experiment']:
            df.loc[counter,'experiment'] = entry['experiment']
            df.loc[counter,'username'] = entry['username']
            df.loc[counter, 'utc time'] = datetime.utcfromtimestamp(entry['creation_time'])
            counter = counter + 1
    #TODO: filter by unique experiment name
    #df.drop_duplicates(cols = 'experiment', inplace = True)
    df.groupby('experiment', group_keys=False)
    options = qgrid.grid.defaults.grid_options
    grid = qgrid.grid.QGridWidget(df=df,
                                  precision=6,
                                  grid_options=options,
                                  remote_js=True)

    def handle_msg(widget, content, buffers=None):
        if content['type'] == 'cell_change':
            obj = objects[content['row']]
            try:
                setattr(obj, content['column'], content['value'])
            except Exception:
                pass

    grid.on_msg(handle_msg)
    return grid    
    #mygrid = qgrid.show_grid(df, remote_js=True,)
    #print "Enter the experiment name here"
    #my_experiment = raw_input()
    #files =  get_files_for_experiment(my_experiment)
    #files = qgrid.get_selected_rows(mygrid)    
    #return files

def get_method_dropdown():
    methods = ['Not Set',
               '6550_RP_0.5min_25ppm_1000counts',
 '6550_ZIC-HILIC_0.5min_25ppm_1000counts',
 '6550_pHILIC-noguard_0.5min_25ppm_1000counts',
 '6550_pHILIC_0.5min_25ppm_1000counts',
 'QE119_RP_0.5min_25ppm_1000counts',
 'QE119_ZIC-HILIC_0.5min_25ppm_1000counts',
 'QE139_RP_0.5min_25ppm_1000counts',
 'QE139_ZIC-HILIC_0.5min_25ppm_1000counts',
 'QE139_pHILIC_0.5min_25ppm_1000counts',
 'QE144_RP_0.5min_25ppm_1000counts',
 'QE144_ZIC-HILIC_0.5min_25ppm_1000counts']
    method_widget = widgets.Dropdown(
        options=methods,
        value=methods[0],
        description='LC-MS Method:'
    )
    display(method_widget)

#    methods_widget.on_trait_change(filter_istd_qc_by_method,'value')

    return method_widget

def get_ms_monitor_reference_data(notebook_name = "20160203 ms-monitor reference data", token='/project/projectdirs/metatlas/projects/google_sheets_auth/ipython to sheets demo-9140f8697062.json', sheet_name = 'QC and ISTD'):
    """
    Returns a pandas data frame from the google sheet containing the reference data.
    Feeds the 
    """
    json_key = json.load(open(token))
    scope = ['https://spreadsheets.google.com/feeds']

    #this is deprecated as of january, but we have pinned the version of oauth2.
    #see https://github.com/google/oauth2client/issues/401
    credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'].encode(), scope)
    
    #here is the new way incase the version pin is removed
    #credentials = ServiceAccountCredentials.from_json_keyfile_name(token, scope)
    
    gc = gspread.authorize(credentials)

    wks = gc.open(notebook_name)
    istd_qc_data = wks.worksheet(sheet_name).get_all_values()
#     blank_data = wks.worksheet('BLANK').get_all_values()
    headers = istd_qc_data.pop(0)
    df = pd.DataFrame(istd_qc_data,columns=headers)
    print('keys',df.keys())
    print('shape',df.shape)
    df = df[(df['mz_POS'] != '') | (df['mz_NEG'] != '')]

    return df#, blank_data

def filter_istd_qc_by_method(method):
    rt_minutes_tolerance = float(method.split('_')[2].replace('min',''))
    mz_ppm_tolerance = float(method.split('_')[3].replace('ppm',''))
    peak_height_minimum = float(method.split('_')[4].replace('counts',''))

    reference_data = make_dict_of_vals(method,get_ms_monitor_reference_data(),rt_minutes_tolerance)

    reference_data['parameters'] = {}
    reference_data['parameters']['rt_minutes_tolerance'] = rt_minutes_tolerance
    reference_data['parameters']['mz_ppm_tolerance'] = mz_ppm_tolerance
    reference_data['parameters']['peak_height_minimum'] = peak_height_minimum

    return reference_data

def convert_float(val):
    try:
        return float(val)
    except ValueError:
        return np.nan

def make_dict_of_vals(method, df,rt_minutes_tolerance,my_fields = ['COMMON-HILIC', u'ISTD-HILIC', u'QC-HILIC'],base_keys = ['label','inchi_key','mz_POS', 'mz_NEG']):
    pat = '_'.join(method.split('_')[:2]).lower()
    method_specific_cols = [col for col in df.columns if col.lower().endswith(pat)]
    pos_method_specific_cols = [col for col in df.columns if col.lower().endswith(pat+'_pos')]
    neg_method_specific_cols = [col for col in df.columns if col.lower().endswith(pat+'_neg')]
    my_dict = {}
    float_fields = ['rt_min','rt_max','rt_peak','pos_mz','neg_mz']#,'peak-height_pos','peak-height_neg']
    for field in my_fields:
        renamed_field = field.split('-')[0].lower()
        #The first field tells you if you are an ISTD or QC
        #The second field is the join of strings to equal column headings
        my_dict[renamed_field] =  df[df[field] == '1'][base_keys + method_specific_cols + pos_method_specific_cols + neg_method_specific_cols]
        my_dict[renamed_field].rename(columns={'mz_POS':'pos_mz','mz_NEG':'neg_mz'},inplace=True)
        my_dict[renamed_field].rename(columns=lambda x: x.lower(), inplace=True)
        my_dict[renamed_field].rename(columns=lambda x: x.replace('_'+pat,''), inplace=True)
        my_dict[renamed_field].rename(columns=lambda x: x.replace('_'+pat+'_pos',''), inplace=True)
        my_dict[renamed_field].rename(columns=lambda x: x.replace('_'+pat+'_neg',''), inplace=True)
        #my_dict[renamed_field] = df[renamed_field].apply(lambda x: pd.to_numeric(x, errors = 'ignore'))
        my_dict[renamed_field][float_fields].replace(r'\s+', '0', regex=True,inplace=True)
        for ff in float_fields:
            my_dict[renamed_field][ff] = my_dict[renamed_field][ff].apply(convert_float)
        my_dict[renamed_field].loc[:,'rt_min'] -= rt_minutes_tolerance
        my_dict[renamed_field].loc[:,'rt_max'] += rt_minutes_tolerance
        my_dict[renamed_field].drop([col for col in my_dict[renamed_field].columns if col.lower().startswith('file')], axis=1, inplace=True)

    return my_dict

def construct_result_table_for_files(files,qc_str,blank_str,neg_str,pos_str,method,reference_data,experiment):
    df = pd.DataFrame()
    counter = 0
    for my_file in files:
        #determine if the file is a blank, qc or istd:
        filetype = 'istd'
        if blank_str.value.upper() in my_file.name.upper():
            filetype = 'istd'
        elif (qc_str.value.upper() in my_file.name.upper()):
            filetype = 'qc'

        is_pos = pos_str.value.upper() in my_file.name.upper()

        atlas = reference_data[filetype]
        for cidx in atlas.index:
            rt_ref = metob.RtReference()
            rt_ref.rt_peak = atlas.loc[cidx,'rt_peak']
            rt_ref.rt_min = atlas.loc[cidx,'rt_min']
            rt_ref.rt_max = atlas.loc[cidx,'rt_max']

            mz_ref = metob.MzReference()
            mz_ref.mz_tolerance = reference_data['parameters']['mz_ppm_tolerance']
            mz_ref.mz_tolerance_units = 'ppm'

            if is_pos:
                mz_ref.mz = atlas.loc[cidx,'pos_mz']
                mz_ref.detected_polarity = 'positive'
                ref_intensity = 1000.0 #float(atlas.loc[cidx,'peak-height_pos'])
            else:
                mz_ref.mz = atlas.loc[cidx,'neg_mz']
                mz_ref.detected_polarity = 'negative'
                ref_intensity = 1000.0 #float(atlas.loc[cidx,'peak-height_neg'])
            if not np.isfinite(ref_intensity):
                ref_intensity = 1
                
            result = ma_data.get_data_for_a_compound(mz_ref,
                                    rt_ref,[ 'ms1_summary' ],
                                    my_file.hdf5_file,0.3) #extra_time is not used by ms1_summary
            if result['ms1_summary']['rt_peak']:
                if result['ms1_summary']['peak_height'] > reference_data['parameters']['peak_height_minimum']:
                    df.loc[counter,'expected_rt'] = 1
                    df.loc[counter,'expected_rt'] = rt_ref.rt_peak
                    df.loc[counter,'expected_mz'] = mz_ref.mz
                    df.loc[counter,'expected_intensity'] = ref_intensity
                    df.loc[counter,'delta_rt'] = result['ms1_summary']['rt_peak'] - rt_ref.rt_peak
                    df.loc[counter,'delta_mz'] = (result['ms1_summary']['mz_centroid'] - mz_ref.mz)/mz_ref.mz*1e6
                    df.loc[counter,'delta_intensity'] = (result['ms1_summary']['peak_height'] - ref_intensity) / ref_intensity
                    df.loc[counter,'measured_rt'] = result['ms1_summary']['rt_peak']
                    df.loc[counter,'measured_mz'] = result['ms1_summary']['mz_centroid']
                    df.loc[counter,'measured_intensity'] = result['ms1_summary']['peak_height']
                    df.loc[counter,'filetype'] = filetype
                    df.loc[counter, 'name has blank'] = blank_str.value.upper() in my_file.name.upper()
                    df.loc[counter, 'name has QC'] = qc_str.value.upper() in my_file.name.upper()
                    df.loc[counter, 'name has pos'] = pos_str.value.upper() in my_file.name.upper()
                    df.loc[counter, 'name has neg'] = neg_str.value.upper() in my_file.name.upper()
                    df.loc[counter, 'experiment'] = my_file.experiment
                    df.loc[counter, 'filename'] = my_file.name
                    df.loc[counter, 'datestamp'] = my_file.creation_time
                    df.loc[counter, 'utc time'] = datetime.utcfromtimestamp(my_file.creation_time)
                    df.loc[counter, 'lcms method'] = my_file.method #TODO: get instrument and lcms from the method object
                    df.loc[counter, 'sample'] = my_file.sample
                    df.loc[counter, 'username'] = my_file.username
                    df.loc[counter, 'method'] = method.value
                    df.loc[counter,'Compound name'] = atlas.loc[cidx,'label']
                    counter = counter + 1
    timestr = time.strftime("%Y%m%d-%H%M%S")
    df.to_excel('%s_%s_%s.xls'%( timestr, clean_string(experiment.value), clean_string(method.value) ) )
    df.to_excel('%s/%s_%s_%s.xls'%( '/project/projectdirs/metatlas/projects/ms_monitor_tools/ms_monitor_logs', timestr,clean_string(experiment.value), clean_string(method.value) ) )

    f = make_compound_plots(df,'QC',pos_str.value,experiment,method)
    f = make_compound_plots(df,'QC',neg_str.value,experiment,method)
    f = make_compound_plots(df,'ISTD',pos_str.value,experiment,method)
    f = make_compound_plots(df,'ISTD',neg_str.value,experiment,method)

    return df


def make_compound_plots(df,plot_type,polarity,experiment,method):
    if plot_type == 'QC':
        compound_names = df[df['name has QC'] == 1]['Compound name'].unique()
    else:
        compound_names = df[df['name has QC'] == 0]['Compound name'].unique()
    counter = 1
    nRows = len(compound_names)
    nCols = 3
    if nRows > 0:
        f, ax = plt.subplots(nRows, nCols, figsize=(5*nCols,nRows * 2)) #original 8 x 6
        for i,cname in enumerate(compound_names):
        #     sdf = df[df['Compound name'].str.contains(cname, regex=False) & df['filename'].str.contains("POS") & (df['name has blank'] == 0)  & (df['name has QC'] == 1)]
            if plot_type == 'QC':        
                sdf = df[(df['Compound name'] == cname) & df['filename'].str.contains(polarity, case=False, regex=False) & (df['name has blank'] == 0)  & (df['name has QC'] == 1)]
            else:
                sdf = df[(df['Compound name'] == cname) & df['filename'].str.contains(polarity, case=False, regex=False) & (df['name has blank'] == 0)  & (df['name has QC'] == 0)]

            #     ax[i,0].scatter(sdf['measured_rt'].tolist(), sdf['measured_mz'].tolist())#,c=sdf['delta_intensity'])
            ax[i,0].plot(sdf['measured_rt'].tolist(),'.-')#,c=sdf['delta_intensity'])
            ev = sdf['expected_rt'].tolist()
            if not ev:
                ev = [0,0,0]
            l = ax[i,0].axhline(float(ev[0]),color='black',linewidth=2)
            p = ax[i,0].axhspan(float(ev[0]) - 0.5, float(ev[0]) + 0.5, facecolor='0.5', alpha=0.5)
            
            ax[i,1].plot(sdf['measured_mz'].tolist(),'.-')#,c=sdf['delta_intensity'])
            ev = sdf['expected_mz'].tolist()
            if not ev:
                ev = [0,0,0]
            l = ax[i,1].axhline(float(ev[0]),color='black',linewidth=2)
            p = ax[i,1].axhspan(float(ev[0]) - float(ev[0])*5/1e6, float(ev[0]) + float(ev[0])*5/1e6, facecolor='0.5', alpha=0.5)

            
            ax[i,2].plot(sdf['measured_intensity'].tolist(),'.-')#,c=sdf['delta_intensity'])
            ev = sdf['expected_intensity'].tolist()
            if not ev:
                ev = [0,0,0]
            l = ax[i,2].axhline(float(ev[0]),color='black',linewidth=2)
            p = ax[i,2].axhspan(float(ev[0]) - float(ev[0])*0.25, float(ev[0]) + float(ev[0])*0.25, facecolor='0.5', alpha=0.5)

            # plot measured m/z and rt.  Show a line with shading to indicate expected values.  Mouseover to print name of run.
            # run order plots with trendline.  for now use upload time as a proxy for run order.
            ax[i,0].set_title(cname)
            ax[i,1].set_title(cname)
            ax[i,2].set_title(cname)
            
            ax[i,0].get_yaxis().get_major_formatter().set_useOffset(False)
            ax[i,1].get_yaxis().get_major_formatter().set_useOffset(False)
            ax[i,2].get_yaxis().get_major_formatter().set_useOffset(False)
            
            ax[i,0].set_ylabel('RT (min)')
            ax[i,1].set_ylabel('mz')
            ax[i,2].set_ylabel('peak height')
            
        #     ax[i,0].set_ylabel('mz (ppm)')
        #     ax[i,0].set_title(cname)
        #     ax[i,0].get_xaxis().get_major_formatter().set_useOffset(False)
        #     ax[i,0].get_yaxis().get_major_formatter().set_useOffset(False)
        # http://matplotlib.org/examples/pylab_examples/scatter_hist.html
        plt.tight_layout()
#         plt.show()
        timestr = time.strftime("%Y%m%d-%H%M%S")
        f.savefig('plot_summary_%s_%s_%s_%s_%s.png'%( timestr, plot_type, polarity, clean_string(experiment.value), clean_string(method.value) ) )
        f.savefig('/project/projectdirs/metatlas/projects/ms_monitor_tools/ms_monitor_logs/plot_summary_%s_%s_%s_%s_%s.png'%( timestr, plot_type, polarity, clean_string(experiment.value), clean_string(method.value) ) )
        plt.close(f)

# it takes too long to write to a sheet this way.  need to redo it with getting all cells, updating their values and then sending the data over as a large transfer
#import json
#import gspread
#from oauth2client.client import SignedJwtAssertionCredentials
## def append_result_to_google_sheet(df):
#json_key = json.load(open('/project/projectdirs/metatlas/projects/google_sheets_auth/ipython to sheets demo-9140f8697062.json'))
#scope = ['https://spreadsheets.google.com/feeds']
#credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'].encode(), scope)
#gc = gspread.authorize(credentials)
#wks = gc.open("lcms_run_log")
#sheet_data = wks.worksheet('active').get_all_values()
##     blank_data = wks.worksheet('BLANK').get_all_values()
#print sheet_data
#keys = df.keys()
#for index,row in df.iterrows():
#    vals = []
#    for i,k in enumerate(keys):
#        vals.append(row[k])
#    wks.worksheet('active').insert_row(vals, index=1)
