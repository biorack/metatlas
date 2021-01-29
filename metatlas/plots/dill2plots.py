from __future__ import absolute_import
from __future__ import print_function
import sys
import os
import os.path
import multiprocessing as mp
# os.environ['R_LIBS_USER'] = '/project/projectdirs/metatlas/r_pkgs/'
#curr_ld_lib_path = ''


from metatlas.datastructures import metatlas_objects as metob
from metatlas.io import h5_query as h5q
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.tools import spectralprocessing as sp
from metatlas.plots import chromplotplus as cpp
# from metatlas import gui

from textwrap import fill, TextWrapper
# import qgrid
import pandas as pd
import os
import tables
import pickle
import dill
import numpy as np
import re
import json
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem, Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D, IPythonConsole
from itertools import cycle
from collections import defaultdict
from IPython.display import SVG,display,clear_output


from ipywidgets import interact, interactive, fixed, IntProgress
import ipywidgets as widgets
from IPython.display import display

import getpass

from ast import literal_eval
# from datetime import datetime

from matplotlib.widgets import Slider, Button, RadioButtons

from matplotlib.widgets import AxesWidget

import gspread
# from oauth2client.client import SignedJwtAssertionCredentials
from oauth2client.service_account import ServiceAccountCredentials

import sys
import six
from six.moves import range
from six.moves import zip
from functools import reduce
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO

ADDUCT_INFO = {'[2M+H]': {'charge': '1',
              'color': '#fabebe',
              'common': True,
              'comp_num': '2',
              'mass': '1.0073'},
             '[2M-H]': {'charge': '-1',
              'color': '#008080',
              'common': True,
              'comp_num': '2',
              'mass': '-1.0073'},
             '[M+2H]': {'charge': '2',
              'color': '#ffe119',
              'common': True,
              'comp_num': '1',
              'mass': '2.0146'},
             '[M+2Na]': {'charge': '2',
              'color': '#fffac8',
              'common': False,
              'comp_num': '1',
              'mass': '45.9784'},
             '[M+Cl]': {'charge': '-1',
              'color': '#d2f53c',
              'common': True,
              'comp_num': '1',
              'mass': '34.9694'},
             '[M+H-H2O]': {'charge': '1',
              'color': '#911eb4',
              'common': True,
              'comp_num': '1',
              'mass': '-17.0033'},
             '[M+H]': {'charge': '1',
              'color': '#3cb44b',
              'common': True,
              'comp_num': '1',
              'mass': '1.0073'},
             '[M+K]': {'charge': '1',
              'color': '#aa6e28',
              'common': False,
              'comp_num': '1',
              'mass': '38.963158'},
             '[M+NH4]': {'charge': '1',
              'color': '#0082c8',
              'common': True,
              'comp_num': '1',
              'mass': '18.0338'},
             '[M+Na]': {'charge': '1',
              'color': '#f58231',
              'common': True,
              'comp_num': '1',
              'mass': '22.9892'},
             '[M+acetate]': {'charge': '-1',
              'color': '#808000',
              'common': False,
              'comp_num': '1',
              'mass': '59.0139'},
             '[M+formate]': {'charge': '-1',
              'color': '#5500FF',
              'common': False,
              'comp_num': '1',
              'mass': '44.998201'},
             '[M-2H]': {'charge': '-2',
              'color': '#f032e6',
              'common': True,
              'comp_num': '1',
              'mass': '-2.014552904'},
             '[M-H+2Na]': {'charge': '1',
              'color': '#000080',
              'common': False,
              'comp_num': '1',
              'mass': '44.9711'},
             '[M-H+Cl]': {'charge': '-2',
              'color': '#ffd8b1',
              'common': False,
              'comp_num': '1',
              'mass': '33.9621'},
             '[M-H+Na]': {'charge': '0',
              'color': '#e6beff',
              'common': False,
              'comp_num': '1',
              'mass': '21.98194425'},
             '[M-H]': {'charge': '-1',
              'color': '#46f0f0',
              'common': True,
              'comp_num': '1',
              'mass': '-1.0073'},
             '[M-e]': {'charge': '1',
              'color': '#aaffc3',
              'common': False,
              'comp_num': '1',
              'mass': '-0.0005'},
             '[M]': {'charge': '0',
              'color': '#e6194b',
              'common': True,
              'comp_num': '1',
              'mass': '0'}}

def get_google_sheet(notebook_name = "Sheet name",
                     token='/project/projectdirs/metatlas/projects/google_sheets_auth/ipython to sheets demo-9140f8697062.json',
                     sheet_name = 'Sheet1',
                    literal_cols=None):
    """
    Returns a pandas data frame from the google sheet.
    Assumes header row is first row.

    To use the token hard coded in the token field,
    the sheet must be shared with:
    metatlas-ipython-nersc@ipython-to-sheets-demo.iam.gserviceaccount.com
    Unique sheet names are a requirement of this approach.

    """
#     scope = ['https://spreadsheets.google.com/feeds']
#     scope = ['https://www.googleapis.com/auth/spreadsheets']
    scope = ['https://spreadsheets.google.com/feeds', 'https://www.googleapis.com/auth/drive']
    #this is deprecated as of january, but we have pinned the version of oauth2.
    #see https://github.com/google/oauth2client/issues/401
#     json_key = json.load(open(token))
#     credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'].encode(), scope)
    credentials = ServiceAccountCredentials.from_json_keyfile_name(token, scope)
    #here is the new way incase the version pin is removed
    #credentials = ServiceAccountCredentials.from_json_keyfile_name(token, scope)

    gc = gspread.authorize(credentials)
    wks = gc.open(notebook_name)
    istd_qc_data = wks.worksheet(sheet_name).get_all_values()
    headers = istd_qc_data.pop(0)
    df = pd.DataFrame(istd_qc_data,columns=headers)

    # Use round trip through read_csv to infer dtypes
    s = StringIO()
    df.to_csv(s)
    df2 = pd.read_csv(StringIO(s.getvalue()))
    if 'Unnamed: 0' in df2.columns:
        df2.drop(columns=['Unnamed: 0'],inplace=True)

    #turn list elements into lists instead of strings
    if literal_cols is not None:
        for col in literal_cols:
            df2[col] = df2[col].apply(literal_eval)
    df2 = df2.fillna('')

    return df2

class VertSlider(AxesWidget):
    """
    A slider representing a floating point range

    The following attributes are defined
      *ax*        : the slider :class:`matplotlib.axes.Axes` instance

      *val*       : the current slider value

      *vline*     : a :class:`matplotlib.lines.Line2D` instance
                     representing the initial value of the slider

      *poly*      : A :class:`matplotlib.patches.Polygon` instance
                     which is the slider knob

      *valfmt*    : the format string for formatting the slider text

      *label*     : a :class:`matplotlib.text.Text` instance
                     for the slider label

      *closedmin* : whether the slider is closed on the minimum

      *closedmax* : whether the slider is closed on the maximum

      *slidermin* : another slider - if not *None*, this slider must be
                     greater than *slidermin*

      *slidermax* : another slider - if not *None*, this slider must be
                     less than *slidermax*

      *dragging*  : allow for mouse dragging on slider

    Call :meth:`on_changed` to connect to the slider event
    """
    def __init__(self, ax, label, valmin, valmax, valinit=0.5, valfmt='%.1e',
                 closedmin=True, closedmax=True, slidermin=None,
                 slidermax=None, dragging=True, **kwargs):
        """
        Create a slider from *valmin* to *valmax* in axes *ax*

        *valinit*
            The slider initial position

        *label*
            The slider label

        *valfmt*
            Used to format the slider value

        *closedmin* and *closedmax*
            Indicate whether the slider interval is closed

        *slidermin* and *slidermax*
            Used to constrain the value of this slider to the values
            of other sliders.

        additional kwargs are passed on to ``self.poly`` which is the
        :class:`matplotlib.patches.Rectangle` which draws the slider
        knob.  See the :class:`matplotlib.patches.Rectangle` documentation
        valid property names (e.g., *facecolor*, *edgecolor*, *alpha*, ...)
        """
        AxesWidget.__init__(self, ax)

        self.valmin = valmin
        self.valmax = valmax
        self.val = valinit
        self.valinit = valinit
        self.poly = ax.axhspan(valmin, valinit, 0, 1, **kwargs)

        self.vline = ax.axhline(valinit, 0, 1, color='r', lw=1)

        self.valfmt = valfmt
        ax.set_xticks([])
        ax.set_ylim((valmin, valmax))
        ax.set_yticks([])
        ax.set_navigate(False)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        if dragging:
            self.connect_event('motion_notify_event', self._update)
        self.label = ax.text(0.5, 1.03, label, transform=ax.transAxes,
                             verticalalignment='center',
                             horizontalalignment='center')

        self.valtext = ax.text(0.5, -0.03, valfmt % valinit,
                               transform=ax.transAxes,
                               verticalalignment='center',
                               horizontalalignment='center')

        self.cnt = 0
        self.observers = {}

        self.closedmin = closedmin
        self.closedmax = closedmax
        self.slidermin = slidermin
        self.slidermax = slidermax
        self.drag_active = False

    def _update(self, event):
        """update the slider position"""
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes == self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif ((event.name == 'button_release_event') or
              (event.name == 'button_press_event' and
               event.inaxes != self.ax)):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return

        val = event.ydata
        if val <= self.valmin:
            if not self.closedmin:
                return
            val = self.valmin
        elif val >= self.valmax:
            if not self.closedmax:
                return
            val = self.valmax

        if self.slidermin is not None and val <= self.slidermin.val:
            if not self.closedmin:
                return
            val = self.slidermin.val

        if self.slidermax is not None and val >= self.slidermax.val:
            if not self.closedmax:
                return
            val = self.slidermax.val

        self.set_val(val)

    def set_val(self, val):
        xy = self.poly.xy
        xy[1] = 0, val
        xy[2] = 1, val
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % val)
        if self.drawon:
            self.ax.figure.canvas.draw()
        self.val = val
        if not self.eventson:
            return
        for cid, func in six.iteritems(self.observers):
            func(val)

    def on_changed(self, func):
        """
        When the slider value is changed, call *func* with the new
        slider position

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        try:
            del self.observers[cid]
        except KeyError:
            pass

    def reset(self):
        """reset the slider to the initial value if needed"""
        if (self.val != self.valinit):
            self.set_val(self.valinit)


class adjust_rt_for_selected_compound(object):
    def __init__(self,
                 data,
                 include_lcmsruns = None,
                 exclude_lcmsruns = None,
                 include_groups = None,
                 exclude_groups = None,
                 msms_hits = None,
                 color_me = '',
                 compound_idx = 0,
                 width = 10,
                 height = 4,
                 y_scale='linear',
                 alpha = 0.5,
                 min_max_color = 'green',
                 peak_color = 'darkviolet',
                 slider_color = 'ghostwhite',
                 y_max = 'auto',
                 y_min = 0,
                 peak_flags = ('keep', 'remove', 'unresolvable isomers','poor peak shape'),
                 msms_flags = ('no selection', '-1, bad match - should remove compound', '0, no ref match available or no MSMS collected', '0.5, co-isolated precursor, partial match', '0.5, partial match of fragments', '1, perfect match to internal reference library', '1, perfect match to external reference library', '1 co-isolated precursor but all reference ions are in sample spectrum'),
                 adjustable_rt_peak = False):
        """
        data: a metatlas_dataset where files and compounds are stored.
        for example,
        self.metatlas_dataset[file_idx][compound_idx]['identification'].rt_references[-1].unique_id
        is the unique id to the retention time reference for a compound in a file.
        msms_hits: output of get_msms_hits routine
        width: specify a width value in inches for the plots and slides
        height: specify a width value in inches for the plots and slides
        min_max_color & peak_color: specify a valid matplotlib color string for the slider and vertical bars
        slider_color: background color for the sliders. Must be a valid matplotlib color

        Press Left and Right arrow keys to move to the next or previous compound
        """

        self.msms_hits = msms_hits
        self.color_me = color_me
        self.compound_idx = compound_idx
        self.width = width
        self.height = height
        self.y_scale = y_scale
        self.alpha = alpha
        self.min_max_color = min_max_color
        self.peak_color = peak_color
        self.slider_color = slider_color
        self.y_max = y_max
        self.y_min = y_min
        self.peak_flags = peak_flags
        self.msms_flags = msms_flags
        self.adjustable_rt_peak = adjustable_rt_peak

        # filter runs from the metatlas dataset
        if include_lcmsruns:
            data = filter_lcmsruns_in_dataset_by_include_list(data,'lcmsrun',include_lcmsruns)

        if include_groups:
            data = filter_lcmsruns_in_dataset_by_include_list(data,'group',include_groups)
        if exclude_lcmsruns:
            data = filter_lcmsruns_in_dataset_by_exclude_list(data,'lcmsrun',exclude_lcmsruns)
        if exclude_groups:
            data = filter_lcmsruns_in_dataset_by_exclude_list(data,'group',exclude_groups)
        self.data = data

        if self.peak_flags == '':
            self.peak_flags = ('keep', 'remove', 'unresolvable isomers','poor peak shape')
        if self.msms_flags == '':
            self.msms_flags = ('no selection','-1, bad match - should remove compound', '0, no ref match available or no MSMS collected', '0.5, co-isolated precursor, partial match', '0.5, partial match of fragments', '1, perfect match to internal reference library', '1, perfect match to external reference library', '1 co-isolated precursor but all reference ions are in sample spectrum')

        #Turn On interactive plot
        plt.ion()
        # create figure and first axes
        self.fig,(self.ax2, self.ax) = plt.subplots(2, 1, figsize=(width, height*2.2))
        plt.subplots_adjust(left=0.09, bottom=0.275, hspace=0.4)
#         plt.ticklabel_format(style='plain', axis='x')
#         plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        # warn the user if they do not own the atlas; and can not edit its values
        self.enable_edit = True
        self.atlas = metob.retrieve('Atlas',unique_id = self.data[0][0]['atlas_unique_id'],username='*')[-1]
        print(("loaded file for username = ", self.atlas.username))
        if getpass.getuser() != self.atlas.username:
            self.ax.set_title("YOUR ARE %s YOU ARE NOT ALLOWED TO EDIT VALUES THE RT CORRECTOR. USERNAMES ARE NOT THE SAME"%getpass.getuser())
            self.enable_edit = False

        #create all event handlers
        self.fig.canvas.callbacks.connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect('key_press_event', self.press)

        #create the plot
        self.hit_ctr = 0
        self.set_plot_data()


    def set_plot_data(self):
        #set y-scale and bounds if provided
        file_names = ma_data.get_file_names(self.data)
        self.ax.set_yscale(self.y_scale)
        if self.y_max != 'auto':
            self.ax.set_ylim(self.y_min,self.y_max)

        self.ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        self.ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        default_data = self.data[0][self.compound_idx]
        if default_data['identification'].name:
            compound_str = default_data['identification'].name.split('///')[0]
        elif default_data['identification'].compound[-1].name:
            compound_str = default_data['identification'].compound[-1].name
        else:
            compound_str = 'nameless compound'

        cname = compound_str
        if len(default_data['identification'].compound) > 0 and hasattr(default_data['identification'].compound[-1],"inchi_key"):
            inchi_key = default_data['identification'].compound[-1].inchi_key
        else:
            inchi_key = ""
        compound_str = '%d, %s'%(self.compound_idx, compound_str)

        try:
            adduct = default_data['identification'].mz_references[0].adduct
        except (KeyError, AttributeError):
            adduct = None

        mz_theoretical = default_data['identification'].mz_references[0].mz
        mz_measured = default_data['data']['ms1_summary']['mz_centroid']
        if not mz_measured:
            mz_measured = 0
        
        delta_mz = abs(mz_theoretical - mz_measured)
        delta_ppm = delta_mz / mz_theoretical * 1e6

        mz_header = "m/z theoretical = %5.4f, m/z measured = %5.4f, ppm diff = %5.4f" % (mz_theoretical, mz_measured, delta_ppm)

        self.ax.set_title('')
        if adduct != None:
            self.ax.set_ylabel('%s\n%s'%(compound_str,adduct))
        else:
            self.ax.set_ylabel('%s'%compound_str)

        self.ax.set_xlabel('Retention Time')
        self.my_rt = metob.retrieve('RTReference',
                               unique_id = default_data['identification'].rt_references[-1].unique_id, username='*')[-1]
        for d in self.data: #this loops through the files
            if d[self.compound_idx]['data']['eic']:
                if len(d[self.compound_idx]['data']['eic']['rt']) > 0:
                    x = d[self.compound_idx]['data']['eic']['rt']
                    y = d[self.compound_idx]['data']['eic']['intensity']
                    x = np.asarray(x)
                    y = np.asarray(y)
                    x = x[y>0]
                    y = y[y>0]#y[y<0.0] = 0.0
                    if self.color_me != '':
                        for i, cl in enumerate(self.color_me):
                            zorder = len(self.color_me)+ 2 - i
                            if cl[1] in d[self.compound_idx]['lcmsrun'].name:
                                self.ax.plot(x,y,'k-',zorder=zorder,linewidth=2.0,alpha=self.alpha, picker=5, color=cl[0], label = d[self.compound_idx]['lcmsrun'].name.replace('.mzML',''))
                            else:
                                self.ax.plot(x,y,'k-',zorder=1,linewidth=2.0,alpha=self.alpha, picker=5, label = d[self.compound_idx]['lcmsrun'].name.replace('.mzML',''))
                    else:
                        self.ax.plot(x,y,'k-',zorder=1,linewidth=2.0,alpha=self.alpha, picker=5, label = d[self.compound_idx]['lcmsrun'].name.replace('.mzML',''))

        self.min_line = self.ax.axvline(self.my_rt.rt_min, color=self.min_max_color,linewidth=4.0)
        self.max_line = self.ax.axvline(self.my_rt.rt_max, color=self.min_max_color,linewidth=4.0)
        self.peak_line = self.ax.axvline(self.my_rt.rt_peak, color=self.peak_color,linewidth=4.0)

        self.rt_peak_ax = plt.axes([0.09, 0.11, 0.81, 0.02], facecolor=self.slider_color)
        self.rt_max_ax = plt.axes([0.09, 0.14, 0.81, 0.02], facecolor=self.slider_color)
        self.rt_min_ax = plt.axes([0.09, 0.17, 0.81, 0.02], facecolor=self.slider_color)

        self.y_scale_ax = plt.axes([0.925, 0.2755555, 0.02, 0.26], facecolor=self.slider_color)

        min_x = self.ax.get_xlim()[0]
        max_x = self.ax.get_xlim()[1]
#
        self.rt_min_slider = Slider(self.rt_min_ax, 'RT min', min_x, max_x, valinit=self.my_rt.rt_min,color=self.min_max_color)
        self.rt_min_slider.vline.set_color('black')
        self.rt_min_slider.vline.set_linewidth(4)
        self.rt_max_slider = Slider(self.rt_max_ax, 'RT max', min_x, max_x, valinit=self.my_rt.rt_max,color=self.min_max_color)
        self.rt_max_slider.vline.set_color('black')
        self.rt_max_slider.vline.set_linewidth(4)
        self.rt_peak_slider = Slider(self.rt_peak_ax,'RT peak', min_x, max_x, valinit=self.my_rt.rt_peak,color=self.peak_color)
        self.rt_peak_slider.vline.set_color('black')
        self.rt_peak_slider.vline.set_linewidth(4)
        if self.enable_edit:
            self.rt_min_slider.on_changed(self.update_rt)
            self.rt_max_slider.on_changed(self.update_rt)
            self.rt_peak_slider.on_changed(self.update_rt)

        (self.slider_y_min,self.slider_y_max) = self.ax.get_ylim()
        self.slider_val = self.slider_y_max
        self.y_scale_slider = VertSlider(self.y_scale_ax,'',self.slider_y_min,self.slider_y_max, valfmt = '', valinit=self.slider_y_max,color=self.peak_color)
        self.y_scale_slider.vline.set_color('black')
        self.y_scale_slider.vline.set_linewidth(8)
        self.y_scale_slider.on_changed(self.update_yscale)

        self.lin_log_ax = plt.axes([0.1, 0.38, 0.1, 0.15])#, axisbg=axcolor)
        self.lin_log_ax.axis('off')
        self.lin_log_radio = RadioButtons(self.lin_log_ax, ('linear', 'log'))
        self.lin_log_radio.on_clicked(self.set_lin_log)

        self.peak_flag_ax = plt.axes([0.76, 0.38, 0.1, 0.15])#, axisbg=axcolor)
        self.peak_flag_ax.axis('off')
        peak_flags = self.peak_flags
        my_id = metob.retrieve('CompoundIdentification',
                               unique_id = self.data[0][self.compound_idx]['identification'].unique_id, username='*')[-1]
        if my_id.ms1_notes in peak_flags:
            peak_flag_index = peak_flags.index(my_id.ms1_notes)
        else:
            peak_flag_index = 0
        self.peak_flag_radio = RadioButtons(self.peak_flag_ax, peak_flags)
        self.peak_flag_radio.on_clicked(self.set_peak_flag)
        self.peak_flag_radio.set_active(peak_flag_index)
        
        self.msms_flag_ax = plt.axes([0.76, 0.68, 0.1, 0.15])#, axisbg=axcolor)
        self.msms_flag_ax.axis('off')
        msms_flags = self.msms_flags
        if my_id.ms2_notes in msms_flags:
            msms_flag_index = msms_flags.index(my_id.ms2_notes)
        else:
            msms_flag_index = 0
        self.msms_flag_radio = RadioButtons(self.msms_flag_ax, msms_flags)
        self.msms_flag_radio.on_clicked(self.set_msms_flag)
        self.msms_flag_radio.set_active(msms_flag_index)


        #self.fig2,self.ax2 = plt.subplots(figsize=(14, 6))
        my_scan_rt = self.msms_hits.index.get_level_values('msms_scan')
        my_file_name = self.msms_hits.index.get_level_values('file_name')
        #hits_mz_tolerance = 0.005
        hits_mz_tolerance = default_data['identification'].mz_references[-1].mz_tolerance*1e-6
        
        hits = self.msms_hits[(my_scan_rt >= float(self.my_rt.rt_min)) & (my_scan_rt <= float(self.my_rt.rt_max)) & (self.msms_hits['inchi_key'] == inchi_key) \
                & (abs(self.msms_hits['measured_precursor_mz'] - mz_theoretical)/mz_theoretical <= hits_mz_tolerance)]
        self.hits = hits.sort_values('score', ascending=False)

        if len(self.hits) > 0:
            hit_file_name = self.hits.index.get_level_values('file_name')[self.hit_ctr]
            hit_ref_id = self.hits.index.get_level_values('id')[self.hit_ctr]
            hit_score = self.hits['score'].iloc[self.hit_ctr]
            rt_theoretical = self.data[int(file_names.index(hit_file_name))][self.compound_idx]['identification'].rt_references[0].rt_peak
            rt_ms1 = self.data[int(file_names.index(hit_file_name))][self.compound_idx]['data']['ms1_summary']['rt_peak']
            rt_ms2 = self.hits.index.get_level_values('msms_scan')[self.hit_ctr]
            hit_query = self.hits['msv_query_aligned'][self.hit_ctr:self.hit_ctr+1].iloc[0]
            hit_ref = self.hits['msv_ref_aligned'][self.hit_ctr:self.hit_ctr+1].iloc[0]
            mz_precursor = self.hits['measured_precursor_mz'].iloc[self.hit_ctr]
            mz_theoretical = self.data[int(file_names.index(hit_file_name))][self.compound_idx]['identification'].mz_references[0].mz
            mz_measured = self.data[int(file_names.index(hit_file_name))][self.compound_idx]['data']['ms1_summary']['mz_centroid']
            delta_mz = abs(mz_theoretical - mz_measured)
            delta_ppm = delta_mz / mz_theoretical * 1e6
            rt_header = "RT theoretical = %3.2f, RT MS1 measured = %3.2f, RT MS2 measured = %3.2f" % (rt_theoretical, rt_ms2, rt_ms2)
            mz_header = "precursor m/z = %5.4f, m/z theoretical = %5.4f, m/z measured = %5.4f, ppm diff = %3.2f" % (mz_precursor, mz_theoretical, mz_measured, delta_ppm)
            plot_msms_comparison2(0, mz_header, rt_header, hit_ref_id, hit_file_name, hit_score, self.ax2, hit_query, hit_ref)
        
        else:
            hit_query = np.empty((2,2,))
            hit_query[:] = np.nan
            hit_ref = hit_query
            file_idx = file_with_max_ms1_intensity(self.data,self.compound_idx)[0]
            rt_theoretical = default_data['identification'].rt_references[0].rt_peak
            mz_theoretical = default_data['identification'].mz_references[0].mz
            if file_idx != None:
                rt_ms1 = self.data[file_idx][self.compound_idx]['data']['ms1_summary']['rt_peak']
                mz_measured = self.data[file_idx][self.compound_idx]['data']['ms1_summary']['mz_centroid']
            else:
                rt_ms1 = np.nan
                mz_measured = np.nan
            delta_mz = abs(mz_theoretical - mz_measured)
            delta_ppm = delta_mz / mz_theoretical * 1e6
            rt_header = "RT theoretical = %3.2f, RT MS1 measured = %3.2f" % (rt_theoretical, rt_ms1)
            mz_header = "m/z theoretical = %5.4f, m/z measured = %5.4f, ppm diff = %3.2f" % (mz_theoretical, mz_measured, delta_ppm)
            hit_ref_id = "N/A"
            hit_file_name= file_names[file_idx]
            hit_score = np.nan
            plot_msms_comparison2(0, mz_header, rt_header, hit_ref_id, hit_file_name, hit_score, self.ax2, hit_query, hit_ref)

    def set_lin_log(self,label):
        self.ax.set_yscale(label)
        self.ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        self.fig.canvas.draw_idle()

    def set_peak_flag(self,label):
        my_id = metob.retrieve('CompoundIdentification',
                               unique_id = self.data[0][self.compound_idx]['identification'].unique_id, username='*')[-1]
        my_id.ms1_notes = label
        metob.store(my_id)

    def set_msms_flag(self,label):
        my_id = metob.retrieve('CompoundIdentification',
                               unique_id = self.data[0][self.compound_idx]['identification'].unique_id, username='*')[-1]
        my_id.ms2_notes = label
        metob.store(my_id)

    def on_pick(self,event):
        thisline = event.artist
        thisline.set_color('cyan')
        self.ax.set_title(thisline.get_label(), fontsize=7)

    def press(self,event):
        if event.key == 'right':
            if self.compound_idx + 1 < len(self.data[0]):
                self.compound_idx += 1
                self.hit_ctr = 0
                self.ax.cla()
                self.ax2.cla()
                self.rt_peak_ax.cla()
                self.rt_min_ax.cla()
                self.rt_max_ax.cla()
                self.y_scale_ax.cla()
                self.set_plot_data()
        if event.key == 'left':
            if self.compound_idx > 0:
                self.compound_idx -= 1
                self.hit_ctr = 0
                self.ax.cla()
                self.ax2.cla()
                self.rt_peak_ax.cla()
                self.rt_min_ax.cla()
                self.rt_max_ax.cla()
                self.y_scale_ax.cla()
                self.set_plot_data()
        if event.key == 'up':
            if self.hit_ctr > 0:
                self.hit_ctr -= 1
            self.ax.cla()
            self.ax2.cla()
            self.rt_peak_ax.cla()
            self.rt_min_ax.cla()
            self.rt_max_ax.cla()
            self.y_scale_ax.cla()
            self.set_plot_data()
        if event.key == 'down':
            if self.hit_ctr < len(self.hits) - 1:
                self.hit_ctr += 1
            self.ax.cla()
            self.ax2.cla()
            self.rt_peak_ax.cla()
            self.rt_min_ax.cla()
            self.rt_max_ax.cla()
            self.y_scale_ax.cla()
            self.set_plot_data()
        if event.key == 'x':
            self.peak_flag_radio.set_active(1)
            #This is really hacky, but using set_peak_flag function above didn't work.
            my_id = metob.retrieve('CompoundIdentification',
                               unique_id = self.data[0][self.compound_idx]['identification'].unique_id, username='*')[-1]
            my_id.description = 'remove'
            metob.store(my_id)

    def update_yscale(self,val):
        self.y_scale_slider.valinit = self.slider_val
        self.slider_val = self.y_scale_slider.val
        if self.slider_y_min < 0:
            self.slider_y_min = -0.2*self.slider_val
        else:
            self.slider_y_min = 0.02
        self.ax.set_ylim(self.slider_y_min,self.slider_val)
        self.fig.canvas.draw_idle()

    def update_rt(self,val):
        self.my_rt.rt_min = self.rt_min_slider.val
        self.my_rt.rt_max = self.rt_max_slider.val
        if self.adjustable_rt_peak:
            self.my_rt.rt_peak = self.rt_peak_slider.val

        self.rt_min_slider.valinit = self.my_rt.rt_min
        self.rt_max_slider.valinit = self.my_rt.rt_max
        self.rt_peak_slider.valinit = self.my_rt.rt_peak

        metob.store(self.my_rt)
        self.min_line.set_xdata((self.my_rt.rt_min,self.my_rt.rt_min))
        self.max_line.set_xdata((self.my_rt.rt_max,self.my_rt.rt_max))
        self.peak_line.set_xdata((self.my_rt.rt_peak,self.my_rt.rt_peak))
        self.fig.canvas.draw_idle()
        #self.ax.cla()
        #self.ax2.cla()
        #self.rt_peak_ax.cla()
        #self.rt_min_ax.cla()
        #self.rt_max_ax.cla()
        #self.y_scale_ax.cla()
        #self.set_plot_data()
        #self.fig.canvas.draw_idle()

class adjust_mz_for_selected_compound(object):
    def __init__(self,
                 data,
                 include_lcmsruns = None,
                 exclude_lcmsruns = None,
                 include_groups = None,
                 exclude_groups = None,
                 compound_idx = 0,
                 width = 12,
                 height = 6,
                 y_scale='linear',
                 alpha = 0.5,
                 min_max_color = 'sage',
                 peak_color = 'darkviolet',
                 slider_color = 'ghostwhite',
                 y_max = 'auto',
                 y_min = 0):
        """
        data: a metatlas_dataset where files and compounds are stored.
        for example,
        self.metatlas_dataset[file_idx][compound_idx]['identification'].rt_references[-1].unique_id
        is the unique id to the retention time reference for a compound in a file.

        width: specify a width value in inches for the plots and slides
        height: specify a width value in inches for the plots and slides
        min_max_color & peak_color: specify a valid matplotlib color string for the slider and vertical bars
        slider_color: background color for the sliders. Must be a valid matplotlib color

        Press Left and Right arrow keys to move to the next or previous compound
        """

        self.compound_idx = compound_idx
        self.width = width
        self.height = height
        self.y_scale = y_scale
        self.alpha = alpha
        self.min_max_color = min_max_color
        self.peak_color = peak_color
        self.slider_color = slider_color
        self.y_max = y_max
        self.y_min = y_min

        # filter runs from the metatlas dataset
        if include_lcmsruns:
            data = filter_lcmsruns_in_dataset_by_include_list(data,'lcmsrun',include_lcmsruns)

        if include_groups:
            data = filter_lcmsruns_in_dataset_by_include_list(data,'group',include_groups)
        if exclude_lcmsruns:
            data = filter_lcmsruns_in_dataset_by_exclude_list(data,'lcmsrun',exclude_lcmsruns)
        if exclude_groups:
            data = filter_lcmsruns_in_dataset_by_exclude_list(data,'group',exclude_groups)
        self.data = data

        # create figure and first axes
        self.fig,self.ax = plt.subplots(figsize=(width, height))
        plt.subplots_adjust(left=0.09, bottom=0.275)
#         plt.ticklabel_format(style='plain', axis='x')
#         plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        # warn the user if they do not own the atlas; and can not edit its values
        self.enable_edit = True
        self.atlas = metob.retrieve('Atlas',unique_id = self.data[0][0]['atlas_unique_id'],username='*')[-1]
        print(("loaded file for username = ", self.atlas.username))
        if getpass.getuser() != self.atlas.username:
            self.ax.set_title("YOUR ARE %s YOU ARE NOT ALLOWED TO EDIT VALUES THE RT CORRECTOR. USERNAMES ARE NOT THE SAME"%getpass.getuser())
            self.enable_edit = False

        #create all event handlers
        self.fig.canvas.callbacks.connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect('key_press_event', self.press)

        #create the plot
        self.set_plot_data()


    def set_plot_data(self):
        self.ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        self.ax.ticklabel_format(useOffset=False, style='plain', axis='x')

        default_data = self.data[0][self.compound_idx]
        if default_data['identification'].name:
            compound_str = default_data['identification'].name.split('///')[0]
        elif default_data['identification'].compound[-1].name:
            compound_str = default_data['identification'].compound[-1].name
        else:
            compound_str = 'nameless compound'

        compound_str = '%d, %s'%(self.compound_idx, compound_str)

        self.ax.set_title('')
        self.ax.set_ylabel('%s'%compound_str)
        self.ax.set_xlabel('Retention Time')
        self.my_mz = metob.retrieve('MZReference',
                               unique_id = default_data['identification'].mz_references[-1].unique_id, username='*')[-1]
        for i,d in enumerate(self.data): #this loops through the files
            if d[self.compound_idx]['data']['ms1_summary']:
#                 if len(d[self.compound_idx]['data']['ms1_summary']['rt']) > 0:
                    x = d[self.compound_idx]['data']['ms1_summary']['mz_centroid']
                    y = d[self.compound_idx]['data']['ms1_summary']['peak_height']
                    x = np.asarray(x)
                    y = np.asarray(y)
                    self.ax.plot(x,y,'k.',linewidth=2.0,alpha=self.alpha, picker=5, label = d[self.compound_idx]['lcmsrun'].name.replace('.mzML',''))


        mz_delta = self.my_mz.mz_tolerance*self.my_mz.mz/1e6

        self.min_line = self.ax.axvline(self.my_mz.mz-mz_delta, color=self.min_max_color,linewidth=4.0)
        self.max_line = self.ax.axvline(self.my_mz.mz+mz_delta, color=self.min_max_color,linewidth=4.0)
        self.peak_line = self.ax.axvline(self.my_mz.mz, color=self.peak_color,linewidth=4.0)

        min_x = self.ax.get_xlim()[0]
        max_x = self.ax.get_xlim()[1]
        print((min_x,max_x))

        self.mz_peak_ax = plt.axes([0.09, 0.05, 0.81, 0.03], axisbg=self.slider_color)
        self.mz_max_ax = plt.axes([0.09, 0.1, 0.81, 0.03], axisbg=self.slider_color)
        self.mz_min_ax = plt.axes([0.09, 0.15, 0.81, 0.03], axisbg=self.slider_color)

        self.mz_min_slider = Slider(self.mz_min_ax, 'mz min', min_x, max_x, valinit=self.my_mz.mz-mz_delta,color=self.min_max_color,valfmt='%1.4f')
        self.mz_min_slider.vline.set_color('black')
        self.mz_min_slider.vline.set_linewidth(4)

        self.mz_max_slider = Slider(self.mz_max_ax, 'mz max', min_x, max_x, valinit=self.my_mz.mz+mz_delta,color=self.min_max_color,valfmt='%1.4f')
        self.mz_max_slider.vline.set_color('black')
        self.mz_max_slider.vline.set_linewidth(4)

        self.mz_peak_slider = Slider(self.mz_peak_ax,'mz peak', min_x, max_x, valinit=self.my_mz.mz,color=self.peak_color,valfmt='%1.4f')
        self.mz_peak_slider.vline.set_color('black')
        self.mz_peak_slider.vline.set_linewidth(4)
#         if self.enable_edit:
#             self.rt_min_slider.on_changed(self.update_rt)
#             self.rt_max_slider.on_changed(self.update_rt)
#             self.rt_peak_slider.on_changed(self.update_rt)

        self.lin_log_ax = plt.axes([0.1, 0.75, 0.1, 0.15])#, axisbg=axcolor)
        self.lin_log_ax.axis('off')
        self.lin_log_radio = RadioButtons(self.lin_log_ax, ('linear', 'log'))
        self.lin_log_radio.on_clicked(self.set_lin_log)

    def set_lin_log(self,label):
        self.ax.set_yscale(label)
        self.ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        self.fig.canvas.draw_idle()

    def on_pick(self,event):
        thisline = event.artist
        thisline.set_color('red')
        self.ax.set_title(thisline.get_label(), fontsize=7)

    def press(self,event):
        if event.key == 'right':
            if self.compound_idx + 1 < len(self.data[0]):
                self.compound_idx += 1
                self.ax.cla()
                self.mz_peak_ax.cla()
                self.mz_min_ax.cla()
                self.mz_max_ax.cla()
                self.set_plot_data()
        if event.key == 'left':
            if self.compound_idx > 0:
                self.compound_idx -= 1
                self.ax.cla()
                self.mz_peak_ax.cla()
                self.mz_min_ax.cla()
                self.mz_max_ax.cla()
                self.set_plot_data()

#     def update_rt(self,val):
#         self.my_rt.rt_min = self.rt_min_slider.val
#         self.my_rt.rt_max = self.rt_max_slider.val
#         self.my_rt.rt_peak = self.rt_peak_slider.val

#         self.rt_min_slider.valinit = self.my_rt.rt_min
#         self.rt_max_slider.valinit = self.my_rt.rt_max
#         self.rt_peak_slider.valinit = self.my_rt.rt_peak

#         metob.store(self.my_rt)
#         self.min_line.set_xdata((self.my_rt.rt_min,self.my_rt.rt_min))
#         self.max_line.set_xdata((self.my_rt.rt_max,self.my_rt.rt_max))
#         self.peak_line.set_xdata((self.my_rt.rt_peak,self.my_rt.rt_peak))
#         self.fig.canvas.draw_idle()


def replace_compound_id_with_name(x):
    id_list = literal_eval(x)
    if id_list:
        found_compound = metob.retrieve('Compounds',unique_id=id_list[0],username='*')
        return found_compound[-1].name
    else:
        return ''

def make_compound_id_df(data):
    ids = []
    for d in data[0]:
        ids.append(d['identification'])
    df = metob.to_dataframe(ids)
    df['compound'] = df['compound'].apply(replace_compound_id_with_name).astype('str')
    df['rt_unique_id'] = df['rt_references'].apply(lambda x: literal_eval(x))
#     df['mz_unique_id'] = df['mz_references'].apply(lambda x: literal_eval(x))
#     df['frag_unique_id'] = df['frag_references'].apply(lambda x: literal_eval(x))
    df = df[['compound','name','username','rt_unique_id']]#,'mz_unique_id','frag_unique_id']]
    return df

def show_compound_grid(input_fname = '',input_dataset=[]):
    """
    Provide a valid path to data in or a dataset
    """
    if not input_dataset:
        print("loading...")
        data = ma_data.get_dill_data(input_fname)
    else:
        data = input_dataset
    atlas_in_data = metob.retrieve('Atlas',unique_id = data[0][0]['atlas_unique_id'],username='*')
    print(("loaded file for username = ", atlas_in_data[0].username))
    username = getpass.getuser()
    if username != atlas_in_data[0].username:
        print(("YOUR ARE", username, "YOU ARE NOT ALLOWED TO EDIT VALUES THE RT CORRECTOR. USERNAMES ARE NOT THE SAME"))
        #return
    compound_df = make_compound_id_df(data)
    #compound_grid = gui.create_qgrid([])
    #compound_grid.df = compound_df
    compound_grid = qgrid.QGridWidget(df=compound_df)#,set_grid_option={'show_toolbar',True})
    #qgrid.show_grid(compound_df,show_toolbar=True)
    compound_grid.export()
    #display(compound_grid)
    return data,compound_grid



def getcommonletters(strlist):
    """
    Parameters
    ----------
    strlist

    Returns
    -------

    """
    return ''.join([x[0] for x in zip(*strlist) if reduce(lambda a,b:(a == b) and a or None,x)])


def findcommonstart(strlist):
    """
    Parameters
    ----------
    strlist

    Returns
    -------

    """
    strlist = strlist[:]
    prev = None
    while True:
        common = getcommonletters(strlist)
        if common == prev:
            break
        strlist.append(common)
        prev = common

    return getcommonletters(strlist)



def plot_all_compounds_for_each_file(input_dataset = [], input_fname = '', include_lcmsruns = [],exclude_lcmsruns = [], nCols = 8, scale_y=True , output_loc=''):

    """
    Parameters
    ----------
    kwargs

    Returns
    -------

    """

    if not input_dataset:
        data = ma_data.get_dill_data(os.path.expandvars(input_fname))
    else:
        data = input_dataset

    # filter runs from the metatlas dataset
    if include_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_include_list(data,'lcmsrun',include_lcmsruns)
        data = filter_lcmsruns_in_dataset_by_include_list(data,'group',include_lcmsruns)

    if exclude_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_exclude_list(data,'lcmsrun',exclude_lcmsruns)
        data = filter_lcmsruns_in_dataset_by_exclude_list(data,'group',exclude_lcmsruns)

    compound_names = ma_data.get_compound_names(data)[0]
    file_names = ma_data.get_file_names(data)


    output_loc = os.path.expandvars('output_loc')

    nRows = int(np.ceil(len(compound_names)/float(nCols)))


    xmin = 0
    xmax = 210
    subrange = float(xmax-xmin)/float(nCols) # scale factor for the x-axis

    y_max = list()
    if scale_y:
        for file_idx,my_file in enumerate(file_names):
            temp = -1
            counter = 0
            for compound_idx,compound in enumerate(compound_names):
                d = data[file_idx][compound_idx]
                if len(d['data']['eic']['rt']) > 0:
                    counter += 1
                    y = max(d['data']['eic']['intensity'])
                    if y > temp:
                        temp = y
            #y_max.append(temp)
            y_max += [temp] * counter
    else:
        for file_idx,my_file in enumerate(file_names):
            for compound_idx,compound in enumerate(compound_names):
                d = data[file_idx][compound_idx]
                if len(d['data']['eic']['rt']) > 0:
                    y_max.append(max(d['data']['eic']['intensity']))
    y_max = cycle(y_max)

    # create ouput dir
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)


    for file_idx,my_file in enumerate(file_names):
        ax = plt.subplot(111)#, aspect='equal')
        plt.setp(ax, 'frame_on', False)
        ax.set_ylim([0, nRows+7])

        col = 0
        row = nRows+6
        counter = 1

        for compound_idx,compound in enumerate(compound_names):
            if col == nCols:
                row -= 1.3
                col = 0

            d = data[file_idx][compound_idx]

            rt_min = d['identification'].rt_references[0].rt_min
            rt_max = d['identification'].rt_references[0].rt_max
            rt_peak = d['identification'].rt_references[0].rt_peak

            if len(d['data']['eic']['rt']) > 0:
                x = d['data']['eic']['rt']
                y = d['data']['eic']['intensity']
                y = y/next(y_max)
                new_x = (x-x[0])*subrange/float(x[-1]-x[0])+col*(subrange+2) ## remapping the x-range
                xlbl = np.array_str(np.linspace(min(x), max(x), 8), precision=2)
                rt_min_ = (rt_min-x[0])*subrange/float(x[-1]-x[0])+col*(subrange+2)
                rt_max_ = (rt_max-x[0])*subrange/float(x[-1]-x[0])+col*(subrange+2)
                rt_peak_ = (rt_peak-x[0])*subrange/float(x[-1]-x[0])+col*(subrange+2)
                ax.plot(new_x, y+row,'k-')#,ms=1, mew=0, mfc='b', alpha=1.0)]
                #ax.annotate('plot={}'.format(col+1),(max(new_x)/2+col*subrange,row-0.1), size=5,ha='center')
                ax.annotate(xlbl,(min(new_x),row-0.1), size=2)
                ax.annotate('{0},{1},{2},{3}'.format(compound,rt_min, rt_peak, rt_max),(min(new_x),row-0.2), size=2)#,ha='center')
                myWhere = np.logical_and(new_x>=rt_min_, new_x<=rt_max_ )
                ax.fill_between(new_x,min(y)+row,y+row,myWhere, facecolor='c', alpha=0.3)
                col += 1
            else:
                new_x = np.asarray([0,1])#(x-x[0])*subrange/float(x[-1]-x[0])+col*(subrange+2) ## remapping the x-range
                ax.plot(new_x, new_x-new_x+row,'r-')#,ms=1, mew=0, mfc='b', alpha=1.0)]
                ax.annotate(compound,(min(new_x),row-0.1), size=2)
                col += 1
            counter += 1

        plt.title(my_file)
        fig = plt.gcf()
        fig.set_size_inches(nRows*1.0, nCols*4.0)
        fig.savefig(os.path.join(output_loc, my_file + '-' + str(counter) + '.pdf'))
        plt.clf()


def plot_all_files_for_each_compound(input_dataset = [], input_fname = '', include_lcmsruns = [],exclude_lcmsruns = [], nCols = 8, scale_y=True , output_loc=''):

    """
    Parameters
    ----------
    kwargs

    Returns
    -------

    """

    if not input_dataset:
        data = ma_data.get_dill_data(os.path.expandvars(input_fname))
    else:
        data = input_dataset

    # filter runs from the metatlas dataset
    if include_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_include_list(data,'lcmsrun',include_lcmsruns)
        data = filter_lcmsruns_in_dataset_by_include_list(data,'group',include_lcmsruns)

    if exclude_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_exclude_list(data,'lcmsrun',exclude_lcmsruns)
        data = filter_lcmsruns_in_dataset_by_exclude_list(data,'group',exclude_lcmsruns)

    compound_names = ma_data.get_compound_names(data)[0]
    file_names = ma_data.get_file_names(data)
    output_loc = os.path.expandvars(output_loc)

    nRows = int(np.ceil(len(file_names)/float(nCols)))
    print(('nrows = ', nRows))

    xmin = 0
    xmax = 210
    subrange = float(xmax-xmin)/float(nCols) # scale factor for the x-axis


    y_max = list()
    if scale_y:
        for compound_idx,compound in enumerate(compound_names):
            temp = -1
            counter = 0
            for file_idx,my_file in enumerate(file_names):
                d = data[file_idx][compound_idx]
                if len(d['data']['eic']['rt']) > 0:
                    counter += 1
                    y = max(d['data']['eic']['intensity'])
                    if y > temp:
                        temp = y
            y_max += [temp] * counter
    else:
        for compound_idx,compound in enumerate(compound_names):
            for file_idx,my_file in enumerate(file_names):
                d = data[file_idx][compound_idx]
                if len(d['data']['eic']['rt']) > 0:
                    y_max.append(max(d['data']['eic']['intensity']))

    print(("length of ymax is ", len(y_max)))
    y_max = cycle(y_max)



    # create ouput dir
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)
    plt.ioff()
    for compound_idx,compound in enumerate(compound_names):
        ax = plt.subplot(111)#, aspect='equal')
        plt.setp(ax, 'frame_on', False)
        ax.set_ylim([0, nRows+7])

        col = 0
        row = nRows+6
        counter = 1

        for file_idx,my_file in enumerate(file_names):
            if col == nCols:
                row -= 1.3
                col = 0

            d = data[file_idx][compound_idx]
            #file_name = compound_names[compound_idx]

            rt_min = d['identification'].rt_references[0].rt_min
            rt_max = d['identification'].rt_references[0].rt_max
            rt_peak = d['identification'].rt_references[0].rt_peak

            if len(d['data']['eic']['rt']) > 0:
                x = d['data']['eic']['rt']
                y = d['data']['eic']['intensity']
                y = y/next(y_max)
                new_x = (x-x[0])*subrange/float(x[-1]-x[0])+col*(subrange+2) ## remapping the x-range
                xlbl = np.array_str(np.linspace(min(x), max(x), 8), precision=2)
                rt_min_ = (rt_min-x[0])*subrange/float(x[-1]-x[0])+col*(subrange+2)
                rt_max_ = (rt_max-x[0])*subrange/float(x[-1]-x[0])+col*(subrange+2)
                rt_peak_ = (rt_peak-x[0])*subrange/float(x[-1]-x[0])+col*(subrange+2)
                ax.plot(new_x, y+row,'k-')#,ms=1, mew=0, mfc='b', alpha=1.0)]
                #ax.annotate('plot={}'.format(col+1),(max(new_x)/2+col*subrange,row-0.1), size=5,ha='center')
                ax.annotate(xlbl,(min(new_x),row-0.1), size=2)
                ax.annotate('{0},{1},{2},{3}'.format(my_file,rt_min, rt_peak, rt_max),(min(new_x),row-0.2), size=2)#,ha='center')
                myWhere = np.logical_and(new_x>=rt_min_, new_x<=rt_max_ )
                #ax.fill_between(new_x,min(y)+row,y+row,myWhere, facecolor='c', alpha=0.3)
                col += 1
            else:
                new_x = np.asarray([0,1])
                ax.plot(new_x, new_x-new_x+row,'r-')#,ms=1, mew=0, mfc='b', alpha=1.0)]
#                 y = [0,1]#(x-x[0])*subrange/float(x[-1]-x[0])+col*(subrange+2) ## remapping the x-range
#                 ax.plot(new_x, y-y+row,'r-')#,ms=1, mew=0, mfc='b', alpha=1.0)]
                ax.annotate(my_file,(min(new_x),row-0.1), size=1)
                col += 1
            counter += 1

        plt.title(compound)
        fig = plt.gcf()
        fig.set_size_inches(nRows*1.0,nCols*4.0)

        fig.savefig(os.path.join(output_loc, compound + '-' + str(counter) + '.pdf'))
        plt.close(fig)






""" contribution from Hans de Winter """
def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]


def desalt(mol):
    #input is an rdkit mol
    #returns an rdkit mol keeping the biggest component
    #returns original mol if only one component
    #returns a boolean indicated if cleaning was necessary
    d = Chem.rdmolops.GetMolFrags(mol) #these are atom indices
    if len(d) == 1: #If there are fragments or multiple molecules this will be greater than 1
        return mol,False
    my_smiles=Chem.MolToSmiles(mol)
    parent_atom_count=0;
    disconnected=my_smiles.split('.')
    #With GetMolFrags, we've already established that there is more than one disconnected structure
    for s in disconnected:
        little_mol=Chem.MolFromSmiles(s)
        count = little_mol.GetNumAtoms()
        if count > parent_atom_count:
            parent_atom_count = count
            parent_mol = little_mol
    return parent_mol,True

""" contribution from Hans de Winter """
def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

def NeutraliseCharges(mol, reactions=None):
    reactions=_InitialiseNeutralisationReactions()
    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = Chem.AllChem.ReplaceSubstructs(mol, reactant, product)
            rms_smiles = Chem.MolToSmiles(rms[0])
            mol = Chem.MolFromSmiles(rms_smiles)
    if replaced:
        return (mol, True) #Chem.MolToSmiles(mol,True)
    else:
        return (mol, False)


def drawStructure_Fragment(pactolus_tree,fragment_idx,myMol,myMol_w_Hs):
    from copy import deepcopy
    fragment_atoms = np.where(pactolus_tree[fragment_idx]['atom_bool_arr'])[0]
    depth_of_hit = np.sum(pactolus_tree[fragment_idx]['bond_bool_arr'])
    mol2 = deepcopy(myMol_w_Hs)
    # Now set the atoms you'd like to remove to dummy atoms with atomic number 0
    fragment_atoms = np.where(pactolus_tree[fragment_idx]['atom_bool_arr']==False)[0]
    for f in fragment_atoms:
        mol2.GetAtomWithIdx(f).SetAtomicNum(0)

    # Now remove dummy atoms using a query
    mol3 = Chem.DeleteSubstructs(mol2, Chem.MolFromSmarts('[#0]'))
    mol3 = Chem.RemoveHs(mol3)
    # You get what you are looking for
    return moltosvg(mol3),depth_of_hit


def moltosvg(mol,molSize=(450,150),kekulize=True):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # It seems that the svg renderer used doesn't quite hit the spec.
    # Here are some fixes to make it work in the notebook, although I think
    # the underlying issue needs to be resolved at the generation step
    return svg.replace('svg:','')

def get_ion_from_fragment(frag_info,spectrum):
    hit_indices = np.where(np.sum(frag_info,axis=1))
    hit = spectrum[hit_indices,:][0]
    return hit,hit_indices


def calculate_median_of_internal_standards(dataset_for_median,atlas,include_lcmsruns = [],exclude_lcmsruns = [], include_groups = [],exclude_groups = []):
    """

    """

    # filter runs from the metatlas dataset

    # dataset_for_median = copy.deepcopy(dataset_for_median)
    if include_lcmsruns:
        dataset_for_median = filter_lcmsruns_in_dataset_by_include_list(dataset_for_median,'lcmsrun',include_lcmsruns)
    if include_groups:
        dataset_for_median = filter_lcmsruns_in_dataset_by_include_list(dataset_for_median,'group',include_groups)
    if exclude_lcmsruns:
        dataset_for_median = filter_lcmsruns_in_dataset_by_exclude_list(dataset_for_median,'lcmsrun',exclude_lcmsruns)
    if exclude_groups:
        dataset_for_median = filter_lcmsruns_in_dataset_by_exclude_list(dataset_for_median,'group',exclude_groups)
    internal_standard_vals = []
    for i,dd in enumerate(dataset_for_median): #loop through files
        for j,d in enumerate(dd): #loop through compounds
            if atlas.compound_identifications[j].internal_standard_id != 'nan':
                save_dict = {'file_name':d['lcmsrun'].name,'internal_standard_id':atlas.compound_identifications[j].internal_standard_id}
                for fieldname in ['peak_height','peak_area']:
                    if (not d['data']['ms1_summary']) or (not d['data']['ms1_summary'][fieldname]):
                        v = 0
                    else:
                        v = d['data']['ms1_summary'][fieldname]
                    save_dict[fieldname] = v
                internal_standard_vals.append(save_dict)
    return internal_standard_vals


def normalize_peaks_by_internal_standard(metatlas_dataset,atlas,include_lcmsruns = [],exclude_lcmsruns = [], include_groups = [],exclude_groups = []):
    """
    Takes in a metatlas dataset and an atlas. Returns a metatlas dataset with
    ms1_summary peak_height and peak_area normalized by internal standard where
    user selected in their atlas.

    The compound_identification in the atlas has the followign fields:
        internal_standard_id = MetUnicode(help='Freetext identifier for an internal standard')
        do_normalization = MetBool(False)
        internal_standard_to_use = MetUnicode(help='identifier of which internal standard to normalize by')

    Peaks are normalized by:

    I_normalized = I_molecule_in_file / I_standard_in_file * MEDIAN(I_standard_in_good_files)

    "good files" for calculating the median intensity of the standard are identified
    by exclude_lcmsruns=[]

    The patterns in exclude_lcmsruns will remove files that you don't want to use for calculating the median intensity

    """
    internal_standard_vals = calculate_median_of_internal_standards(metatlas_dataset,atlas,include_lcmsruns = include_lcmsruns,
        exclude_lcmsruns =exclude_lcmsruns, include_groups =include_groups,exclude_groups =exclude_groups)

    median_vals = pd.DataFrame(internal_standard_vals).drop('file_name',axis=1).groupby('internal_standard_id').median()
    df = pd.DataFrame(internal_standard_vals)#.drop('peak_height',axis=1)
    norm_dfs = {}
    norm_dfs['peak_area'] = df.pivot(index='internal_standard_id', columns='file_name', values='peak_area')
    norm_dfs['peak_height'] = df.pivot(index='internal_standard_id', columns='file_name', values='peak_height')

    for i,dd in enumerate(metatlas_dataset): #loop through files
        if dd[0]['lcmsrun'].name in norm_dfs['peak_area'].columns: #make sure the file name is in the normalization dataframe
            for j,d in enumerate(dd): #loop through compounds
                if atlas.compound_identifications[j].do_normalization == True:
                    for fieldname in ['peak_height','peak_area']:
                        if (not d['data']['ms1_summary']) or (not d['data']['ms1_summary'][fieldname]):
                            v = 0
                        else:
                            norm_val = norm_dfs[fieldname].loc[atlas.compound_identifications[j].internal_standard_to_use,d['lcmsrun'].name]
                            median_val = median_vals.loc[atlas.compound_identifications[j].internal_standard_to_use,fieldname]
                            metatlas_dataset[i][j]['data']['ms1_summary'][fieldname] = d['data']['ms1_summary'][fieldname] / norm_val * median_val

    return metatlas_dataset

#plot msms and annotate
#compound name
#formula
#adduct
#theoretical m/z
#histogram of retention times
#scatter plot of retention time with peak area
#retention time
#print all chromatograms
#structure

def make_output_dataframe(input_fname = '',input_dataset = [],include_lcmsruns = [],exclude_lcmsruns = [], include_groups = [],exclude_groups = [], output_loc = [], fieldname = 'peak_height', use_labels=False, short_names_df=pd.DataFrame()):
    """
    fieldname can be: peak_height, peak_area, mz_centroid, rt_centroid, mz_peak, rt_peak
    """
    if not input_dataset:
        data = ma_data.get_dill_data(os.path.expandvars(input_fname))
    else:
        data = input_dataset

    # filter runs from the metatlas dataset
    if include_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_include_list(data,'lcmsrun',include_lcmsruns)
    if include_groups:
        data = filter_lcmsruns_in_dataset_by_include_list(data,'group',include_groups)

    if exclude_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_exclude_list(data,'lcmsrun',exclude_lcmsruns)
    if exclude_groups:
        data = filter_lcmsruns_in_dataset_by_exclude_list(data,'group',exclude_groups)

    compound_names = ma_data.get_compound_names(data,use_labels=use_labels)[0]
    file_names = ma_data.get_file_names(data)
    group_names = ma_data.get_group_names(data)
    group_shortnames = ma_data.get_group_shortnames(data)
    output_loc = os.path.expandvars(output_loc)
    fieldname = fieldname

    df = pd.DataFrame( index=compound_names, columns=file_names, dtype=float)

    # peak_height['compound'] = compound_list
    # peak_height.set_index('compound',drop=True)
    for i,dd in enumerate(data):
        for j,d in enumerate(dd):
            if (not d['data']['ms1_summary']) or (not d['data']['ms1_summary'][fieldname]):
                df.loc[compound_names[j],file_names[i]] = 0
            else:
                df.loc[compound_names[j],file_names[i]] = d['data']['ms1_summary'][fieldname]
    columns = []
    if short_names_df.empty:
        for i,f in enumerate(file_names):
            columns.append((group_names[i],f))
        df.columns = pd.MultiIndex.from_tuples(columns,names=['group', 'file'])
    else:
        for i,f in enumerate(file_names):
            temp = [group_names[i],f, group_shortnames[i]]
            temp.extend(short_names_df.loc[f.split('.')[0]].values.tolist())
            columns.append(tuple(temp))
        print(columns)
        df.columns = pd.MultiIndex.from_tuples(columns,names=['group', 'file', 'short groupname', 'sample treatment', 'short filename','short samplename'])

    if output_loc:
        if not os.path.exists(output_loc):
            os.makedirs(output_loc)
        df.to_csv(os.path.join(output_loc, fieldname + '.tab'),sep='\t')

    return df

def file_with_max_precursor_intensity(data,compound_idx):
    idx = None
    my_max = 0
    for i,d in enumerate(data):
        if 'data' in list(d[compound_idx]['data']['msms'].keys()):
            if type(d[compound_idx]['data']['msms']['data']) != list:#.has_key('precursor_intensity'):
                temp = d[compound_idx]['data']['msms']['data']['precursor_intensity']
                if len(temp)>0:
                    m = max(temp)
                    if m > my_max:
                        my_max = m
                        idx = i
    return idx,my_max

def file_with_max_ms1_intensity(data,compound_idx):
    idx = None
    my_max = 0
    for i,d in enumerate(data):
        if 'intensity' in list(d[compound_idx]['data']['eic'].keys()) and d[compound_idx]['data']['eic']['intensity'] != []:
            temp = max(d[compound_idx]['data']['eic']['intensity'])
            if temp > my_max:
                my_max = temp
                idx = i
    return idx,my_max

def file_with_max_score(data, frag_refs, compound_idx, filter_by):
    idx = []
    max_score = np.nan
    best_ref_spec = []

    for file_idx in range(len(data)):
        #empty can look like this:
        # {'eic': {'rt': [], 'intensity': [], 'mz': []}, 'ms1_summary': {'num_ms1_datapoints': 0.0, 'rt_centroid': nan, 'mz_peak': nan, 'peak_height': nan, 'rt_peak': nan, 'peak_area': nan, 'mz_centroid': nan}, 'msms': {'data': {'rt': array([], dtype=float64), 'collision_energy': array([], dtype=float64), 'i': array([], dtype=float64), 'precursor_intensity': array([], dtype=float64), 'precursor_MZ': array([], dtype=float64), 'mz': array([], dtype=float64)}}}
        #or empty can look like this:
        # {'eic': None, 'ms1_summary': None, 'msms': {'data': []}}

        if ('data' in list(data[file_idx][compound_idx]['data']['msms'].keys())) and \
            (isinstance(data[file_idx][compound_idx]['data']['msms']['data'],dict)) and \
            ('rt' in list(data[file_idx][compound_idx]['data']['msms']['data'].keys())) and \
            (len(data[file_idx][compound_idx]['data']['msms']['data']['rt'])>0):

            msv_sample_scans = np.array([data[file_idx][compound_idx]['data']['msms']['data']['mz'], data[file_idx][compound_idx]['data']['msms']['data']['i']])
            rt_of_msv_sample = np.array(data[file_idx][compound_idx]['data']['msms']['data']['rt'])

            scan_idxs = [i+1
                         for i in range(rt_of_msv_sample.size-1)
                         if rt_of_msv_sample[i] != rt_of_msv_sample[i+1]]

            for i, msv_sample in enumerate(np.split(msv_sample_scans, scan_idxs, axis=1)):

                for f, frag in sp.filter_frag_refs(data, frag_refs, compound_idx, file_idx, filter_by).iterrows():
                    msv_ref = sp.sort_ms_vector_by_mz(np.array(frag['mz_intensities']).T)

                    score = sp.score_ms_vectors_composite(*sp.pairwise_align_ms_vectors(msv_sample, msv_ref, .005, 'shape'))

                    if score > max_score or np.isnan(max_score):
                        max_score = score
                        idx = file_idx
                        best_ref_spec = [frag['mz_intensities']]

    return idx, max_score, best_ref_spec

def plot_errorbar_plots(df,output_loc='', use_shortnames=True, ylabel=""):

    output_loc = os.path.expandvars(output_loc)
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)

    plt.ioff()
    for compound in df.index:
        if 'short groupname' in df.columns.names and use_shortnames:
            m = df.loc[compound].groupby(level='short groupname').mean()
            e = df.loc[compound].groupby(level='short groupname').std()
            c = df.loc[compound].groupby(level='short groupname').count()
        else:
            m = df.loc[compound].groupby(level='group').mean()
            e = df.loc[compound].groupby(level='group').std()
            c = df.loc[compound].groupby(level='group').count()

        for i in range(len(e)):
            if c[i]>0:
                e[i] = e[i] / c[i]**0.5

        f, ax = plt.subplots(1, 1,figsize=(12,12))
        m.plot(yerr=e, kind='bar',ax=ax)
        ax.set_title(compound,fontsize=12,weight='bold')
        if ylabel != "":
            plt.ylabel(ylabel)
        plt.tight_layout()
        f.savefig(os.path.join(output_loc, compound + '_errorbar.pdf'))

        #f.clear()
        plt.close(f)#f.clear()

def make_boxplot_plots(df,output_loc='', use_shortnames=True, ylabel=""):

    output_loc = os.path.expandvars(output_loc)
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)

    plt.ioff()
    for compound in df.index:
        f, ax = plt.subplots(1, 1,figsize=(12,12))
        if 'short groupname' in df.columns.names and use_shortnames:
            df.loc[compound].groupby(level='short groupname').apply(pd.DataFrame).plot(kind='box',ax=ax)
        else:
            df.loc[compound].groupby(level='group').apply(pd.DataFrame).plot(kind='box',ax=ax)
        
        for i, (n, grp) in enumerate(df.loc[compound].groupby(level='short groupname')):
            x = [i+1] *len(grp)
            x = np.random.normal(x, 0.04, size=len(x))
            plt.scatter(x, grp)
        ax.set_title(compound,fontsize=12,weight='bold')
        plt.xticks(rotation=90)
        if ylabel != "":
            plt.ylabel(ylabel)
        plt.tight_layout()
        f.savefig(os.path.join(output_loc, compound + '_boxplot.pdf'))
        #f.clear()
        plt.close(f)#f.clear()

def frag_refs_to_json(json_dir = '/project/projectdirs/metatlas/projects/sharepoint/', name = 'frag_refs', save = True):
    ids = metob.retrieve('CompoundIdentification',username='*')
    frag_refs = [cid for cid in ids if cid.frag_references]

    data = {'head_id': [],
            'inchi_key': [],
            'neutralized_inchi_key': [],
            'neutralized_2d_inchi_key': [],
            'polarity': [],
            'collision_energy': [],
            'technique': [],
            'precursor_mz': [],
            'mz_intensities': []}

    for fr in frag_refs:
        data['head_id'].append(fr.frag_references[0].head_id),
        data['inchi_key'].append(fr.compound[0].inchi_key)
        data['neutralized_inchi_key'].append(fr.compound[0].neutralized_inchi_key)
        data['neutralized_2d_inchi_key'].append(fr.compound[0].neutralized_2d_inchi_key)
        data['polarity'].append(fr.frag_references[0].polarity)
        data['precursor_mz'].append(fr.frag_references[0].precursor_mz)
        data['mz_intensities'].append([(m.mz, m.intensity) for m in fr.frag_references[0].mz_intensities])
        data['collision_energy'].append(fr.frag_references[0].collision_energy)
        data['technique'].append(fr.frag_references[0].technique)

    if save:
        with open(os.path.join(json_dir, name + '.json'), 'w') as text_file:
            text_file.write(json.dumps(data))
    else:
        return json.dumps(data)

# def get_idenficications_with_fragrefs():
#     """
#     Select all CompoundIdentifications that have a fragmentation reference
#     """


def make_identification_figure(frag_json_dir = '/project/projectdirs/metatlas/projects/sharepoint/', frag_json_name = 'frag_refs',
    input_fname = '', input_dataset = [], include_lcmsruns = [],
    exclude_lcmsruns = [], include_groups = [], exclude_groups = [], output_loc = [], use_labels=False):
    output_loc = os.path.expandvars(output_loc)
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)

    if not input_dataset:
        data = ma_data.get_dill_data(os.path.expandvars(input_fname))
    else:
        data = input_dataset

    # filter runs from the metatlas dataset
    if include_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_include_list(data,'lcmsrun',include_lcmsruns)
    if include_groups:
        data = filter_lcmsruns_in_dataset_by_include_list(data,'group',include_groups)

    if exclude_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_exclude_list(data,'lcmsrun',exclude_lcmsruns)
    if exclude_groups:
        data = filter_lcmsruns_in_dataset_by_exclude_list(data,'lcmsrun',exclude_groups)
        #data = filter_lcmsruns_in_dataset_by_exclude_list(data,'group',exclude_lcmsruns)


    compound_names = ma_data.get_compound_names(data,use_labels=use_labels)[0]
    file_names = ma_data.get_file_names(data)
    # print(len(data),len(data[0]),len(compound_names))


    frag_refs = pd.read_json(os.path.join(frag_json_dir, frag_json_name + ".json"))


    for compound_idx in range(len(compound_names)):
        file_idx = None
        file_precursor_intensity = 0
        score = None
        ref_spec = []

        if any([(len(data[i][compound_idx]['identification'].compound)!=0) and (data[i][compound_idx]['identification'].compound is not None) for i in range(len(file_names))]):
            # print('checking for compound ids')
            file_idx, score, ref_spec = file_with_max_score(data, frag_refs, compound_idx, 'inchi_key and rt and polarity')
            if ~isinstance(file_idx,int): #There is not a reference for that compound
                file_idx = file_with_max_precursor_intensity(data,compound_idx)[0]
            # print('found one',file_idx)
        else:
            file_idx = file_with_max_precursor_intensity(data,compound_idx)[0]
        # print(file_idx,compound_idx, compound_names[compound_idx])
        if isinstance(file_idx,int):
            # print('printing')
            # print(file_idx,compound_idx)
            fig = plt.figure(figsize=(20,20))
        #     fig = plt.figure()
            ax = fig.add_subplot(211)
            ax.set_title(compound_names[compound_idx],fontsize=12,weight='bold')
            ax.set_xlabel('m/z',fontsize=12,weight='bold')
            ax.set_ylabel('intensity',fontsize=12,weight='bold')

            #TODO: iterate across all collision energies
            precursor_intensity = data[file_idx][compound_idx]['data']['msms']['data']['precursor_intensity']
            idx_max = np.argwhere(precursor_intensity == np.max(precursor_intensity)).flatten()

            mz = data[file_idx][compound_idx]['data']['msms']['data']['mz'][idx_max]
            zeros = np.zeros(data[file_idx][compound_idx]['data']['msms']['data']['mz'][idx_max].shape)
            intensity = data[file_idx][compound_idx]['data']['msms']['data']['i'][idx_max]

            ax.vlines(mz,zeros,intensity,colors='r',linewidth = 2)
            sx = np.argsort(intensity)[::-1]
            labels = [1.001e9]
            for i in sx:
                if np.min(np.abs(mz[i] - labels)) > 0.1 and intensity[i] > 0.02 * np.max(intensity):
                    ax.annotate('%5.4f'%mz[i], xy=(mz[i], 1.01*intensity[i]),rotation = 90, horizontalalignment = 'center', verticalalignment = 'left')
                    labels.append(mz[i])


#                                        precursor_mz = data[file_idx][compound_idx]['data']['msms']['precursor_mz'])
#             print data[file_idx][compound_idx]['data']['msms']['polarity']
            if ref_spec:
                ref_mz = []
                ref_intensity = []
                ref_zeros = []
                for s in ref_spec[0]:
                    ref_mz.append(s[0])
                    ref_intensity.append(s[1]*-1)
                    ref_zeros.append(0)
                s = -1* intensity[sx[0]] / min(ref_intensity)

#                 L = plt.ylim()
#                 print data[file_idx][compound_idx]['identification'].compound[0].name, float(intensity[sx[0]]), float(min(ref_intensity))
                ax.vlines(ref_mz,ref_zeros,[r*s for r in ref_intensity],colors='r',linewidth = 2)
#                 print "we have reference spectra", len(ref_spec[0])
            plt.ioff()
            plt.axhline()
            plt.tight_layout()
            L = plt.ylim()
            plt.ylim(L[0],L[1]*1.12)
            if data[file_idx][compound_idx]['identification'].compound:
                inchi =  data[file_idx][compound_idx]['identification'].compound[0].inchi
                myMol = Chem.MolFromInchi(inchi.encode('utf-8'))
                # myMol,neutralised = NeutraliseCharges(myMol)
                if myMol:
                    image = Draw.MolToImage(myMol, size = (300,300) )
                    ax2 = fig.add_subplot(223)
                    ax2.imshow(image)
                    ax2.axis('off')
            #     SVG(moltosvg(myMol))

            ax3 = fig.add_subplot(224)
            ax3.set_xlim(0,1)
            mz_theoretical = data[file_idx][compound_idx]['identification'].mz_references[0].mz
            mz_measured = data[file_idx][compound_idx]['data']['ms1_summary']['mz_centroid']
            if not mz_measured:
                mz_measured = 0

            delta_mz = abs(mz_theoretical - mz_measured)
            delta_ppm = delta_mz / mz_theoretical * 1e6

            rt_theoretical = data[file_idx][compound_idx]['identification'].rt_references[0].rt_peak
            rt_measured = data[file_idx][compound_idx]['data']['ms1_summary']['rt_peak']
            if not rt_measured:
                rt_measured = 0
            ax3.text(0,1,'%s'%os.path.basename(data[file_idx][compound_idx]['lcmsrun'].hdf5_file),fontsize=12)
            ax3.text(0,0.95,'%s %s'%(compound_names[compound_idx], data[file_idx][compound_idx]['identification'].mz_references[0].adduct),fontsize=12)
            ax3.text(0,0.9,'m/z theoretical = %5.4f, measured = %5.4f, %5.4f ppm difference'%(mz_theoretical, mz_measured, delta_ppm),fontsize=12)
            ax3.text(0,0.85,'Expected Elution of %5.2f minutes, %5.2f min actual'%(rt_theoretical,rt_measured),fontsize=12)
            if score != None:
                ax3.text(0,0.80,'Score: %f'%(score),fontsize=12)
            ax3.set_ylim(0.2,1.01)
            ax3.axis('off')
        #     plt.show()
            fig.savefig(os.path.join(output_loc, compound_names[compound_idx] + '.pdf'))
            plt.close()


def top_five_scoring_files(data, frag_refs, compound_idx, filter_by):
    file_idxs = []
    ref_idxs = []
    scores = []
    msv_sample_list = []
    msv_ref_list = []
    rt_list = []

    for file_idx in range(len(data)):
        try:
            assert(isinstance(data[file_idx][compound_idx]['data']['msms']['data'], dict))
        except AssertionError:
            continue
        except IndexError:
            continue
        except KeyError:
            continue

        msv_sample_scans = np.array([data[file_idx][compound_idx]['data']['msms']['data']['mz'], data[file_idx][compound_idx]['data']['msms']['data']['i']])
        rt_of_msv_sample = np.array(data[file_idx][compound_idx]['data']['msms']['data']['rt'])

        scan_idxs = [i+1
                     for i in range(rt_of_msv_sample.size-1)
                     if rt_of_msv_sample[i] != rt_of_msv_sample[i+1]]

        for i, msv_sample in enumerate(np.split(msv_sample_scans, scan_idxs, axis=1)):
            current_best_score = None
            current_best_ref_idx = None
            current_best_msv_sample = None
            current_best_msv_ref = None
            current_best_rt = None

            for ref_idx, frag in sp.filter_frag_refs(data, frag_refs, compound_idx, file_idx, filter_by).iterrows():
                msv_ref = np.array(frag['mz_intensities']).T

                msv_sample_aligned, msv_ref_aligned = sp.pairwise_align_ms_vectors(msv_sample, msv_ref, .005, 'shape')

                score = sp.score_ms_vectors_composite(msv_sample_aligned, msv_ref_aligned)

                if current_best_score == None or score > current_best_score:
                    current_best_score = score
                    current_best_ref_idx = ref_idx
                    current_best_msv_sample = msv_sample_aligned
                    current_best_msv_ref = msv_ref_aligned
                    current_best_rt = np.split(rt_of_msv_sample, scan_idxs)[i][0]


            if current_best_score:
                scores.append(current_best_score)
                file_idxs.append(file_idx)
                ref_idxs.append(current_best_ref_idx)
                msv_sample_list.append(current_best_msv_sample)
                msv_ref_list.append(current_best_msv_ref)
                rt_list.append(current_best_rt)

    return list(zip(*sorted(zip(file_idxs, ref_idxs, scores, msv_sample_list, msv_ref_list, rt_list), key=lambda l: l[2], reverse=True)[:5]))

def plot_msms_comparison(i, score, ax, msv_sample, msv_ref):

    msv_sample_matches, msv_ref_matches, msv_sample_nonmatches, msv_ref_nonmatches = sp.partition_aligned_ms_vectors(msv_sample, msv_ref)

    msv_sample_unaligned = np.concatenate((msv_sample_matches, msv_sample_nonmatches), axis=1)
    msv_ref_unaligned = np.concatenate((msv_ref_matches, msv_ref_nonmatches), axis=1)

    sample_mz = msv_sample_nonmatches[0]
    sample_zeros = np.zeros(msv_sample_nonmatches[0].shape)
    sample_intensity = msv_sample_nonmatches[1]

    ax.vlines(sample_mz, sample_zeros, sample_intensity, colors='r', linewidth=1)

    shared_mz = msv_sample_matches[0]
    shared_zeros = np.zeros(msv_sample_matches[0].shape)
    shared_sample_intensity = msv_sample_matches[1]

    ax.vlines(shared_mz, shared_zeros, shared_sample_intensity, colors='g', linewidth=1)

    most_intense_idxs = np.argsort(msv_sample_unaligned[1])[::-1]

    if i == 0:
        ax.set_title('%.4f' % score, fontsize=8, weight='bold')
        ax.set_xlabel('m/z', fontsize=8, weight='bold')
        ax.set_ylabel('intensity', fontsize=8, weight='bold')
        ax.tick_params(axis='both', which='major', labelsize=6)

        labels = [1.001e9]

        intensity_requirement = [m for m in most_intense_idxs
                                 if
                                 np.min(np.abs(msv_sample_unaligned[0][m] - labels)) > 0.1
                                 and msv_sample_unaligned[1][m] > 0.2 * np.max(msv_sample_unaligned[1])]

        for m in max([most_intense_idxs[:6], intensity_requirement], key=len):
            if np.min(np.abs(msv_sample_unaligned[0][m] - labels)) > 0.1 and msv_sample_unaligned[1][m] > 0.02 * np.max(msv_sample_unaligned[1]):
                ax.annotate('%5.4f' % msv_sample_unaligned[0][m],
                            xy=(msv_sample_unaligned[0][m], 1.01 * msv_sample_unaligned[1][m]),
                            rotation=90,
                            horizontalalignment='left', verticalalignment='center',
                            size=4)
                labels.append(msv_sample_unaligned[0][m])

    if msv_ref_unaligned[0].size > 0:
        ref_scale = -1 * np.max(msv_sample_unaligned[1]) / np.max(msv_ref_unaligned[1])

        ref_mz = msv_ref_nonmatches[0]
        ref_zeros = np.zeros(msv_ref_nonmatches[0].shape)
        ref_intensity = ref_scale * msv_ref_nonmatches[1]
        shared_ref_intensity = ref_scale * msv_ref_matches[1]

        ax.vlines(ref_mz, ref_zeros, ref_intensity, colors='r', linewidth=1)

        ax.vlines(shared_mz, shared_zeros, shared_ref_intensity, colors='g', linewidth=1)

        ax.axhline()

    ylim = ax.get_ylim()
    ax.set_ylim(ylim[0], ylim[1] * 1.33)

def plot_msms_comparison2(i, mz_header, rt, ref_id, filename, score, ax, msv_sample, msv_ref):

    msv_sample_matches, msv_ref_matches, msv_sample_nonmatches, msv_ref_nonmatches = sp.partition_aligned_ms_vectors(msv_sample, msv_ref)

    msv_sample_unaligned = np.concatenate((msv_sample_matches, msv_sample_nonmatches), axis=1)
    msv_ref_unaligned = np.concatenate((msv_ref_matches, msv_ref_nonmatches), axis=1)

    sample_mz = msv_sample_nonmatches[0]
    sample_zeros = np.zeros(msv_sample_nonmatches[0].shape)
    sample_intensity = msv_sample_nonmatches[1]

    ax.vlines(sample_mz, sample_zeros, sample_intensity, colors='r', linewidth=1)

    shared_mz = msv_sample_matches[0]
    shared_zeros = np.zeros(msv_sample_matches[0].shape)
    shared_sample_intensity = msv_sample_matches[1]

    ax.vlines(shared_mz, shared_zeros, shared_sample_intensity, colors='g', linewidth=1)

    most_intense_idxs = np.argsort(msv_sample_unaligned[1])[::-1]

    if i == 0:
        ax.set_title('MSMS ref ID = %s\n%s' % (ref_id, filename), fontsize=7)
        ax.set_xlabel('m/z\nScore = %.4f, %s\n%s' % (score, rt, mz_header), fontsize=10, weight='bold')
        ax.set_ylabel('intensity', fontsize=10)
        #ax.tick_params(axis='both', which='major', labelsize=8)

        labels = [1.001e9]

        intensity_requirement = [m for m in most_intense_idxs
                                 if
                                 np.min(np.abs(msv_sample_unaligned[0][m] - labels)) > 0.1
                                 and msv_sample_unaligned[1][m] > 0.2 * np.max(msv_sample_unaligned[1])]

        for m in max([most_intense_idxs[:6], intensity_requirement], key=len):
            if np.min(np.abs(msv_sample_unaligned[0][m] - labels)) > 0.1 and msv_sample_unaligned[1][m] > 0.02 * np.max(msv_sample_unaligned[1]):
                ax.annotate('%5.4f' % msv_sample_unaligned[0][m],
                            xy=(msv_sample_unaligned[0][m], 1.01 * msv_sample_unaligned[1][m]),
                            rotation=90,
                            horizontalalignment='left', verticalalignment='center',
                            size=6)
                labels.append(msv_sample_unaligned[0][m])

    if msv_ref_unaligned[0].size > 0:
        ref_scale = -1 * np.max(msv_sample_unaligned[1]) / np.max(msv_ref_unaligned[1])

        ref_mz = msv_ref_nonmatches[0]
        ref_zeros = np.zeros(msv_ref_nonmatches[0].shape)
        ref_intensity = ref_scale * msv_ref_nonmatches[1]
        shared_ref_intensity = ref_scale * msv_ref_matches[1]

        ax.vlines(ref_mz, ref_zeros, ref_intensity, colors='r', linewidth=1)

        ax.vlines(shared_mz, shared_zeros, shared_ref_intensity, colors='g', linewidth=1)

        ax.axhline()

    ylim = ax.get_ylim()
    ax.set_ylim(ylim[0], ylim[1] * 1.33)


def plot_structure(ax, compound, dimensions):
    if compound:
        inchi =  compound[0].inchi
        myMol = Chem.MolFromInchi(inchi.encode('utf-8'))

        if myMol:
            image = Draw.MolToImage(myMol, size=(dimensions, dimensions))
            ax.imshow(image)

    ax.axis('off')


def plot_ema_compound_info(ax, compound_info, label=''):
    wrapper = TextWrapper(width=28, break_on_hyphens=True)

    if compound_info.compound:
        name = ['Name:', wrapper.fill(compound_info.compound[0].name)]
        label = ['Label:', wrapper.fill(label)]
        formula = ['Formula:', compound_info.compound[0].formula]
        polarity = ['Polarity:', compound_info.mz_references[0].detected_polarity]
        neutral_mass = ['Monoisotopic Mass:', compound_info.compound[0].mono_isotopic_molecular_weight]
        theoretical_mz = ['Theoretical M/Z:', compound_info.mz_references[0].mz]
        adduct = ['Adduct:', compound_info.mz_references[0].adduct]

        cell_text = [name, label, formula, polarity, neutral_mass, theoretical_mz, adduct]

        ema_compound_info_table = ax.table(cellText=cell_text,
                                           colLabels=['', 'EMA Compound Info'],
                                           bbox=[0.0, 0.0, 1, 1], loc='top left')
        ema_compound_info_table.scale(1, .7)
        ema_compound_info_table.auto_set_font_size(False)
        ema_compound_info_table.set_fontsize(4)

        cellDict = ema_compound_info_table.get_celld()
        for i in range(len(cell_text)+1):
            cellDict[(i,0)].set_width(0.3)
            cellDict[(i,1)]._loc = 'center'

    ax.axis('off')


def plot_eic(ax, data, compound_idx):
    for file_idx in range(len(data)):

        rt_min = data[file_idx][compound_idx]['identification'].rt_references[0].rt_min
        rt_max = data[file_idx][compound_idx]['identification'].rt_references[0].rt_max
        rt_peak = data[file_idx][compound_idx]['identification'].rt_references[0].rt_peak

        try:
            assert len(data[file_idx][compound_idx]['data']['eic']['rt']) > 1
            x = np.asarray(data[file_idx][compound_idx]['data']['eic']['rt'])
            y = np.asarray(data[file_idx][compound_idx]['data']['eic']['intensity'])

            ax.plot(x, y, 'k-', linewidth=.1, alpha=min(1, 10*(1./len(data))))
            myWhere = np.logical_and(x>=rt_min, x<=rt_max )
            ax.fill_between(x,0,y,myWhere, facecolor='c', alpha=min(1, 2*(1./len(data))))
        except (AssertionError, TypeError):
            pass

    # ax.tick_params(labelbottom='off')
    ax.xaxis.set_tick_params(labelsize=5)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().set_visible(False)
    ax.axvline(rt_min, color='k', linewidth=1.0)
    ax.axvline(rt_max, color='k', linewidth=1.0)
    ax.axvline(rt_peak, color='r', linewidth=1.0)


def plot_score_and_ref_file(ax, score, rt, ref):
    ax.text(0.5, 1, '%.4f'%score,
        weight='bold',
        horizontalalignment='center',
        verticalalignment='top',
        fontsize=4,
        transform=ax.transAxes)

    ax.text(0, .45, fill(str(ref) + ' RT=%5.3f'%rt, width=26),
        horizontalalignment='left',
        verticalalignment='center',
        rotation='vertical',
        fontsize=2,
        transform=ax.transAxes)

def get_msms_hits(metatlas_dataset, use_labels=False, extra_time=False, keep_nonmatches=False,
                  pre_query='database == "metatlas"',
                  # pre_query = 'index == index or index == @pd.NaT',
                  query='(@inchi_key == inchi_key) and (@polarity == polarity) and ((@precursor_mz - .5*(((.5*(@pre_mz_ppm**-decimal)/(decimal+1)) + .005 + ((.5*(@pre_mz_ppm**-decimal)/(decimal+1)) - .005)**2)**.5)) <= precursor_mz <= (@precursor_mz + .5*(((.5*(@pre_mz_ppm**-decimal)/(decimal+1)) + .005 + ((.5*(@pre_mz_ppm**-decimal)/(decimal+1)) - .005)**2)**.5)))',
                  # query='(@inchi_key == inchi_key) and (@polarity == polarity) and ((@precursor_mz - (.5*(@pre_mz_ppm**-decimal)/(decimal+1)) - @pre_mz_ppm*(@precursor_mz*1e-6)) <= precursor_mz <= (@precursor_mz + (.5*(@pre_mz_ppm**-decimal)/(decimal+1)) + @pre_mz_ppm*(@precursor_mz*1e-6)))',
                  # query='(@inchi_key == inchi_key) and (@polarity == polarity) and (@rt-.1 < rt < @rt+.1)  and ((@precursor_mz - (.5*(@pre_mz_ppm**-decimal)/(decimal+1)) - @pre_mz_ppm*(@precursor_mz*1e-6)) <= precursor_mz <= (@precursor_mz + (.5*(@pre_mz_ppm**-decimal)/(decimal+1)) + @pre_mz_ppm*(@precursor_mz*1e-6)))',
                  **kwargs):
    kwargs = dict(locals(), **kwargs)

    resolve_by = kwargs.pop('resolve_by', 'shape')
    frag_mz_tolerance = kwargs.pop('frag_mz_tolerance', .005)

    # Reference parameters
    ref_loc = kwargs.pop('ref_loc', '/global/project/projectdirs/metatlas/projects/spectral_libraries/msms_refs_v2.tab')
    ref_dtypes = kwargs.pop('ref_dtypes', {'database':str, 'id':str, 'name':str,
                                           'spectrum':object,'decimal':int, 'precursor_mz':float,
                                           'polarity':str, 'adduct':str, 'fragmentation_method':str,
                                           'collision_energy':str, 'instrument':str, 'instrument_type':str,
                                           'formula':str, 'exact_mass':float,
                                           'inchi_key':str, 'inchi':str, 'smiles':str})

    ref_index = kwargs.pop('ref_index', ['database', 'id'])
    if 'do_centroid' in kwargs:
        do_centroid = kwargs.pop('do_centroid')
    else:
        do_centroid = False
        
    if 'ref_df' in kwargs:
        ref_df = kwargs.pop('ref_df')
    else:
        ref_df = pd.read_csv(ref_loc,
                             sep='\t',
                             dtype=ref_dtypes
                            ).set_index(ref_index)

    ref_df = ref_df.query(pre_query, local_dict=dict(locals(), **kwargs))

    if ref_df['spectrum'].apply(type).eq(str).all():
        ref_df['spectrum'] = ref_df['spectrum'].apply(lambda s: eval(s)).apply(np.array)

    file_names = ma_data.get_file_names(metatlas_dataset)
    compound_names = ma_data.get_compound_names(metatlas_dataset)[0]

    msms_hits = []

    for compound_idx,compound_name in enumerate(compound_names):
        sys.stdout.write('\r'+'Processing: {} / {} compounds.'.format(compound_idx+1,len(compound_names)))
        sys.stdout.flush()

        #Code below is commented out to make get_msms_hits work when there isn't a compound in identification - VS, Nov 2019
        #if len(metatlas_dataset[0][compound_idx]['identification'].compound) == 0:
            # exit here if there isn't a compound in the identification
        #    continue

        if metatlas_dataset[0][compound_idx]['identification'].name:
            name = metatlas_dataset[0][compound_idx]['identification'].name.split('///')[0]
        elif metatlas_dataset[0][compound_idx]['identification'].compound[-1].name:
            name = metatlas_dataset[0][compound_idx]['identification'].compound[-1].name
        else:
            name = None

        try:
            adduct = metatlas_dataset[0][compound_idx]['identification'].mz_references[0].adduct
        except (KeyError, AttributeError):
            adduct = None

        if len(metatlas_dataset[0][compound_idx]['identification'].compound) > 0:
            inchi_key = metatlas_dataset[0][compound_idx]['identification'].compound[0].inchi_key
        else:
            inchi_key = ''
        pre_mz_ppm = metatlas_dataset[0][compound_idx]['identification'].mz_references[0].mz_tolerance
        precursor_mz = metatlas_dataset[0][compound_idx]['identification'].mz_references[0].mz

        rt_min = metatlas_dataset[0][compound_idx]['identification'].rt_references[0].rt_min
        rt_max = metatlas_dataset[0][compound_idx]['identification'].rt_references[0].rt_max

        compound_hits = []

        for file_idx,file_name in enumerate(file_names):

            polarity = metatlas_dataset[file_idx][compound_idx]['identification'].mz_references[0].detected_polarity

            try:
                assert set(['rt', 'i', 'precursor_MZ', 'mz']).issubset(set(metatlas_dataset[file_idx][compound_idx]['data']['msms']['data'].keys()))
            except (KeyError, AssertionError, AttributeError):
                continue

            rt_mz_i_df = pd.DataFrame({k:metatlas_dataset[file_idx][compound_idx]['data']['msms']['data'][k]
                                      for k in ['rt', 'mz', 'i', 'precursor_MZ', 'precursor_intensity']}
                                      ).sort_values(['rt', 'mz'])

            for rt in rt_mz_i_df.rt.unique():
                if not extra_time:
                    if not rt_min <= rt <= rt_max:
                        continue

                msv_sample = rt_mz_i_df.loc[rt_mz_i_df['rt'] == rt,['mz', 'i','rt','precursor_MZ','precursor_intensity']]
                precursor_mz_sample = msv_sample['precursor_MZ'].values[0]
                precursor_intensity_sample = msv_sample['precursor_intensity'].values[0]
                msv_sample = msv_sample[['mz','i']].values.T
                
                if do_centroid:
                    max_peaks, min_peaks = sp.peakdet(msv_sample[1], 1000.0)
                    if max_peaks.shape[0]>0:
                        idx = max_peaks[:,0].astype(int).flatten()
                        msv_sample = msv_sample[:,idx]
                    else:
                        msv_sample = np.zeros((0,0))
#                 msv_sample.sort_values('mz',inplace=True)
#                 msv_sample = msv_sample
 
                #Filter ions greater than 2.5 + precursor M/Z 
#                 msv_sample[1] = msv_sample[1] / msv_sample[1].sum()
#                 print(msv_sample)
                if msv_sample.size > 0:
                    msv_sample = msv_sample[:,msv_sample[0] < precursor_mz_sample + 2.5]
                    scan_df = sp.search_ms_refs(msv_sample, **dict(locals(), **kwargs))
                else:
                    scan_df = {}

                if len(scan_df) > 0:
                    scan_df['file_name'] = file_name
                    scan_df['msms_scan'] = rt
                    scan_df['name'] = name
                    scan_df['adduct'] = adduct
                    scan_df['inchi_key'] = inchi_key
                    scan_df['precursor_mz'] = precursor_mz
                    scan_df['measured_precursor_mz'] = precursor_mz_sample
                    scan_df['measured_precursor_intensity'] = precursor_intensity_sample

                    scan_df.set_index('file_name', append=True, inplace=True)
                    scan_df.set_index('msms_scan', append=True, inplace=True)

                    msms_hits.append(scan_df)

                elif keep_nonmatches:
                    scan_df = {}

                    scan_df['file_name'] = file_name
                    scan_df['msms_scan'] = rt
                    scan_df['name'] = name
                    scan_df['adduct'] = adduct
                    scan_df['inchi_key'] = inchi_key
                    scan_df['precursor_mz'] = precursor_mz
                    scan_df['num_matches'] = 0
                    scan_df['measured_precursor_mz'] = precursor_mz_sample
                    scan_df['measured_precursor_intensity'] = precursor_intensity_sample

                    scan_df['score'] = precursor_intensity_sample
                    scan_df['msv_query_aligned'] = msv_sample
                    scan_df['msv_ref_aligned'] = np.full_like(msv_sample, np.nan)

                    scan_df = pd.DataFrame([scan_df])

                    for idx in ref_df.index.names:
                        scan_df[idx] = None
                    scan_df.set_index('database', append=False, inplace=True)
                    scan_df.set_index('id', append=True, inplace=True)
                    scan_df.set_index('file_name', append=True, inplace=True)
                    scan_df.set_index('msms_scan', append=True, inplace=True)

                    msms_hits.append(scan_df)

    sys.stdout.write('\n'+'Done!!!')
    if len(msms_hits)>0:
        hits = pd.concat(msms_hits)
        return hits
        #Check if number of matches for a compound across all files is 1 or less and set the score to its maximum intensity.
        #This will identify MSMS with single ion / no fragmentation
#         print(hits.groupby(['inchi_key', 'adduct'])['num_matches'])
#         if keep_nonmatches==True:
#             idxs = hits.groupby(['inchi_key', 'adduct'])['num_matches'].transform(max) <= 1 
#             hits['score'][idxs] = hits['measured_precursor_intensity'][idxs]
        g = hits.groupby(['inchi_key', 'adduct'])['num_matches'].transform(max)
        idxs =  g<= 1
        proc_idx = g[idxs].index
        for i in proc_idx:
            hits.loc[i,'score'] = hits.loc[i,'measured_precursor_intensity']
        return hits
    else:
        return pd.DataFrame(columns=ref_df.index.names+['file_name', 'msms_scan', 'score', 'num_matches','inchi_key','precursor_mz','adduct','score']
                           ).set_index(ref_df.index.names+['file_name', 'msms_scan'])

def make_chromatograms(
    input_dataset = [], group='index', share_y = True, save=True, output_loc=[], short_names_df=pd.DataFrame(), short_names_header=None):
    
    file_names = ma_data.get_file_names(input_dataset)
    
    if short_names_df.empty:
        if short_names_header != None:
            sys.stdout.write('short_names_df not provided. Using full_filename for the plots!')
            short_names_df = pd.DataFrame()
    elif short_names_header == None:
            sys.stdout.write('short_names_header not provided. Using full_filename for the plots!')
            short_names_df = pd.DataFrame()
    elif short_names_header not in short_names_df.columns:
            sys.stdout.write('short_names_header not found in short_names_df. Using full_filename for the plots!')
            short_names_df = pd.DataFrame()
    else:
        short_names_df = short_names_df[[short_names_header]]
        short_names_df.columns=['shortname']


    
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)
    compound_names = ma_data.get_compound_names(input_dataset,use_labels=True)[0]
    args_list = []

    chromatogram_str = 'compound_EIC_chromatograms'

    if not os.path.exists(os.path.join(output_loc,chromatogram_str)):
        os.makedirs(os.path.join(output_loc,chromatogram_str))

    for compound_idx, my_compound in enumerate(compound_names):
        my_data = list()
        for file_idx, my_file in enumerate(file_names):
            my_data.append(input_dataset[file_idx][compound_idx])
        kwargs = {'data': my_data,
                 'file_name': os.path.join(output_loc, chromatogram_str, my_compound+'.pdf'),
                 'group': group,
                 'save': save,
                 'share_y': share_y,
                 'names': file_names,
                 #'shortname':findcommonstart(file_names)}
                 'shortname':short_names_df}
        args_list.append(kwargs)
    max_processes = 4
    pool = mp.Pool(processes=min(max_processes, len(input_dataset[0])))
    pool.map(cpp.chromplotplus, args_list)
    pool.close()
    pool.terminate()

def make_identification_figure_v2(
    input_fname = '', input_dataset = [], include_lcmsruns = [], exclude_lcmsruns = [], include_groups = [],
    exclude_groups = [], output_loc = [], msms_hits = None, use_labels=False,intensity_sorted_matches=False, short_names_df=pd.DataFrame()):
    #empty can look like this:
    # {'eic': {'rt': [], 'intensity': [], 'mz': []}, 'ms1_summary': {'num_ms1_datapoints': 0.0, 'rt_centroid': nan, 'mz_peak': nan, 'peak_height': nan, 'rt_peak': nan, 'peak_area': nan, 'mz_centroid': nan},
    #'msms': {'data': {'rt': array([], dtype=float64), 'collision_energy': array([], dtype=float64), 'i': array([], dtype=float64), 'precursor_intensity': array([], dtype=float64), 'precursor_MZ': array([], dtype=float64), 'mz': array([], dtype=float64)}}}
    #or empty can look like this:
    # {'eic': None, 'ms1_summary': None, 'msms': {'data': []}}
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)
    if not input_dataset:
        data = ma_data.get_dill_data(os.path.expandvars(input_fname))
    else:
        data = input_dataset

    #Filter runs from the metatlas dataset
    if include_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_include_list(data, 'lcmsrun', include_lcmsruns)
    if include_groups:
        data = filter_lcmsruns_in_dataset_by_include_list(data, 'group', include_groups)

    if exclude_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_exclude_list(data, 'lcmsrun', exclude_lcmsruns)
    if exclude_groups:
        data = filter_lcmsruns_in_dataset_by_exclude_list(data, 'group', exclude_groups)

    #msms_hits_df = get_msms_hits(data, use_labels, ref_index=['database', 'id', 'inchi_key', 'precursor_mz'])
    msms_hits_df = msms_hits.copy()

    if msms_hits_df is not None:
        #if 'inchi_key' in msms_hits_df.columns:
        #    msms_hits_df.rename(columns={'inchi_key':'inchi_key_2'},inplace=True)
        #msms_hits_df.reset_index(['inchi_key', 'precursor_mz'], inplace=True)
        msms_hits_df.reset_index(inplace = True)
        msms_hits_df.sort_values('score', ascending=False, inplace=True)
        # msms_hits_df.drop_duplicates(['inchi_key', 'file_name'], keep='first', inplace=True)
        # msms_hits_df = msms_hits_df.groupby(['inchi_key']).head(5).sort_values(['inchi_key'], kind='mergesort')

    #Obtain compound and file names
    compound_names = ma_data.get_compound_names(data,use_labels)[0]
    file_names = ma_data.get_file_names(data)

    df = pd.DataFrame()
    #Turn off interactive plotting
    plt.ioff()
    plt.clf()
    #Iterate over compounds
    for compound_idx in range(len(compound_names)):
        sys.stdout.write('\r'+'Making Identification Figure for: {} / {} compounds.'.format(compound_idx+1,len(compound_names)))
        sys.stdout.flush()
        file_idxs, scores, msv_sample_list, msv_ref_list, rt_list = [], [], [], [], []
        
        if len(data[0][compound_idx]['identification'].compound) > 0 and hasattr(data[0][compound_idx]['identification'].compound[0],"inchi_key"):
            inchi_key = data[0][compound_idx]['identification'].compound[0].inchi_key
        else:
            inchi_key = ""

        #Find 5 best file and reference pairs by score
        try:
            comp_msms_hits = msms_hits_df[(msms_hits_df['inchi_key'] == inchi_key) \
                                          & (msms_hits_df['msms_scan'] >= data[0][compound_idx]['identification'].rt_references[0].rt_min) \
                                          & (msms_hits_df['msms_scan'] <= data[0][compound_idx]['identification'].rt_references[0].rt_max) \
                                          & ((abs(msms_hits_df['precursor_mz'].values.astype(float) - data[0][compound_idx]['identification'].mz_references[0].mz)/data[0][compound_idx]['identification'].mz_references[0].mz) \
                                             <= data[0][compound_idx]['identification'].mz_references[0].mz_tolerance*1e-6)].drop_duplicates('file_name').head(5)
            # Dont need assert anymore, keep_nonmatch in get_msms_hits should replace the assert
            #assert len(comp_msms_hits) > 0
            file_idxs = [file_names.index(f) for f in comp_msms_hits['file_name']]
            scores = comp_msms_hits['score'].values.tolist()
            msv_sample_list = comp_msms_hits['msv_query_aligned'].values.tolist()
            msv_ref_list = comp_msms_hits['msv_ref_aligned'].values.tolist()
            rt_list = comp_msms_hits['msms_scan'].values.tolist()

        #except (IndexError, AssertionError, TypeError) as e:
        except (IndexError, TypeError) as e:
            file_idx = None
            max_intensity = 0

            for fi in range(len(data)):
                try:
                    temp = max(data[fi][compound_idx]['data']['eic']['intensity'])
                    if temp > max_intensity:
                        file_idx = fi
                        max_intensity = temp
                except (ValueError,TypeError):
                    continue

            file_idxs = [file_idx]
            msv_sample_list = [np.array([0, np.nan]).T]
            msv_ref_list = [np.array([0, np.nan]).T]
            scores = [np.nan]


        #Plot if compound yields any scores
        if file_idxs and file_idxs[0] is not None:
            #Top 5 MSMS Spectra
            ax1 = plt.subplot2grid((24, 24), (0, 0), rowspan=12, colspan=12)
            ax2a = plt.subplot2grid((24, 24), (0, 12), rowspan=3, colspan=3)
            ax2a.tick_params(axis='both', length=2)
            ax2a.set_xticklabels([])
            ax2a.set_yticklabels([])
            ax2b = plt.subplot2grid((24, 24), (3, 12), rowspan=3, colspan=3)
            ax2b.tick_params(axis='both', length=2)
            ax2b.set_xticklabels([])
            ax2b.set_yticklabels([])
            ax2c = plt.subplot2grid((24, 24), (6, 12), rowspan=3, colspan=3)
            ax2c.tick_params(axis='both', length=2)
            ax2c.set_xticklabels([])
            ax2c.set_yticklabels([])
            ax2d = plt.subplot2grid((24, 24), (9, 12), rowspan=3, colspan=3)
            ax2d.tick_params(axis='both', length=2)
            ax2d.set_xticklabels([])
            ax2d.set_yticklabels([])


            for i,(score,ax) in enumerate(zip(scores,[ax1, ax2a, ax2b, ax2c, ax2d])):
                plot_msms_comparison(i, score, ax, msv_sample_list[i], msv_ref_list[i])



            #Next Best Scores and Filenames
            ax4a = plt.subplot2grid((24, 24), (0, 15), rowspan=3, colspan=1)
            ax4a.axis('off')
            ax4b = plt.subplot2grid((24, 24), (3, 15), rowspan=3, colspan=1)
            ax4b.axis('off')
            ax4c = plt.subplot2grid((24, 24), (6, 15), rowspan=3, colspan=1)
            ax4c.axis('off')
            ax4d = plt.subplot2grid((24, 24), (9, 15), rowspan=3, colspan=1)
            ax4d.axis('off')


            if short_names_df.empty:
                for i,(score,ax) in enumerate(zip(scores[1:],[ax4a, ax4b, ax4c, ax4d])):
                    plot_score_and_ref_file(ax, score, rt_list[i+1], os.path.basename(data[file_idxs[i+1]][compound_idx]['lcmsrun'].hdf5_file))
            else:
                for i,(score,ax) in enumerate(zip(scores[1:],[ax4a, ax4b, ax4c, ax4d])):
                    short_samplename  = short_names_df.loc[os.path.basename(data[file_idxs[i+1]][compound_idx]['lcmsrun'].hdf5_file).split('.')[0], 'short_samplename'][0]
                    plot_score_and_ref_file(ax, score, rt_list[i+1], short_samplename)

        
        #EMA Compound Info
        if file_idxs and file_idxs[0] is not None:
            ax3 = plt.subplot2grid((24, 24), (0, 16), rowspan=6, colspan=8)
            plot_ema_compound_info(ax3, data[file_idxs[0]][compound_idx]['identification'])#,
                               # ma_data.get_compound_names(data,use_labels=True)[0][compound_idx])
        else:
            ax3 = plt.subplot2grid((24, 24), (0, 0), rowspan=6, colspan=8)
            plot_ema_compound_info(ax3, data[0][compound_idx]['identification'])#,


        #Structure
        if file_idxs and file_idxs[0] is not None:
            ax5 = plt.subplot2grid((24, 24), (13, 0), rowspan=6, colspan=6)
            plot_structure(ax5, data[file_idxs[0]][compound_idx]['identification'].compound, 100)
        else:
            ax5 = plt.subplot2grid((24, 24), (13, 0), rowspan=6, colspan=6)
            plot_structure(ax5, data[0][compound_idx]['identification'].compound, 100)

        #EIC
        if file_idxs and file_idxs[0] is not None:
            ax6 = plt.subplot2grid((24, 24), (6, 16), rowspan=6, colspan=6)
            plot_eic(ax6, data, compound_idx)
        else:
            ax6 = plt.subplot2grid((24, 24), (6, 0), rowspan=6, colspan=6)
            plot_eic(ax6, data, compound_idx)

#             #Reference and Sample Info
#             ax10 = plt.subplot2grid((24, 24), (14, 6), rowspan=10, colspan=20)
#             plot_ref_sample_info(ax10, 1, 1)

        #Old code
        if file_idxs and file_idxs[0] is not None:
            ax7 = plt.subplot2grid((24, 24), (15, 6), rowspan=9, colspan=20)
            mz_theoretical = data[file_idxs[0]][compound_idx]['identification'].mz_references[0].mz
            mz_measured = data[file_idxs[0]][compound_idx]['data']['ms1_summary']['mz_centroid']
            if not mz_measured:
                mz_measured = 0

            delta_mz = abs(mz_theoretical - mz_measured)
            delta_ppm = delta_mz / mz_theoretical * 1e6

            rt_theoretical = data[file_idxs[0]][compound_idx]['identification'].rt_references[0].rt_peak
            rt_measured = data[file_idxs[0]][compound_idx]['data']['ms1_summary']['rt_peak']
            if not rt_measured:
                rt_measured = 0
            ax7.text(0,1,'%s'%fill(os.path.basename(data[file_idxs[0]][compound_idx]['lcmsrun'].hdf5_file), width=54),fontsize=8)
            ax7.text(0,0.9,'%s %s'%(compound_names[compound_idx], data[file_idxs[0]][compound_idx]['identification'].mz_references[0].adduct),fontsize=8)
            ax7.text(0,0.85,'Measured M/Z = %5.4f, %5.4f ppm difference'%(mz_measured, delta_ppm),fontsize=8)
            ax7.text(0,0.8,'Expected Elution of %5.2f minutes, %5.2f min actual'%(rt_theoretical,rt_measured),fontsize=8)
            if len(rt_list) > 0:
                ax7.text(0,0.7,'MSMS Scan at %5.3f minutes'%rt_list[0],fontsize=8)
                msv_sample_matches = sp.partition_aligned_ms_vectors(msv_sample_list[0], msv_ref_list[0])[0]
                if intensity_sorted_matches:
                    msv_sample_matches = msv_sample_matches[:, msv_sample_matches[1].argsort()[::-1]]
                if len(msv_sample_matches[0]) > 0:
                    mz_sample_matches = msv_sample_matches[0].tolist()
                    threshold_mz_sample_matches = sp.remove_ms_vector_noise(msv_sample_matches)[0].tolist()
                else:
                    mz_sample_matches = [np.nan]
                    threshold_mz_sample_matches = [np.nan]
                ax7.text(0,0.6,
                         fill('Matching M/Zs above 1E-3*max: ' + ', '.join(['%5.3f'%m for m in threshold_mz_sample_matches]), width=90) + '\n\n' +
                         fill('All Matching M/Zs: ' + ', '.join(['%5.3f'%m for m in mz_sample_matches]), width=90),
                         fontsize=6, verticalalignment='top')
                df.loc[compound_idx, 'label'] = compound_names[compound_idx]
                df.loc[compound_idx, 'file name'] = file_names[file_idxs[0]]
                df.loc[compound_idx, 'RT'] = rt_list[0]
                df.loc[compound_idx, 'score'] = scores[0]
                df.loc[compound_idx, 'Matching M/Zs above 1E-3*max'] =', '.join(['%5.3f'%m for m in threshold_mz_sample_matches])
                df.loc[compound_idx, 'All matching M/Zs'] = ','.join(['%5.3f'%m for m in mz_sample_matches])

            ax7.set_ylim(.5,1.1)
            ax7.axis('off')

        plt.savefig(os.path.join(output_loc, compound_names[compound_idx] + '.pdf'))
        plt.close()
    df.to_csv(os.path.join(output_loc, 'MatchingMZs.tab'),sep='\t')


def plot_ms1_spectra(polarity = None, mz_min = 5, mz_max = 5, input_fname = '', input_dataset = [], compound_names = [],  include_lcmsruns = [], exclude_lcmsruns = [], include_groups = [], exclude_groups = [], output_loc = []):
    """
    Plot three views of ms1 spectra for compounds in input_dataset using file with highest RT peak of a polarity:
    Unscaled: plots ms1 spectra within window of mz_min and mz_max
    Scaled: plots ms1 spectra within window of mz_min and mz_max scaling mz of compound to 70%
    Full Range: plots ms1 spectra without window (unscaled)
    """
    print('here I am')
    if not input_dataset:
        data = ma_data.get_dill_data(os.path.expandvars(input_fname))
    else:
        data = input_dataset

    if include_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_include_list(data, 'lcmsrun', include_lcmsruns)
    if include_groups:
        data = filter_lcmsruns_in_dataset_by_include_list(data, 'group', include_groups)

    if exclude_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_exclude_list(data, 'lcmsrun', exclude_lcmsruns)
    if exclude_groups:
        data = filter_lcmsruns_in_dataset_by_exclude_list(data, 'group', exclude_groups)

    #Make sure there is data
    assert(len(data) != 0)

    all_compound_names = ma_data.get_compound_names(data)[0]

    #Set default compound list to all compounds in input_dataset
    if not compound_names:
        compound_names = all_compound_names

    #Find implicit polarity and make sure there is not more than one
    if 'POS' in include_lcmsruns or 'NEG' in exclude_lcmsruns:
        assert(polarity == None or polarity == 'positive')
        polarity = 'positive'
    if 'NEG' in include_lcmsruns or 'POS' in exclude_lcmsruns:
        assert(polarity == None or polarity == 'negative')
        polarity = 'negative'

    if 'POS' in include_groups or 'NEG' in exclude_groups:
        assert(polarity == None or polarity == 'positive')
        polarity = 'positive'
    if 'NEG' in include_groups or 'POS' in exclude_groups:
        assert(polarity == None or polarity == 'negative')
        polarity = 'negative'

    assert(polarity == 'positive' or polarity == 'negative')

    #Additional variables used acorss all compounds
    lcms_polarity = 'ms1_' + polarity[:3]
    titles = ['Unscaled', 'Scaled', 'Full Range']

    for compound_idx in [i for i,c in enumerate(all_compound_names) if c in compound_names]:
        print(('compound is',compound_idx))
        #Find file_idx of with highest RT peak
        highest = 0
        file_idx = None
        for i,d in enumerate(data):
            if d[compound_idx]['identification'].mz_references[0].detected_polarity == polarity:
                if d[compound_idx]['data']['ms1_summary']['peak_height'] > highest:
                    highest = d[compound_idx]['data']['ms1_summary']['peak_height']
                    file_idx = i

        lcms_data = ma_data.df_container_from_metatlas_file(data[file_idx][compound_idx]['lcmsrun'].hdf5_file)

        #Find RT and mz peak for compound in file
        rt_peak = data[file_idx][compound_idx]['data']['ms1_summary']['rt_peak']
        rt_peak_actual = lcms_data[lcms_polarity].iloc[(lcms_data[lcms_polarity].rt - rt_peak).abs().argsort()[0]].rt
        mz_peak_actual = data[file_idx][compound_idx]['data']['ms1_summary']['mz_peak']

        #Create and sort dataframe containing RT peak, mz and intensity
        df_all = lcms_data[lcms_polarity][(lcms_data[lcms_polarity].rt == rt_peak_actual)]
        df_all.sort_values('i',ascending=False,inplace=True)

        #Limit prior dataframe to +/- mz_min, mz_max
        df_window = df_all[(df_all['mz'] > mz_peak_actual - mz_min) &
                           (df_all['mz'] < mz_peak_actual + mz_max) ]

        #Plot compound name, mz, and RT peak
        plt.ioff()
        fig = plt.gcf()
        fig.suptitle('%s, m/z: %5.4f, rt: %f'%(all_compound_names[compound_idx], mz_peak_actual, rt_peak_actual),
                                                fontsize=8,weight='bold')

        #Create axes for different views of ms1 spectra (unscaled, scaled, and full range)
        ax1 = plt.subplot2grid((11, 12), (0, 0), rowspan=5, colspan=5)
        ax2 = plt.subplot2grid((11, 12), (0, 7), rowspan=5, colspan=5)
        ax3 = plt.subplot2grid((11, 12), (6, 0), rowspan=5, colspan=12)

        #Plot ms1 spectra
        for ax_idx,(ax,df) in enumerate(zip([ax1, ax2, ax3], [df_window, df_window, df_all])):

            ax.set_xlabel('m/z',fontsize=8,weight='bold')
            ax.set_ylabel('intensity',fontsize=8,weight='bold')
            ax.tick_params(axis='both', which='major', labelsize=6)
            ax.set_title(titles[ax_idx],fontsize=8,weight='bold')

            mzs = df['mz']
            zeros = np.zeros(len(df['mz']))
            intensities = df['i']

            ax.vlines(mzs, zeros, intensities, colors='r',linewidth = 2)

            labels = [1.001e9]
            for i,row in df.iloc[:6].iterrows():
                ax.annotate('%.4f'%row.mz, xy=(row.mz, 1.03*row.i),rotation = 90, horizontalalignment = 'center', verticalalignment = 'left', fontsize=6)
                labels.append(row.mz)

            ax.axhline(0)

            if ax_idx != 2:
                ax.set_xlim(mz_peak_actual - mz_min,  mz_peak_actual + mz_max)

            ylim = ax.get_ylim()

            if ax_idx == 1:
                ax.set_ylim(ylim[0], df[((mz_peak_actual - .05 < df['mz']) & (df['mz'] < mz_peak_actual + .05))].iloc[0]['i']*1.43)
            else:
                ax.set_ylim(ylim[0], ylim[1]*1.43)

            if not os.path.exists(output_loc):
                os.makedirs(output_loc)

        plt.savefig(os.path.join(output_loc, all_compound_names[compound_idx] + '.pdf'))


def export_atlas_to_spreadsheet(myAtlas, output_filename='', input_type = 'atlas'):
    """
    Return a pandas dataframe containing Atlas info.  Optionally save it.
    This function can also work on a MetAtlas dataset (list of lists returned by get_data_for_atlas_and_groups).
    """
    cols = [c for c in metob.Compound.class_trait_names() if not c.startswith('_')]
    cols = sorted(cols)
    atlas_export = pd.DataFrame( )

    if input_type != 'atlas':
        num_compounds = len(myAtlas[0])
    else:
        num_compounds = len(myAtlas.compound_identifications)

    for i in range(num_compounds):
        if input_type != 'atlas':
            my_id = myAtlas[0][i]['identification']
            n = my_id.name
        else:
            my_id = myAtlas.compound_identifications[i]

        if my_id.compound:
            for c in cols:
                g = getattr(my_id.compound[0],c)
                if g:
                    atlas_export.loc[i,c] = g
                else:
                    atlas_export.loc[i,c] = ''
        atlas_export.loc[i, 'label'] = my_id.name
        atlas_export.loc[i, 'id_notes'] = my_id.description
        atlas_export.loc[i,'rt_min'] = my_id.rt_references[0].rt_min
        atlas_export.loc[i,'rt_max'] = my_id.rt_references[0].rt_max
        atlas_export.loc[i,'rt_peak'] = my_id.rt_references[0].rt_peak
        atlas_export.loc[i,'mz'] = my_id.mz_references[0].mz
        atlas_export.loc[i,'mz_tolerance'] = my_id.mz_references[0].mz_tolerance
        atlas_export.loc[i,'adduct'] = my_id.mz_references[0].adduct
        atlas_export.loc[i,'polarity'] = my_id.mz_references[0].detected_polarity
        # if my_id.frag_references:
        #     atlas_export.loc[i,'has_fragmentation_reference'] = True
        #     # TODO: Gather the frag reference information and export it
        # else:
        #     atlas_export.loc[i,'has_fragmentation_reference'] = False

    if output_filename:
        if not os.path.exists(os.path.dirname(output_filename)):
            os.makedirs(os.path.dirname(output_filename))
        atlas_export.to_csv(output_filename)

    return atlas_export

def get_data_for_groups_and_atlas(group,myAtlas,output_filename,use_set1 = False):
    """
    get and pickle everything This is MSMS, raw MS1 datapoints, compound, group info, and file info
    """
    data = []
    import copy as copy
    for i,treatment_groups in enumerate(group):
        for j in range(len(treatment_groups.items)):
            myFile = treatment_groups.items[j].hdf5_file
    #         try:
    #             rt_reference_index = int(treatment_groups.name[-1]) - 1
    #         except:
    #             rt_reference_index = 3
            print((i, len(group), myFile))
            row = []
            for compound in myAtlas.compound_identifications:
                result = {}
                result['atlas_name'] = myAtlas.name
                result['atlas_unique_id'] = myAtlas.unique_id
                result['lcmsrun'] = treatment_groups.items[j]
                result['group'] = treatment_groups
                temp_compound = copy.deepcopy(compound)
                if use_set1:
                    if '_Set1' in treatment_groups.name:
                        temp_compound.rt_references[0].rt_min -= 0.2
                        temp_compound.rt_references[0].rt_max -= 0.2
                        temp_compound.rt_references[0].rt_peak -= 0.2
                    temp_compound.mz_references[0].mz_tolerance = 20
                result['identification'] = temp_compound
                result['data'] = ma_data.get_data_for_a_compound(temp_compound.mz_references[0],
                                        temp_compound.rt_references[0],
                                        [ 'ms1_summary', 'eic', 'msms' ],
                                        myFile,0.2)
    #                 print result['data']['ms1_summary']
                row.append(result)
            data.append(row)
        with open(output_filename,'w') as f:
            dill.dump(data,f)

def filter_atlas(atlas_df = '', input_dataset = [], num_data_points_passing = 5, peak_height_passing = 1e6):
    metatlas_dataset = input_dataset
    num_data_points_passing = np.array([[((metatlas_dataset[i][j]['identification'].rt_references[0].rt_min-0.75 <= 
                                       np.array(metatlas_dataset[i][j]['data']['eic']['rt'])) &
                                      (np.array(metatlas_dataset[i][j]['data']['eic']['rt']) <=
                                       metatlas_dataset[i][j]['identification'].rt_references[0].rt_max+0.75)).sum()>num_data_points_passing
                                    for i in range(len(metatlas_dataset))]
                                    for j in range(len(metatlas_dataset[0]))]).any(axis=1)
    peak_height_passing = np.array([[('data' in list(metatlas_dataset[i][j].keys())) &
                                 ('ms1_summary' in list(metatlas_dataset[i][j]['data'].keys())) &
                                 (metatlas_dataset[i][j]['data']['ms1_summary']['peak_height']>=peak_height_passing)
                                for i in range(len(metatlas_dataset))]
                                for j in range(len(metatlas_dataset[0]))]).any(axis=1)
    compound_passing = num_data_points_passing & peak_height_passing
    return atlas_df[compound_passing].reset_index(drop=True)

def filter_metatlas_objects_to_most_recent(object_list,field):
    #from datetime import datetime, date
    #remove from list if another copy exists that is newer
    unique_values = []
    for i,a in enumerate(object_list):
        unique_values.append( getattr(a,field) )
    unique_values = list(set(unique_values))
    keep_object_list = []
    for u in unique_values:
        old_last_modified = 0
        for i,a in enumerate(object_list):
            if getattr(a,field) == u:
                last_modified = getattr(a,'last_modified')
                if last_modified > old_last_modified:
                    keep_object = a
                    old_last_modified = last_modified
        keep_object_list.append(keep_object)
    return keep_object_list
#        print i, a.name,  datetime.utcfromtimestamp(a.last_modified)

def get_metatlas_atlas(name = '%%',username = '*', most_recent = True,do_print = True):
    from datetime import datetime, date
    atlas = metob.retrieve('Atlas',name = name,username=username)
    if most_recent:
        atlas = filter_metatlas_objects_to_most_recent(atlas,'name')
    if do_print:
        for i,a in enumerate(atlas):
            print((i, len(a.compound_identifications),a.name,  datetime.utcfromtimestamp(a.last_modified)))

    return atlas

class interact_get_metatlas_files():
    def __init__(self, experiment = '%violacein%', name = '%_%', most_recent = True):
        self.experiment = experiment
        self.name = name
        self.most_recent = most_recent
#         http://ipywidgets.readthedocs.io/en/latest/examples/Using%20Interact.html
        self.w = interact(self.Task, experiment=self.experiment, name=self.name, most_recent = self.most_recent,__manual=True)#continuous_update=False)#

    def Task(self,experiment,name,most_recent):
        self.experiment = experiment
        self.name = name
        self.most_recent = most_recent
        self.files = get_metatlas_files(experiment = experiment,name = name,most_recent = most_recent)#self.most_recent)
        txt = widgets.Text()
        txt.value = '%d Files were found matching that pattern'%len(self.files)
        display(txt)



def get_metatlas_files(experiment = '%%',name = '%%',most_recent = True):
    """
    experiment is the folder name
    name is the filename
    """
    files = metob.retrieve('LcmsRun',experiment=experiment,name=name, username='*')
    if most_recent:
        files = filter_metatlas_objects_to_most_recent(files,'mzml_file')
    return files

def make_prefilled_fileinfo_sheet(groups, filename):
    #make a prefilled fileinfo sheet for editing groups manually and reimport to workflow
    with open(filename,'w') as fid:
        fid.write('mzml_file\tgroup\tdescription\tshort_name\n')
        for g in groups:
            for f in g.items:
                fid.write('%s\t%s\t%s\t%s\n'% (f.mzml_file, g.name, f.description, g.short_name))


def make_empty_fileinfo_sheet(filename,flist):
    #dump all the files to a spreadheet, download it, and make a "filled in" one.
    with open(filename,'w') as fid:
        fid.write('mzml_file\tgroup\tdescription\n')
        for f in flist:
            fid.write('%s\t\t\n'%f.mzml_file)

def make_groups_from_fileinfo_sheet(filename,filetype='tab',store=False):
    '''

    '''
    if filetype == 'tab':
        df = pd.read_csv(filename,sep='\t')
    elif filetype == 'csv':
        df = pd.read_csv(filename,sep=',')
    elif filetype == 'df':
        df = filename
    else:
        df = pd.read_excel(filename)
    grouped = df.groupby(by='group')
    return_groups = []
    for g in grouped.groups.keys():
        indices = grouped.groups[g]
        myGroup = metob.Group()
        myGroup.name = '%s'%g
        myGroup.description = df.loc[indices[0],'description']
        file_set = []
        for i in indices:
            file_set.append(metob.retrieve('LcmsRun',mzml_file='%%%s'%df.loc[i,'mzml_file'],username='*')[0])
        myGroup.items = file_set
        return_groups.append(myGroup)
        if store:
            metob.store(myGroup)
    return return_groups

def check_compound_names(df):
    # compounds that have the wrong compound name will be listed
    # Keep running this until no more compounds are listed
    bad_names = []
    for i,row in df.iterrows():
        #if type(df.name[x]) != float or type(df.label[x]) != float:
            #if type(df.name[x]) != float:
        if (not pd.isnull(row.inchi_key)) and (len(row.inchi_key)>0):# or type(df.inchi_key[x]) != float:
            if not metob.retrieve('Compounds',inchi_key=row.inchi_key, username = '*'):
                print((row.inchi_key, "compound is not in database. Exiting Without Completing Task!"))
                bad_names.append(row.inchi_key)
    return bad_names


def check_file_names(df,field):
    bad_files = []
    for i,row in df.iterrows():
        if row[field] != '':
            if not metob.retrieve('Lcmsruns',name = '%%%s%%'%row[field],username = '*'):
                print((row[field], "file is not in the database. Exiting Without Completing Task!"))
                bad_files.append(row[field])
    return bad_files


# def get_formatted_atlas_from_google_sheet(polarity='POS',
#                                           method='QE_HILIC',
#                                           mz_tolerance=10):
#     import metatlas.ms_monitor_util as mmu
#     df2 = mmu.get_ms_monitor_reference_data()
#     #print df.head()
#     #df2 = pd.DataFrame(df[1:],columns=df[0])

#     fields_to_keep = [ 'name',
#                     'label',
#                       'inchi_key',
#                     'mz_%s'%polarity,
#                     'rt_min_%s'%method,
#                     'rt_max_%s'%method,
#                     'rt_peak_%s'%method,
#                     'file_mz_%s_%s'%(method,polarity),
#                     'file_rt_%s_%s'%(method,polarity),
#                     'file_msms_%s_%s'%(method,polarity)]

#     fields_there = []
#     for f in fields_to_keep:
#          if f in df2.keys():
#                 fields_there.append(f)

#     df3 = df2.loc[:,fields_there]

#     df3['mz_tolerance'] = mz_tolerance

#     if polarity == 'POS':
#         df3['polarity'] = 'positive'
#     else:
#         df3['polarity'] = 'negative'

#     renamed_columns = [c.replace('_%s'%method,'').replace('_%s'%polarity,'') for c in df3.columns]
#     for i,c in enumerate(df3.columns):
#         df3 = df3.rename(columns = {c:renamed_columns[i]})
#     df3 = df3[df3['mz'] != '']

#     return df3


def make_atlas_from_spreadsheet(filename='valid atlas file.csv',
                                atlas_name='20161007_MP3umZHILIC_BPB_NEG_ExampleAtlasName',
                                filetype=('excel','csv','tab','dataframe'),
                                sheetname='only for excel type input',
                                polarity = ('positive','negative'),
                                store=False,
                                mz_tolerance=None):
    '''
    specify polarity as 'positive' or 'negative'

    '''
    if isinstance(filename,pd.DataFrame):
        df = filename
    else:
        if ( filetype=='excel' ) and sheetname:
            df = pd.read_excel(filename,sheetname=sheetname)
        elif ( filetype=='excel' ):
            df = pd.read_excel(filename)
        elif filetype == 'tab':
            df = pd.read_csv(filename,sep='\t')
        else:
            df = pd.read_csv(filename,sep=',')
    df.dropna(how="all", inplace=True)
    df.columns = [x.lower() for x in df.columns]

    if 'inchi_key' not in df.columns:
        df['inchi_key'] = ""
    if 'adduct' not in df.columns:
        df['adduct'] = ""

    bad_names = check_compound_names(df)
    if bad_names:
        return bad_names
    #Make sure all the files specified for references are actually there
    #if 'file_rt' in df.keys():
        #strip '.mzmL' from cells
        #df.file_rt = df.file_rt.str.replace('\..+', '')
        #bad_files = check_file_names(df,'file_rt')
        #if bad_files:
        #     return bad_files
    #if 'file_mz' in df.keys():
    #    #strip '.mzmL' from cells
    #    df.file_mz = df.file_mz.str.replace('\..+', '')
    #    bad_files = check_file_names(df,'file_mz')
    #    if bad_files:
    #         return bad_files
    if 'file_msms' in list(df.keys()):
        #strip '.mzmL' from cells
        df.file_msms = df.file_msms.str.replace('\..+', '')
        bad_files = check_file_names(df,'file_msms')
        if bad_files:
             return bad_files



    all_identifications = []

#     for i,row in df.iterrows():
    for i,row in df.iterrows():
        if type(row.inchi_key) != float or type(row.label) != float: #this logic is to skip empty rows

            myID = metob.CompoundIdentification()

            if (not pd.isnull(row.inchi_key)) and (len(row.inchi_key)>0): # this logic is where an identified metabolite has been specified
                c = metob.retrieve('Compounds',inchi_key=row.inchi_key,username = '*') #currently, all copies of the molecule are returned.  The 0 is the most recent one.
                if c:
                    c = c[-1]
            else:
                c = 'use_label'
            if type(row.label) != float:
                compound_label = row.label #if no name, then use label as descriptor
            else:
                compound_label = 'no label'

            if c:
                if c != 'use_label':
                    myID.compound = [c]
                myID.name = compound_label

                try:
                    myID.do_normalization = row.do_normalization
                    myID.internal_standard_id = row.internal_standard_id
                    myID.internal_standard_to_use = row.internal_standard_to_use
                except:
                    # no internal standard information was provided
                    pass

                try:
                    myID.identification_notes = row.identification_notes
                except:
                    # no identification_notes were provided
                    pass



                mzRef = metob.MzReference()
                # take the mz value from the spreadsheet
                mzRef.mz = row.mz
                #TODO: calculate the mz from theoretical adduct and modification if provided.
                #     mzRef.mz = c.MonoIso topic_molecular_weight + 1.007276
                if mz_tolerance:
                    mzRef.mz_tolerance = mz_tolerance
                else:
                    try:
                        mzRef.mz_tolerance = row.mz_tolerance
                    except:
                        if 'mz_threshold' in df.columns:
                            mzRef.mz_tolerance = row.mz_threshold
                        else:
                            sys.exit("mz_tolerance or mz_threshold not provided. Can't make atlas.")
                
                mzRef.mz_tolerance_units = 'ppm'
                mzRef.detected_polarity = polarity
                #if 'file_mz' in df.keys():
                #    f = metob.retrieve('Lcmsruns',name = '%%%s%%'%df.file_mz[x],username = '*')[0]
                #    mzRef.lcms_run = f
                if 'adduct' in row:
                    if ~pd.isnull(row.adduct):
                        mzRef.adduct = row.adduct

                myID.mz_references = [mzRef]

                rtRef = metob.RtReference()
                rtRef.rt_units = 'min'
                rtRef.rt_min = row.rt_min
                rtRef.rt_max = row.rt_max
                rtRef.rt_peak = row.rt_peak
                #if 'file_rt' in df.keys():
                #    f = metob.retrieve('Lcmsruns',name = '%%%s%%'%df.file_rt[x],username = '*')[0]
                #    rtRef.lcms_run = f
                myID.rt_references = [rtRef]

                if ('file_msms' in list(df.keys())) and (c != 'use_label'):
                    if (type(row.file_msms) != float) and (row.file_msms != ''):
                        frag_ref = metob.FragmentationReference()
                        f = metob.retrieve('Lcmsruns',name = '%%%s%%'%row.file_msms,username = '*')[0]
                        frag_ref.lcms_run = f
                        frag_ref.polarity = polarity
                        frag_ref.precursor_mz = row.mz

                        data = ma_data.get_data_for_a_compound(mzRef, rtRef, [ 'msms' ],f.hdf5_file,0.3)
                        if isinstance(data['msms']['data'], np.ndarray):
                            precursor_intensity = data['msms']['data']['precursor_intensity']
                            idx_max = np.argwhere(precursor_intensity == np.max(precursor_intensity)).flatten()
                            mz = data['msms']['data']['mz'][idx_max]
                            intensity = data['msms']['data']['i'][idx_max]
                            spectrum = []
                            for i in range(len(mz)):
                                mzp = metob.MzIntensityPair()
                                mzp.mz = mz[i]
                                mzp.intensity = intensity[i]
                                spectrum.append(mzp)
                            frag_ref.mz_intensities = spectrum
                            myID.frag_references = [frag_ref]
                            print('')
                            print(('found reference msms spectrum for ',myID.compound[0].name, 'in file',row.file_msms))

                all_identifications.append(myID)

    myAtlas = metob.Atlas()
    myAtlas.name = atlas_name
    myAtlas.compound_identifications = all_identifications
    if store:
        metob.store(myAtlas)
    return myAtlas

def filter_empty_metatlas_objects(object_list,field):
    filtered_list = []
    for i,g in enumerate(object_list):
        try:
            #This bare try/accept is to handle the invalid groups left over in the database from the original objects.
            #These groups don't conform to the current schema and will throw an error when you query their attributes.
            if (len(getattr(g,field)) > 0):
                filtered_list.append(g)
        except:
            pass
    return filtered_list

def filter_metatlas_objects_by_list(object_list,field,filter_list):
    filtered_list = []
    for i,g in enumerate(object_list):
        if any(ext in getattr(g,field) for ext in filter_list):
            filtered_list.append(g)
    return filtered_list

def remove_metatlas_objects_by_list(object_list,field,filter_list):
    filtered_list = []
    for i,g in enumerate(object_list):
        if not any(ext in getattr(g,field) for ext in filter_list):
            filtered_list.append(g)
    return filtered_list

def filter_lcmsruns_in_dataset_by_include_list(metatlas_dataset,selector,include_list):
    """
    Returns a metatlas dataset containing LCMS runs or groups (denoted by selector) that have substrings listed in the include list
    selector can be 'lcmsrun' or 'group'
    include_list will look something like this: ['QC','Blank']
    """
    filtered_dataset = []
    for d in metatlas_dataset:
        if any(ext in d[0][selector].name for ext in include_list):
            filtered_dataset.append(d)
    return filtered_dataset

def filter_lcmsruns_in_dataset_by_exclude_list(metatlas_dataset,selector,exclude_list):
    """
    Returns a metatlas dataset containing LCMS runs or groups (denoted by selector) that have substrings not listed in the include list
    selector can be 'lcmsrun' or 'group'
    exclude_list will look something like this: ['QC','Blank']
    """
    filtered_dataset = []
    for d in metatlas_dataset:
        if not any(ext in d[0][selector].name for ext in exclude_list):
            filtered_dataset.append(d)
    return filtered_dataset


def filter_compounds_in_dataset_by_exclude_list(metatlas_dataset,exclude_list):
    """
    Since the rows of the dataset are expected to line up with an atlas export, this is probably not a good idea to use.
    """
    filtered_dataset = []
    for d_row in metatlas_dataset:
        filtered_row = []
        for d in d_row:
            if not any(ext in d['identification'].name for ext in exclude_list):
                if not any(ext in d['identification'].compound[0].name for ext in exclude_list):
                    filtered_row.append(d)
        filtered_dataset.append(filtered_row)
    return filtered_dataset

def filter_compounds_in_dataset_by_include_list(metatlas_dataset,include_list):
    """
    Since the rows of the dataset are expected to line up with an atlas export, this is probably not a good idea to use.
    """
    filtered_dataset = []
    for d_row in metatlas_dataset:
        filtered_row = []
        for d in d_row:
            if any(ext in d['identification'].name for ext in include_list):
                if any(ext in d['identification'].compound[0].name for ext in include_list):
                    filtered_row.append(d)
        filtered_dataset.append(filtered_row)
    return filtered_dataset

def select_groups_for_analysis(name = '%', description = [], username = '*', do_print = True, most_recent = True, remove_empty = True, include_list = [], exclude_list = []):
    if description:
        groups = metob.retrieve('Groups', name = name, description = description, username=username)
    else:
        groups = metob.retrieve('Groups', name = name, username=username)
    if most_recent:
        groups = filter_metatlas_objects_to_most_recent(groups,'name')

    if include_list:
        groups = filter_metatlas_objects_by_list(groups,'name',include_list)

    if exclude_list:
        groups = remove_metatlas_objects_by_list(groups,'name',exclude_list)

    print((len(groups)))

    if remove_empty:
        groups = filter_empty_metatlas_objects(groups,'items')
    if do_print:
        from datetime import datetime, date
        for i,a in enumerate(groups):
            print((i, a.name,  datetime.utcfromtimestamp(a.last_modified)))

    return groups



#
