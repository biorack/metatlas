import sys
import os
import os.path
# os.environ['R_LIBS_USER'] = '/project/projectdirs/metatlas/r_pkgs/'
#curr_ld_lib_path = ''


from metatlas import metatlas_objects as metob
from metatlas import h5_query as h5q
from metatlas.helpers import metatlas_get_data_helper_fun as ma_data
from metatlas import gui

import qgrid
import pandas as pd
import os
import tables
import pickle
import dill
import numpy as np
import re
from matplotlib import pyplot as plt

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem, Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D, IPythonConsole
from itertools import cycle
from collections import defaultdict
from IPython.display import SVG,display


from ipywidgets import interact, interactive, fixed
import ipywidgets as widgets
from IPython.display import display

import getpass

from ast import literal_eval
# from datetime import datetime

from matplotlib.widgets import Slider, Button, RadioButtons

from matplotlib.widgets import AxesWidget

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
        for cid, func in self.observers.iteritems():
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
        print("loaded file for username = ", self.atlas.username)
        if getpass.getuser() != self.atlas.username:
            self.ax.set_title("YOUR ARE %s YOU ARE NOT ALLOWED TO EDIT VALUES THE RT CORRECTOR. USERNAMES ARE NOT THE SAME"%getpass.getuser())
            self.enable_edit = False
            
        #create all event handlers
        self.fig.canvas.callbacks.connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect('key_press_event', self.press)

        #create the plot
        self.set_plot_data()
        

    def set_plot_data(self):
        #set y-scale and bounds if provided
        self.ax.set_yscale(self.y_scale)
        if self.y_max != 'auto':
            self.ax.set_ylim(self.y_min,self.y_max)
            
        self.ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        
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
        self.my_rt = metob.retrieve('RTReference',
                               unique_id = default_data['identification'].rt_references[-1].unique_id, username='*')[-1]
        for d in self.data: #this loops through the files
            if d[self.compound_idx]['data']['eic']:
                if len(d[self.compound_idx]['data']['eic']['rt']) > 0:
                    x = d[self.compound_idx]['data']['eic']['rt']
                    y = d[self.compound_idx]['data']['eic']['intensity']
                    x = np.asarray(x)
                    y = np.asarray(y)
                    minval = np.min(y[np.nonzero(y)])
                    y = y - minval
                    x = x[y>0]
                    y = y[y>0]#y[y<0.0] = 0.0
                    self.ax.plot(x,y,'k-',linewidth=2.0,alpha=self.alpha, picker=5, label = d[self.compound_idx]['lcmsrun'].name.replace('.mzML',''))
                    

        self.min_line = self.ax.axvline(self.my_rt.rt_min, color=self.min_max_color,linewidth=4.0)
        self.max_line = self.ax.axvline(self.my_rt.rt_max, color=self.min_max_color,linewidth=4.0)
        self.peak_line = self.ax.axvline(self.my_rt.rt_peak, color=self.peak_color,linewidth=4.0)
      
        self.rt_peak_ax = plt.axes([0.09, 0.05, 0.81, 0.03], axisbg=self.slider_color)
        self.rt_max_ax = plt.axes([0.09, 0.1, 0.81, 0.03], axisbg=self.slider_color)
        self.rt_min_ax = plt.axes([0.09, 0.15, 0.81, 0.03], axisbg=self.slider_color)

        self.y_scale_ax = plt.axes([0.925, 0.275, 0.02, 0.63], axisbg=self.slider_color)

        min_x = self.ax.get_xlim()[0]
        max_x = self.ax.get_xlim()[1]
        
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
        self.ax.set_title(thisline.get_label())        
    
    def press(self,event):
        if event.key == 'right':
            if self.compound_idx + 1 < len(self.data[0]):
                self.compound_idx += 1
                self.ax.cla()
                self.rt_peak_ax.cla()
                self.rt_min_ax.cla()
                self.rt_max_ax.cla()
                self.y_scale_ax.cla()
                self.set_plot_data()
        if event.key == 'left':
            if self.compound_idx > 0:
                self.compound_idx -= 1
                self.ax.cla()
                self.rt_peak_ax.cla()
                self.rt_min_ax.cla()
                self.rt_max_ax.cla()
                self.y_scale_ax.cla()
                self.set_plot_data()
    
    def update_yscale(self,val):
        self.y_scale_slider.valinit = self.slider_val
        self.slider_val = self.y_scale_slider.val
        self.ax.set_ylim(self.slider_y_min,self.slider_val)
        self.fig.canvas.draw_idle()
        
    def update_rt(self,val):
        self.my_rt.rt_min = self.rt_min_slider.val
        self.my_rt.rt_max = self.rt_max_slider.val
        self.my_rt.rt_peak = self.rt_peak_slider.val
        
        self.rt_min_slider.valinit = self.my_rt.rt_min
        self.rt_max_slider.valinit = self.my_rt.rt_max
        self.rt_peak_slider.valinit = self.my_rt.rt_peak
        
        metob.store(self.my_rt)
        self.min_line.set_xdata((self.my_rt.rt_min,self.my_rt.rt_min))
        self.max_line.set_xdata((self.my_rt.rt_max,self.my_rt.rt_max))
        self.peak_line.set_xdata((self.my_rt.rt_peak,self.my_rt.rt_peak))
        self.fig.canvas.draw_idle()
        
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
        print("loaded file for username = ", self.atlas.username)
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
        print(min_x,max_x)
        
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
        self.ax.set_title(thisline.get_label())        
    
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
    print("loaded file for username = ", atlas_in_data[0].username)
    username = getpass.getuser()
    if username != atlas_in_data[0].username:
        print("YOUR ARE", username, "YOU ARE NOT ALLOWED TO EDIT VALUES THE RT CORRECTOR. USERNAMES ARE NOT THE SAME")
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
                y = y/y_max.next()
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
    print('nrows = ', nRows)
    
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

    print("length of ymax is ", len(y_max))
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
                y = y/y_max.next()
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

def make_output_dataframe(input_fname = '',input_dataset = [],include_lcmsruns = [],exclude_lcmsruns = [], include_groups = [],exclude_groups = [], output_loc = [], fieldname = 'peak_height'):
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

    compound_names = ma_data.get_compound_names(data)[0]
    file_names = ma_data.get_file_names(data)
    group_names = ma_data.get_group_names(data)
    output_loc = os.path.expandvars(output_loc)
    fieldname = fieldname
    
    df = pd.DataFrame( index=compound_names, columns=file_names, dtype=float)

    # peak_height['compound'] = compound_list
    # peak_height.set_index('compound',drop=True)
    for i,dd in enumerate(data):
        for j,d in enumerate(dd):
            if (not d['data']['ms1_summary']) or (not d['data']['ms1_summary'][fieldname]):
                df.ix[compound_names[j],file_names[i]] = 0
            else:
                df.ix[compound_names[j],file_names[i]] = d['data']['ms1_summary'][fieldname]  
    columns = []
    for i,f in enumerate(file_names):
        columns.append((group_names[i],f))
    df.columns = pd.MultiIndex.from_tuples(columns,names=['group', 'file'])

    if output_loc:    
        if not os.path.exists(output_loc):
            os.makedirs(output_loc)
        df.to_csv(os.path.join(output_loc, fieldname + '.tab'),sep='\t')

    return df

def file_with_max_precursor_intensity(data,compound_idx):
    idx = []
    my_max = 0
    for i,d in enumerate(data):
        if 'data' in d[compound_idx]['data']['msms'].keys():
            if type(d[compound_idx]['data']['msms']['data']) != list:#.has_key('precursor_intensity'):
                temp = d[compound_idx]['data']['msms']['data']['precursor_intensity']
                if len(temp)>0:
                    m = max(temp)
                    if m > my_max:
                        my_max = m
                        idx = i
    return idx,my_max

def plot_errorbar_plots(df,output_loc=''):
    
    output_loc = os.path.expandvars(output_loc)
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)
        
    plt.ioff()
    for compound in df.index:
        m = df.ix[compound].groupby(level='group').mean()
        e = df.ix[compound].groupby(level='group').std()
        c = df.ix[compound].groupby(level='group').count()

        for i in range(len(e)):
            if c[i]>0:
                e[i] = e[i] / c[i]**0.5
        
        f, ax = plt.subplots(1, 1,figsize=(12,12))
        m.plot(yerr=e, kind='bar',ax=ax)
        ax.set_title(compound,fontsize=12,weight='bold')
        plt.tight_layout()
        f.savefig(os.path.join(output_loc, compound + '_errorbar.pdf'))

        #f.clear()
        plt.close(f)#f.clear()

def get_reference_msms_spectra(frag_refs, compound_name = '', polarity = '', precursor_mz = 0.0):
    spectra = []
    if polarity ==0:
        polarity = 'negative'
    else:
        polarity = 'positive'
    for fr in frag_refs:
        if (fr.compound[0].name == compound_name) and (fr.frag_references[0].polarity == polarity ):# and (abs(fr.frag_references[0].precursor_mz - precursor_mz)<0.2):
            spectra.append( [(m.mz, m.intensity) for m in fr.frag_references[0].mz_intensities] )
    return spectra

# def get_idenficications_with_fragrefs():
#     """
#     Select all CompoundIdentifications that have a fragmentation reference
#     """
    

def make_identification_figure(input_fname = '',input_dataset = [],include_lcmsruns = [],exclude_lcmsruns = [],include_groups = [],exclude_groups = [], output_loc = []):
    output_loc = os.path.expandvars(output_loc)    
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)
    
    
    if not input_dataset:
        data = ma_data.get_dill_data(os.path.expandvars(input_fname))
    else:
        data = input_dataset
    print(len(data))
    print(len(data[0]))
    # filter runs from the metatlas dataset
    if include_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_include_list(data,'lcmsrun',include_lcmsruns)
    if include_groups:
        data = filter_lcmsruns_in_dataset_by_include_list(data,'group',include_groups)
        
    print(len(data))
    print(len(data[0]))
    if exclude_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_exclude_list(data,'lcmsrun',exclude_lcmsruns)
    if exclude_groups:
        data = filter_lcmsruns_in_dataset_by_exclude_list(data,'lcmsrun',exclude_groups)
        #data = filter_lcmsruns_in_dataset_by_exclude_list(data,'group',exclude_lcmsruns)
    
    print(len(data))
    print(len(data[0]))
    compound_names = ma_data.get_compound_names(data)[0]
    file_names = ma_data.get_file_names(data)

    print('loading preexisting compound identifications')

    ids = metob.retrieve('CompoundIdentification',username='*')
    frag_refs = [cid for cid in ids if cid.frag_references]
    
    print('getting spectra from files')
    
    
    for compound_idx in range(len(compound_names)):
        file_idx, m = file_with_max_precursor_intensity(data,compound_idx)
        if m:
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

            
#             data[file_idx][compound_idx]['identification'].mz_references[0].polarity
#             print data[file_idx][compound_idx]['data']['msms']
            if data[file_idx][compound_idx]['identification'].compound:
                if data[file_idx][compound_idx]['identification'].mz_references[0].detected_polarity == 'positive':
                    my_polarity = 1
                else:
                    my_polarity = 0
                ref_spec = get_reference_msms_spectra(frag_refs, 
                                       compound_name = data[file_idx][compound_idx]['identification'].compound[0].name, 
                                       polarity = my_polarity)
            else:
                ref_spec = []
    #TODO: get the precursor_mz sorted out
            
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
            plt.axhline()
            plt.tight_layout()
            L = plt.ylim()
            plt.ylim(L[0],L[1]*1.12)
            if data[file_idx][compound_idx]['identification'].compound:
                inchi =  data[file_idx][compound_idx]['identification'].compound[0].inchi
                myMol = Chem.MolFromInchi(inchi.encode('utf-8'))
                myMol,neutralised = NeutraliseCharges(myMol)
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
            ax3.set_ylim(0.2,1.01)
            ax3.axis('off')
        #     plt.show()
            fig.savefig(os.path.join(output_loc, compound_names[compound_idx] + '.pdf'))
            #fig.clear()
            plt.cla()
            del fig
            plt.close('all')#f.clear()


    

    
            
def export_atlas_to_spreadsheet(myAtlas, output_filename='', input_type = 'atlas'):
    """
    Return a pandas dataframe containing Atlas info.  Optionally save it.
    This function can also work on a MetAtlas dataset (list of lists returned by get_data_for_atlas_and_groups).
    """
    cols = [c for c in metob.Compound.class_trait_names() if not c.startswith('_')]
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
        atlas_export.loc[i, 'label'] = my_id.name
        atlas_export.loc[i,'rt_min'] = my_id.rt_references[0].rt_min
        atlas_export.loc[i,'rt_max'] = my_id.rt_references[0].rt_max
        atlas_export.loc[i,'rt_peak'] = my_id.rt_references[0].rt_peak
        atlas_export.loc[i,'mz'] = my_id.mz_references[0].mz
        atlas_export.loc[i,'mz_tolerance'] = my_id.mz_references[0].mz_tolerance
        atlas_export.loc[i,'polarity'] = my_id.mz_references[0].detected_polarity
        if my_id.frag_references:
            atlas_export.loc[i,'has_fragmentation_reference'] = True
            # TODO: Gather the frag reference information and export it
        else:
            atlas_export.loc[i,'has_fragmentation_reference'] = False
    
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
            print(i, len(group), myFile)
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
            print(i, len(a.compound_identifications),a.name,  datetime.utcfromtimestamp(a.last_modified))

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
        if not pd.isnull(row.inchi_key):# or type(df.inchi_key[x]) != float:
            if not metob.retrieve('Compounds',inchi_key=row.inchi_key, username = '*'):
                print(row.inchi_key, "compound is not in database. Exiting Without Completing Task!")
                bad_names.append(row.inchi_key)
    return bad_names


def check_file_names(df,field):
    bad_files = []
    for i,row in df.iterrows():
        if row[field] != '':
            if not metob.retrieve('Lcmsruns',name = '%%%s%%'%row[field],username = '*'):
                print(row[field], "file is not in the database. Exiting Without Completing Task!")
                bad_files.append(row[field])
    return bad_files


def get_formatted_atlas_from_google_sheet(polarity='POS',
                                          method='QE_HILIC',
                                          mz_tolerance=10):
    import metatlas.ms_monitor_util as mmu
    df2 = mmu.get_ms_monitor_reference_data()
    #print df.head()
    #df2 = pd.DataFrame(df[1:],columns=df[0])

    fields_to_keep = [ 'name',
                    'label',
                      'inchi_key',
                    'mz_%s'%polarity,
                    'rt_min_%s'%method,
                    'rt_max_%s'%method,
                    'rt_peak_%s'%method,
                    'file_mz_%s_%s'%(method,polarity),
                    'file_rt_%s_%s'%(method,polarity),
                    'file_msms_%s_%s'%(method,polarity)]
    
    fields_there = []
    for f in fields_to_keep:
         if f in df2.keys():
                fields_there.append(f)
    
    df3 = df2.loc[:,fields_there]
    
    df3['mz_tolerance'] = mz_tolerance

    if polarity == 'POS':
        df3['polarity'] = 'positive'
    else:
        df3['polarity'] = 'negative'

    renamed_columns = [c.replace('_%s'%method,'').replace('_%s'%polarity,'') for c in df3.columns]
    for i,c in enumerate(df3.columns):
        df3 = df3.rename(columns = {c:renamed_columns[i]})
    df3 = df3[df3['mz'] != '']

    return df3


def make_atlas_from_spreadsheet(filename='valid atlas file.csv',
                                atlas_name='20161007_MP3umZHILIC_BPB_NEG_ExampleAtlasName',
                                filetype=('excel','csv','tab','dataframe'),
                                sheetname='only for excel type input',
                                polarity = ('positive','negative'),
                                store=False,
                                mz_tolerance=10):
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
    if 'file_msms' in df.keys():
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
            
            if not pd.isnull(row.inchi_key): # this logic is where an identified metabolite has been specified
                c = metob.retrieve('Compounds',neutralized_inchi_key=row.inchi_key,username = '*') #currently, all copies of the molecule are returned.  The 0 is the most recent one. 
                if c:
                    c = c[0]
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
                        mzRef.mz_tolerance = row.mz_threshold    
                
                mzRef.mz_tolerance_units = 'ppm'
                mzRef.detected_polarity = polarity
                #if 'file_mz' in df.keys():
                #    f = metob.retrieve('Lcmsruns',name = '%%%s%%'%df.file_mz[x],username = '*')[0]
                #    mzRef.lcms_run = f
                #     mzRef.adduct = '[M-H]'   
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
                    
                if ('file_msms' in df.keys()) and (c != 'use_label'):
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
                            print('found reference msms spectrum for ',myID.compound[0].name, 'in file',row.file_msms)

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

    print(len(groups))

    if remove_empty:
        groups = filter_empty_metatlas_objects(groups,'items')
    if do_print:
        from datetime import datetime, date
        for i,a in enumerate(groups):
            print(i, a.name,  datetime.utcfromtimestamp(a.last_modified))

    return groups



