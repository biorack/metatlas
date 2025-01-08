import itertools
import logging
import os
import os.path
import warnings
# os.environ['R_LIBS_USER'] = '/project/projectdirs/metatlas/r_pkgs/'
# curr_ld_lib_path = ''

from copy import deepcopy
from pathlib import Path
from typing import Callable, List, Optional, Sequence, Sized, Union, Dict
from enum import Enum

from metatlas.datastructures import metatlas_objects as metob
from metatlas.tools.logging import log_errors
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.io import write_utils
from metatlas.io.metatlas_get_data_helper_fun import extract
from metatlas.plots.compound_eic import save_compound_eic_pdf
from metatlas.plots.tic import save_sample_tic_pdf
from metatlas.plots.utils import pdf_with_text
from metatlas.tools import fastanalysis
from metatlas.tools import parallel
from metatlas.tools import spectralprocessing as sp
from metatlas.tools.notebook import in_papermill

from tqdm.notebook import tqdm
from textwrap import fill, TextWrapper
# import qgrid
import pandas as pd
import dill
import numpy as np
import numpy.typing as npt
import json
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from ipywidgets import interact
import ipywidgets as widgets
from IPython.display import display, HTML
from IPython import get_ipython

import getpass
from typing import Type, Any

from ast import literal_eval
from datetime import datetime

from matplotlib.widgets import Slider, RadioButtons

import matplotlib.patches

import gspread
# from oauth2client.client import SignedJwtAssertionCredentials
from oauth2client.service_account import ServiceAccountCredentials

from metatlas.tools.util import or_default

from functools import reduce
from io import StringIO

from matchms.similarity import CosineHungarian
from matchms.typing import SpectrumType
from matchms import Spectrum

logger = logging.getLogger(__name__)

MetatlasDataset = List[List[Any]]  # avoiding a circular import

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

GUI_FIG_LABEL = 'Annotation GUI'

LOGGING_WIDGET = widgets.Output()

# INSTRUCTIONS_PATH should have a header row and the following columns:
# inchi_key, adduct, chromatography, polarity, note
# any field can be left blank. The notes fields will be concatonated togther
# for all rows where the non-note, non-blank fields match the current context.

INSTRUCTIONS_PATH = '/global/cfs/cdirs/m2650/targeted_analysis/instructions_for_analysts.csv'

class InstructionSet(object):
    def __init__(self, instructions_path):
        try:
            self.data = pd.read_csv(instructions_path, dtype=str, na_values=[], keep_default_na=False)
        except FileNotFoundError:
            logger.warning('Could not find instructions file %s.', instructions_path)
            self.data = pd.DataFrame()

    def query(self, inchi_key, adduct, chromatography, polarity):
        inputs = {"inchi_key": inchi_key, "adduct": adduct, "chromatography": chromatography, "polarity": polarity}
        assert any(v is not None for v in inputs.values())
        conditions = [f"({key}.str.startswith('{value}').values | {key}=='')"
                      for key, value in inputs.items() if value != '']
        out = self.data.query(' & '.join(conditions))['note'].tolist()
        return out if out else ['No instructions for this data']


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

class LayoutPosition(Enum):
    """Define vertical ordering of element in GUI"""

    GUI = 0
    ID_NOTE = 1
    INFO = 2

class adjust_rt_for_selected_compound(object):
    def __init__(self,
                 data,
                 include_lcmsruns=None,
                 exclude_lcmsruns=None,
                 include_groups=None,
                 exclude_groups=None,
                 msms_hits=None,
                 color_me=None,
                 compound_idx=0,
                 width=10,
                 height=4,
                 y_scale='linear',
                 alpha=0.5,
                 min_max_color='green',
                 peak_color='darkviolet',
                 slider_color='ghostwhite',
                 y_max='auto',
                 y_min=0,
                 peak_flags=None,
                 msms_flags=None,
                 adjustable_rt_peak=False):
        """
        INPUTS:
            data: a metatlas_dataset where info on files and compounds are stored
            include/exclude lcmsruns/groups: list of substrings to filter on
            msms_hits: output of get_msms_hits routine
            color_me: '' or list of tuples with (color string, filename substring)
                      that define line colors in the EIC plot
            compound_idx: atlas index number of the first compound to display
            width: width value in inches for the plots and sliders
            height: height value in inches for each of the plots
            y_scale: 'linear' or 'log' for y-axis of EIC plot
            alpha: transparency factor for lines in EIC plot
            min_max_color: matplotlib color string for the sliders and vertical lines
            peak_color: matplotlib color string for the slider and vertical line
            slider_color: matplotlib color string for the background of the sliders
            y_max: y-axis maximum on EIC plot or 'auto' to fix to data
            y_min: y-axis minimum on EIC plot
            peak_flags: list of strings that define radio buttons for EIC plot
                        None or '' gets a default set of buttons
            msms_flags: list of strings that define radio buttons for MSMS plot
                        None or '' gets a default set of buttons
            adjustable_rt_peak: Boolean - is the RT peak slider movable?
        OUTPUTS:
            Writes RT min/max/peak changes to database
            Writes changes to the 'Flag' radio buttons to database

        Key Bindings:
           Next compound: 'l' or right arrow
           Previous compound: 'h' or left arrow
           Next MSMS hit: 'j' or down arrow
           Previous MSMS hit: 'k' or up arrow
           Cycle zoom on MSMS plot: 'z'
           Flag for removal: 'x'
           Toggle display of similar compounds list and highlighting: 's'
        """
        logger.debug("Initializing new instance of %s.", self.__class__.__name__)
        self.data = data
        self.msms_hits = msms_hits.sort_values('score', ascending=False)
        self.color_me = or_default(color_me, [('black', '')])
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

        self.file_names = ma_data.get_file_names(self.data)
        self.configure_flags()
        self.data = filter_runs(self.data, include_lcmsruns, include_groups, exclude_lcmsruns, exclude_groups)

        self.similar_rects = []
        # only the atlas owner can change RT limits or flags
        self.enable_edit = getpass.getuser() == self.data.atlas.username
        self.hit_ctr = 0
        self.msms_zoom_factor = 1
        self.match_idx = None
        self.enable_similar_compounds = True
        self.similar_compounds = []
        self.in_switch_event = True
        self.instruction_set = InstructionSet(INSTRUCTIONS_PATH)
        # native matplotlib key bindings that we want to override
        disable_keyboard_shortcuts({'keymap.yscale': ['l'],
                                    'keymap.xscale': ['k'],
                                    'keymap.save': ['s'],
                                    'keymap.home': ['h']})
        adjust_rt_for_selected_compound.disable()
        self.create_notes_widgets()
        # Turn On interactive plot
        ipy = get_ipython()
        if ipy:  # test suite does not run ipython, so need to bypass
            ipy.magic('matplotlib widget')
        self.layout_figure()
        # create all event handlers
        self.fig.canvas.callbacks.connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect('key_press_event', self.press)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.id_note.observe(self.on_id_note_change, names='value')
        self.set_plot_data()
        # Turn On interactive plot
        plt.ion()

    def set_plot_data(self):
        logger.debug('Starting replot')
        self.enable_similar_compounds = True
        self.similar_compounds = self.get_similar_compounds()
        self.eic_plot()
        self.filter_hits()
        self.msms_plot()
        self.flag_radio_buttons()
        self.notes()
        self.in_switch_event = False
        logger.debug('Finished replot')

    def notes(self):
        self.copy_button_area.clear_output()
        self.id_note.value = self.current_id.identification_notes or ''
        inchi_key = self.current_inchi_key
        if inchi_key is None:
            return
        adduct = self.current_adduct
        polarity = self.data.ids.polarity
        chromatography = self.data.ids.chromatography
        instruction_list = self.instruction_set.query(inchi_key, adduct, chromatography, polarity)
        self.instructions.value = '; '.join(instruction_list)
        clipboard_text = f"{inchi_key}\t{adduct}\t{chromatography}\t{polarity}"
        with self.copy_button_area:
            make_copy_to_clipboard_button(clipboard_text, 'Copy Index')

    def on_id_note_change(self, change):
        if not self.in_switch_event:
            logger.debug('ID_NOTE change handler got value: %s', change['new'])
            self.data.set_note(self.compound_idx, "identification_notes", change["new"])

    def eic_plot(self):
        logger.debug('Starting eic_plot')
        self.ax.set_title('')
        self.ax.set_xlabel('Retention Time')
        # set y-scale and bounds if provided
        self.ax.set_yscale(self.y_scale)
        if self.y_max != 'auto':
            self.ax.set_ylim(self.y_min, self.y_max)
        self.ax.set_ylabel(self.get_ms1_y_axis_label())
        if self.y_scale == 'linear':
            self.ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        self.display_eic_data()
        self.y_max_slider()
        idx = 0 if self.y_scale == 'linear' else 1
        self.lin_log_radio = self.create_radio_buttons(self.lin_log_ax, ('linear', 'log'),
                                                       self.set_lin_log, 0.07, active_idx=idx)
        self.rt_bounds()
        self.highlight_similar_compounds()
        logger.debug('Finished eic_plot')

    def flag_radio_buttons(self):
        if self.current_id.ms1_notes in self.peak_flags:
            peak_flag_index = self.peak_flags.index(self.current_id.ms1_notes)
        else:
            peak_flag_index = 0
        radius = 0.02 * self.gui_scale_factor
        logger.debug('Setting peak flag radio button with index %d', peak_flag_index)
        self.peak_flag_radio = self.create_radio_buttons(self.peak_flag_ax, self.peak_flags,
                                                         self.set_peak_flag,
                                                         radius,
                                                         active_idx=peak_flag_index)
        self.peak_flag_radio.active = self.enable_edit
        if self.current_id.ms2_notes in self.msms_flags:
            msms_flag_index = self.msms_flags.index(self.current_id.ms2_notes)
        else:
            msms_flag_index = 0
        logger.debug('Setting msms flag radio button with index %d', msms_flag_index)
        self.msms_flag_radio = self.create_radio_buttons(self.msms_flag_ax, self.msms_flags,
                                                         self.set_msms_flag,
                                                         radius,
                                                         active_idx=msms_flag_index)
        self.msms_flag_radio.active = self.enable_edit

    def y_max_slider(self):
        (self.slider_y_min, self.slider_y_max) = self.ax.get_ylim()
        self.y_scale_slider = Slider(self.y_scale_ax, '', self.slider_y_min, self.slider_y_max,
                                     valfmt='', valinit=self.slider_y_max, color=self.peak_color,
                                     orientation='vertical')
        self.y_scale_slider.hline.set_linewidth(0)
        self.y_scale_slider.on_changed(self.update_y_scale)

    def rt_bounds(self):
        # put vlines on plot before creating sliders, as adding the vlines may increase plot
        # width, as the vline could occur outside of the data points
        rt = self.data.rts[self.compound_idx]
        self.min_line = self.ax.axvline(rt.rt_min, color=self.min_max_color, linewidth=4.0)
        self.max_line = self.ax.axvline(rt.rt_max, color=self.min_max_color, linewidth=4.0)
        self.peak_line = self.ax.axvline(rt.rt_peak, color=self.peak_color, linewidth=4.0)
        self.rt_min_slider = self.rt_slider(self.rt_min_ax, 'RT min', rt.rt_min,
                                            self.min_max_color, self.update_rt_min)
        self.rt_max_slider = self.rt_slider(self.rt_max_ax, 'RT max', rt.rt_max,
                                            self.min_max_color, self.update_rt_max)
        self.rt_peak_slider = self.rt_slider(self.rt_peak_ax, 'RT peak', rt.rt_peak,
                                             self.peak_color, self.update_rt_peak)
        self.rt_peak_slider.active = self.adjustable_rt_peak and self.enable_edit
        self.rt_min_slider.active = self.enable_edit
        self.rt_max_slider.active = self.enable_edit

    def rt_slider(self, axes, label, valinit, color, on_changed):
        min_x, max_x = self.ax.get_xlim()
        slider = Slider(axes, label, min_x, max_x, valinit=valinit, color=color)
        slider.vline.set_linewidth(0)
        slider.on_changed(on_changed)
        return slider

    def unhighlight_similar_compounds(self):
        [i.remove() for i in self.similar_rects]
        self.similar_rects = []

    def highlight_similar_compounds(self):
        self.unhighlight_similar_compounds()
        min_y, max_y = self.ax.get_ylim()
        min_x, max_x = self.ax.get_xlim()
        height = max_y - min_y
        for compound in self.similar_compounds:
            if compound['index'] == self.compound_idx:
                continue
            width = abs(compound['rt'].rt_max - compound['rt'].rt_min)
            if compound['rt'].rt_min+width < min_x or compound['rt'].rt_min > max_x:
                continue  # would be off the plot
            color = 'red' if compound['overlaps'] else 'blue'
            rect = matplotlib.patches.Rectangle((compound['rt'].rt_min, min_y), width, height,
                                                linewidth=0, alpha=0.12, facecolor=color)
            text_x = max(compound['rt'].rt_min, min_x)
            text_y = max_y - (max_y - min_y)*0.05  # drop just below top of plot
            text = self.ax.text(text_x, text_y, compound['label'], fontsize='small')
            self.ax.add_patch(rect)
            self.similar_rects.append(text)
            self.similar_rects.append(rect)

    def display_eic_data(self):
        for sample in self.data:  # loop through the files
            eic = sample[self.compound_idx]['data']['eic']
            if eic and len(eic['rt']) > 0:
                x = np.asarray(eic['rt'])
                y = np.asarray(eic['intensity'])
                x = x[y > 0]
                y = y[y > 0]  # y[y<0.0] = 0.0
                label = sample[self.compound_idx]['lcmsrun'].name.replace('.mzML', '')
                for i, (color, label_filter) in enumerate(self.color_me):
                    if label_filter in label:
                        zorder = len(self.color_me) + 2 - i
                    else:
                        zorder = 1
                        color = 'black'
                    self.ax.plot(x, y, '-', zorder=zorder, linewidth=2, alpha=self.alpha,
                                 picker=True, pickradius=5, color=color, label=label)

    def configure_flags(self):
        default_peak = ['keep', 'remove', 'unresolvable isomers', 'poor peak shape']
        default_msms = ['no selection',
                        '-1, bad match - should remove compound',
                        '0, no ref match available or no MSMS collected',
                        '0.5, co-isolated precursor, partial match',
                        '0.5, partial match of fragments',
                        '1, perfect match to internal reference library',
                        '1, perfect match to external reference library',
                        '1, co-isolated precursor but all reference ions are in sample spectrum']
        if self.peak_flags is None or self.peak_flags == '':
            self.peak_flags = default_peak
        if self.msms_flags is None or self.msms_flags == '':
            self.msms_flags = default_msms

    def get_ms1_y_axis_label(self):
        if self.current_id.name:
            compound_name = self.current_id.name.split('///')[0]
        elif self.current_id.compound[-1].name:
            compound_name = self.current_id.compound[-1].name
        else:
            compound_name = 'nameless compound'
        try:
            adduct = self.current_adduct
        except (KeyError, AttributeError):
            return '%d, %s' % (self.compound_idx, compound_name)
        return '%d, %s\n%s' % (self.compound_idx, compound_name, adduct)

    def filter_hits(self):
        inchi_key = extract(self.current_id, ['compound', -1, 'inchi_key'], None)
        hits_mz_tolerance = self.current_id.mz_references[-1].mz_tolerance*1e-6
        mz_theoretical = self.current_id.mz_references[0].mz
        my_scan_rt = self.msms_hits.index.get_level_values('msms_scan')
        filtered = self.msms_hits[(my_scan_rt >= float(self.data.rts[self.compound_idx].rt_min)) &
                                  (my_scan_rt <= float(self.data.rts[self.compound_idx].rt_max)) &
                                  within_tolerance(self.msms_hits['measured_precursor_mz'],
                                                   mz_theoretical, hits_mz_tolerance)]

        self.hits = filtered if inchi_key is None else filtered[(filtered['inchi_key'] == inchi_key)]

    def msms_plot(self, font_scale=10.0):
        logger.debug('Starting msms_plot')
        compound = None
        hit_file_name = None
        if not self.hits.empty:
            hit_file_name, compound = get_hit_metadata(self.data, self.hits, self.file_names,
                                                       self.hit_ctr, self.compound_idx)
        mz_header, rt_header, cpd_header = get_msms_plot_headers(self.data, self.hits, self.hit_ctr,
                                                                 self.compound_idx, self.similar_compounds)
        cpd_header_wrap = fill(cpd_header, width=int(self.width * font_scale))  # text wrap
        hit_ref_id, hit_score, hit_query, hit_ref = get_msms_plot_data(self.hits, self.hit_ctr)
        self.ax2.cla()
        self.ax2.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        self.mz_lines = plot_msms_comparison2(0, mz_header, rt_header, cpd_header_wrap, hit_ref_id,
                                              hit_file_name, hit_score, self.ax2, hit_query,
                                              hit_ref, self.msms_zoom_factor)
        min_x = self.ax2.get_xlim()[0]  # fails if original location is not within plot
        self.mz_annot = self.ax2.annotate('', xy=(min_x, 0), visible=False)
        logger.debug('Finished msms_plot')

    def create_notes_widgets(self):
        wide_layout = widgets.Layout(width="85%")
        self.instructions = widgets.HTML(value="Compound Info will go here", layout=wide_layout)
        self.copy_button_area = widgets.Output()
        self.id_note = widgets.Textarea(
            description="ID Notes", value="", placeholder="No note entered", continuous_update=False,
            layout=wide_layout
        )
        display(widgets.HBox([self.instructions, self.copy_button_area],
                layout=widgets.Layout(justify_content='space-between')))
        display(self.id_note)

    def layout_figure(self):
        self.gui_scale_factor = self.height/3.25 if self.height < 3.25 else 1
        base_font_size = 10
        y_slider_width = 0.01

        plt.rcParams.update({'font.size': base_font_size * self.gui_scale_factor})

        max_radio_label_len = max([len(x) for x in self.peak_flags+self.msms_flags])

        self.plot_hspace = 2.2/(self.height*2)
        self.plot_left_pos = 1.0/self.width
        # not using a fixed width font, so this assume a even distribution of character widths
        # in the radio button labels
        self.plot_right_pos = 1.0-(0.7/self.width +
                                   max_radio_label_len*0.075/self.width*self.gui_scale_factor)
        self.plot_top_pos = 1.0-(0.396/(self.height*(2+self.plot_hspace)))
        self.plot_bottom_pos = 1.32/(self.height*(2+self.plot_hspace))

        # create figure and first axes
        combined_plot_height = self.height * (2 + self.plot_hspace)
        self.fig, (self.ax2, self.ax) = plt.subplots(2, 1, num=GUI_FIG_LABEL,
                                                     figsize=(self.width, combined_plot_height))
        plt.subplots_adjust(left=self.plot_left_pos, right=self.plot_right_pos,
                            bottom=self.plot_bottom_pos, top=self.plot_top_pos,
                            hspace=self.plot_hspace)

        y_axis_height = (self.plot_top_pos - self.plot_bottom_pos) * \
                        (1-self.plot_hspace/(2+self.plot_hspace))/2
        self.layout_rt_sliders()
        self.layout_y_scale_slider(y_slider_width, y_axis_height)
        self.layout_radio_buttons(y_slider_width, y_axis_height)

    def layout_y_scale_slider(self, y_slider_width, y_axis_height):
        self.y_scale_ax = plt.axes([self.plot_right_pos, self.plot_bottom_pos,
                                    y_slider_width, y_axis_height],
                                   facecolor=self.slider_color)

    def layout_radio_buttons(self, y_slider_width, y_axis_height):
        self.radio_button_radius = 0.02 * self.gui_scale_factor
        radio_button_axes_width = 1-self.plot_right_pos
        lin_log_height_fraction = 0.3
        self.lin_log_ax = layout_radio_button_set([self.plot_left_pos,
                                                   self.plot_bottom_pos + y_axis_height*(1-lin_log_height_fraction),
                                                   radio_button_axes_width,
                                                   y_axis_height*lin_log_height_fraction],
                                                  anchor='NW')
        self.peak_flag_ax = layout_radio_button_set([self.plot_right_pos + y_slider_width,
                                                     self.plot_bottom_pos,
                                                     radio_button_axes_width,
                                                     y_axis_height])
        self.msms_flag_ax = layout_radio_button_set([self.plot_right_pos,
                                                     self.plot_top_pos - y_axis_height,
                                                     radio_button_axes_width,
                                                     y_axis_height])

    def layout_rt_sliders(self):
        x_axis_label_height = 0.5
        combined_plot_height = self.height * (2 + self.plot_hspace)
        rt_slider_height = (self.plot_bottom_pos-x_axis_label_height/combined_plot_height)/4.0
        rt_slider_width = self.plot_right_pos - self.plot_left_pos

        self.rt_peak_ax = plt.axes([self.plot_left_pos, 0, rt_slider_width, rt_slider_height],
                                   facecolor=self.slider_color)
        self.rt_max_ax = plt.axes([self.plot_left_pos, rt_slider_height*1.5,
                                   rt_slider_width, rt_slider_height],
                                  facecolor=self.slider_color)
        self.rt_min_ax = plt.axes([self.plot_left_pos, rt_slider_height*3.0,
                                   rt_slider_width, rt_slider_height],
                                  facecolor=self.slider_color)

    def create_radio_buttons(self, axes, labels, on_click_handler, radius, active_idx=0):
        buttons = RadioButtons(axes, labels, active=active_idx)
        buttons.set_radio_props({'s':1600*radius, 'facecolor':'blue'})
        buttons.on_clicked(on_click_handler)
        return buttons

    def set_lin_log(self, label):
        logger.debug('Y-scale of EIC plot set to %s scale.', label)
        self.ax.set_yscale(label)
        if label == 'linear':
            self.ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        self.fig.canvas.draw_idle()

    def set_flag(self, name, value):
        logger.debug('Setting flag "%s" to "%s".', name, value)
        self.data.set_note(self.compound_idx, name, value)

    def set_peak_flag(self, label):
        old_label = ma_data.extract(self.current_id, ["ms1_notes"])
        self.set_flag('ms1_notes', label)
        if not self.enable_similar_compounds:
            return
        was_remove = is_remove(old_label)
        now_remove = is_remove(label)
        if (not was_remove) and now_remove:
            self.unhighlight_similar_compounds()
        if was_remove and (not now_remove):
            self.highlight_similar_compounds()

    def set_msms_flag(self, label):
        self.set_flag('ms2_notes', label)

    @log_errors(output_context=LOGGING_WIDGET)
    def on_pick(self, event):
        thisline = event.artist
        thisline.set_color('cyan')
        label = thisline.get_label()
        self.ax.set_title(label, fontsize=7)
        logger.debug("Sample %s selected on EIC plot via mouse click event.", label)

    @log_errors(output_context=LOGGING_WIDGET)
    def on_motion(self, event):
        if event.inaxes == self.ax2:  # in msms mirror plot
            for collection in self.mz_lines:
                found, ind = collection.contains(event)
                if found:
                    segments = collection.get_segments()
                    vertice = segments[ind["ind"][0]]
                    mz = vertice[0][0]
                    self.mz_annot.set_text(f"{mz:.5f}")
                    self.mz_annot.xyann = (mz, event.ydata)
                    self.mz_annot.set_visible(True)
                    self.fig.canvas.draw_idle()
                    return
            if self.mz_annot.get_visible():
                self.mz_annot.set_visible(False)
                self.fig.canvas.draw_idle()

    def hide_radio_buttons(self):
        self.peak_flag_radio.set_radio_props({'facecolor':'white'})
        self.lin_log_radio.set_radio_props({'facecolor':'white'})
        self.msms_flag_radio.set_radio_props({'facecolor':'white'})
        for label in self.peak_flag_radio.labels:
            label.set_color('white')
        for label in self.lin_log_radio.labels:
            label.set_color('white')
        for label in self.msms_flag_radio.labels:
            label.set_color('white')

    def update_plots(self):
        self.msms_zoom_factor = 1
        self.ax.cla()
        self.ax2.cla()
        self.rt_peak_ax.cla()
        self.rt_min_ax.cla()
        self.rt_max_ax.cla()
        self.y_scale_ax.cla()
        self.hide_radio_buttons()
        self.set_plot_data()

    @log_errors(output_context=LOGGING_WIDGET)
    def press(self, event):
        if event.key in ['right', 'l']:
            if self.compound_idx + 1 < len(self.data[0]):
                self.in_switch_event = True
                self.compound_idx += 1
                logger.debug("Increasing compound_idx to %d (inchi_key:%s adduct:%s).",
                             self.compound_idx,
                             self.current_inchi_key,
                             self.current_adduct
                             )
                self.hit_ctr = 0
                self.match_idx = None
                self.update_plots()
                self.in_switch_event = False
        elif event.key in ['left', 'h']:
            if self.compound_idx > 0:
                self.in_switch_event = True
                self.compound_idx -= 1
                logger.debug("Decreasing compound_idx to %d (inchi_key:%s adduct:%s).",
                             self.compound_idx,
                             self.current_inchi_key,
                             self.current_adduct
                             )
                self.hit_ctr = 0
                self.match_idx = None
                self.update_plots()
                self.in_switch_event = False
        elif event.key in ['up', 'k']:
            if self.hit_ctr > 0:
                self.hit_ctr -= 1
                logger.debug("Decreasing hit_ctr to %d.", self.hit_ctr)
                self.update_plots()
        elif event.key in ['down', 'j']:
            if self.hit_ctr < len(self.hits) - 1:
                logger.debug("Increasing hit_ctr to %d.", self.hit_ctr)
                self.hit_ctr += 1
                self.update_plots()
        elif event.key == 'x':
            if not self.enable_edit:
                self.warn_if_not_atlas_owner()
                return
            logger.debug("Removing compound %d via 'x' key event.", self.compound_idx)
            self.peak_flag_radio.set_active(1)
        elif event.key == 'z':
            self.msms_zoom_factor = 1 if self.msms_zoom_factor == 25 else self.msms_zoom_factor * 5
            logger.debug("Setting msms zoom factor to %d.", self.msms_zoom_factor)
            self.msms_plot()
        elif event.key == 's':
            self.enable_similar_compounds = not self.enable_similar_compounds
            self.similar_compounds = self.get_similar_compounds()
            if self.enable_similar_compounds:
                logger.debug("Enabling similar compounds list and EIC plot highlighting.")
                self.highlight_similar_compounds()
            else:
                logger.debug("Removing similar compounds list and EIC plot hightlighting.")
                self.unhighlight_similar_compounds()
            self.msms_plot()
        elif event.key == 'm':
            num_sim = len(self.similar_compounds)
            if num_sim > 0:
                self.match_idx = 0 if self.match_idx is None else (self.match_idx + 1) % num_sim
                self.match_rts()

    def match_rts(self):
        """Sets RT min and max to match similar compound referenced by match_idx"""
        source = self.similar_compounds[self.match_idx]['rt']
        logger.debug("Matching RT bounds to index %s with min %d, max %d.",
                     self.match_idx, source.rt_min, source.rt_max)
        self.update_rt('rt_min', source.rt_min)
        self.rt_min_slider.set_val(source.rt_min)
        self.update_rt('rt_max', source.rt_max)
        self.rt_max_slider.set_val(source.rt_max)

    def update_y_scale(self, val):
        if self.slider_y_min < 0:
            self.slider_y_min = -0.2 * val
        else:
            self.slider_y_min = 0.02
        self.ax.set_ylim(self.slider_y_min, val)
        self.fig.canvas.draw_idle()

    def warn_if_not_atlas_owner(self):
        user = getpass.getuser()
        if user != self.data.atlas.username:
            text = ("YOU ARE %s. YOU ARE NOT THE ATLAS OWNER."
                    "YOU ARE NOT ALLOWED TO EDIT VALUES WITH THE RT CORRECTOR.")
            self.ax.set_title(text % user)
            logger.warning(text, user)

    def update_rt(self, which, val):
        """
        inputs:
            which: 'rt_min', 'rt_max', or 'rt_peak'
            val: new RT value
        """
        logger.debug("Updating %s to %0.4f", which, val)
        slider = {'rt_min': self.rt_min_slider, 'rt_peak': self.rt_peak_slider,
                  'rt_max': self.rt_max_slider}
        line = {'rt_min': self.min_line, 'rt_peak': self.peak_line, 'rt_max': self.max_line}
        self.data.set_rt(self.compound_idx, which, val)
        slider[which].valinit = val
        line[which].set_xdata((val, val))
        if which != 'rt_peak':
            self.msms_zoom_factor = 1
            self.filter_hits()
            self.similar_compounds = self.get_similar_compounds()
            self.highlight_similar_compounds()
            self.msms_plot()
        self.fig.canvas.draw_idle()

    def update_rt_min(self, val):
        self.update_rt('rt_min', val)

    def update_rt_max(self, val):
        self.update_rt('rt_max', val)

    def update_rt_peak(self, val):
        self.update_rt('rt_peak', val)

    def get_similar_compounds(self, use_labels=True):
        """
        inputs:
            use_labels: if True use compound labels in output instead of compound names
        returns:
            list of dicts containing information on compounds with similar mz values or mono
                 isotopic MW when compared to self.compound_idx
            each dict contains:
                index: position in self.data[0]
                label: compound name or label string
                rt: a metatlas.datastructures.metatlas_objects.RtReference
                overlaps: True if compound has RT bounds overlapping with those of self.compound_idx
        """
        if not self.enable_similar_compounds:
            return []
        cid = self.current_id
        if len(cid.compound) == 0 or is_remove(ma_data.extract(cid, ["ms1_notes"])):
            return []
        out = []
        cid_mz_ref = cid.mz_references[0].mz
        cid_mass = cid.compound[0].mono_isotopic_molecular_weight
        if cid.compound[0].inchi_key != '':
            try:
                cid_inchikey_prefix = cid.compound[0].inchi_key.split('-')[0]
            except:
                logger.warning("Inchi key is anomalous for compound %d.", self.compound_idx)
                cid_inchikey_prefix = None
        else:
            logger.warning("Compound %d does not have inchi key.", self.compound_idx)
            cid_inchikey_prefix = None
        for compound_iter_idx, _ in enumerate(self.data[0]):
            cpd_iter_id = self.data[0][compound_iter_idx]['identification']
            if len(cpd_iter_id.compound) == 0 or is_remove(ma_data.extract(cpd_iter_id, ["ms1_notes"])):
                continue
            mass = cpd_iter_id.compound[0].mono_isotopic_molecular_weight
            mz_ref = cpd_iter_id.mz_references[0].mz
            try:
                inchikey_prefix = cpd_iter_id.compound[0].inchi_key.split('-')[0]
            except:
                logger.warning("Inchi key is anomalous for compound %d.", cpd_iter_id.compound[0].name)
                cid_inchikey_prefix = ''
            if (mz_ref-0.005 <= cid_mz_ref <= mz_ref+0.005) or \
                (mass-0.005 <= cid_mass <= mass+0.005) or \
                (cid_inchikey_prefix == inchikey_prefix):
                out.append({'index': compound_iter_idx,
                            'label': cpd_iter_id.name if use_labels else cpd_iter_id.compound[0].name,
                            'rt': self.data.rts[compound_iter_idx],
                            'overlaps': rt_range_overlaps(self.data.rts[self.compound_idx],
                                                          self.data.rts[compound_iter_idx])})
                logger.debug("Adding similar compound with index %d and min %d, max %d.",
                             out[-1]["index"], out[-1]["rt"].rt_min, out[-1]["rt"].rt_max)
        return out

    @property
    def current_id(self):
        return self.data[0][self.compound_idx]['identification']

    @property
    def current_inchi_key(self):
        return extract(self.current_id, ['compound', 0, 'inchi_key'], None)

    @property
    def current_adduct(self):
        return self.current_id.mz_references[0].adduct

    @staticmethod
    def disable():
        """Stops the GUI from being updated and interfering with the generation of output figures"""
        plt.close(GUI_FIG_LABEL)


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
        self.data = filter_runs(data, include_lcmsruns, include_groups, exclude_lcmsruns, exclude_groups)
        self.match_idx = None

        # create figure and first axes
        self.fig,self.ax = plt.subplots(figsize=(width, height))
        plt.subplots_adjust(left=0.09, bottom=0.275)
#         plt.ticklabel_format(style='plain', axis='x')
#         plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        # warn the user if they do not own the atlas; and can not edit its values
        self.enable_edit = True
        self.atlas = metob.retrieve('Atlas',unique_id = self.data[0][0]['atlas_unique_id'],username='*')[-1]
        logger.info("loaded file for username = %s", self.atlas.username)
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
    data = filter_runs(data, include_lcmsruns, include_lcmsruns, exclude_lcmsruns, exclude_lcmsruns)
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
    y_max = itertools.cycle(y_max)

    # create ouput dir
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)

    disable_interactive_plots()
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
    data = filter_runs(data, include_lcmsruns, include_lcmsruns, exclude_lcmsruns, exclude_lcmsruns)

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
    y_max = itertools.cycle(y_max)

    # create ouput dir
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)
    disable_interactive_plots()
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
    dataset_for_median = filter_runs(dataset_for_median, include_lcmsruns, include_groups,
                                     exclude_lcmsruns, exclude_groups)
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


def filter_runs(data, include_lcmsruns=None, include_groups=None, exclude_lcmsruns=None, exclude_groups=None):
    """filter runs from the metatlas dataset"""
    if include_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_include_list(data, 'lcmsrun', include_lcmsruns)
    if include_groups:
        data = filter_lcmsruns_in_dataset_by_include_list(data, 'group', include_groups)
    if exclude_lcmsruns:
        data = filter_lcmsruns_in_dataset_by_exclude_list(data, 'lcmsrun', exclude_lcmsruns)
    if exclude_groups:
        data = filter_lcmsruns_in_dataset_by_exclude_list(data, 'group', exclude_groups)
    return data


def make_output_dataframe(input_fname: Optional[Path] = None, input_dataset=None, include_lcmsruns=None, exclude_lcmsruns=None, include_groups=None, exclude_groups=None, output_loc: Optional[Path] = None, fieldname='peak_height', use_labels=False, short_names_df=None, summarize=False, polarity='', overwrite=True):
    """
    fieldname can be: peak_height, peak_area, mz_centroid, rt_centroid, mz_peak, rt_peak
    """
    assert input_dataset or input_fname
    if input_dataset:
        full_data = input_dataset
    elif input_fname:
        full_data = ma_data.get_dill_data(os.path.expandvars(input_fname))
    data = filter_runs(full_data, include_lcmsruns, include_groups, exclude_lcmsruns, exclude_groups)
    compound_names = ma_data.get_compound_names(data, use_labels=use_labels)[0]
    file_names = ma_data.get_file_names(data)
    group_names = ma_data.get_group_names(data)
    group_shortnames = ma_data.get_group_shortnames(data)
    out = pd.DataFrame(index=compound_names, columns=file_names, dtype=float)

    for i, sample in enumerate(data):
        for j, compound in enumerate(sample):
            ids = ['data', 'ms1_summary', fieldname]
            out.loc[compound_names[j], file_names[i]] = ma_data.extract(compound, ids, 0)
    columns = []
    if short_names_df is None:
        short_names_df = pd.DataFrame()
        for i, name in enumerate(file_names):
            columns.append((group_names[i], name))
        out.columns = pd.MultiIndex.from_tuples(columns, names=['group', 'file'])
    else:
        for i, name in enumerate(file_names):
            temp = [group_names[i], name, group_shortnames[i]]
            temp.extend(short_names_df.loc[name.split('.')[0]].values.tolist())
            columns.append(tuple(temp))
        out.columns = pd.MultiIndex.from_tuples(columns, names=['group', 'file', 'short groupname', 'sample treatment', 'short filename', 'short samplename'])
    out = out.reindex(sorted(out.columns), axis=1)
    if summarize:
        out.columns = out.columns.droplevel()
        out = append_stats_columns(out)
    if output_loc:
        output_loc = Path(os.path.expandvars(output_loc))
        prefix = f"{polarity}_" if polarity != '' else ''
        df_path = output_loc / f"{prefix}{fieldname}.tab"
        write_utils.export_dataframe_die_on_diff(out, df_path, fieldname, overwrite=overwrite, sep="\t", float_format="%.9e")
    return out


def append_stats_columns(in_df):
    stats = pd.DataFrame(dtype=float)
    stats['mean'] = in_df.mean(numeric_only=True, axis=1)
    stats['median'] = in_df.median(numeric_only=True, axis=1)
    stats['min'] = in_df.min(numeric_only=True, axis=1)
    stats['max'] = in_df.max(numeric_only=True, axis=1)
    stats['standard deviation'] = in_df.std(numeric_only=True, axis=1)
    stats['standard error'] = in_df.sem(numeric_only=True, axis=1)
    stats['#NaNs'] = in_df.isin(['NaN']).sum(axis=1)
    return pd.concat([in_df, stats], axis=1)


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
    return idx, my_max

def file_with_max_ms1_intensity(data, compound_idx, limit_to_rt_range=False):
    file_idx_max = None
    value_max = 0
    for file_idx, sample in enumerate(data):
        try:
            temp = get_ms1_df(sample[compound_idx], limit_to_rt_range)['intensity'].max()
            if temp > value_max:
                value_max = temp
                file_idx_max = file_idx
        except KeyError:
            pass
    return file_idx_max, value_max

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


            cos = CosineHungarian(tolerance=0.005)

            for i, msv_sample in enumerate(np.split(msv_sample_scans, scan_idxs, axis=1)):

                mms_sample = Spectrum(mz=msv_sample[0], intensities=np.sqrt(msv_sample[1]), metadata={'precursor_mz':np.nan})

                for f, frag in sp.filter_frag_refs(data, frag_refs, compound_idx, file_idx, filter_by).iterrows():
                    msv_ref = sp.sort_ms_vector_by_mz(np.array(frag['mz_intensities']).T)

                    mms_ref = Spectrum(mz=msv_ref[0], intensities=np.sqrt(msv_ref[1]), metadata={'precursor_mz':np.nan})

                    score = cos.pair(mms_sample, mms_ref)['score'].item()

                    if score > max_score or np.isnan(max_score):
                        max_score = score
                        idx = file_idx
                        best_ref_spec = [frag['mz_intensities']]

    return idx, max_score, best_ref_spec

def plot_errorbar_plots(df,output_loc='', use_shortnames=True, ylabel=""):

    output_loc = os.path.expandvars(output_loc)
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)

    disable_interactive_plots()
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


def make_boxplot_plots(df: pd.DataFrame, output_loc: Path, use_shortnames: bool = True, ylabel: str = "",
                       overwrite: bool = True, max_cpus: int = 1, logy: bool = False) -> None:
    output_loc = Path(os.path.expandvars(output_loc))
    logger.info('Exporting box plots of %s to %s.', ylabel, output_loc)
    disable_interactive_plots()
    args = [(compound, df, output_loc, use_shortnames, ylabel, overwrite, logy) for compound in df.index]
    parallel.parallel_process(_make_boxplot_single_arg, args, max_cpus, unit='plot')

def _make_boxplot_single_arg(arg_list):
    """ this is a hack, but multiprocessing constrains the functions that can be passed """
    make_boxplot(*arg_list)

def make_boxplot(compound: int, df: pd.DataFrame, output_loc: Path, use_shortnames: bool, ylabel: str, overwrite: bool, logy: bool) -> None:
    fig_path = output_loc / f"{compound}{'_log' if logy else ''}_boxplot.pdf"
    write_utils.check_existing_file(fig_path, overwrite)
    level = 'short groupname' if use_shortnames and 'short groupname' in df.columns.names else 'group'
    num_points = 0

    # Reorder columns for boxplots
    istd_cols = [col for col in df.columns if 'istd' in col[-1].lower()]
    exctrl_cols = [col for col in df.columns if 'exctrl' in col[-1].lower() or 'txctrl' in col[-1].lower()]
    refstd_cols = [col for col in df.columns if 'refstd' in col[-1].lower()]
    sample_cols = sorted([col for col in df.columns if col not in istd_cols and col not in exctrl_cols and col not in refstd_cols])
    df_sorted = df[istd_cols + exctrl_cols + sample_cols + refstd_cols]
    g_sorted = df_sorted.loc[compound].groupby(level=level, sort=False)

    plt.rcParams.update({'font.size': 12})
    f, ax = plt.subplots(1, 1, figsize=(max(len(g_sorted)*0.5, 12), 12))

    g_sorted.apply(pd.DataFrame).plot(kind='box', ax=ax)
    for i, (n, grp) in enumerate(g_sorted):
        x = [i+1] *len(grp)
        x = np.random.normal(x, 0.04, size=len(x))
        plt.scatter(x, grp)
        num_points += np.sum(~np.isnan(grp))
    if num_points == 0:
        logger.warning('Zero data points in box plot of %s for %s.', ylabel, compound)
        pdf_with_text("Molecule not detected", fig_path)
        return
    ax.set_title(compound,fontsize=12,weight='bold')
    plt.xticks(rotation=90)
    if logy:
        plt.yscale('log')
    if ylabel != "":
        plt.ylabel(ylabel)
    plt.tight_layout()
    f.savefig(fig_path)
    plt.close(f)
    plt.rcdefaults()  # reset font size to default value
    logger.debug('Exported box plot of %s for %s at %s.', ylabel, compound, fig_path)


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
    data = filter_runs(data, include_lcmsruns, include_groups, exclude_lcmsruns, exclude_groups)
    compound_names = ma_data.get_compound_names(data,use_labels=use_labels)[0]
    file_names = ma_data.get_file_names(data)
    # print(len(data),len(data[0]),len(compound_names))


    frag_refs = pd.read_json(os.path.join(frag_json_dir, frag_json_name + ".json"))
    disable_interactive_plots()
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

        cos = CosineHungarian(tolerance=0.005)

        for i, msv_sample in enumerate(np.split(msv_sample_scans, scan_idxs, axis=1)):
            current_best_score = None
            current_best_ref_idx = None
            current_best_msv_sample = None
            current_best_msv_ref = None
            current_best_rt = None

            mms_sample = Spectrum(mz=msv_sample[0], intensities=np.sqrt(msv_sample[1]), metadata={'precursor_mz':np.nan})

            for ref_idx, frag in sp.filter_frag_refs(data, frag_refs, compound_idx, file_idx, filter_by).iterrows():
                msv_ref = np.array(frag['mz_intensities']).T

                msv_sample_aligned, msv_ref_aligned = sp.pairwise_align_ms_vectors(msv_sample, msv_ref, 0.005)

                mms_ref = Spectrum(mz=msv_ref[0], intensities=np.sqrt(msv_ref[1]), metadata={'precursor_mz':np.nan})

                score = cos.pair(mms_sample, mms_ref)['score'].item()

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

def plot_msms_comparison2(i, mz_header, rt, cpd_header, ref_id, filename, score, ax, msv_sample, msv_ref, zoom_factor=1):
    pickradius = 10
    msv_sample_matches, msv_ref_matches, msv_sample_nonmatches, msv_ref_nonmatches = sp.partition_aligned_ms_vectors(msv_sample, msv_ref)

    msv_sample_unaligned = np.concatenate((msv_sample_matches, msv_sample_nonmatches), axis=1)
    msv_ref_unaligned = np.concatenate((msv_ref_matches, msv_ref_nonmatches), axis=1)

    sample_mz = msv_sample_nonmatches[0]
    sample_zeros = np.zeros(msv_sample_nonmatches[0].shape)
    sample_intensity = msv_sample_nonmatches[1]

    lines = [ax.vlines(sample_mz, sample_zeros, sample_intensity, colors='r', linewidth=1, pickradius=pickradius)]

    shared_mz = msv_sample_matches[0]
    shared_zeros = np.zeros(msv_sample_matches[0].shape)
    shared_sample_intensity = msv_sample_matches[1]

    lines.append(ax.vlines(shared_mz, shared_zeros, shared_sample_intensity, colors='g', linewidth=1, pickradius=pickradius))

    most_intense_idxs = np.argsort(msv_sample_unaligned[1])[::-1]

    if i == 0:
        ax.set_title('MSMS ref ID = %s\n%s' % (ref_id, filename), fontsize='small', fontstretch='condensed')
        if cpd_header == "":
            ax.set_xlabel('m/z\nScore = %.4f, %s\n%s' % (score, rt, mz_header), weight='bold', fontsize=7)
        else:
            ax.set_xlabel('m/z\nScore = %.4f, %s\n%s\n%s' % (score, rt,  mz_header, cpd_header), weight='bold', fontsize=7)

        ax.set_ylabel('intensity')

        labels = [1.001e9]

        intensity_requirement = [m for m in most_intense_idxs
                                 if
                                 np.min(np.abs(msv_sample_unaligned[0][m] - labels)) > 0.1
                                 and msv_sample_unaligned[1][m] > 0.2 * np.max(msv_sample_unaligned[1])]
        for m in max([most_intense_idxs[:6], intensity_requirement], key=len):
            if np.min(np.abs(msv_sample_unaligned[0][m] - labels)) > 0.1 and msv_sample_unaligned[1][m] > 0.02 * np.max(msv_sample_unaligned[1]):
                ax.annotate('%5.4f' % msv_sample_unaligned[0][m],
                            xy=(msv_sample_unaligned[0][m], msv_sample_unaligned[1][m]),
                            size=6)
                labels.append(msv_sample_unaligned[0][m])
    if msv_ref_unaligned[0].size > 0:
        ref_scale = -1 * np.max(msv_sample_unaligned[1]) / np.max(msv_ref_unaligned[1])
        ref_mz = msv_ref_nonmatches[0]
        ref_zeros = np.zeros(msv_ref_nonmatches[0].shape)
        ref_intensity = ref_scale * msv_ref_nonmatches[1]
        shared_ref_intensity = ref_scale * msv_ref_matches[1]
        lines.append(ax.vlines(ref_mz, ref_zeros, ref_intensity, colors='r', linewidth=1, pickradius=pickradius))
        lines.append(ax.vlines(shared_mz, shared_zeros, shared_ref_intensity, colors='g', linewidth=1, pickradius=pickradius))
        ax.axhline()
    ylim = ax.get_ylim()
    new_ylim = 0 if ylim[0] == 0 else ylim[0]/zoom_factor
    ax.set_ylim(new_ylim, ylim[1]/zoom_factor)
    return lines


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
            ax.fill_between(x,0,y,myWhere, facecolor='c', alpha=.25)
        except (AssertionError, TypeError):
            pass

    # ax.tick_params(labelbottom='off')
    ax.xaxis.set_tick_params(labelsize=5)
    ax.yaxis.set_tick_params(labelsize=5)
    ax.tick_params(axis='y', labelsize=5)
    ax.get_yaxis().get_major_formatter().set_useOffset(True)
    #ax.get_yaxis().set_visible(False)
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

def build_msms_refs_spectra(msms_refs: pd.DataFrame) -> pd.DataFrame:
    """Converts MS2 spectral arrays from strings to numpy arrays and creates MatchMS spectra."""

    msms_refs['spectrum'] = msms_refs['spectrum'].apply(lambda s: np.array(json.loads(s)))
    msms_refs['matchms_spectrum'] = msms_refs.apply(lambda s: Spectrum(s.spectrum[0], np.sqrt(s.spectrum[1]), metadata={'precursor_mz':s.precursor_mz}), axis=1)

    return msms_refs

def load_msms_refs_file(refs_path: str, pre_query: str, query: str,
                        ref_dtypes: Dict[str, Type[Any]], ref_index: str | List[str]) -> pd.DataFrame:
    """Load MSMS refs from file path."""

    if ref_dtypes is None:
        ref_dtypes = {'database': str, 'id': str, 'name': str,
                          'spectrum': object, 'decimal': int, 'precursor_mz': float,
                          'polarity': str, 'adduct': str, 'fragmentation_method': str,
                          'collision_energy': str, 'instrument': str, 'instrument_type': str,
                          'formula': str, 'exact_mass': float,
                          'inchi_key': str, 'inchi': str, 'smiles': str}

    msms_refs = pd.read_csv(refs_path, sep='\t', dtype=ref_dtypes)
    msms_refs = msms_refs.query(pre_query)
    if query is not None:
        msms_refs = msms_refs.query(query)
    if ref_index is not None:
        msms_refs.set_index(ref_index)

    return msms_refs

def convert_to_centroid(sample_df):

    max_peaks, _ = sp.peakdet(sample_df[1], 1000.0)
    if max_peaks.shape[0] > 0:
        idx = max_peaks[:, 0].astype(int).flatten()
        return sample_df[:, idx]
    return np.zeros((0, 0))

def create_nonmatched_msms_hits(msms_data: pd.DataFrame, inchi_key: str, compound_idx: int) -> pd.DataFrame:
    """
    Create MSMS hits DataFrame with no corresponding match in MSMS refs.

    If there is no InChi Key in the atlas (compound identifications) for a particular compound,
    Setting inchi_key to an empty string is necessary for correct plotting.
    """

    compound_msms_hits = msms_data[msms_data['compound_idx'] == compound_idx].reset_index(drop=True).copy()

    if inchi_key is None:
        compound_msms_hits['inchi_key'] = ''

    compound_msms_hits['database'] = np.nan
    compound_msms_hits['id'] = np.nan
    compound_msms_hits['adduct'] = ''

    compound_msms_hits['score'] = compound_msms_hits['measured_precursor_intensity']
    compound_msms_hits['num_matches'] = 0
    compound_msms_hits['spectrum'] = [np.array([np.array([]), np.array([])])] * len(compound_msms_hits)

    return compound_msms_hits


def score_matrix(cos: Type[CosineHungarian], references: List[SpectrumType], queries: List[SpectrumType]) -> npt.NDArray:
    """
    Calculate matrix of scores and matching ion counts using MatchMS spectrum objects.

    This is a replacement for (and is derived from) the native MatchMS BaseSimilarity.matrix method.
    Source for is here: https://github.com/matchms/matchms/blob/master/matchms/similarity/BaseSimilarity.py

    This additional function is necessary to fix numpy errors that arise when the score matrix is filled with all 0 values.
    """

    score_datatype = [("score", "float"), ("matches", "int")]

    n_rows = len(references)
    n_cols = len(queries)
    idx_row = []
    idx_col = []
    scores = []
    for i_ref, reference in enumerate(references[:n_rows]):
        for i_query, query in enumerate(queries[:n_cols]):
            score = cos.pair(reference, query)

            idx_row.append(i_ref)
            idx_col.append(i_query)
            scores.append(score)

    idx_row = np.array(idx_row)
    idx_col = np.array(idx_col)
    scores_data = np.array(scores, dtype=score_datatype)

    scores_array = np.zeros(shape=(n_rows, n_cols), dtype=score_datatype)
    scores_array[idx_row, idx_col] = scores_data.reshape(-1)

    return scores_array


def get_hits_per_compound(cos: Type[CosineHungarian], compound_idx: int,
                          msms_data: pd.DataFrame, msms_refs: pd.DataFrame) -> pd.DataFrame:
    """
    Get MSMS hits for an individual compound index and InChiKey.

    If no key is given, or it isn't present in the MSMS refs, create nonmatched dummy hits DataFrame.
    """
    unique_inchi_keys = msms_data[msms_data['compound_idx']==compound_idx]['inchi_key'].unique()

    if len(unique_inchi_keys) > 1:
        raise Exception('Only 1 inchi_key is allowed per compound name.')
    elif len(unique_inchi_keys) == 0:
        inchi_key = None
    else:
        inchi_key = unique_inchi_keys[0]

    if inchi_key not in msms_refs['inchi_key'].tolist() or inchi_key is None:
        nonmatched_msms_hits = create_nonmatched_msms_hits(msms_data, inchi_key, compound_idx)
        return nonmatched_msms_hits

    filtered_msms_refs = msms_refs[msms_refs['inchi_key']==inchi_key].reset_index(drop=True).copy()
    filtered_msms_refs = build_msms_refs_spectra(filtered_msms_refs)

    filtered_msms_data = msms_data[msms_data['inchi_key']==inchi_key].reset_index(drop=True).drop(columns=['inchi_key', 'precursor_mz']).copy()

    scores_matches = score_matrix(cos, filtered_msms_data.matchms_spectrum.tolist(), filtered_msms_refs.matchms_spectrum.tolist())

    inchi_msms_hits = pd.merge(filtered_msms_data, filtered_msms_refs.drop(columns=['name', 'adduct']), how='cross')
    inchi_msms_hits['score'] = scores_matches['score'].flatten()
    inchi_msms_hits['num_matches'] = scores_matches['matches'].flatten()

    inchi_msms_hits['precursor_ppm_error'] = (abs(inchi_msms_hits['measured_precursor_mz'] - inchi_msms_hits['precursor_mz']) / inchi_msms_hits['precursor_mz']) * 1000000
    inchi_msms_hits = inchi_msms_hits[inchi_msms_hits['precursor_ppm_error']<=inchi_msms_hits['cid_pmz_tolerance']]

    if inchi_msms_hits.empty:
        nonmatched_msms_hits = create_nonmatched_msms_hits(msms_data, inchi_key, compound_idx)
        return nonmatched_msms_hits

    return inchi_msms_hits

def get_msms_hits(metatlas_dataset: MetatlasDataset, extra_time: bool | float = False, keep_nonmatches: bool = False,
                  pre_query: str= 'database == "metatlas"', query: str | None = None, ref_dtypes: Dict[str, Type[Any]] | None = None,
                  ref_loc: str | None = None, ref_df: pd.DataFrame | None = None, frag_mz_tolerance: float = 0.005,
                  ref_index: List[str] | str | None = None, do_centroid: bool = False, resolve_by: str | None = None) -> pd.DataFrame:
    """
    Get MSMS Hits from metatlas dataset and MSMS refs.

    Input:
        Metatlas Dataset
        MSMS refs path or DataFrame
    Parameters:
        extra_time: If float instead of False, extra time will be added to the rt_min and max set in the atlas
        keep_nonmatches: Whether or not to keep MSMS hits that haven't matched any spectra in the MSMS refs file
        pre_query: DataFrame query applied after reading in the MSMS refs file
        query: Optional additional DataFrame query applied
        ref_dtypes: Dtypes used when reading in MSMS refs file
        ref_loc: path to MSMS refs file
        ref_df: DataFrame of MSMS refs
        frag_mz_tolerance: m/z tolerance used when matching peaks during scoring
        ref_index: index to use in MSMS refs DataFrame
        do_centroid: centroid data before scoring (currently not implemented)
        resolve_by: how to resolve multiple matching peaks within tolerance (currently not implemented)

    Output:
        MSMS hits DataFrame with scores, aligned spectra for plotting, and metadata
    """

    msms_hits_cols = ['database', 'id', 'file_name', 'msms_scan', 'score', 'num_matches',
                     'msv_query_aligned', 'msv_ref_aligned', 'name', 'adduct', 'inchi_key',
                     'precursor_mz', 'measured_precursor_mz',
                     'measured_precursor_intensity']

    if ref_df is None:
        msms_refs = load_msms_refs_file(ref_loc, pre_query, query, ref_dtypes, ref_index)
    else:
        msms_refs = ref_df

    msms_data = ma_data.arrange_ms2_data(metatlas_dataset, do_centroid)

    if msms_data.empty:
        return pd.DataFrame(columns=msms_hits_cols).set_index(['database', 'id', 'file_name', 'msms_scan'])

    if not extra_time:
        cid_rt_min = msms_data['cid_rt_min']
        cid_rt_max = msms_data['cid_rt_max']
    else:
        cid_rt_min = msms_data['cid_rt_min'] - extra_time
        cid_rt_max = msms_data['cid_rt_max'] + extra_time

    msms_data = msms_data[(msms_data['msms_scan']>=cid_rt_min) & (msms_data['msms_scan']<=cid_rt_max)]

    if len(msms_data) == 0:
        return pd.DataFrame(columns=msms_hits_cols).set_index(['database', 'id', 'file_name', 'msms_scan'])

    compound_names = ma_data.get_compound_names(metatlas_dataset)[0]

    cos = CosineHungarian(tolerance=frag_mz_tolerance)

    msms_hits = []
    for compound_idx in tqdm(range(len(compound_names)), unit='compound', disable=in_papermill()):
        compound_msms_hits = get_hits_per_compound(cos, compound_idx, msms_data, msms_refs)
        msms_hits.append(compound_msms_hits)

    msms_hits = pd.concat(msms_hits)

    msv_queries_aligned, msv_refs_aligned = sp.align_all_ms_vectors(zip(msms_hits['query_spectrum'].tolist(), msms_hits['spectrum'].tolist()), frag_tolerance=frag_mz_tolerance)
    msms_hits['msv_query_aligned'] = msv_queries_aligned
    msms_hits['msv_ref_aligned'] = msv_refs_aligned

    if not keep_nonmatches:
        msms_hits.dropna(subset=['database', 'id'], how='all', inplace=True)

    return msms_hits[msms_hits_cols].set_index(['database', 'id', 'file_name', 'msms_scan'])

def make_chromatograms(input_dataset, include_lcmsruns=None, exclude_lcmsruns=None, include_groups=None, exclude_groups=None, group='index', share_y=True, save=True, output_loc=None, short_names_df=None, short_names_header=None, polarity='', overwrite=False, max_cpus=1, suffix='', max_plots_per_page=30):
    bad_parameters = {"group": group != "index", "save": not save,
                      "short_names_df": short_names_df is not None}
    if any(bad_parameters.values()):
        warnings.warn((f"Parameters {', '.join([k for k, v in bad_parameters.items() if v])} of "
                       "make_chromatograms() are no longer utilized and will be removed in an "
                       "up coming release"), FutureWarning, stacklevel=2)
    data = filter_runs(input_dataset, include_lcmsruns, include_groups, exclude_lcmsruns, exclude_groups)
    prefix = f"{polarity}_" if polarity != '' else ''
    out_dir = Path(output_loc) / f"{prefix}compound_EIC_chromatograms{suffix}"
    logger.info('Saving chromatograms to %s.', out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    disable_interactive_plots()
    compound_names = ma_data.get_compound_names(data, use_labels=True)[0]
    args = [(data, i, out_dir / f"{name}.pdf", overwrite, share_y, max_plots_per_page)
            for i, name in enumerate(compound_names)]
    parallel.parallel_process(_save_eic_pdf, args, max_cpus, unit='plot')


def _save_eic_pdf(multi_args):
    save_compound_eic_pdf(*multi_args)


def make_identification_figure_v2(input_fname: Optional[Path] = None, input_dataset=[], include_lcmsruns=[], exclude_lcmsruns=[],
                                  include_groups=[], exclude_groups=[], output_loc: Path = None, msms_hits=None,
                                  use_labels=False, intensity_sorted_matches=False,
                                  short_names_df=pd.DataFrame(), polarity='', overwrite=True):
    assert output_loc
    assert input_fname or input_dataset
    prefix = '' if polarity == '' else f"{polarity}_"
    output_loc = output_loc / f"{prefix}msms_mirror_plots"
    logger.info("Exporting indentification figures to %s", output_loc)
    data = input_dataset if input_dataset else ma_data.get_dill_data(os.path.expandvars(input_fname))
    data = filter_runs(data, include_lcmsruns, include_groups, exclude_lcmsruns, exclude_groups)

    if msms_hits is not None:
        msms_hits_df = msms_hits.reset_index().sort_values('score', ascending=False)
    compound_names = ma_data.get_compound_names(data, use_labels)[0]
    file_names = ma_data.get_file_names(data)

    match_dtypes = {'label': str,
                    'file name': str,
                    'Matching M/Zs above 1E-3*max': str,
                    'All matching M/Zs': str}

    match = pd.DataFrame(columns=match_dtypes).astype(match_dtypes)
    disable_interactive_plots()
    plt.clf()
    for compound_idx, _ in enumerate(compound_names):
        file_idxs, scores, msv_sample_list, msv_ref_list, rt_list = [], [], [], [], []
        inchi_key = extract(data, [0, compound_idx, "identification", "compound", 0, "inchi_key"], "")
        #  Find 5 best file and reference pairs by score
        try:
            rt_ref = data[0][compound_idx]['identification'].rt_references[0]
            mz_ref = data[0][compound_idx]['identification'].mz_references[0]
            comp_msms_hits = msms_hits_df[(msms_hits_df['inchi_key'] == inchi_key)
                                          & (msms_hits_df['msms_scan'] >= rt_ref.rt_min)
                                          & (msms_hits_df['msms_scan'] <= rt_ref.rt_max)
                                          & within_tolerance(
                                                msms_hits_df['precursor_mz'].values.astype(float),
                                                mz_ref.mz,
                                                mz_ref.mz_tolerance*1e-6
                                            )
                                          ].drop_duplicates('file_name').head(5)
            comp_msms_hits = comp_msms_hits[comp_msms_hits['file_name'].isin(file_names)]
            file_idxs = [file_names.index(f) for f in comp_msms_hits['file_name']]
            scores = comp_msms_hits['score'].values.tolist()
            msv_sample_list = comp_msms_hits['msv_query_aligned'].values.tolist()
            msv_ref_list = comp_msms_hits['msv_ref_aligned'].values.tolist()
            rt_list = comp_msms_hits['msms_scan'].values.tolist()
        except (IndexError, TypeError):
            file_idx = None
            max_intensity = 0
            for file_idx, _ in enumerate(data):
                try:
                    temp = max(data[file_idx][compound_idx]['data']['eic']['intensity'])
                    if temp > max_intensity:
                        max_file_idx = file_idx
                        max_intensity = temp
                except (ValueError, TypeError):
                    continue

            file_idxs = [max_file_idx]
            msv_sample_list = [np.array([0, np.nan]).T]
            msv_ref_list = [np.array([0, np.nan]).T]
            scores = [np.nan]

        #  Plot if compound yields any scores
        if file_idxs and file_idxs[0] is not None:
            #  Top 5 MSMS Spectra
            top_5_axis = [plt.subplot2grid((24, 24), (0, 0), rowspan=12, colspan=12)]
            for i in [0, 3, 6, 9]:
                top_5_axis.append(plt.subplot2grid((24, 24), (i, 12), rowspan=3, colspan=3))
                top_5_axis[-1].tick_params(axis='both', length=2)
                top_5_axis[-1].set_xticklabels([])
                top_5_axis[-1].set_yticklabels([])
            for i, (score, axis) in enumerate(zip(scores, top_5_axis)):
                plot_msms_comparison(i, score, axis, msv_sample_list[i], msv_ref_list[i])

            def no_axis_plot(i):
                axis = plt.subplot2grid((24, 24), (i, 15), rowspan=3, colspan=1)
                axis.axis('off')
                return axis

            # Next Best Scores and Filenames
            next_best = [no_axis_plot(i) for i in [0, 3, 6, 9]]

            if short_names_df.empty:
                for i, (score, axis) in enumerate(zip(scores[1:], next_best)):
                    plot_score_and_ref_file(axis, score, rt_list[i+1], os.path.basename(data[file_idxs[i+1]][compound_idx]['lcmsrun'].hdf5_file))
            else:
                for i, (score, ax) in enumerate(zip(scores[1:], next_best)):
                    short_samplename = short_names_df.loc[os.path.basename(data[file_idxs[i+1]][compound_idx]['lcmsrun'].hdf5_file).split('.')[0], 'short_samplename'][0]
                    plot_score_and_ref_file(ax, score, rt_list[i+1], short_samplename)

        # EMA Compound Info
        if file_idxs and file_idxs[0] is not None:
            ax3 = plt.subplot2grid((24, 24), (0, 16), rowspan=6, colspan=8)
            plot_ema_compound_info(ax3, data[file_idxs[0]][compound_idx]['identification'])
        else:
            ax3 = plt.subplot2grid((24, 24), (0, 0), rowspan=6, colspan=8)
            plot_ema_compound_info(ax3, data[0][compound_idx]['identification'])

        # Structure
        if file_idxs and file_idxs[0] is not None:
            ax5 = plt.subplot2grid((24, 24), (13, 0), rowspan=6, colspan=6)
            plot_structure(ax5, data[file_idxs[0]][compound_idx]['identification'].compound, 100)
        else:
            ax5 = plt.subplot2grid((24, 24), (13, 0), rowspan=6, colspan=6)
            plot_structure(ax5, data[0][compound_idx]['identification'].compound, 100)

        # EIC
        if file_idxs and file_idxs[0] is not None:
            ax6 = plt.subplot2grid((21, 21), (6, 15), rowspan=5, colspan=6)
            plot_eic(ax6, data, compound_idx)
        else:
            ax6 = plt.subplot2grid((21, 21), (6, 0), rowspan=5, colspan=6)
            plot_eic(ax6, data, compound_idx)

        # Old code
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
                match.loc[compound_idx, 'label'] = compound_names[compound_idx]
                match.loc[compound_idx, 'file name'] = file_names[file_idxs[0]]
                match.loc[compound_idx, 'RT'] = rt_list[0]
                match.loc[compound_idx, 'score'] = scores[0]
                match.loc[compound_idx, 'Matching M/Zs above 1E-3*max'] = ', '.join(['%5.3f' % m for m in threshold_mz_sample_matches])
                match.loc[compound_idx, 'All matching M/Zs'] = ','.join(['%5.3f' % m for m in mz_sample_matches])

            ax7.set_ylim(.5,1.1)
            ax7.axis('off')

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="tight_layout not applied: number of rows in subplot specifications must be multiples of one another.")
            plt.tight_layout()
        fig_path = output_loc / f"{compound_names[compound_idx]}.pdf"
        write_utils.check_existing_file(fig_path, overwrite)
        plt.savefig(fig_path)
        plt.close()
        logger.debug('Exported identification figures for %s to %s.', compound_names[compound_idx], fig_path)
    match_path = output_loc / "MatchingMZs.tab"
    write_utils.export_dataframe(match, match_path, 'matching MZs', overwrite, sep='\t', float_format="%.12e")


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
    data = filter_runs(data, include_lcmsruns, include_groups, exclude_lcmsruns, exclude_groups)

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

    disable_interactive_plots()
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


def export_atlas_to_spreadsheet(atlas, output_filename=None):
    """
    inputs:
        atlas: metatlas.datastructures.metatlas_objects.Atlas or metatlas_dataset
        output_filename: location to save csv
    output:
        returns a pandas DataFrame containing atlas
        Saves output DataFrame to output_filename in csv format.
    """
    # cols is a list of tuples, with column name as first value, and extract() ids list as second value
    cols = [(c, ['compound', 0, c]) for c in metob.Compound.class_trait_names() if not c.startswith('_')]
    cols = sorted(cols, key=lambda x: x[0])
    cols.extend([('label', ['name']), ('id_notes', ['description'])])
    cols.extend([(c, [c]) for c in ['ms1_notes', 'ms2_notes', 'identification_notes']])
    cols.extend([(c, ['rt_references', 0, c]) for c in ['rt_min', 'rt_max', 'rt_peak']])
    cols.extend([(c, ['mz_references', 0, c]) for c in ['mz', 'mz_tolerance', 'adduct']])
    cols.append(('polarity', ['mz_references', 0, 'detected_polarity']))

    cols_dtypes = {'chebi_id': str,
                   'chebi_url': str,
                   'creation_time': float,
                   'description': str,
                   'formula': str,
                   'head_id': str,
                   'hmdb_id': str,
                   'hmdb_url': str,
                   'img_abc_id': str,
                   'inchi': str,
                   'inchi_key': str,
                   'iupac_name': str,
                   'kegg_id': str,
                   'kegg_url': str,
                   'last_modified': float,
                   'lipidmaps_id': str,
                   'lipidmaps_url': str,
                   'metacyc_id': str,
                   'mono_isotopic_molecular_weight': float,
                   'name': str,
                   'neutralized_2d_inchi': str,
                   'neutralized_2d_inchi_key': str,
                   'neutralized_inchi': str,
                   'neutralized_inchi_key': str,
                   'num_free_radicals': float,
                   'number_components': float,
                   'permanent_charge': float,
                   'prev_uid': str,
                   'pubchem_compound_id': str,
                   'pubchem_url': str,
                   'source': str,
                   'synonyms': str,
                   'unique_id': str,
                   'username': str,
                   'wikipedia_url': str,
                   'label': str,
                   'id_notes': str,
                   'ms1_notes': str,
                   'ms2_notes': str,
                   'identification_notes': str,
                   'rt_min': float,
                   'rt_max': float,
                   'rt_peak': float,
                   'mz': float,
                   'mz_tolerance': float,
                   'adduct': str,
                   'polarity': str}

    out = pd.DataFrame(columns=cols_dtypes).astype(cols_dtypes)
    is_atlas = isinstance(atlas, metob.Atlas)
    compound_ids = atlas.compound_identifications if is_atlas else [i['identification'] for i in atlas[0]]
    for i, my_id in enumerate(compound_ids):
        for column_name, ids in cols:
            out.loc[i, column_name] = extract(my_id, ids)
    if output_filename:
        os.makedirs(os.path.dirname(output_filename), exist_ok=True)
        out.to_csv(output_filename)
    return out

def get_data_for_groups_and_atlas(group,myAtlas,output_filename,use_set1 = False):
    """
    get and pickle everything This is MSMS, raw MS1 datapoints, compound, group info, and file info
    """
    data = []
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
                temp_compound = deepcopy(compound)
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


def compound_indices_marked_remove(data):
    """
    inputs:
        data: metatlas_dataset
    outputs:
        list of compound_idx of the compound identifications with ms1_notes to remove
    """
    return [i for i, j in enumerate(data[0]) if is_remove(extract(j, ['identification', 'ms1_notes']))]


def is_remove(obj):
    """ is obj a string that starts with 'remove' (case insensitive)? """
    return isinstance(obj, str) and obj.lower().startswith('remove')


def first_not_none(obj, default):
    """ returns obj if it is not None, otherwise returns default """
    return default if obj is None else obj


def filter_by_remove(atlas_df, data):
    """
    inputs:
        atlas_df: pandas DataFrame containing an atlas
        data: metatlas_dataset
    outputs:
        a tuple containing 2 pandas DataFrames:
            atlas_df where ms1_notes begins with 'remove'
            atlas_df where ms1_notes does not begin with 'remove'

    ms1_notes comparison with 'remove' is case insensitive.
    """
    rm_idxs = compound_indices_marked_remove(data)
    keep_idxs = atlas_df.index.difference(rm_idxs)
    return(atlas_df.iloc[keep_idxs].copy(), atlas_df.iloc[rm_idxs].copy())


def get_intensity(compound):
    """
    inputs:
        compound: a CompoundIdentification object
    returns a list of intensity values or an empty list if the intensity attribute does not exist
    """
    return ma_data.extract(compound, ['data', 'eic', 'intensity'], [])


def filter_atlas(atlas_df, data, num_data_points_passing=5, peak_height_passing=1e6):
    """
    inputs:
        atlas_df: panda DataFrame containing an atlas
        data: metatlas_dataset
        num_data_points_passing: number of points in EIC that must be exceeded in one or more samples
                                 in order for the compound to remain in the atlas
        peak_height_passing: max intensity in EIC that must be exceeded in one or more samples
                             in order for the compound to remain in the atlas
    returns a pandas DataFrame containing the updated atlas
    """
    keep_idxs = strong_signal_compound_idxs(data, num_data_points_passing, peak_height_passing)
    return atlas_df.iloc[keep_idxs].reset_index(drop=True)


def strong_signal_compound_idxs(
        data,
        num_points_passing: Optional[int] = None,
        peak_height_passing: Optional[float] = None,
        msms_score_passing: Optional[float] = None) -> List[int]:
    """
    inputs:
        data: metatlas_dataset
        num_data_points_passing: number of points in EIC that must be exceeded in one or more samples
                                 in order for the compound to remain in the atlas
        peak_height_passing: max intensity in EIC that must be exceeded in one or more samples
                             in order for the compound to remain in the atlas
        msms_score_passing: max msms spectra similarity score tha tmust be exceeded in one or
                             more samples in order for the compound to remain in the atlas
    returns list of indices that are above the thresholds
    """
    num_points_idxs = num_points_passing_idxs(data, num_points_passing)
    height_idxs = peak_height_passing_idxs(data, peak_height_passing)
    score_idxs = msms_score_passing_idxs(data, msms_score_passing)
    return list(set.intersection(set(num_points_idxs), set(height_idxs), set(score_idxs)))


def all_idxs(in_list: Sized) -> List[int]:
    """Return a list of all the indices of in_list"""
    return list(range(len(in_list)))


def peak_height_passing_idxs(data, height_threshold: Optional[float]) -> List[int]:
    """" Returns a list of compound indices that have intensity exceeding height_threshold"""
    def eval_peak_height(compound):
        return (max(get_intensity(compound)+[0])) > height_threshold
    if len(data) == 0:
        return []
    return all_idxs(data[0]) if height_threshold is None else get_passing_idxs(data, eval_peak_height)


def num_points_passing_idxs(data, min_points: Optional[int]) -> List[int]:
    """"
    Returns a list of compound indices that have more than min_points of intensity data
    for atleast one sample
    """
    def eval_num_points(compound):
        return len(get_intensity(compound)) > min_points
    if len(data) == 0:
        return []
    return all_idxs(data[0]) if min_points is None else get_passing_idxs(data, eval_num_points)


def get_passing_idxs(data, eval_func: Callable) -> List[int]:
    per_sample_compound = np.array([[eval_func(compound) for compound in sample] for sample in data]
                                   ).any(axis=0)
    return np.flatnonzero(per_sample_compound).tolist()


def msms_score_passing_idxs(data, msms_score: Optional[float]) -> List[int]:
    """
    inputs:
        msms_score: spectra similarity score must be exceeded to pass
    """
    if len(data) == 0:
        return []
    if msms_score is None:
        return all_idxs(data[0])
    scores_df = fastanalysis.make_scores_df(data.data, data.hits)
    return [i for i, row in scores_df.iterrows() if row["max_msms_score"] > msms_score]


def filter_metatlas_objects_to_most_recent(object_list, field):
    # remove from list if another copy exists that is newer
    unique_values = []
    for a in object_list:
        unique_values.append(getattr(a, field))
    unique_values = list(set(unique_values))
    keep_object_list = []
    for u in unique_values:
        old_last_modified = 0
        for a in object_list:
            if getattr(a, field) == u:
                last_modified = getattr(a, 'last_modified')
                if last_modified > old_last_modified:
                    keep_object = a
                    old_last_modified = last_modified
        keep_object_list.append(keep_object)
    return keep_object_list


def get_metatlas_atlas(name = '%%',username = '*', most_recent = True,do_print = True):
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


def get_metatlas_files(experiment: Union[str, Sequence[str]] = '%', name: str = '%', most_recent: bool = True):
    """
    experiment is the folder name or a list of folder names
    name is the filename
    """
    batches = [experiment] if isinstance(experiment, str) else experiment
    files = list(itertools.chain.from_iterable(
        [metob.retrieve('LcmsRun', experiment=batch.rstrip('%'), name=name, username='*') for batch in batches]
    ))
    if most_recent:
        files = filter_metatlas_objects_to_most_recent(files, 'mzml_file')
    return sorted(files, key=lambda x: x.name)


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


def check_compound_names(atlas_df):
    """
    inputs:
        atlas_df: pandas dataframe representation of an atlas
    throws ValueError if some compounds are not found in the database
    """
    inchi_keys = {
        row.inchi_key
        for _, row in atlas_df.iterrows()
        if (not pd.isnull(row.inchi_key)) and (len(row.inchi_key) > 0) and row.inchi_key != 'None'
    }
    if len(inchi_keys) == 0:
        return
    found_keys = {result.inchi_key for result in
                  metob.retrieve('Compounds', inchi_key=list(inchi_keys), username='*')}
    bad_names = inchi_keys.difference(found_keys)
    if bad_names:
        raise ValueError(f"Compound not found in database: {', '.join(bad_names)}.")


def check_filenames(atlas_df, field):
    """
    inputs:
        atlas_df: pandas dataframe representation of an atlas
        field: column name in atlas_df to test for valid lcmsruns
    throws ValueError if values in atlas_df[field] are not in database as lcmsruns
    """
    try:
        to_test = set(atlas_df[field].to_list())
    except KeyError:
        return  # nothing to check against
    found = {row['name'] for row in metob.retrieve('Lcmsruns', name=list(to_test), username='*')
             if 'name' in row}
    bad_files = set(to_test).difference(found)
    if bad_files:
        raise ValueError(f"LCMS runs not found in database: {', '.join(bad_files)}.")


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


def _clean_dataframe(dataframe, required_columns=None, lower_case_col_names=True):
    """
    inputs:
        dataframe: pandas dataframe
        required_columns: list of column names that must have a non-NA values
        lower_case_col_names: should column names be modified to lower case
    Modifies dataframe in place. The following rows removed:
        fully empty (all fields have NA values)
        containing required_columns with 1 or more NA values
    """
    dataframe.dropna(how="all", inplace=True)
    if required_columns is not None and len(required_columns) > 0:
        dataframe.dropna(how="any", subset=required_columns, inplace=True)
    if lower_case_col_names:
        dataframe.columns = [x.lower() for x in dataframe.columns]


def _add_columns(dataframe, column_names, default_values=None):
    """
    inputs:
        dataframe: pandas dataframe
        column_names: list of column names to add to dataframe if they do not already exist
        default_values: a single default value for all columns or a list of default values
                        the same length as column_names
    Modifies the dataframe in place
    """
    assert isinstance(column_names, list)
    num_col = len(column_names)
    if isinstance(default_values, str):
        default_values = [default_values]
    num_default = 1 if default_values is None else len(default_values)
    assert num_default in [1, num_col]
    default_values = [default_values]*num_col if num_default == 1 else default_values
    for name, default in zip(column_names, default_values):
        if name not in dataframe.columns:
            dataframe[name] = default


def _get_dataframe(filename_or_df=None, filetype=None, sheetname=None):
    """
    inputs:
        filename_or_df: a filename to an excel, tsv or csv file, or a pandas DataFrame
        filetype: a string in dataframe, excel, tab, csv
        sheetname: name of a sheet in an excel file, or get first sheet if None
    returns a pandas Dataframe
    """
    assert filetype in ['dataframe', 'excel', 'tab', 'csv']
    if filetype == 'dataframe':
        return filename_or_df.copy()
    if filetype == 'excel':
        return pd.read_excel(filename_or_df, sheetname=0 if sheetname is None else sheetname)
    return pd.read_csv(filename_or_df, sep='\t' if filetype == 'tab' else ',')


def get_compounds(row):
    # currently, all copies of the molecule are returned.  The 0 is the most recent one.
    results = (None if pd.isnull(row.inchi_key) else
               metob.retrieve('Compounds', inchi_key=row.inchi_key, username='*'))
    compound = results[0] if results else metob.Compound()
    compound.name = row.label if isinstance(row.label, str) and len(row.label) > 0 else 'no label'
    return [compound]


def get_compound_identification(row, polarity, mz_tolerance):
    my_id = metob.CompoundIdentification()
    my_id.compound = get_compounds(row)
    my_id.name = my_id.compound[0].name
    _copy_attributes(row, my_id, ['do_normalization', 'internal_standard_id', 'internal_standard_to_use',
                                  'identification_notes', 'ms1_notes', 'ms2_notes'])
    my_id.mz_references = get_mz_references(row, polarity, mz_tolerance)
    my_id.rt_references = get_rt_references(row)
    my_id.frag_references = get_frag_references(row, my_id.name, polarity,
                                                my_id.mz_references[0], my_id.rt_references[0])
    my_id.intensity_references = []
    return my_id


def get_mz_references(row, polarity, mz_tolerance=None):
    assert polarity in ['positive', 'negative']
    mzr = metob.MzReference()
    mzr.mz = row.mz
    # TODO: calculate the mz from theoretical adduct and modification if provided.
    #     my_id.mz_references[0].mz = c.MonoIso topic_molecular_weight + 1.007276
    if mz_tolerance is not None:
        mzr.mz_tolerance = mz_tolerance
    else:
        try:
            mzr.mz_tolerance = row.mz_tolerance
        except AttributeError:
            mzr.mz_tolerance = row.mz_threshold
    mzr.mz_tolerance_units = 'ppm'
    mzr.detected_polarity = polarity
    # if 'file_mz' in atlas_df.keys():
    #     f = metob.retrieve('Lcmsruns',name = '%%%s%%'%atlas_df.file_mz[x],username = '*')[0]
    #     mzRef.lcms_run = f
    if pd.notna(row.adduct):
        mzr.adduct = row.adduct
    return [mzr]


def get_rt_references(row):
    rtr = metob.RtReference()
    rtr.rt_units = 'min'
    _copy_attributes(row, rtr, ['rt_min', 'rt_max', 'rt_peak'], error_on_missing=True)
    # if 'file_rt' in atlas_df.keys():
    #     f = metob.retrieve('Lcmsruns',name = '%%%s%%'%atlas_df.file_rt[x],username = '*')[0]
    #     rtr.lcms_run = f
    return [rtr]


def get_frag_references(row, name, polarity, mz_ref, rt_ref):
    """
    inputs:
        row: atlas_df row for the compound identification of interest
        name: compound name
        polarity: positive or negative
        mz_ref: MzReference object
        rt_ref: RtReference object
    returns an array of FragmentationReferences or empty array if no msms data is found
    """
    assert polarity in ['positive', 'negative']
    try:
        run_name = row.file_msms.replace('.mzmL', '')
        run = metob.retrieve('Lcmsruns', name=f"%{run_name}%", username='*')[0]
    except (AttributeError, IndexError):
        return []
    data = ma_data.get_data_for_a_compound(mz_ref, rt_ref, ['msms'], run.hdf5_file, extra_time=0.3)
    if not isinstance(data['msms']['data'], np.ndarray):
        return []
    frag_ref = metob.FragmentationReference()
    frag_ref.lcms_run = run
    frag_ref.polarity = polarity
    frag_ref.precursor_mz = row.mz
    precursor_intensity = data['msms']['data']['precursor_intensity']
    idx_max = np.argwhere(precursor_intensity == np.max(precursor_intensity)).flatten()
    mz_list = data['msms']['data']['mz'][idx_max]
    intensity_list = data['msms']['data']['i'][idx_max]
    frag_ref.mz_intensities = get_spectrum(mz_list, intensity_list)
    logger.info('Found reference msms spectrum for %s in file %s.', name, row.file_msms)
    return [frag_ref]


def get_spectrum(mz_list, intensity_list):
    """
    inputs:
        mz_list: list of mz values
        intensity_list: list of intensities values
    returns a list of MzIntensityPairs()
    """
    assert len(mz_list) == len(intensity_list)
    spectrum = []
    for msms_mz, intensity in zip(mz_list, intensity_list):
        spectrum.append(metob.MzIntensityPair())
        spectrum[-1].mz = msms_mz
        spectrum[-1].intensity = intensity
    return spectrum


def get_atlas(name, atlas_df, polarity, mz_tolerance):
    """
    inputs:
        name: string with name of atlas
        atlas_df: pandas DataFrame with atlas definition
        polarity: positive or negative
        mz_tolerance: float to set for all mz_tolerance values
    returns an Atlas object

    atlas_df should not contain empty strings, use np.NaN instead
    """
    atlas = metob.Atlas()
    atlas.name = name
    atlas.compound_identifications = []
    for _, row in atlas_df.iterrows():
        my_id = get_compound_identification(row, polarity, mz_tolerance)
        if my_id is None:
            logger.warning(('get_atlas() dropping compound %s '
                            '(inchi_key %s) because it is not in the database.'), row.label, row.inchi_key)
        else:
            atlas.compound_identifications.append(my_id)
    return atlas


def make_atlas_from_spreadsheet(filename, atlas_name, filetype, sheetname=None,
                                polarity=None, store=False, mz_tolerance=None):
    '''
    specify polarity as 'positive' or 'negative'

    '''
    required_columns = ['rt_min', 'rt_peak', 'rt_max', 'mz']
    logger.debug('Generating atlas named %s from %s source.', atlas_name, filetype)
    atlas_df = _get_dataframe(filename, filetype, sheetname)
    _clean_dataframe(atlas_df, required_columns=[])  # only remove fully empty
    initial_row_count = len(atlas_df)
    _clean_dataframe(atlas_df, required_columns)
    _add_columns(
        atlas_df,
        column_names=['label', 'inchi_key', 'adduct', 'mz_tolerance', 'polarity'],
        default_values=[np.NaN, np.NaN, np.NaN, mz_tolerance, polarity]
    )
    check_compound_names(atlas_df)
    check_filenames(atlas_df, 'file_msms')
    atlas = get_atlas(atlas_name, atlas_df, polarity, mz_tolerance)
    rows_removed = initial_row_count - len(atlas_df)
    if rows_removed > 0:
        raise ValueError(
            f"Required columns ({', '.join(required_columns)}) missing in {rows_removed} rows."
        )
    if store:
        logger.info('Saving atlas named %s to DB.', atlas_name)
        metob.store(atlas)
    return atlas


def _copy_attributes(source, dest, attribute_list, default_list=None, error_on_missing=False):
    """
    inputs:
        source: object to copy attributes from
        dest: object to copy attributes to
        attribute_list: list of string containing attribute names
        default_list: list of default values corresponding to same positions in attribute_list
    Modifies dest in place to have all attributes from attribute_list with values coming from
    source or default_list. If source does not contain the attribute and default_list is None,
    then do not add the attribute to dest if it does not already exist.
    """
    if default_list is None:
        for attribute in attribute_list:
            if error_on_missing:
                setattr(dest, attribute, getattr(source, attribute))
            else:
                try:
                    setattr(dest, attribute, getattr(source, attribute))
                except AttributeError:
                    pass
    else:
        for attribute, default in zip(attribute_list, default_list):
            setattr(dest, attribute, getattr(source, attribute, default))


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


def filter_metatlas_objects_by_list(object_list, field, filter_list):
    """
    inputs:
        object_list: list to be filtered by its attribute values
        field: name of attribute to filter on
        filter_list: strings that are tested to see if they are substrings of the attribute value
    returns filtered list of objects that have a match in filter_list
    if filter_list is empty, then return object_list
    """
    if filter_list:
        return filter_by_list(object_list, lambda x: getattr(x, field), filter_list)
    return object_list


def remove_metatlas_objects_by_list(object_list, field, filter_list):
    """
    inputs:
        object_list: iterable to be filtered by its attribute values
        field: name of attribute to filter on
        filter_list: strings that are tested to see if they are substrings of the attribute value
    returns filtered list of objects that do not have matches to filter_list
    """
    return filter_by_list(object_list, lambda x: getattr(x, field), filter_list, include=False)


def filter_by_list(data, key_func, term_list, include=True):
    """
    inputs:
        data: iterable to be filtered
        key_func: function that takes a member of d and returns string to compare with term_list
        term_list: strings that are tested to see if they are substrings of key_func return value
        include: if True, then matches are included in output, else matches are excluded
    """
    allow = any if include else lambda x: not any(x)
    return [d for d in data if allow(ext in key_func(d) for ext in term_list)]


def filter_lcmsruns_in_dataset_by_include_list(metatlas_dataset, selector, include_list):
    """
    Returns a metatlas dataset containing LCMS runs or groups (denoted by selector) that have substrings
    listed in the include list.
    selector can be 'lcmsrun' or 'group'
    include_list will look something like this: ['QC','Blank']
    """
    return filter_by_list(metatlas_dataset, lambda x: x[0][selector].name, include_list)


def filter_lcmsruns_in_dataset_by_exclude_list(metatlas_dataset, selector, exclude_list):
    """
    Returns a metatlas dataset containing LCMS runs or groups (denoted by selector) that have substrings
    not listed in the include list.
    selector can be 'lcmsrun' or 'group'
    exclude_list will look something like this: ['QC','Blank']
    """
    return filter_by_list(metatlas_dataset, lambda x: x[0][selector].name, exclude_list, include=False)


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


def select_groups_for_analysis(name='%', description=[], username='*', do_print=True, most_recent=True,
                               remove_empty=True, include_list=[], exclude_list=[]):
    if description:
        groups = metob.retrieve('Groups', name=name, description=description, username=username)
    else:
        groups = metob.retrieve('Groups', name=name, username=username)
    if most_recent:
        groups = filter_metatlas_objects_to_most_recent(groups, 'name')
    if include_list:
        groups = filter_metatlas_objects_by_list(groups, 'name', include_list)
    if exclude_list:
        groups = remove_metatlas_objects_by_list(groups, 'name', exclude_list)
    if remove_empty:
        groups = filter_empty_metatlas_objects(groups, 'items')
    if do_print:
        for i, group in enumerate(groups):
            print((i, group.name,  datetime.utcfromtimestamp(group.last_modified)))
    return groups


def disable_keyboard_shortcuts(mapping):
    """
    Takes a dictionary with a subset of keys from plot.rcParams and values
    are arrays of strings, which are keyboard short cuts to be removed
    """
    for action, remove_keys_list in mapping.items():
        for key_combo in remove_keys_list:
            if action in plt.rcParams and key_combo in plt.rcParams[action]:
                plt.rcParams[action].remove(key_combo)


def get_ms1_df(sample, limit_to_rt_range=True):
    ms1_df = pd.DataFrame(data=sample['data']['eic'])
    if not limit_to_rt_range:
        return ms1_df
    rt_min = sample['identification'].rt_references[0].rt_min
    rt_max = sample['identification'].rt_references[0].rt_max
    keep = np.logical_and(rt_min < ms1_df['rt'], ms1_df['rt'] < rt_max)
    return ms1_df[keep]


def get_mz_centroid(ms1_df):
    if ms1_df.empty:
        return np.nan
    return sum(ms1_df.mz * ms1_df.intensity) / sum(ms1_df.intensity)


def get_rt_peak(ms1_df):
    if ms1_df.empty:
        return np.nan
    ms1_peak_df = ms1_df.loc[ms1_df['intensity'].idxmax()]
    return ms1_peak_df.rt


def get_msms_plot_headers(data, hits, hit_ctr, compound_idx, similar_compounds):
    """
    inputs:
        data: metatlas_dataset-like object
        hits: dataframe
        hit_ctr: the index in hits of the current hit
        compound_idx: index of current compound in 2nd dim of data
        compound: object for current compound
    returns:
        tuple of strings
            (mz_header, rt_header, cpd_header)
    """
    if not hits.empty:
        rt_ms2 = hits.index.get_level_values('msms_scan')[hit_ctr]
        mz_precursor = hits['measured_precursor_mz'].iloc[hit_ctr]

    file_idx = file_with_max_ms1_intensity(data, compound_idx, limit_to_rt_range=True)[0]
    if file_idx is None:
        return ('', '', '')
    rt_theoretical = data[file_idx][compound_idx]['identification'].rt_references[0].rt_peak
    mz_theoretical = data[file_idx][compound_idx]['identification'].mz_references[0].mz
    ms1_df = get_ms1_df(data[file_idx][compound_idx])
    mz_measured = get_mz_centroid(ms1_df)
    rt_ms1 = get_rt_peak(ms1_df)

    delta_mz = abs(mz_theoretical - mz_measured)
    delta_ppm = delta_mz / mz_theoretical * 1e6
    mz_header = ["m/z theoretical = %5.4f" % mz_theoretical,
                 "m/z measured = %5.4f" % mz_measured,
                 "ppm diff = %3.2f" % delta_ppm]
    rt_header = ["RT theoretical = %3.2f" % rt_theoretical,
                 "RT MS1 measured = %3.2f" % rt_ms1]
    if not hits.empty:
        mz_header.insert(0, "precursor m/z = %5.4f" % mz_precursor)
        rt_header.append("RT MS2 measured = %3.2f" % rt_ms2)
    return (', '.join(mz_header), ', '.join(rt_header),
            get_similar_compounds_header(similar_compounds, compound_idx))


def get_similar_compounds_header(similar_compounds, compound_idx):
    """
    inputs:
        similar_compounds: the output from get_similar_compounds()
        compound_index: index of current compound being considered
    returns:
    """
    if len(similar_compounds) < 2:
        return ''
    joined = '; '.join([_similar_compound_to_str(compound, compound_idx) for compound in similar_compounds])
    return f"Similar Compounds = {joined}"


def _similar_compound_to_str(sdict, compound_idx):
    """
    inputs:
        sdict: a dict returned from get_similar_compounds()
        compound_index: index of current compound being considered
    returns:
        string with only non-breaking spaces or '*' if dict represents current compound
    """
    if sdict['index'] == compound_idx:
        return '*'
    return f"{sdict['index']}, {sdict['label']} {{RT-{sdict['rt'].rt_peak:.2f}}}".replace(' ', '\xa0')


def get_msms_plot_data(hits, hit_ctr):
    """
    inputs:
        hits: dataframe of msms hits
        hit_ctr: index of current hit in hits
    returns:
        tuple: (hit_ref_id, hit_score, hit_query, hit_ref)
    """
    if hits.empty:
        hit_ref_id = "N/A"
        hit_score = np.nan
        hit_query = np.full((2, 2, ), np.nan)
        hit_ref = hit_query
    else:
        hit_ref_id = hits.index.get_level_values('id')[hit_ctr]
        hit_score = hits['score'].iloc[hit_ctr]
        hit_query = hits['msv_query_aligned'][hit_ctr:hit_ctr+1].iloc[0]
        hit_ref = hits['msv_ref_aligned'][hit_ctr:hit_ctr+1].iloc[0]
    return (hit_ref_id, hit_score, hit_query, hit_ref)


def get_hit_metadata(data, hits, file_names, hit_ctr, compound_idx):
    """
    returns a tuple containing:
        file name (without path) of the hit or the string 'None'
        compound object for the hit or None
    """
    if not hits.empty:
        hit_file_name = hits.index.get_level_values('file_name')[hit_ctr]
        return (hit_file_name, data[int(file_names.index(hit_file_name))][compound_idx])
    file_idx = file_with_max_ms1_intensity(data, compound_idx, limit_to_rt_range=True)[0]
    if file_idx:
        return (os.path.basename(data[file_idx][compound_idx]['lcmsrun'].hdf5_file),
                data[file_idx][compound_idx])
    return ('None', None)


def within_tolerance(measured, theoretical, tolerance):
    """ Returns True if normalized, absolute difference is with tolerance """
    return abs(measured - theoretical)/theoretical <= tolerance


def layout_radio_button_set(area, anchor='SW'):
    """
    inputs:
        area: [left, bottom, width, height]
        anchor: string for anchor direction
    returns:
        an axes for radiobuttons at area with axis off and equal aspect ratio
    """
    axes = plt.axes(area, anchor=anchor, aspect='equal')
    axes.axis('off')
    return axes


def rt_range_overlaps(rt1, rt2):
    """
    inputs:
        rt1: metatlas.datastructures.metatlas_objects.RtReference
        rt2: metatlas.datastructures.metatlas_objects.RtReference
    returns:
        True if there is overlap in the RT min-max regions of rt1 and and rt2
    """
    return ((rt2.rt_min <= rt1.rt_min <= rt2.rt_max) or (rt2.rt_min <= rt1.rt_max <= rt2.rt_max) or
            (rt1.rt_min <= rt2.rt_min <= rt1.rt_max) or (rt1.rt_min <= rt2.rt_max <= rt1.rt_max))


def disable_interactive_plots():
    """Close interactive figures and turn off interactive plotting"""
    adjust_rt_for_selected_compound.disable()
    plt.ioff()


def tic_pdf(data, polarity, file_name, overwrite=False, sharey=True,
            x_min=1.5, x_max=None, y_min=0, y_max=None, max_plots_per_page=30):
    save_sample_tic_pdf(
        data, polarity, file_name, overwrite, sharey, x_min, x_max, y_min, y_max, max_plots_per_page
    )

def make_copy_to_clipboard_button(text: str, button_text: str) -> None:
    display(HTML(f"""
        <button type="button" onclick="copy_to_clipboard()">{button_text}</button>
        <script>
            function copy_to_clipboard() {{
                navigator.clipboard.writeText("{text}").then(function() {{
                    console.log('Async: Copying to clipboard was successful!');
                }}, function(err) {{
                    console.error('Async: Could not copy text: ', err);
                }});
            }}
        </script>
    """))
