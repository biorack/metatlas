from datetime import datetime
from metatlas import metatlas_objects as metob
from metatlas.helpers import metatlas_get_data_helper_fun as ma_data
from metatlas import gui


import pandas as pd
import os
import os.path

from IPython.display import display
import matplotlib.pyplot as plt
try:
    import ipywidgets as widgets
except ImportError:
    from IPython.html import widgets
try:
    import traitlets
except ImportError:
    from IPython.utils import traitlets

from ipywidgets import interact, fixed, FloatSlider


try:
    foo = widgets.Select()
except Exception as e:
    print(e)
    sys.exit(0)



data = []
groups = []
file_names = []
compound_names = []
compound_objects = []
files_idx = dict()
compound_idx = dict()
groups_idx = dict()
grid = gui.create_qgrid([])
grid2 = gui.create_qgrid([])

compounds_list_dict = dict()


atlas_header = ['Atlas Name', 'No. Compounds', 'Last Modified']
compound_header = ['Compound', 'rt_max', 'rt_min', 'rt_peak', 'rt_units', 'mz', 'mz_tolerance', 'mz_tolerance_units',
                   'lcms_run']

# --------------------------------------------------------------------------------------------------------------------
# --------------------- W I D G E T S ---------------- W I D G E T S -------------------- W I D G E T S --------------
# --------------------------------------------------------------------------------------------------------------------

# single select widget for atlas
atlases_select = widgets.Select(description='Atlases:', options=[])

# text widget to select for atlases in a database. accepts wild cards
search_string = widgets.Text(description="")
search_button = widgets.Button(description="Search for Atlas")
display_compounds_and_files = widgets.Button(description="Display Compounds & Files")


# single select widget for the compound
wcompounds = widgets.Select(
    description="compounds",
    options=[]
)

# multiple select widget for the files
wfiles = widgets.Select(
    description="files",
    options=[]
)

# dill file name
wfname = widgets.Text(
    description='Atlas Name',
    value='myAtlas',
)

# button that plots the intensity
plot_button = widgets.Button(description='Plot me')
save_atlas_button = widgets.Button(description="Save Atlas")
save_as_button = widgets.Button(description="Save Atlas As")
save_atlas_as_txt = widgets.Text()

# radio buttons to control how atlas is to be written
choose_output = widgets.RadioButtons(description='choose output:',
                                     options=['new atlas name', 'old atlas name'],
                                     )

# button that creates the atlas based on the selection from the radio buttons
create_atlas_btn = widgets.Button(description="Create Atlas from Dill")
modify_atlas_btn = widgets.Button(description="Modify Selected Atlas")

# text box that holds the selected atlas' RT values for the selected compound. displays
atlas_ref_vals = widgets.Text(description="RT values for compound in Atlas")
atlas_ref_vals.value = "RT values for compound in Atlas go here"
atlas_ref_vals.color = 'red'

# sliders for the user to change the RT values
rtmin_widget = FloatSlider()
rtpeak_widget = FloatSlider()
rtmax_widget = FloatSlider()


# --------------------------------------------------------------------------------------------------------------------
# --------------- F U N C T I O N S ------------- F U N C T I O N S ----------- F U N C T I O N S --------------------
# --------------------------------------------------------------------------------------------------------------------

###########################################################################
###
def mod_atlas_compound_RT_values(**kwargs):
    """
    Parameters
    ----------
    kwargs: dictionary that holds atlas (object or name), compound, rt_min, rt_max, and rt_peak

    Returns a modified atlas object
    -------
    """
    atlas = kwargs['atlas']
    compound = kwargs['compound']
    rt_min = kwargs['rt_min']
    rt_max = kwargs['rt_max']
    rt_peak = kwargs['rt_peak']

    if isinstance(atlas, str):
        atlas = metob.retrieve('Atlas', name=atlas, username='*')

    num_compounds = len(atlas[0].compound_identifications)

    for x in range(num_compounds):
        cpd_name = atlas[0].compound_identifications[x].compound[0].name
        if compound == cpd_name:  # adjust the rt_values
            atlas[0].compound_identifications[x].rt_references[0].rt_min = rt_min
            atlas[0].compound_identifications[x].rt_references[0].rt_max = rt_max
            atlas[0].compound_identifications[x].rt_references[0].rt_peak = rt_peak

            break

    return atlas



###########################################################################
###
def atlas_grid(sender):
    atlas_dict = dict()
    for i in atlas_header:
        atlas_dict[i] = list()

    wild_card = search_string.value
    atlas = metob.retrieve('Atlas', name=wild_card, username='*')

    for i, a in enumerate(atlas):
        atlas_dict['Atlas Name'].append(a.name)
        atlas_dict['No. Compounds'].append(str(len(a.compound_identifications)))
        atlas_dict['Last Modified'].append(str(datetime.utcfromtimestamp(a.last_modified)))

    grid.df = pd.DataFrame.from_dict(atlas_dict)
    grid.width = "100%"



###########################################################################
###
def plot_intensity(cval, fvals, rt_min, rt_max, rt_peak):
    global data

    for idx, _fs in enumerate(fvals):
        # d = data[idx][compound_idx]
        d = data[idx][cval]

        if len(d['data']['eic']['rt']) > 0:
            x = d['data']['eic']['rt']
            y = d['data']['eic']['intensity']
            plt.plot(x, y, 'k-', ms=1, mew=0, mfc='b', alpha=1.0)

    plt.axvline(rt_min, color='b', linewidth=2.0)
    plt.axvline(rt_max, color='g', linewidth=2.0)
    plt.axvline(rt_peak, color='r', linewidth=2.0)


###########################################################################
##
def plot_button_clicked(sender):
    global data, rtmin_widget, rtmax_widget, rtpeak_widget

    plt.cla()
    plt.clf()
    plt.close()

    # get the pkl file name from the selection box
    pkl_fname = wfiles.value

    # get data and compound names from pickled file
    data = ma_data.get_dill_data(pkl_fname)
    file_names = ma_data.get_file_names(data)
    (compound_names, compound_objects) = ma_data.get_compound_names(data)
    
    # get the name of the compound as selected from the grid
    print(grid2.get_selected_rows())
    n = grid2.get_selected_rows()[0]
    atlas_compound = grid2.df.loc[n]['Compound']

    min_x = list()
    max_x = list()

    # see if selected atlas compound is in the pickle file
    if atlas_compound not in compound_names:
        print("Compound not found")
        return

    compound_idx = compound_names.index(atlas_compound)

    for idx, _fs in enumerate(file_names):
        # d = data[idx][compound_idx]
        d = data[idx][0]
        rt_min = d['identification'].rt_references[0].rt_min
        rt_max = d['identification'].rt_references[0].rt_max
        rt_peak = d['identification'].rt_references[0].rt_peak

        if len(d['data']['eic']['rt']) > 0:
            x = d['data']['eic']['rt']
            y = d['data']['eic']['intensity']
            min_x.append(min(x))
            max_x.append(max(x))
            plt.plot(x, y, 'k-', ms=1, mew=0, mfc='b', alpha=1.0)

    plt.axvline(rt_min, color='b', linewidth=2.0)
    plt.axvline(rt_max, color='g', linewidth=2.0)
    plt.axvline(rt_peak, color='r', linewidth=2.0)

    rtmin_widget.close()
    rtpeak_widget.close()
    rtmax_widget.close()
    rtmin_widget = FloatSlider(min=min(min_x), max=max(max_x), step=0.01, value=rt_min, color='blue')
    rtpeak_widget = FloatSlider(min=min(min_x), max=max(max_x), step=0.01, value=rt_peak, color='red')
    rtmax_widget = FloatSlider(min=min(min_x), max=max(max_x), step=0.01, value=rt_max, color='green')

    interact(plot_intensity,
             cval=fixed(compound_idx),
             fvals=fixed(file_names),
             rt_min=rtmin_widget,
             rt_peak=rtpeak_widget,
             rt_max=rtmax_widget)


###########################################################################
##
def display_atlases():
    search_string.width = '100%'
    atlases_select.width = '100%'
    hbox = widgets.HBox((search_string, search_button))
    display(hbox)
    search_button.on_click(atlas_grid)
    grid.on_msg(my_handler)
    display(grid)
    display(grid2)


###########################################################################
###
def display_pkl_files_and_plot_data(pkl_path='$HOME',filter_str = '', extension = '*.pkl'):
    import subprocess
    import fnmatch
    # if location == '$HOME':
    #     pkl_path = os.path.expandvars(location)
    # else:
    #     pkl_path = os.path.expandvars(os.path.join('$HOME', location))
    # print pkl_path
    # pkl_file_list = subprocess.Popen(["find", pkl_path,"-iname",filter_str, "-iname", extension],
                                     # stdout=subprocess.PIPE).communicate()[0]
# 
    # pkl_file_list = pkl_file_list.split('\n')
    pkl_file_list = []
    for root, dirnames, filenames in os.walk(pkl_path):
        for filename in fnmatch.filter(filenames, extension):
            fname = os.path.join(root, filename)
            if filter_str:
                if filter_str.lower() in os.path.basename(fname).lower():
                    pkl_file_list.append(fname)
            else:
                pkl_file_list.append(fname)
    wfiles.options = pkl_file_list

    display(widgets.HBox((plot_button, save_atlas_button, save_as_button, save_atlas_as_txt)))
    display(wfiles)
    wfiles.width = "100%"
    plot_button.on_click(plot_button_clicked)
    save_atlas_button.on_click(save_atlas)
    save_as_button.on_click(save_atlas_as)


###########################################################################
###
def save_atlas(sender):
    n = grid.get_selected_rows()
    m = grid2.get_selected_rows()
    if len(m) == 0 or len(n) == 0:
        ok_to_save_atlas = False
    else:
        ok_to_save_atlas = True
        m = m[0]
        n = n[0]

    if ok_to_save_atlas:
        kwargs = dict()
        kwargs['atlas'] = grid.df.loc[n]['Atlas']
        kwargs['compound'] = grid2.df.loc[m]['Compound']
        kwargs['rt_min'] = rtmin_widget.value
        kwargs['rt_max'] = rtmax_widget.value
        kwargs['rt_peak'] = rtpeak_widget.value

        atlas = mod_atlas_compound_RT_values(**kwargs)

        metob.store(atlas)


    else:
        print("cannot save atlas")


###########################################################################
###
def save_atlas_as(sender):
    n = grid.get_selected_rows()
    m = grid2.get_selected_rows()
    if len(m) == 0 or len(n) == 0 or len(save_atlas_as_txt.value) == 0:
        ok_to_save_atlas = False
    else:
        ok_to_save_atlas = True
        m = m[0]
        n = n[0]

    if ok_to_save_atlas:
        kwargs = dict()
        kwargs['atlas'] = grid.df.loc[n]['Atlas']
        kwargs['compound'] = grid2.df.loc[m]['Compound']
        kwargs['rt_min'] = rtmin_widget.value
        kwargs['rt_max'] = rtmax_widget.value
        kwargs['rt_peak'] = rtpeak_widget.value

        atlas = mod_atlas_compound_RT_values(**kwargs)

        if len(save_atlas_as_txt.value) > 1:
            atlas.name = save_atlas_as_txt.value

        metob.store(atlas)


    else:
        print("cannot save atlas")


###########################################################################
###
def my_handler(widget, content, buffers=None):
    if content['type'] == 'selection_change':
        row = content['rows'][0]
        # get the compounds in that atlas and display their content
        atlas_name = grid.df['Atlas Name'][row]

        atlas = metob.retrieve('Atlas', name=atlas_name, username="*")
        compound_vals_dict = dict()
        for i in compound_header:
            compound_vals_dict[i] = list()

        for x in atlas[0].compound_identifications:
            
            if x.compound:
                compound_vals_dict['Compound'].append(str(x.compound[0].name))
            else:
                compound_vals_dict['Compound'].append(str(x.name))

            compound_vals_dict['rt_max'].append(str(x.rt_references[0].rt_max))
            compound_vals_dict['rt_min'].append(str(x.rt_references[0].rt_min))
            compound_vals_dict['rt_peak'].append(str(x.rt_references[0].rt_peak))
            compound_vals_dict['rt_units'].append(str(x.rt_references[0].rt_units))
            compound_vals_dict['mz'].append(str(x.mz_references[0].mz))
            compound_vals_dict['mz_tolerance'].append(str(x.mz_references[0].mz_tolerance))
            compound_vals_dict['mz_tolerance_units'].append(str(x.mz_references[0].mz_tolerance_units))
            compound_vals_dict['lcms_run'].append(str(x.rt_references[0].lcms_run))

        grid2.df = pd.DataFrame.from_dict(compound_vals_dict)
        grid2.width = '100%'







