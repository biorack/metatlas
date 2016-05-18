import sys
sys.path.insert(0,'/home/jimmy/dev/metatlas' )


from metatlas import metatlas_objects as metob
import metatlas.helpers.metatlas_get_data_helper_fun as ma_data
import os


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

from ipywidgets import interact, interactive, fixed, FloatSlider

# one select for the compound
wcompounds = widgets.Select(
    description="compounds",
    options=[]
)

# have a multiple select for the files
wfiles = widgets.SelectMultiple(
    description="files",
    options=[]
)

wfname = widgets.Text(
    description='Atlas Name',
    value='myAtlas',
)

all_files = widgets.Checkbox(
    description='Select All Files',
    value=False,
)

plot_button = widgets.Button(description='Plot me')
create_atlas_btn = widgets.Button(description="Create Atlas")

rtmin_widget = FloatSlider()
rtpeak_widget = FloatSlider()
rtmax_widget = FloatSlider()


###########################################################################
###
def plot_intensity(cval, fvals, rt_min, rt_max, rt_peak):
    for i in range(len(fvals)):
        d = data[files_idx[fvals[i]]][compound_idx[cval]]

        if len(d['data']['eic']['rt']) > 0:
            x = d['data']['eic']['rt']
            y = d['data']['eic']['intensity']
            plt.plot(x, y, 'k-', ms=1, mew=0, mfc='b', alpha=1.0)
    plt.axvline(rt_min, color='b', linewidth=2.0)
    plt.axvline(rt_max, color='g', linewidth=2.0)
    plt.axvline(rt_peak, color='r', linewidth=2.0)


###########################################################################
###
def create_atlas(b):
    identifications = list()
    file_names = wfiles.value
    compound_name = wcompounds.value
    idx2 = compound_idx[compound_name]

    atlas = metob.Atlas()
    atlas.name = wfname.value

    # create an empty rt reference
    rt_ref = metob.RtReference()

    rt_ref.rt_min = rtmin_widget.value
    rt_ref.rt_max = rtmax_widget.value
    rt_ref.rt_peak = rtpeak_widget.value
    rt_ref.rt_units = data[0][idx2]['identification'].rt_references[0].rt_units

    # create an empty mz_reference
    mz_ref = metob.MzReference()

    mz_ref.mz = data[0][idx2]['identification'].mz_references[0].mz
    mz_ref.mz_tolerance = data[0][idx2]['identification'].mz_references[0].mz_tolerance
    mz_ref.mz_tolerance_units = data[0][idx2]['identification'].mz_references[0].mz_tolerance_units
    mz_ref.detected_polarity = data[0][idx2]['identification'].mz_references[0].detected_polarity

    identification = metob.CompoundIdentification()
    identification.compoud = compound_name
    identification.name = compound_name
    identification.rt_references = [rt_ref]
    identification.mz_references = [mz_ref]

    identifications.append(identification)

    atlas.compound_identifications = identifications
    metob.store(atlas)


def select_files(b):
    all_files.value = not all_files.value


###########################################################################
##
def plot_button_clicked(b):
    plt.cla()
    plt.clf()
    plt.close()

    fvals = wfiles.value
    cval = wcompounds.value
    global rtmin_widget, rtmax_widget, rtpeak_widget

    min_x = list()
    max_x = list()

    if all_files.value == True:
        fvals = file_names
    else:
        fvals = wfiles.value

    for i in range(len(fvals)):
        d = data[files_idx[fvals[i]]][compound_idx[cval]]
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
             cval=fixed(cval),
             fvals=fixed(fvals),
             rt_min=rtmin_widget,
             rt_peak=rtpeak_widget,
             rt_max=rtmax_widget)


def dill2atlas(fname):
    data = ma_data.get_dill_data(fname)
    groups = ma_data.get_group_names(data)
    file_names = ma_data.get_file_names(data)
    (compound_names, compound_objects) = ma_data.get_compound_names(data)

    files_idx = dict()
    for f_idx, f_name in enumerate(file_names):
        files_idx[f_name] = f_idx

    compound_idx = dict()
    for cpd_idx, cpd_name in enumerate(compound_names):
        compound_idx[cpd_name] = cpd_idx

    groups_idx = dict()
    for grp_idx, grp_name in enumerate(groups):
        groups_idx[grp_name] = grp_idx





    wcompounds.options=compound_names

    wfiles.options=file_names



    display(widgets.HBox((wfname, create_atlas_btn)))
    display(widgets.HBox((wcompounds, widgets.VBox((all_files, wfiles)))))
    display(plot_button)



    plot_button.on_click(plot_button_clicked)
    create_atlas_btn.on_click(create_atlas)
    all_files.observe(select_files)