from collections import defaultdict
import getpass
import numpy as np
import tables
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
from . import (
    retrieve, LcmsRun, to_dataframe, database, get_chromatogram,
    get_spectrogram, get_info
)
from .object_helpers import format_timestamp


def create_qgrid(objects, options=None):
    """Create a qgrid from a list of metatlas objects.
    """
    import qgrid
    qgrid.nbinstall(overwrite=False)

    dataframe = to_dataframe(objects)
    if options:
        defaults = qgrid.grid.defaults.grid_options
        options = defaults.update(options)
    else:
        options = qgrid.grid.defaults.grid_options
    grid = qgrid.grid.QGridWidget(df=dataframe,
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


def show_experiments(username=None):
    """Create a gui to browse the available experiments.
    """
    query = 'SELECT DISTINCT username, experiment FROM lcmsruns'
    entries = [e for e in database.query(query)]
    experiments = defaultdict(set)
    for entry in entries:
        if entry['experiment']:
            experiments[entry['username']].add(entry['experiment'])
            experiments['all'].add(entry['experiment'] or '')

    user = getpass.getuser()
    if user not in experiments:
        user = 'all'
    user_experiments = list(experiments.get(user, {}))

    if user_experiments:
        objects = retrieve('lcmsruns', username=user,
                           experiment=user_experiments[0])
    else:
        objects = [LcmsRun()]

    grid = create_qgrid(objects)

    user_widget = widgets.Dropdown(
        options=list(experiments.keys()),
        value=user,
        description='Username:'
    )
    ex_widget = widgets.Dropdown(
        options=list(experiments[user]),
        value=list(experiments[user])[1],
        description='Experiment:'
    )

    def user_change(trait, value):
        user_widget.value = value
        ex_widget.options = list(experiments[value])
        ex_widget.value = list(experiments[value])[0]

    def experiment_change(trait, value):
        objects = retrieve('lcmsruns', username=user_widget.value,
                           experiment=value)
        df = to_dataframe(objects)
        if df is not None:
            grid.df = df

    user_widget.on_trait_change(user_change, 'value')
    ex_widget.on_trait_change(experiment_change, 'value')

    gui = widgets.VBox([widgets.HBox((user_widget, ex_widget)), grid])
    display(gui)


def show_lcms_run(run, min_mz=None, max_mz=None, polarity=None, ms_level=1):
    """Interact with LCMS data - XIC linked to a Spectrogram plot.

    Parameters
    ----------
    run: LcmsRun
        The object to interact with.
    min_mz: float
        Minimum m/z (defaults to file min)
    max_mz: float
        Maximum m/z (defaults to file max)
    polarity: {0, 1}
        Polarity (defaults to neg. if present in file, else pos.)
    ms_level: {0, 1}
        The ms level.
    """
    fid = tables.open_file(run.hdf5_file)

    info = get_info(fid)
    if polarity is None:
        if info['ms%s_neg' % ms_level]['nrows']:
            polarity = 0
        else:
            polarity = 1
    if polarity == 0:
        table_name = 'ms%s_neg' % ms_level
    else:
        table_name = 'ms%s_pos' % ms_level
    if min_mz is None:
        min_mz = info[table_name]['min_mz']
    if max_mz is None:
        max_mz = info[table_name]['max_mz']

    rt, irt = get_chromatogram(fid, min_mz, max_mz, 1, polarity)
    mz, imz = get_spectrogram(fid, rt[0], rt[1], 1, polarity)

    fig, (ax1, ax2) = plt.subplots(ncols=2)
    ax1.plot(rt, irt)
    ax1.set_title('XIC: %0.1f - %0.1f m/z' % (min_mz, max_mz))
    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('Intensity')
    ax1._vline = ax1.axvline(rt[0])

    ax2.vlines(mz, 0, imz)
    ax2.set_xlabel('Mass (m/z)')
    ax2.set_title('MS%s Spectrogram at %0.1f min' % (ms_level, rt.min()))

    def callback(event):
        if event.inaxes == ax1:
            rt_event = event.xdata
            # get the closest actual RT
            idx = (np.abs(rt - rt_event)).argmin()
            mz, imz = get_spectrogram(fid, rt[idx], rt[idx], 1, polarity)

            ax1._vline.remove()
            ax1._vline = ax1.axvline(rt_event, color='k')

            ax2.clear()
            ax2.vlines(mz, 0, imz)
            ax2.set_xlabel('Mass (m/z)')
            ax2.set_title('Spectrogram at %0.1f min' % rt_event)
            fig.canvas.draw()
    fig.canvas.mpl_connect('button_press_event', callback)
    fig.canvas.draw()


def _edit_object(obj):
    """Create an IPython widget editor for a Metatlas object"""
    try:
        from ipywidgets import Text, Dropdown, HBox, VBox
    except ImportError:
        from IPython.html.widgets import Text, Dropdown, HBox, VBox
    from IPython.display import display
    names = sorted(obj.trait_names())
    names.remove('name')
    names = ['name'] + [n for n in names if not n.startswith('_')]
    items = [Text('', disabled=True)]
    for name in names:
        if name.startswith('_'):
            continue
        try:
            value = getattr(obj, name)
        except traitlets.TraitError:
            value = None
        trait = obj.traits()[name]
        if name in ['created', 'last_modified']:
            value = format_timestamp(value)
        if (trait.get_metadata('readonly') or
                isinstance(trait, (traitlets.Instance, traitlets.List)) or
                value is None):
            if isinstance(trait, traitlets.Instance):
                value = value.unique_id
            elif isinstance(trait, traitlets.List):
                value = [o.unique_id for o in value]
            items.append(Text(str(value), disabled=True))

        elif isinstance(trait, traitlets.Enum):
            # create a closure around "name" for the on_trait_change
            # callback
            def create_dropdown(name):
                dd = Dropdown(value=value, options=trait.values)

                def callback(dummy, value):
                    setattr(obj, name, value)
                dd.on_trait_change(callback, 'value')
                items.append(dd)

            create_dropdown(name)
        else:
            # create a closure around "name" for the on_trait_change
            # callback
            def create_textbox(name):
                textbox = Text(str(value))

                def callback(dummy, value):
                    try:
                        setattr(obj, name, value)
                    except Exception:
                        textbox.color = 'red'
                textbox.on_trait_change(callback, 'value')
                items.append(textbox)

            create_textbox(name)

    labels = [Text(name, disabled=True) for name in names]
    labels = [Text(obj.__class__.__name__, disabled=True)] + labels
    display(HBox(children=[VBox(children=labels), VBox(children=items)]))


def edit_objects(objects):
    """Edit an object or set of Metatlas in a QGrid.

    If more than one object is given, a grid will be shown, otherwise
    a single editor will be shown.
    """
    if isinstance(object, traitlets.HasTraits):
        return _edit_object(objects)
    if len(objects) == 1:
        return _edit_object(objects[0])

    grid = create_qgrid(objects)
    add_row = widgets.Button(description="Add Row")
    add_row.on_click(grid.add_row)

    rem_row = widgets.Button(description="Remove Row")
    rem_row.on_click(grid.remove_row)

    export = widgets.Button(description="Export")
    export.on_click(grid.export)

    display(widgets.HBox((add_row, rem_row, export)), grid)
