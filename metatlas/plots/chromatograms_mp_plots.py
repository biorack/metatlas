import matplotlib 
#matplotlib.use('Agg')
import sys
#import yaml
import os
#import multiprocessing as mp
from matplotlib import pyplot as plt
import numpy as np
import warnings
from textwrap import wrap
#warnings.filterwarnings("ignore")

from metatlas.io import metatlas_get_data_helper_fun as ma_data

def plot_chromatogram(d,file_name, ax=None):
    """
    
    """
    if ax is None:
        ax = plt.gca()

    rt_min = d['identification'].rt_references[0].rt_min
    rt_max = d['identification'].rt_references[0].rt_max
    rt_peak = d['identification'].rt_references[0].rt_peak
        
    if len(d['data']['eic']['rt']) > 1:
        x = np.asarray(d['data']['eic']['rt'])
        y = np.asarray(d['data']['eic']['intensity'])

        ax.plot(x,y,'k-',linewidth=2.0,alpha=1.0)  
        myWhere = np.logical_and(x>=rt_min, x<=rt_max )
        ax.fill_between(x,0,y,myWhere, facecolor='c', alpha=0.3)

    ax.axvline(rt_min, color='k',linewidth=2.0)
    ax.axvline(rt_max, color='k',linewidth=2.0)
    ax.axvline(rt_peak, color='r',linewidth=2.0)
    ax.set_title("\n".join(wrap(file_name,54)),fontsize=12,weight='bold')


def plot_compounds_and_files_mp(kwargs):
    #print(mp.current_process())

    my_data= kwargs['data'] # data for all compounds for one file
    file_name  = kwargs['file_name'] # full path of output file name
    nRows, nCols = kwargs['rowscols']
    names = kwargs['names']
    share_y = kwargs['share_y']
    # plt.ioff()

    f,ax = plt.subplots(nRows, nCols, sharey=share_y,figsize=(8*nCols,nRows*6))
    ax = ax.flatten()
    plt.rcParams['pdf.fonttype']=42
    plt.rcParams['pdf.use14corefonts'] = True
#    matplotlib.rc('font', family='sans-serif') 
#    matplotlib.rc('font', serif='Helvetica') 
    plt.rcParams['text.usetex'] = False
    plt.rcParams.update({'font.size': 12})
    plt.rcParams.update({'font.weight': 'bold'})
    plt.rcParams['axes.linewidth'] = 2 # set the value globally
    
    for i,name in enumerate(names):
        plot_chromatogram(my_data[i], name, ax=ax[i])

    f.savefig(file_name)    
    plt.close(f)


def plot_compounds_and_files(output_dir,
                             data,
                             nCols = 8,
                             share_y = False,
                             pool=None,
                             plot_types='both'):
    '''
    Parameters
    ----------
    output_dir location of saved pdf plots
    nCols number of columns per pdf file
    share_y subplots share/not share they y axis
    processes number of cores to use
    plot_types compounds per file or files per compound or both

    Returns
    -------
    nothing
    '''

    file_names = ma_data.get_file_names(data)
    compound_names = ma_data.get_compound_names(data)[0]

    # create directory if necessary
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # setup the parameters according to the request
    if 'files' in plot_types.lower():
        nRows = int(np.ceil(len(compound_names)/float(nCols)))
        args_list = []
        for file_idx, my_file in enumerate(file_names):
            kwargs = {'data': data[file_idx],
                      'file_name': os.path.join(output_dir, my_file +'.pdf'),
                      'rowscols': (nRows, nCols),
                      'share_y': share_y,
                      'names': compound_names}
            args_list.append(kwargs)

    if 'compounds' in plot_types.lower():
        nRows = int(np.ceil(len(file_names)/float(nCols)))
        args_list = []
        for compound_idx, my_compound in enumerate(compound_names):
            my_data = list()
            for file_idx, my_file in enumerate(file_names):
                my_data.append(data[file_idx][compound_idx])

            kwargs = {'data': my_data,
                      'file_name': os.path.join(output_dir, my_compound+'.pdf'),
                      'rowscols': (nRows, nCols),
                      'share_y': share_y,
                      'names': file_names}
            args_list.append(kwargs)

    pool.map(plot_compounds_and_files_mp, args_list)


#if __name__ == '__main__':
#    #sys.path.insert(0, '/global/homes/j/jtouma/metatlas')
#    import pickle
#
#    # load pickled data
#    info = pickle.load(open(sys.argv[1], "rb"))
#    sys.path.insert(info['path_idx'], info['path'])
#    from metatlas.helpers import metatlas_get_data_helper_fun as ma_data
#    data = ma_data.get_dill_data(info['pkl_file'])
#    file_names = ma_data.get_file_names(data)
#    compound_names = ma_data.get_compound_names(data)[0]
#
#    print("\n")
#    print(50*'-')
#    print("Number of file: " + str(len(file_names)))
#    print("Number of compounds: " + str(len(compound_names)))
#    if info['plot_types'].lower() == 'both':
#        print("Processing both files and compounds")
#    else:
#        print("processing " + info['plot_types'].lower() + " only")
#    print("Using " + str(info['processes']) + " out of " + str(mp.cpu_count()) + " available cores")
#    print(50*'-')
#    print("\n")
#    plot_compounds_and_files(output_dir=info['output_dir'],
#                             data=data,
#                             compound_names=compound_names,
#                             file_names=file_names,
#                             nCols=info['nCols'],
#                             share_y=info['share_y'],
#                             processes=info['processes'],
#                             plot_types=info['plot_types'])
#
