#import sys
#from metatlas import metatlas_objects as metob
#from metatlas import h5_query as h5q
#import qgrid
from matplotlib import pyplot as plt
#import pandas as pd
import re
import os
#import tables
#import pickle
import dill
import numpy as np
#from rdkit import Chem
#from rdkit.Chem import AllChem
#from rdkit.Chem import Draw
## from rdkit.Chem.rdMolDescriptors import ExactMolWt
#from rdkit.Chem import Descriptors
#from rdkit.Chem import rdMolDescriptors
#from rdkit.Chem import AllChem
#from rdkit.Chem import Draw
#from rdkit.Chem import rdDepictor
#from rdkit.Chem.Draw import rdMolDraw2D
#from rdkit.Chem.Draw import IPythonConsole
#from IPython.display import SVG,display
from collections import defaultdict
#import time
#from textwrap import wrap
#from matplotlib.backends.backend_pdf import PdfPages
import os.path
from itertools import cycle


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


def get_data(fname):
    """
    Parameters
    ----------
    fname

    Returns
    -------

    """
    with open(my_file,'r') as f:
        data = dill.load(f)

    return data


def get_group_names(data):
    """
    Parameters
    ----------
    data

    Returns
    -------

    """
    group_names = []
    for i,d in enumerate(data):
        newstr = d[0]['group'].name
        group_names.append(newstr)

    return group_names


def get_file_names(data):
    """
    Parameters
    ----------
    data

    Returns
    -------

    """
    file_names = []
    for i,d in enumerate(data):
        newstr = os.path.basename(d[0]['lcmsrun'].hdf5_file)
        file_names.append(newstr)
   
    return file_names


def get_compound_names(data):
    """
    Parameters
    ----------
    data

    Returns
    -------

    """
    compound_names = []
    compound_objects = []
    for i,d in enumerate(data[0]):
        # if label: use label
        # else if compound: use compound name
        # else no name
        compound_objects.append(d['identification'])
        if len(d['identification'].compound) > 0:
            _str = d['identification'].compound[0].name
        else:
            _str = d['identification'].name
        newstr = '%s_%s_%s_%5.2f'%(_str,d['identification'].mz_references[0].detected_polarity,
                d['identification'].mz_references[0].adduct,d['identification'].rt_references[0].rt_peak)
        newstr = re.sub('\.', 'p', newstr) #2 or more in regexp

        newstr = re.sub('[\[\]]','',newstr)
        newstr = re.sub('[^A-Za-z0-9+-]+', '_', newstr)
        newstr = re.sub('i_[A-Za-z]+_i_', '', newstr)
        if newstr[0] == '_':
            newstr = newstr[1:]
        if newstr[0] == '-':
            newstr = newstr[1:]
        if newstr[-1] == '_':
            newstr = newstr[:-1]

        newstr = re.sub('[^A-Za-z0-9]{2,}', '', newstr) #2 or more in regexp
        compound_names.append(newstr)

    #If duplicate compound names exist, then append them with a number
    D = defaultdict(list)
    for i,item in enumerate(compound_names):
        D[item].append(i)
    D = {k:v for k,v in D.items() if len(v)>1}
    for k in D.keys():
        for i,f in enumerate(D[k]):
            compound_names[f] = '%s%d'%(compound_names[f],i)
   
    return (compound_names, compound_objects)


def plot_all_compounds_for_each_file(**kwargs):
    """
    Parameters
    ----------
    kwargs

    Returns
    -------

    """
    data = get_data(kwargs['input_fname'])
    compound_names = get_compound_names(data)[0]
    file_names = get_file_names(data)

    nCols = kwargs['nCols']
    scale_y = kwargs['scale_y']
    output_loc = kwargs['output_loc']
#
#    data = kwargs['data']
#    nCols = kwargs['nCols']
#    file_names = kwargs['file_names']
#    compound_names = kwargs['compound_names']
#    scale_y = kwargs['scale_y']
#    output_fname = kwargs['output_fname']

    nRows = int(np.ceil(len(compound_names)/float(nCols)))
    print nRows
    print len(compound_names) 
    
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
        print my_file
        
        ax = plt.subplot(111)#, aspect='equal')
        plt.setp(ax, 'frame_on', False)
        ax.set_ylim([0, nRows+4])
      
        col = 0
        row = nRows+3
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
                new_x = (x-x[0])*subrange/float(x[-1]-x[0])+col*(subrange+2) ## remapping the x-range
                ax.plot(new_x, new_x-new_x+row,'r-')#,ms=1, mew=0, mfc='b', alpha=1.0)]
                ax.annotate(compound,(min(new_x),row-0.1), size=2)
                col += 1
            counter += 1
        
        plt.title(my_file)
        fig = plt.gcf()
        fig.set_size_inches(11, 8.5)
        #fig.savefig('/home/jimmy/ben2/neg/' + my_file + '-' + str(counter) + '.pdf')
        fig.savefig(os.path.join(output_loc, my_file + '-' + str(counter) + '.pdf'))
        plt.clf()


def plot_all_files_for_each_compound(**kwargs):
    """
    Parameters
    ----------
    kwargs

    Returns
    -------

    """

    data = get_data(kwargs['input_fname'])
    compound_names = get_compound_names(data)[0]
    file_names = get_file_names(data)
    nCols = kwargs['nCols']
    scale_y = kwargs['scale_y']
    output_loc = kwargs['output_loc']

    nRows = int(np.ceil(len(file_names)/float(nCols)))
    print 'nrows = ', nRows 
    
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

    print "length of ymax is ", len(y_max)
    y_max = cycle(y_max)



    # create ouput dir
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)

    for compound_idx,compound in enumerate(compound_names):
        print 10*'*'  
        print compound
        print 10*'*'  

        
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
                ax.annotate('{0},{1},{2},{3}'.format(my_file,rt_min, rt_peak, rt_max),(min(new_x),row-0.2), size=1)#,ha='center')
                myWhere = np.logical_and(new_x>=rt_min_, new_x<=rt_max_ )
                ax.fill_between(new_x,min(y)+row,y+row,myWhere, facecolor='c', alpha=0.3)
                col += 1
            else:
                new_x = (x-x[0])*subrange/float(x[-1]-x[0])+col*(subrange+2) ## remapping the x-range
                ax.plot(new_x, y-y+row,'r-')#,ms=1, mew=0, mfc='b', alpha=1.0)]
                ax.annotate(my_file,(min(new_x),row-0.1), size=1)
                col += 1
            counter += 1
        
        plt.title(compound)
        fig = plt.gcf()
        fig.set_size_inches(11, 8.5)
        #fig.savefig('/tmp/' + compound + '-' + str(counter) + '.pdf')
        fig.savefig(os.path.join(output_loc, my_file + '-' + str(counter) + '.pdf'))
        plt.clf()


if __name__ == '__main__':
    import sys

    input_fname = os.path.expandvars(sys.argv[1])
    output_loc = os.path.expandvars(sys.argv[2])


    nCols = 10
    argument = {'input_fname':input_fname,
                'nCols': nCols,
                'scale_y' : True,
                'output_loc': output_loc
               }

    plot_all_compounds_for_each_file(**argument)
    #plot_all_files_for_each_compound(**argument)




