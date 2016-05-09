import sys

# os.environ['R_LIBS_USER'] = '/project/projectdirs/metatlas/r_pkgs/'
curr_ld_lib_path = ''

# os.environ['LD_LIBRARY_PATH'] = curr_ld_lib_path + ':/project/projectdirs/openmsi/jupyterhub_libs/boost_1_55_0/lib' + ':/project/projectdirs/openmsi/jupyterhub_libs/lib'

# sys.path.insert(0, '/project/projectdirs/metatlas/python_pkgs/')
sys.path.insert(0,'/global/project/projectdirs/metatlas/anaconda/lib/python2.7/site-packages' )

from metatlas import metatlas_objects as metob
from metatlas import h5_query as h5q
sys.path.append('/global/project/projectdirs/openmsi/jupyterhub_libs/anaconda/lib/python2.7/site-packages')

import metatlas_get_data_helper_fun as ma_data
ma_data = reload(ma_data)

import qgrid

from matplotlib import pyplot as plt
import pandas as pd
import os
import tables
import pickle


import dill

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
# from rdkit.Chem.rdMolDescriptors import ExactMolWt
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG,display




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
    with open(fname,'r') as f:
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
    data = get_data(os.path.expandvars(kwargs['input_fname']))
    compound_names = get_compound_names(data)[0]
    file_names = get_file_names(data)

    nCols = kwargs['nCols']
    scale_y = kwargs['scale_y']
    output_loc = os.path.expandvars(kwargs['output_loc'])
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
                new_x = np.asarray([0,1])#(x-x[0])*subrange/float(x[-1]-x[0])+col*(subrange+2) ## remapping the x-range
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

    data = get_data(os.path.expandvars(kwargs['input_fname']))
    compound_names = get_compound_names(data)[0]
    file_names = get_file_names(data)
    nCols = kwargs['nCols']
    scale_y = kwargs['scale_y']
    output_loc = os.path.expandvars(kwargs['output_loc'])

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
                ax.fill_between(new_x,min(y)+row,y+row,myWhere, facecolor='c', alpha=0.3)
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
        fig.set_size_inches(11, 8.5)
        #fig.savefig('/tmp/' + compound + '-' + str(counter) + '.pdf')
        fig.savefig(os.path.join(output_loc, compound + '-' + str(counter) + '.pdf'))
        plt.clf()


        
  


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

def make_output_dataframe(**kwargs):
    data = get_data(os.path.expandvars(kwargs['input_fname']))
    compound_names = get_compound_names(data)[0]
    file_names = get_file_names(data)
    group_names = get_group_names(data)
    output_loc = os.path.expandvars(kwargs['output_loc'])
    fieldname = kwargs['fieldname']
    
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)
    
    df = pd.DataFrame( index=compound_names, columns=file_names, dtype=float)

    # peak_height['compound'] = compound_list
    # peak_height.set_index('compound',drop=True)
    for i,dd in enumerate(data):
        for j,d in enumerate(dd):
            if not d['data']['ms1_summary'][fieldname]:
                df.ix[compound_names[j],file_names[i]] = 0
            else:
                df.ix[compound_names[j],file_names[i]] = d['data']['ms1_summary'][fieldname]  
    columns = []
    for i,f in enumerate(file_names):
        columns.append((group_names[i],f))
    df.columns = pd.MultiIndex.from_tuples(columns,names=['group', 'file'])

    df.to_csv(os.path.join(output_loc, fieldname + '.tab'),sep='\t')
    return df

def file_with_max_precursor_intensity(data,compound_idx):
    idx = []
    my_max = 0
    for i,d in enumerate(data):
        if type(d[compound_idx]['data']['msms']['data']) != list:#.has_key('precursor_intensity'):
            temp = d[compound_idx]['data']['msms']['data']
            m = np.max(temp['precursor_intensity'])
            if m > my_max:
                my_max = m
                idx = i
    return idx,my_max

def plot_errorbar_plots(df,**kwargs):#df,compound_list,project_label):
    
    data = get_data(os.path.expandvars(kwargs['input_fname']))
    compound_names = get_compound_names(data)[0]
    file_names = get_file_names(data)
    output_loc = os.path.expandvars(kwargs['output_loc'])
    
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)
        

    for compound in compound_names:
        m = df.ix[compound].groupby(level='group').mean()
        e = df.ix[compound].groupby(level='group').std()
        c = df.ix[compound].groupby(level='group').count()

        for i in range(len(e)):
            if c[i]>0:
                e[i] = e[i] / c[i]**0.5

        f, ax = plt.subplots(1, 1,figsize=(20,12))
        m.plot(yerr=e, kind='bar',ax=ax)
        ax.set_title(compound,fontsize=12,weight='bold')
        plt.tight_layout()
        f.savefig(os.path.join(output_loc, compound + '_errorbar.pdf'))

        f.clear()
        plt.close('all')#f.clear()

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
    
        

def make_identification_figure(**kwargs):#data,file_idx,compound_idx,export_name,project_label):
    #  d = 'data/%s/identification/'%project_label
    
    data = get_data(os.path.expandvars(kwargs['input_fname']))
    compound_names = get_compound_names(data)[0]
    file_names = get_file_names(data)
    output_loc = os.path.expandvars(kwargs['output_loc'])

    ids = metob.retrieve('CompoundIdentification')
    frag_refs = [cid for cid in ids if cid.frag_references]
    
    
    
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)
    
    
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
            ref_spec = get_reference_msms_spectra(frag_refs, 
                                       compound_name = data[file_idx][compound_idx]['identification'].compound[0].name, 
                                       polarity = data[file_idx][compound_idx]['data']['msms']['polarity'])
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
                print s
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
            fig.clear()
            plt.close('all')#f.clear()

    
def match_inchi_key_to_lookup_table(df,compound_lookup = '/global/homes/b/bpb/notebooks/thoughts/uniquecompounds.csv'):
    '''
    Until the compound database is updated use this file to check for portable IDs and common names.
    Takes as input a dataframe with inchi_key as a column and adds the columns from the compound lookup table
    '''
    import re
    lookup_df = pd.read_csv(compound_lookup)
    for i, row in df.iterrows():
        if row['neutralized_inchi_key']:
            idx = lookup_df.metatlas_inchi_key == row['neutralized_inchi_key']
            for j,idx_val in enumerate(idx):
                if idx_val:
                    for k in lookup_df.keys():
                        s = str(lookup_df.loc[j,k])
                        s = re.sub('<[^>]*>', '', s)
                        df.loc[i,k] = s
    return df
    
    
            
def export_atlas_to_spreadsheet(myAtlas,output_filename,input_type = 'atlas'):
    # myAtlases = [atlas[0],atlas[1]] #concatenate the atlases you want to use
    # myAtlases = [atlas[0]]

    cols = ['inchi',
     'mono_isotopic_molecular_weight',
     'creation_time',
     'description',
     'formula',
     'functional_sets',
     'last_modified',
     'reference_xrefs',
     'synonyms',
     'unique_id',
     'url',
     'username']

        # print myAtlas[0].compound_identifications[0].compound
    atlas_export = pd.DataFrame( )

#     atlas_export['name'] = compound_list
#     atlas_export.set_index('name',drop=True)
    if input_type != 'atlas':
        num_compounds = len(myAtlas[0])
    else:
        num_compounds = len(myAtlas.compound_identifications)
    for i in range(num_compounds):
        if input_type != 'atlas':
            my_id = myAtlas[0][i]['identification']
        else:
            my_id = myAtlas.compound_identifications[i]
    
        if myAtlas.compound_identifications[i].compound:
            n = my_id.compound[0].name
        else:
            n = my_id.name
        atlas_export.loc[i,'name'] = n
        if my_id.compound:
            for c in cols:
                g = getattr(my_id.compound[0],c)
                if g:
                    atlas_export.ix[i,c] = getattr(my_id.compound[0],c)
        atlas_export.loc[i, 'label'] = my_id.name
        atlas_export.loc[i,'rt_min'] = my_id.rt_references[0].rt_min
        atlas_export.loc[i,'rt_max'] = my_id.rt_references[0].rt_max
        atlas_export.loc[i,'rt_peak'] = my_id.rt_references[0].rt_peak
        atlas_export.loc[i,'mz'] = my_id.mz_references[0].mz
        atlas_export.loc[i,'mz_tolerance'] = my_id.mz_references[0].mz_tolerance
        atlas_export.loc[i,'polarity'] = my_id.mz_references[0].detected_polarity
        if my_id.frag_references:
            atlas_export.loc[i,'has_fragmentation_reference'] = True
        else:
            atlas_export.loc[i,'has_fragmentation_reference'] = False
    
    for i,row in atlas_export.iterrows():
        mol= []
        if row['inchi'] and (type(row['inchi']) != float):
            mol = Chem.MolFromInchi(row['inchi'].encode('utf-8'))
        if mol:
            ds = desalt(mol)
            c = NeutraliseCharges(ds[0])
            mol = c[0]
            atlas_export.loc[i,'permanent_charge'] = Chem.GetFormalCharge(mol)
            
            neutral_string = Chem.MolToInchi(mol)
            atlas_export.loc[i,'neutralized_inchi'] = neutral_string
            
            neutral_inchi_key = Chem.InchiToInchiKey(neutral_string)
            atlas_export.loc[i,'neutralized_inchi_key'] = neutral_inchi_key
            
    atlas_export = match_inchi_key_to_lookup_table(atlas_export)
    
    
    if not os.path.exists(os.path.dirname(output_filename)):
        os.makedirs(os.path.dirname(output_filename))
    
    atlas_export.to_csv(output_filename)
    return atlas_export
    
def get_data_for_groups_and_atlas(group,myAtlas,output_filename,use_set1 = False):
    data = []
    import copy as copy
    for i,treatment_groups in enumerate(group):
        for j in range(len(treatment_groups.items)):
            myFile = treatment_groups.items[j].hdf5_file
    #         try:
    #             rt_reference_index = int(treatment_groups.name[-1]) - 1
    #         except:
    #             rt_reference_index = 3
            print i,len(group),myFile
            row = []
            for compound in myAtlas.compound_identifications:
                result = {}
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

def get_metatlas_atlas(name = '%%',most_recent = True,do_print = True):
    from datetime import datetime, date
    atlas = metob.retrieve('Atlas',name = name,username='*')
    if most_recent:
        atlas = filter_metatlas_objects_to_most_recent(atlas,'name')
    for i,a in enumerate(atlas):
        print i, a.name,  datetime.utcfromtimestamp(a.last_modified)

    return atlas


def get_metatlas_files(experiment = '%%',name = '%%',most_recent = True):
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
    for x in df.index:
        if type(df.name[x]) != float or type(df.label[x]) != float:
            if type(df.name[x]) != float:
                    if not metob.retrieve('Compounds',name=df.name[x], username = '*'):
                        print df.name[x], "is not in database"
                        bad_names.append(df.name[x])
    return bad_names


def check_file_names(df,field):
    bad_files = []
    for i,row in df.iterrows():
        if not metob.retrieve('Lcmsruns',name = '%%%s%%'%row[field],username = '*'):
            print row[field], "is not in the database"
            bad_files.append(row[field])
    return bad_files


def get_formatted_atlas_from_google_sheet(polarity='POS',method='QE_HILIC',mz_tolerance=10):
    sys.path.insert(0,'/project/projectdirs/metatlas/projects/ms_monitor_tools/' )
    import ms_monitor_util as mmu
    df = mmu.get_ms_monitor_reference_data()
    df2 = pd.DataFrame(df[1:],columns=df[0])

    fields_to_keep = [ 'name',
                    'label',
                    'mz_%s'%polarity,
                    'rt_min_%s'%method,
                    'rt_max_%s'%method,
                    'rt_peak_%s'%method,
                    'file_mz_%s_%s'%(method,polarity),
                    'file_rt_%s_%s'%(method,polarity),
                    'file_msms_%s_%s'%(method,polarity)]
    df3 = df2.loc[:,fields_to_keep]

    df3['mz_tolerance'] = mz_tolerance

    if polarity == 'POS':
        df3['polarity'] = 'polarity'
    else:
        df3['polarity'] = 'negative'

    renamed_columns = [c.replace('_%s'%method,'').replace('_%s'%polarity,'') for c in df3.columns]
    for i,c in enumerate(df3.columns):
        df3 = df3.rename(columns = {c:renamed_columns[i]})
    return df3


def make_atlas_from_spreadsheet(filename=False,atlas_name='temp',filetype='excel',sheetname='',polarity = 'positive', store=False,mz_tolerance=False,dataframe=None):
    '''
    specify polarity as 'positive' or 'negative'
    
    '''
    if isinstance(dataframe,pd.DataFrame):
        df = dataframe
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
    if 'file_rt' in df.keys():
        bad_files = check_file_names(df,'file_rt')
        if bad_files:
             return bad_files
    if 'file_mz' in df.keys():
        bad_files = check_file_names(df,'file_mz')
        if bad_files:
             return bad_files
    if 'file_msms' in df.keys():
        bad_files = check_file_names(df,'file_msms')
        if bad_files:
             return bad_files
    

    
    all_identifications = []

#     for i,row in df.iterrows():
    for x in df.index:
        if type(df.name[x]) != float or type(df.label[x]) != float: #this logic is to skip empty rows
            
            myID = metob.CompoundIdentification()
            
            if type(df.name[x]) != float: # this logic is where a name has been specified
                c = metob.retrieve('Compounds',name=df.name[x],username = '*') #currently, all copies of the molecule are returned.  The 0 is the most recent one. 
                if c:
                    c = c[0]
            else:
                c = 'use_label'
            if type(df.label[x]) != float:
                compound_label = df.label[x] #if no name, then use label as descriptor
            else:
                compound_label = 'no label'
            
            if c:
                if c != 'use_label':
                    myID.compound = [c]
                myID.name = compound_label
                
                
                mzRef = metob.MzReference()
                # take the mz value from the spreadsheet
                mzRef.mz = df.mz[x]
                #TODO: calculate the mz from theoretical adduct and modification if provided.
                #     mzRef.mz = c.MonoIso topic_molecular_weight + 1.007276
                if mz_tolerance:
                    mzRef.mz_tolerance = mz_tolerance
                else:
                    try:
                        mzRef.mz_tolerance = df.mz_tolerance[x]
                    except:
                        mzRef.mz_tolerance = df.mz_threshold[x]    
                
                mzRef.mz_tolerance_units = 'ppm'
                mzRef.detected_polarity = polarity
                if 'file_mz' in df.keys():
                    f = metob.retrieve('Lcmsruns',name = '%%%s%%'%df.file_mz[x],username = '*')[0]
                    mzRef.lcms_run = f
                #     mzRef.adduct = '[M-H]'   
                myID.mz_references = [mzRef]

                rtRef = metob.RtReference()
                rtRef.rt_units = 'min'
                rtRef.rt_min = df.rt_min[x]
                rtRef.rt_max = df.rt_max[x]
                rtRef.rt_peak = df.rt_peak[x]
                if 'file_rt' in df.keys():
                    f = metob.retrieve('Lcmsruns',name = '%%%s%%'%df.file_rt[x],username = '*')[0]
                    rtRef.lcms_run = f
                myID.rt_references = [rtRef]
                    
                if 'file_msms' in df.keys():
                    if type(df.file_msms[x]) != float:
                        frag_ref = metob.FragmentationReference()
                        f = metob.retrieve('Lcmsruns',name = '%%%s%%'%df.file_msms[x],username = '*')[0]
                        frag_ref.lcms_run = f
                        frag_ref.polarity = polarity
                        frag_ref.precursor_mz = df.mz[x]
                        
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
        if (len(getattr(g,field)) > 0):
            filtered_list.append(g)
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

      
def select_groups_for_analysis(name = '%%',do_print = True, most_recent = True, remove_empty = True, filter_list = [], exclude_list = []):
    groups = metob.retrieve('Groups', name = name, username='*')
    if most_recent:
        groups = filter_metatlas_objects_to_most_recent(groups,'name')
    
    if filter_list:
        groups = filter_metatlas_objects_by_list(groups,'name',filter_list)
        
    if exclude_list:
        groups = remove_metatlas_objects_by_list(groups,'name',exclude_list)
    
    if remove_empty:
        groups = filter_empty_metatlas_objects(groups,'items')
    if do_print:
        from datetime import datetime, date
        for i,a in enumerate(groups):
            print i, a.name,  datetime.utcfromtimestamp(a.last_modified)

    return groups

if __name__ == '__main__':
    import sys

    input_fname = os.path.expandvars(sys.argv[1])
    output_loc = os.path.expandvars(sys.argv[2])



    nCols = 10
    argument = {'input_fname':input_fname,
                'nCols': nCols,
                'scale_y' : False,
                'output_loc': output_loc
               }

    plot_all_compounds_for_each_file(**argument)
    argument = {'input_fname':input_fname,
                'nCols': 20,
                'scale_y' : False,
                'output_loc': '/home/jimmy/ben/neg/unscaled/allfiles'
                }
    plot_all_files_for_each_compound(**argument)




