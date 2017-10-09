from matplotlib import pyplot as plt
from matplotlib import collections  as mc
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import warnings
from time import time
from textwrap import wrap
from sys import maxsize

#######################
#letter gen
#######################
def letter_gen():
    label = 'A'
    yield label
    
    while True:
        label_noz = label.rstrip('Z')
        label_list = [ord(l) for l in label_noz]  

        if len(label_noz) != 0:
            label_list[-1] += 1
        else:
            label_list.append(65)

        for i in range(len(label) - len(label_noz)):
            label_list.append(65)

        label = ''.join(chr(l) for l in label_list)
        yield label


#######################
#normalize_data
#######################
def normalize_data(data, Names, x_offset, y_offset, sub_x, sub_y):
    
    X_min = maxsize
    X_max = 0
    Y_min = maxsize
    Y_max = 0
    
    for d in data:
        for r in d['data']['eic']['rt']:
            if np.amin(np.asarray(r)) < X_min:
               X_min = r
            if np.amax(np.asarray(r)) > X_max:
               X_max = r
        for i in d['data']['eic']['intensity']:
            if np.amin(np.asarray(i)) < Y_min:
               Y_min = i
            if np.amax(np.asarray(i)) > Y_max:
               Y_max = i
    
    x_mag = np.log10(X_max - X_min)
    if np.isinf(x_mag) or np.isnan(x_mag):
        x_mag = 0
    else:
        x_mag = int(x_mag)
        
    y_mag = np.log10(Y_max - Y_min) - 1
    if np.isinf(y_mag) or np.isnan(y_mag):
        y_mag = 0
    else:
        y_mag = int(y_mag)

    x_range = (max(0, np.round(X_min - .5*np.power(10, x_mag), int(-1*x_mag))), np.round(X_max + .5*np.power(10, x_mag), int(-1*x_mag)))
    y_range = (max(0, np.round(Y_min - .5*np.power(10, y_mag), int(-1*y_mag))), np.round(Y_max + .5*np.power(10, y_mag), int(-1*y_mag)))
    
    def scale_x(x):
        return (x_offset - (x_offset - sub_x))*(x - x_range[0])/(x_range[1] - x_range[0])
    
    def scale_y(y):
        return (y_offset - (y_offset - sub_y))*(y - y_range[0])/(y_range[1] - y_range[0])
       
    Groups = [d['group'].name for d in data]
    X =  [scale_x(np.asarray(d['data']['eic']['rt'])) for d in data]
    Y = [scale_y(np.asarray(d['data']['eic']['intensity'])) for d in data]
    RT_mins = [scale_x(d['identification'].rt_references[0].rt_min) for d in data]
    RT_maxs = [scale_x(d['identification'].rt_references[0].rt_max) for d in data]
    RT_peaks = [scale_x(d['identification'].rt_references[0].rt_peak) for d in data]
        
    return {'x_range': x_range,
            'y_range': y_range,
            'y_mag': y_mag,
            'data': np.array(zip(Groups, Names, X, Y, RT_mins, RT_maxs, RT_peaks), 
                             dtype=[('group', 'object'),('name', 'object'),('x','object'),('y','object'),('rt_min','float'),('rt_max','float'),('rt_peak','float')])
           }
    

#######################
#chromplotplus
#######################
def chromplotplus(kwargs):
    
    file_name = kwargs['file_name']     
    
    ##Options
    warnings.simplefilter('ignore', FutureWarning)
    share_y = kwargs['share_y']
    group = kwargs['group']
    
    #Subplot size and seperations
    sub_x = 8
    sub_y = 6
    x_offset = 9
    y_offset = 7
    
    #Normalize data
    norm = normalize_data(kwargs['data'], kwargs['names'], y_offset, x_offset, sub_x, sub_y)
    
    x_range = norm['x_range']
    y_range = norm['y_range']
    y_mag = norm['y_mag']
    data = norm['data']

    #Groups
    groups = {} # stores "group name: [index of first data of group, # of data belonging to group]"
    if group == 'page' or group == 'index' or group == 'sort':
        data = sorted(data, key=lambda d: d['group'])
        for i, d in enumerate(data):
            if groups.get(d['group']) == None:
                groups[d['group']] = [i, 1]
            else:
                groups[d['group']][1] += 1
        
        
    #Subplot arrangement
    n_plots_list = []
    n_rows_list = []
    n_cols_list = []
    
    if group == 'page':
        for g in sorted(groups.keys()):
            n_plots_list.append(groups[g][1])
            n_rows_list.append(int(np.ceil(np.sqrt(13.0*(groups[g][1])/11))))
            n_cols_list.append(int(np.ceil((groups[g][1])/float(n_rows_list[-1]))))
    elif group == 'index':
        n_plots_list.append(len(data))  
        n_rows_list.append(int(np.ceil(np.sqrt(13.0*(n_plots_list[0]+len(groups))/11))))
        n_cols_list.append(int(np.ceil((n_plots_list[0]+len(groups))/float(n_rows_list[0]))))
    else:
        n_plots_list.append(len(data))  
        n_rows_list.append(int(np.ceil(np.sqrt(13.0*n_plots_list[0]/11))))
        n_cols_list.append(int(np.ceil(n_plots_list[0]/float(n_rows_list[0]))))
    
    #Hashmark variables
    hash_n = 5
    hash_m = 5
    hash_l = .02*min(sub_x, sub_y)
    
    #Axis values
    x_values = np.linspace(x_range[0], x_range[1], num=hash_n)
    y_values = np.linspace(y_range[0]/np.power(10, y_mag), y_range[1]/np.power(10, y_mag), num=hash_m)
    
    def plot():
        #Plot creation
        fig = plt.figure()
        ax = plt.subplot(111)
        plt.setp(ax, 'frame_on', False)
        plt.ioff()
        ax.set_ylim([sub_y - y_offset, (n_cols)*y_offset + (y_offset - sub_y)])
        ax.set_xlim([sub_x - x_offset, (n_rows)*x_offset + (x_offset - sub_x)])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid('off')
        ax.margins(1)
        
        #Group title
        if group == 'page':
            plt.title("\n".join(wrap(g,54)), size = 12., weight='bold')
    
        #Coordinate lists for lines to be drawn
        boxes = []
        hashes = []
        rt_edges = []
        rt_peaks = []      
        
        #Letter generator for labeling groups
        if group == 'index':
            lg = letter_gen()
    
        c = 0 #counter for plots
        for j in np.arange(n_cols - 1, -1, -1):
            for i in np.arange(n_rows):
                
                #Repeating calculations
                x_delta = i*x_offset
                y_delta = j*y_offset
               
                #Break if everything is plotted
                if c >= n_plots:
                    break
                    
                #Label new groups
                if group == 'index' and c in [v[0] for v in groups.values() if v is not None]:
                    ax.annotate("\n".join(wrap(next(lg),54)), 
                            (x_delta + .5*sub_x, y_delta + .5*sub_y - (y_offset - sub_y)), 
                            ha='center', va='center', size = 240./(n_cols+.25), weight='bold')
                    
                    for k,v in groups.items():
                        if v is not None and v[0] == c:
                            groups[k][0] = None
                    continue
                
                #Retention Times
                rt_min = d[c]['rt_min'] + x_delta
                rt_max = d[c]['rt_max'] + x_delta
                rt_peak = d[c]['rt_peak'] + x_delta
                
                rt_edges.append([(rt_min, y_delta), (rt_min, sub_y + y_delta)])
                rt_edges.append([(rt_max, y_delta), (rt_max, sub_y + y_delta)])
                rt_peaks.append([(rt_peak, y_delta), (rt_peak, sub_y + y_delta)])
                
                #Data
                if len(d[c]['x']) > 1:
                    x = d[c]['x'] + x_delta
                    y = d[c]['y'] + y_delta
                    ax.plot(x, y, 'k-',linewidth=2.0/np.sqrt(n_cols*n_rows),alpha=1.0)
                    myWhere = np.logical_and(x>=rt_min, x<=rt_max )
                    ax.fill_between(x, y_delta, y, myWhere, facecolor='c', edgecolor='c', linewidth=0, alpha=0.3)
                
                #Boxes
                boxes.append([(x + x_delta, y + y_delta) for x,y in 
                              [(0, 0), (sub_x, 0), (sub_x, sub_y), (0, sub_y), (0, 0)]])
                
                #Subplot Titles
                ax.annotate("\n".join(wrap(d[c]['name'],48)), 
                            (x_delta + .5*sub_x, y_delta + sub_y + .1*(y_offset - sub_y)), 
                            ha='center', size = 8./(n_cols+.25), weight='bold')
                
                #Hashmarks and associated labels
                for axis in ['bottom', 'left', 'right', 'top']:
                    
                    #Horizontal range
                    for k in range(0, hash_n):
                        if axis == 'bottom':
                            start = (k*(1.0/hash_n)*sub_x + x_delta, y_delta)
                            end = (k*(1.0/hash_n)*sub_x + x_delta, hash_l + y_delta)
                            
                            #X-axis labels
                            ax.annotate(x_values[k], (start[0], start[1] - .15*(y_offset - sub_y)), 
                                        ha='center', size = 8./(n_cols+.25))
                            
                            hashes.append([start, end])
                            
                        if axis == 'top':
                            start = (k*(1.0/hash_n)*sub_x + x_delta, sub_y + y_delta)
                            end = (k*(1.0/hash_n)*sub_x + x_delta, sub_y - hash_l + y_delta)
                            
                            #Y-axis magnitude labels
                            if k == 0 and (share_y == False or i == 0):
                                ax.annotate('1e{}'.format(y_mag), (start[0], start[1] + .1*(y_offset - sub_y)), 
                                            ha='center', va='center', size = 8./(n_cols+.25))
                                   
                            hashes.append([start, end])
                    
                    #Vertical range
                    for l in range(0, hash_m):
                        if axis == 'left':
                            start = (x_delta, l*(1.0/hash_m)*sub_y + y_delta)
                            end = ((hash_l + x_delta, l*(1.0/hash_m)*sub_y + y_delta))
                            
                            #Y-axis labels for leftmost subplots
                            if share_y == False or i == 0:
                                ax.annotate(y_values[l], (start[0] - .15*(x_offset - sub_x), start[1]), 
                                            ha='right', va='center', size = 8./(n_cols+.25))
    
                            hashes.append([start, end])
                                   
                        if axis == 'right':
                            start = (sub_x + x_delta, l*(1.0/hash_m)*sub_y + y_delta)
                            end = (sub_x - hash_l + x_delta, l*(1.0/hash_m)*sub_y + y_delta)
                                   
                            hashes.append([start, end])
                
                c += 1 #increment plot counter
                            
        #Make line colelctions
        bc = mc.LineCollection(boxes, colors=(0,0,0,1), linewidths=1.0/np.sqrt(n_cols*n_rows))
        hc = mc.LineCollection(hashes, colors=(0,0,0,1), linewidths=1.0/np.sqrt(n_cols*n_rows))
        rc = mc.LineCollection(rt_edges, colors=(0,0,0,1), linewidths=2.0/np.sqrt(n_cols*n_rows))
        pc = mc.LineCollection(rt_peaks, colors=(1,0,0,1), linewidths=2.0/np.sqrt(n_cols*n_rows))
        
        #Draw lines
        ax.add_collection(bc)
        ax.add_collection(hc)
        ax.add_collection(rc)
        ax.add_collection(pc)
        
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['pdf.use14corefonts'] = True
        plt.rcParams['text.usetex'] = False
        pdf.savefig()
        plt.close()
        
            
    with PdfPages(file_name) as pdf:
        if group == 'page':
            for i, g in enumerate(sorted(groups.keys())):
                n_rows = n_rows_list[i]
                n_cols = n_cols_list[i]
                n_plots = n_plots_list[i]
                d = data[groups[g][0]:groups[g][0] + groups[g][1]]
                
                plot()
        else:
            n_rows = n_rows_list[0]
            n_cols = n_cols_list[0]
            n_plots = n_plots_list[0]
            d = data
            
            plot()
