from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import collections as mc
import numpy as np
from scipy.interpolate import interp1d
from textwrap import wrap

def chromplotplus(kwargs):
    ChromPlotPlus(**kwargs)

class CompoundFileEIC:
    def __init__(self, compound_file_data, rt_bounds, rt_min, rt_max, shortname):
        self.group_name = compound_file_data['group'].name
        self.file_name = compound_file_data['lcmsrun'].name
        if not shortname.empty:
            self.file_name = shortname.loc[compound_file_data['lcmsrun'].name.split('.')[0], 'shortname'][0]
        try:
            self.eic = np.asarray([compound_file_data['data']['eic']['rt'],
                                   compound_file_data['data']['eic']['intensity']]).astype(float)
            self.eic = self.eic[:, (rt_min <= self.eic[0]) & (self.eic[0] <= rt_max)]

            try:
                self.intensity_max_inbounds = self.eic[1, (rt_bounds[0] <= self.eic[0]) & (self.eic[0] <= rt_bounds[1])].max()
            except ValueError:
                self.intensity_max_inbounds = np.nan

        except (KeyError, AttributeError, TypeError):
            self.eic = np.array([[],[]])
            self.intensity_max_inbounds = np.nan

class ChromPlotPlus:

    def __init__(self, data, shortname,
                 group, file_name, rt_buffer = .5,
                 x_scale = .8, y_scale = .75,
                 x_ratio = 13.0, y_ratio=11.0,
                 num_x_hashes=4, num_y_hashes=4,
                 **kwargs):

        assert len(data) > 0

        self.rt_peak = data[0]['identification'].rt_references[0].rt_peak
        self.rt_bounds = np.array([data[0]['identification'].rt_references[0].rt_min,
                                   data[0]['identification'].rt_references[0].rt_max])

        self.rt_min = min(self.rt_bounds[0], self.rt_peak) - rt_buffer
        self.rt_max = max(self.rt_bounds[1], self.rt_peak) + rt_buffer

        self.compound_eics = [CompoundFileEIC(compound_file_data,
                                              self.rt_bounds,
                                              self.rt_min,
                                              self.rt_max, 
                                              shortname) for compound_file_data in data]
        
        self.compound_eics = sorted(self.compound_eics,
                                    key = lambda c: (c.group_name,
                                                     c.file_name))

        try:
            self.intensity_max = np.nanmax([c.intensity_max_inbounds for c in self.compound_eics])
            if np.isnan(self.intensity_max):
                self.intensity_max = np.concatenate([c.eic[1] for c in self.compound_eics]).max()
            self.intensity_scale = np.floor(np.log10(1.1*self.intensity_max))
            for c in self.compound_eics:
                c.eic[1] = c.eic[1].clip(0, 1.1*self.intensity_max)
        except ValueError:
            self.intensity_max = 0
            self.intensity_scale = 0

        self.group = group
        self.file_name = file_name
        self.x_scale = x_scale
        self.y_scale = y_scale
        self.x_ratio = x_ratio
        self.y_ratio = y_ratio
        self.num_x_hashes = num_x_hashes
        self.num_y_hashes = num_y_hashes

        plt.ioff()

        self.fig, self.ax = plt.subplots()
        plt.setp(self.ax, 'frame_on', False)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.grid('off')
        self.ax.margins(1)

        self.__make_figure()

        plt.ion()

    def __make_figure(self):

        if self.group:
            num_group_labels = len(set([c.group_name for c in self.compound_eics]))

            self.compound_eics = ['^'] + self.compound_eics
            idxs = [i for
                    i,c in enumerate(zip(self.compound_eics,
                                         self.compound_eics[1:]))
                    if c[0] == '^'
                    or c[0].group_name != c[1].group_name]

            label = self.__yield_label()
            labels = [next(label) for _ in range(len(idxs))]

            for i,l in zip(idxs[::-1], labels[::-1]):
                self.compound_eics.insert(i+1, l)
            self.compound_eics = self.compound_eics[1:]

        num_cols = int(np.ceil((self.x_scale*(len(self.compound_eics))/self.y_scale)**.5))
        num_rows = int(np.ceil(float(len(self.compound_eics))/num_cols))

        self.ax.set_xlim((-.1,num_cols+.1))
        self.ax.set_ylim((-.1,num_rows+.1))

        is_subplot = np.array([(np.nan if isinstance(s,str) else 1)
                               for s in self.compound_eics])

        #Affine Transformations
        orgin_transform = np.array([[1, 0, -.5],
                                    [0, 1, -.5],
                                    [0, 0, 1]])
        scale_transform = np.array([[self.x_scale, 0, 0],
                                    [0, self.y_scale, 0],
                                    [0, 0, 1]])
        translate_transform = np.array([[[1, 0, (i%num_cols)+.5],
                                         [0, 1, num_rows-np.ceil(i/num_cols)-.5],
                                         [0, 0, 1]]
                                         for i in range(len(self.compound_eics))])
        tso_transform = np.matmul(translate_transform, np.matmul(scale_transform, orgin_transform))

        #Boxes
        boxes = np.matmul(tso_transform,
                          self.__make_boxes(is_subplot)
                         )[:,0:2].transpose(0,2,1)
        self.ax.add_collection(mc.LineCollection(boxes, colors=(0,0,0,1), linewidths=1.0/(num_cols*num_rows)**.5))

        #EICs
        eic_lines = [np.matmul(tso_transform[i], e)[0:2].T
                     for i,e in enumerate(self.__make_eics())]
        self.ax.add_collection(mc.LineCollection(eic_lines, 2.0/(num_cols*num_rows)**.5))

        #RT bounds
        rt_bound_lines = np.matmul(tso_transform,
                                   self.__make_rt_bounds(is_subplot)
                                  )[:,:,0:2].reshape(-1,2,2,order='F').transpose(0,2,1)
        # rt_bound_lines = rt_bound_lines[~np.isnan(rt_bound_lines[:,0,0])]
        # print np.isnan(rt_bound_lines[:,0,0])
        self.ax.add_collection(mc.LineCollection(rt_bound_lines, colors=(0,0,0,1), linewidths=2.0/(num_cols*num_rows)**.5))

        #Fills
        between = [(rt_bound_lines[2*i,0,0] <= e[:,0]) & (e[:,0] <= rt_bound_lines[(2*i)+1,0,0])
                   for i,e in enumerate(eic_lines)]
        for i,(e,b) in enumerate(zip(eic_lines, between)):
            self.ax.fill_between(e[:,0], e[:,1], np.full_like(e[:,0], np.matmul(tso_transform[i], np.array([0,0,1]))[1]),
            where=b, interpolate=False, facecolor='c', edgecolor='c', linewidth=0, alpha=0.3)
        # self.ax.add_collection(mc.PolyCollection(self.__make_fill(eic_lines, rt_bound_lines, is_subplot)))

        #RT peak
        self.rt_peak_lines = np.matmul(tso_transform,
                                  self.__make_rt_peaks(is_subplot)
                                 )[:,0:2].transpose(0,2,1)
        self.ax.add_collection(mc.LineCollection(self.rt_peak_lines, colors=(1,0,0,1), linewidths=2.0/(num_cols*num_rows)**.5))

        #Hashes
        hashes = np.matmul(tso_transform,
                           self.__make_hashes(is_subplot)
                          )[:,:,0:2].transpose(3,1,0,2).reshape(-1,2,2)
        self.ax.add_collection(mc.LineCollection(hashes, colors=(0,0,0,1), linewidths=1.0/(num_cols*num_rows)**.5))

        # X-axis labels
        x_labels = self.__make_x_labels()
        x_coordinates = np.matmul(tso_transform,
                                  self.__make_x_hash_coordinates(is_subplot)
                                 )[:,0:2].transpose(0,2,1)
        for subplot_xy in x_coordinates:
            for l,xy in zip(x_labels, subplot_xy):
                self.ax.annotate(l, xy, ha='center', va='top', size = 8./(num_cols+.25))

        # Y-axis labels
        y_labels = self.__make_y_labels()
        y_coordinates = np.matmul(tso_transform,
                                  self.__make_y_hash_coordinates(is_subplot)
                                 )[:,0:2].transpose(0,2,1)
        for subplot_xy in y_coordinates:
            for l,xy in zip(y_labels, subplot_xy):
                self.ax.annotate(l, xy, ha='right', va='center', size = 8./(num_rows+.25))

        # Y-scale label
        self.y_scale_coordinates = np.matmul(tso_transform,
                                  self.__make_y_scale_coordinates(is_subplot)
                                 )[:,0:2].transpose(0,2,1)
        for subplot_xy in self.y_scale_coordinates:
            for xy in subplot_xy:
                self.ax.annotate('1E%d'%int(self.intensity_scale), xy, ha='right', va='bottom', size = 8./(num_rows+.25), weight='bold')

        # Titles
        title_coordinates = np.matmul(tso_transform,
                                  self.__make_title_coordinates(is_subplot)
                                 )[:,0:2].transpose(0,2,1)

        for i,subplot_xy in enumerate(title_coordinates):
            if not isinstance(self.compound_eics[i], str):
                self.ax.annotate('\n'.join(wrap(self.compound_eics[i].file_name,48)),
                                 subplot_xy[0], ha='center', va='bottom', size = 8./(num_cols+.25), weight='bold')

        # Groups
        if self.group:
            group_label_coordinates = boxes = np.matmul(tso_transform,
                                                        self.__make_group_label_coordinates(is_subplot)
                                                        )[:,0:2].transpose(0,2,1)
            for i,subplot_xy in enumerate(group_label_coordinates):
                if isinstance(self.compound_eics[i], str):
                    self.ax.annotate(self.compound_eics[i],
                                     subplot_xy[0], ha='center', va='center', size = 100./(num_cols+.25), weight='bold')

        with PdfPages(self.file_name) as pdf:
            plt.rcParams['pdf.fonttype'] = 42
            plt.rcParams['pdf.use14corefonts'] = True
            plt.rcParams['text.usetex'] = False
            pdf.savefig(self.fig)
            plt.close()


    @staticmethod
    def __yield_label():
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


    def __make_boxes(self, is_subplot):
        return np.array([[[0,0,1,1,0],
                          [0,1,1,0,0],
                          [1,1,1,1,1]]]*len(self.compound_eics)
                       )*is_subplot[:,None,None]

    def __make_eics(self):
        return [np.asarray([(s.eic[0]-self.rt_min)/(self.rt_max-self.rt_min),
                             s.eic[1]/(1.1*self.intensity_max),
                             [1]*s.eic.shape[1]]) if not isinstance(s,str)
                            else [[np.nan], [np.nan], [np.nan]]
                           for s in self.compound_eics]

    def __make_rt_bounds(self, is_subplot):
        return (np.array([[[(self.rt_bounds[0]-self.rt_min)/(self.rt_max-self.rt_min)]*2,
                           [0,1],
                           [1,1]],
                          [[(self.rt_bounds[1]-self.rt_min)/(self.rt_max-self.rt_min)]*2,
                           [0,1],
                           [1,1]]]*len(self.compound_eics)).reshape(len(self.compound_eics),2,3,2).transpose(1,0,2,3)
                )*is_subplot[None,:,None,None]

    # def __make_fill(self, eic_lines, rt_bound_lines, is_subplot):
    #     rt_bound_intensities = [interp1d(eic_lines[i][:,0],
    #                                      eic_lines[i][:,1],
    #                                      bounds_error=False, copy=False, assume_sorted=True
    #                                      )(rt_bound_lines[2*i:2*(i+1),0,0])
    #                             if ~np.isnan(is_subplot[i]) and eic_lines[i].shape[0] > 1
    #                             else np.array([np.nan, np.nan])
    #                             for i in range(len(eic_lines))]
    #
    #     return [[[rt_bound_lines[2*i:2*(i+1),0,0].min(),
    #               rt_bound_lines[2*i,0,1]]]+
    #             [[rt_bound_lines[2*i:2*(i+1),0,0].min(),
    #               rt_bound_intensities[i][rt_bound_lines[2*i:2*(i+1),0,0].argmin()]]]+
    #             eic_lines[i][(eic_lines[i][:,0] > rt_bound_lines[2*i:2*(i+1),0,0].min()) &
    #                     (eic_lines[i][:,0] < rt_bound_lines[2*i:2*(i+1),0,0].max())
    #                    ].tolist()+
    #             [[rt_bound_lines[2*i:2*(i+1),0,0].max(),
    #               rt_bound_intensities[i][rt_bound_lines[2*i:2*(i+1),0,0].argmax()]]]+
    #             [[rt_bound_lines[2*i:2*(i+1),0,0].max(),
    #               rt_bound_lines[2*i,0,1]]]+
    #             [[rt_bound_lines[2*i:2*(i+1),0,0].min(),
    #               rt_bound_lines[2*i,0,1]]]
    #             for i in range(len(eic_lines))
    #             if ~(np.isnan(rt_bound_intensities[i]).all())]

    def __make_rt_peaks(self, is_subplot):
        return np.array([[[(self.rt_peak-self.rt_min)/(self.rt_max-self.rt_min)]*2,
                          [0,1],
                          [1,1]]]*len(self.compound_eics))*is_subplot[:,None,None]

    def __make_hashes(self, is_subplot):
        return (np.array([[[[[float(i)/(self.num_x_hashes+1),float(i)/(self.num_x_hashes+1)],
                             [0,.02],
                             [1,1]]]+
                           [[[float(i)/(self.num_x_hashes+1),float(i)/(self.num_x_hashes+1)],
                             [1,.98],
                             [1,1]]] for i in range(1,self.num_x_hashes+1)]+
                          [[[[0,.02],
                             [float(j)/(self.num_y_hashes+1),float(j)/(self.num_y_hashes+1)],
                             [1,1]]]+
                           [[[1,.98],
                             [float(j)/(self.num_y_hashes+1),float(j)/(self.num_y_hashes+1)],
                             [1,1]]] for j in range(1,self.num_y_hashes+1)]]*len(self.compound_eics)
                        ).reshape(len(self.compound_eics),2*(self.num_x_hashes+self.num_y_hashes),3,2
                                 )*is_subplot[:,None,None,None]).transpose(3,0,2,1)

    def __make_x_labels(self):
        return np.linspace(self.rt_min,
                           self.rt_max,
                           self.num_x_hashes+2).round(3)

    def __make_x_hash_coordinates(self, is_subplot):
        return (np.array([[np.linspace(0,1,self.num_x_hashes+2),
                          [-.025]*(self.num_x_hashes+2),
                          [1]*(self.num_x_hashes+2)]]*len(self.compound_eics)))*is_subplot[:,None,None]

    def __make_y_labels(self):
        return np.linspace(0,
                           (1.1*self.intensity_max)/(10**self.intensity_scale),
                           self.num_y_hashes+2).round(2)

    def __make_y_hash_coordinates(self, is_subplot):
        return (np.array([[[0]*(self.num_y_hashes+2),
                           np.linspace(0,1,self.num_y_hashes+2),
                           [1]*(self.num_y_hashes+2)]]*len(self.compound_eics)))*is_subplot[:,None,None]

    def __make_y_scale_coordinates(self, is_subplot):
        return (np.array([[[-.075],
                           [1.05],
                           [1]]]*len(self.compound_eics)))*is_subplot[:,None,None]

    def __make_title_coordinates(self, is_subplot):
        return (np.array([[[.5],
                           [1],
                           [1]]]*len(self.compound_eics)))*is_subplot[:,None,None]

    def __make_group_label_coordinates(self, is_subplot):
        is_not_subplot = np.isnan(is_subplot)
        is_not_subplot[~is_not_subplot] = np.nan
        return (np.array([[[.5],
                           [.5],
                           [1]]]*len(self.compound_eics)))*is_not_subplot[:,None,None]
