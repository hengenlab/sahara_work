from sahara_work.crit_class import Crit
import numpy as np
import criticality as cr
import musclebeachtools as mbt
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import seaborn as sns
import re
import pandas as pd
import os
import glob
import signal
import gc
from datetime import datetime as dt
from datetime import timedelta
import csv
from copy import deepcopy as cdc


class Crit_hlab(Crit):

    def __init__(self, **kwargs):
        super().__init__(
            spikewords = kwargs.get('spikewords'),
            perc = kwargs.get('perc'),
            nfactor_bm = kwargs.get('nfactor_bm'),
            nfactor_tm = kwargs.get('nfactor_tm'),
            nfactor_bm_tail = kwargs.get('nfactor_bm_tail'),
            nfactor_tm_tail = kwargs.get('nfactor_tm_tail'),
            saveloc = kwargs.get('saveloc'),
            pltname = kwargs.get('pltname'),
            plot = kwargs.get('plot'),
            none_fact = kwargs.get('none_fact'),
            exclude = kwargs.get('exclude'),
            exclude_burst = kwargs.get('exclude_burst'),
            exclude_time = kwargs.get('exclude_time'),
            exclude_diff_b = kwargs.get('exclude_diff_b'),
            exclude_diff_t = kwargs.get('exclude_diff_t'),
            bm = kwargs.get('bm'),
            tm = kwargs.get('tm'),
            subsample = kwargs.get('subsample'),
            subsample_factor = kwargs.get('subsample_factor'), 
            subsample_iter = kwargs.get('subsample_iter'),
            subsample_replace = kwargs.get('subsample_replace')
        )
        # all parameters set by run_crit
        self.burst = kwargs.get('burst')
        self.T = kwargs.get('T')
        self.p_value_burst = kwargs.get('p_value_burst')
        self.p_value_t = kwargs.get('p_value_t')
        self.dcc = kwargs.get('dcc')
        self.scaling_plot = kwargs.get('scaling_plot')
        self.burst_cdf_plot = kwargs.get('burst_cdf_plot')
        self.t_cdf_plot = kwargs.get('t_cdf_plot')
        self.pre = kwargs.get('pre')
        self.fit = kwargs.get('fit')
        self.xmin = kwargs.get('xmin')
        self.xmax = kwargs.get('xmax')
        self.tmin = kwargs.get('tmin')
        self.tmax = kwargs.get('tmax')
        self.alpha = kwargs.get('alpha')
        self.EXCLUDED_b = kwargs.get('EXCLUDED_b')
        self.EXCLUDED_t = kwargs.get('EXCLUDED_t')

        # all parameters set by run_branching_ratio
        self.br1 = kwargs.get('br1')
        self.br2 = kwargs.get('br2')
        self.acc1 = kwargs.get('acc1')
        self.acc2 = kwargs.get('acc2')
        self.kmax = kwargs.get('kmax')
        self.br_binsize = kwargs.get('br_binsize')

        # all parameters that can be set by lilo and stitch
        self.time_frame = kwargs.get('time_frame')
        self.block_num = kwargs.get('block_num')
        self.qualities = kwargs.get('qualities')
        self.cell_types = kwargs.get('cell_types')
        self.hour_bins = kwargs.get('hour_bins')
        self.ava_binsize = kwargs.get('ava_binsize')
        self.animal = kwargs.get('animal')
        self.date = kwargs.get('date')
        self.final = kwargs.get('final')
        self.all_cells = kwargs.get('all_cells')
        self.used_cells = kwargs.get('used_cells')
        self.probe = kwargs.get('probe')
        self.scored_by = kwargs.get('scored_by')
        self.pathname = kwargs.get('pathname')
        self.filename = kwargs.get('filename')


    def __repr__(self):
        return str(f'CRIT_HLAB: {self.animal} -- {self.date} -- {self.time_frame}')

    def get_params(self):
        print(f'PARAMETERS:\n'
              f'binsize: {self.ava_binsize}\n'
              f'perc: {self.perc}\n'
              f'hour bins: {self.hour_bins}\n'
              f'qualities: {self.qualities}\n'
              f'cell types: {self.cell_types}\n'
              f'nfactor_bm: {self.nfactor_bm}\n'
              f'nfactor_tm: {self.nfactor_tm}\n'
              f'nfactor_bm_tail: {self.nfactor_bm_tail}\n'
              f'nfactor_tm_tail: {self.nfactor_tm_tail}\n')

    def run_crit_from_start(self, flag = 2, save = False):

        """
        This isn't going to work as a class function, if you need it then fix it
        """
        if self.final:
            print('This crit object is final, there are no cells saved here. If youd like to rerun this block start from lilo_and_stitch')
            return
        total_time = s.__get_totaltime(self.time_frame)
        num_bins = int(total_time / self.hour_bins)
        bin_len = int((self.hour_bins * 3600) / self.ava_binsize)
        good_cells = [cell for cell in self.cells if self.quality in self.qualities and self.cell_type in self.cell_types]
        spikewords = mbt.n_spiketimes_to_spikewords(good_cells, binsz = self.ava_binsize, binarize = 1)
        idx = self.block_num
        if idx == num_bins - 1:
            data = spikewords[:, (idx * bin_len):]
        else:
            data = spikewords[:, (idx * bin_len): ((idx + 1) * bin_len)]
        self.spikewords = data
        param_str = s.__get_paramstr(self.animal, self.probe, self.date, self.time_frame, self.hour_bins, self.perc, self.ava_binsize, self.qualities, self.cell_types, idx)
        self.pltname = param_str
        self.run_crit(flag = flag)
        print(f'BLOCK RESULTS: P_vals - {self.p_value_burst}   {self.p_value_t} \n DCC: {self.dcc}')
        if save:
            to_save = np.array([obj])
            np.save(f'{obj.saveloc}Crit_{param_str}', to_save)

    def run_branching_ratio(self, cells, binsize=0.004, start = 0, end = 1, kmax=500, pltname = None):
        '''
        cells: mbt neuron list
        binsize: ava_binsize length for spikewords
        start: where to start BR (hours)
        end: where to end BR (hours)
        kmax: number of bins to wait between activity essentially
        pltname: what to save the plot as (if None no plot is made)
        '''
        br1, br2, acc1, acc2 = mbt.n_branching_ratio(cells, ava_binsz=binsize, kmax=kmax, start=start, end=end, binarize=1, plotname=pltname)
        self.br1 = br1
        self.br2 = br2
        self.acc1 = acc1 
        self.acc2 = acc2 
        self.br_binsize = binsize
        self.kmax = kmax



    def plot_raster(self, window_start = 200, window_end = 400, saveplot = False, show = True):
        # window of time to look at
        small_spikewords = self.spikewords[:, window_start : window_end]
        n, m = np.shape(small_spikewords)

        # pulling out spikes from the binned spikewords array
        positions = []
        for cell in small_spikewords:
            positions.append(np.where(cell > 0)[0])

        # pulling out threshold to determine avalanches
        network = np.nansum(small_spikewords, axis = 0)
        sortN = np.sort(network)
        threshold = sortN[round(m * self.perc)]
        thresh_max = np.ma.masked_where(network <= threshold, network)

        # determining avalanches
        zdata = cdc(network)
        zdata[~thresh_max.mask] = 1  # avalanches
        zdata[thresh_max.mask] = 0  # intervals

        # pulling out the edges for the plot
        edges = np.diff(zdata)
        ontimes = np.where(edges > 0)[0]
        offtimes = np.where(edges < 0)[0]

        #edge cases

        if zdata[0] == 1:
            ontimes = np.insert(ontimes, 0, 0)
        if zdata[-1] == 1:
            offtimes = np.append(offtimes, len(edges))

        #making the rectagles
        xys = []
        widths = []
        heights = []
        for i, on in enumerate(ontimes):
            xys.append((on, 0))
            widths.append(offtimes[i] - on)

        boxes = [Rectangle(xy, width, n) for xy, width in zip(xys, widths)]
        pc = PatchCollection(boxes, facecolor = 'moccasin', alpha = 0.7, edgecolor = 'firebrick')

        # plotting
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (11, 8))

        ax.eventplot(positions, orientation = "horizontal")
        ax.add_collection(pc)
        ax.set_ylabel("Cell")
        ax.set_xlabel("Bin")
        ax.set_title(f"{self.animal} -- Raster Avalanches")
        sns.despine()

        if show:
            plt.show()

        if saveplot:
            fig.savefig(os.path.join(self.saveloc, 'spike_raster.png'))