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
import sys
import math

class Crit:
    """
    
    Class to organize and run criticality analyses based on an array of binarized spiketimes

    Necessary parameters are passed through the init statement with default values being used
    if not given initially. 

    eg:
    crit = Crit(spikewords, perc = .5)

    the Crit object can then have analyses run on it by executing 
    crit.run_crit()

    following this the Crit object now holds all the values from the analyses 

    eg:
    crit.dcc
    >>>> 0.23

    
    
    """
    def __init__(self, spikewords, perc = 0.35, nfactor_bm = 0, nfactor_tm = 0, nfactor_bm_tail = 1, nfactor_tm_tail = 1,
                bm = None, tm = None, saveloc = '', pltname = '', plot = True, none_fact = 40, 
                exclude = False, exclude_burst=50, exclude_time=20, exclude_diff_b = 20, exclude_diff_t = 10, subsample = False,
                subsample_factor = None, subsample_iter = None, subsample_replace = False):
        # required parameters
        self.perc = perc
        self.spikewords = spikewords
        n, b = np.shape(spikewords)
        self.num_cells = n
        self.num_bins = b
        self.saveloc = saveloc
        self.pltname = pltname
        self.plot = plot
        self.nfactor_bm = nfactor_bm
        self.nfactor_tm = nfactor_tm
        self.nfactor_bm_tail = nfactor_bm_tail
        self.nfactor_tm_tail = nfactor_tm_tail
        self.none_fact = none_fact
        self.exclude = exclude
        self.exclude_burst = exclude_burst
        self.exclude_time = exclude_time
        self.exclude_diff_b = exclude_diff_b
        self.exclude_diff_t = exclude_diff_t
        self.bm = bm
        self.tm = tm
        self.subsample = subsample
        self.subsample_factor = subsample_factor
        self.subsample_iter = subsample_iter
        self.subsample_replace = subsample_replace

        # all parameters set by run_crit
        self.burst = None
        self.T = None
        self.p_value_burst = None
        self.p_value_t = None
        self.dcc = None
        self.scaling_plot = None
        self.burst_cdf_plot = None
        self.t_cdf_plot = None
        self.sigma = None
        self.fit_sigma = None
        self.xmin = None
        self.xmax = None
        self.tmin = None
        self.tmax = None
        self.alpha = None
        self.beta = None
        self.EXCLUDED_b = None
        self.EXCLUDED_t = None

    def __repr__(self):
        """
        mild descriptors of the object - wont work
        """
        return str(f'CRIT: {self.num_cells} cells - {self.num_bins} bins')

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        print('closing file')


    def get_params(self):
        print(f'PARAMETERS:\n'
              f'perc: {self.perc}\n'
              f'nfactor_bm: {self.nfactor_bm}\n'
              f'nfactor_tm: {self.nfactor_tm}\n'
              f'nfactor_bm_tail: {self.nfactor_bm_tail}\n'
              f'nfactor_tm_tail: {self.nfactor_tm_tail}\n')

    def finalize(self):
        self.cells = None
        self.final = True

    def show_plots(self):
        self.t_cdf_plot.show()
        self.burst_cdf_plot.show()
        self.scaling_plot.show()

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
        ax.set_title(f"Raster Avalanches")
        sns.despine()

        if show:
            plt.show()

        if saveplot:
            fig.savefig(os.path.join(self.saveloc, 'spike_raster.png'))

    def run_crit(self, flag = 2, verbose = True):
        """
        So this will do the whole calculation. I think i'm gonna break it up into methods even further for re-running purposes

        need to:
        - run get avalanches
        - run av_analyses
        """
        burst = []
        T = []
        if self.subsample:
            try:
                idxs = np.random.choice(np.arange(self.num_cells), size=[self.subsample_iter, self.subsample_factor], replace = self.subsample_replace)
            except ValueError as err:
                print(err)
                print('Please choose different subsample parameters')
                sys.exit()
            self.subsample_idxs = idxs
            for i in idxs:
                temp_dat = self.spikewords[i]
                temp_R = cr.get_avalanches(temp_dat, perc = self.perc)
                burst.append(temp_R['S'])
                T.append(temp_R['T'])
            burst = np.concatenate(burst)
            T = np.concatenate(T)
        else:
            R = cr.get_avalanches(self.spikewords, perc = self.perc)
            burst = R['S']
            T = R['T']

        self.burst = burst
        self.T = T

        Result = cr.AV_analysis(self.burst, self.T, flag = flag, pltname = self.pltname, saveloc = self.saveloc, plot = self.plot, bm = self.bm, 
                                tm = self.tm, nfactor_bm = self.nfactor_bm, nfactor_tm = self.nfactor_tm, nfactor_bm_tail = self.nfactor_bm_tail,
                                nfactor_tm_tail = self.nfactor_tm_tail, none_fact = self.none_fact, verbose = verbose, exclude = self.exclude, 
                                exclude_burst = self.exclude_burst, exclude_time = self.exclude_time, exclude_diff_b = self.exclude_diff_b, exclude_diff_t=self.exclude_diff_t)

        if flag == 2:
            self.p_value_burst = Result['P_burst']
            self.p_value_t = Result['P_t']

        self.dcc = Result['df']

        if self.plot:
            self.scaling_plot = Result['scaling_relation_plot']
            if flag == 2:
                self.burst_cdf_plot = Result['burst_cdf']
                self.t_cdf_plot = Result['time_cdf']

        self.sigma = Result['pre']
        self.fit_sigma = Result['fit']
        self.xmin = Result['xmin']
        self.xmax = Result['xmax']
        self.tmin = Result['tmin']
        self.tmax = Result['tmax']
        self.alpha = Result['alpha']
        self.beta = Result['beta']
        self.EXCLUDED_b = Result['EX_b']
        self.EXCLUDED_t = Result['EX_t']
        self.__gen_kappa()
        self.__gen_k2()
        self.__gen_kprob()
        self.__gen_k3()

    def __gen_beta(self):
        tm = int(np.max(self.T) / 20)
        _, _, beta = \
            cr.EXCLUDE(self.T[self.T < np.power(np.max(self.T), self.nfactor_tm_tail)], tm,
                       nfactor = self.nfactor_tm)
        self.beta = beta

    def __gen_k2(self, num = 100):
        if self.burst is None:
            print("You must run_crit() before you can run this function")
            return
        if self.xmin is None:
            print('this block failed while running for some reason')
            return
        idx = np.where(np.logical_and(self.burst <= self.xmax, self.burst >= self.xmin))
        b = self.burst[idx]
        n = np.size(b)
        cdf = np.cumsum(np.histogram(b, np.arange(self.xmin, self.xmax + 2))[0] / n)
        s = np.unique(b)
        A = 1 / np.sum(np.power(s, -self.alpha))
        fit = np.cumsum(A * np.power(np.arange(self.xmin, self.xmax + 1), -self.alpha))

        idxs = np.geomspace(1, np.size(cdf) - 1, num = num, dtype = int)
        half = int(num / 2)
        second_half = idxs[half:]
        diffs = [cdf[i] - fit[i] for i in second_half]
        mean_diff = np.mean(diffs)
        k2b = mean_diff

        self.k2b = k2b
        try:
            idx2 = np.where(np.logical_and(self.T <= self.tmax, self.T >= self.tmin))
            t = self.T[idx2]
            n = np.size(t)
            cdf = np.cumsum(np.histogram(t, np.arange(self.tmin, self.tmax + 2))[0] / n)
            s = np.unique(t)
            A = 1 / np.sum(np.power(s, -self.beta))
            fit = np.cumsum(A * np.power(np.arange(self.tmin, self.tmax + 1), -self.beta))

            idxs = np.geomspace(1, np.size(cdf) - 1, num = num, dtype = int)
            half = int(num / 2)
            second_half = idxs[half:]
            diffs = [cdf[i] - fit[i] for i in second_half]
            mean_diff = np.mean(diffs)
            k2t = mean_diff

            self.k2t = k2t
        except AttributeError:
            self.k2t = None
        return self.k2b, self.k2t

    def __gen_kappa(self, num = 100):
        if self.burst is None:
            print("You must run_crit() before you can run this function")
            return
        if self.xmin is None:
            print('block failed')
            return
        idx = np.where(np.logical_and(self.burst <= self.xmax, self.burst >= self.xmin))
        b = self.burst[idx]
        n = np.size(b)
        cdf = np.cumsum(np.histogram(b, np.arange(self.xmin, self.xmax + 2))[0] / n)
        s = np.unique(b)
        A = 1 / np.sum(np.power(s, -self.alpha))
        fit = np.cumsum(A * np.power(np.arange(self.xmin, self.xmax + 1), -self.alpha))

        idxs = np.geomspace(1, np.size(cdf) - 1, num = num, dtype = int)
        diffs = [fit[i] - cdf[i] for i in idxs]
        mean_diff = np.mean(diffs)
        kappa_burst = 1 + mean_diff

        self.kappa_burst = kappa_burst
        try:
            idx2 = np.where(np.logical_and(self.T <= self.tmax, self.T >= self.tmin))
            t = self.T[idx2]
            n = np.size(t)
            cdf = np.cumsum(np.histogram(t, np.arange(self.tmin, self.tmax + 2))[0] / n)
            s = np.unique(t)
            A = 1 / np.sum(np.power(s, -self.beta))
            fit = np.cumsum(A * np.power(np.arange(self.tmin, self.tmax + 1), -self.beta))

            idxs = np.geomspace(1, np.size(cdf) - 1, num = num, dtype = int)
            diffs = [fit[i] - cdf[i] for i in idxs]
            mean_diff = np.mean(diffs)
            kappa_t = 1 + mean_diff

            self.kappa_t = kappa_t
        except AttributeError:
            self.kappa_t = None
        return self.kappa_burst, self.kappa_t

    def __gen_kprob(self):
        if self.burst is None:
            print("You must run_crit() before you can run this function")
            return
        if self.xmin is None:
            print('this block failed while running for some reason')
            return

        pdf = np.histogram(self.burst, bins = np.arange(1, np.max(self.burst) + 2))[0]
        p = pdf / np.sum(pdf)
        x = np.arange(self.xmin, self.xmax + 1)
        y = (np.size(np.where(self.burst == self.xmin + 6)[0]) / np.power(self.xmin + 6, -self.alpha)) * \
            np.power(x, -self.alpha)
        y = y / np.sum(pdf)

        ps = p[self.xmin:self.xmax]
        diffs = np.diff([ps, y[:-1]], axis = 0)
        above_idx = np.where(diffs < 0)[1]
        below_idx = np.where(diffs > 0)[1]
        prob_above = np.sum(ps[above_idx])
        prob_below = np.sum(ps[below_idx])

        kprob_b = prob_above / prob_below
        self.kprob_b = kprob_b
        try:
            tdf = np.histogram(self.T, bins = np.arange(1, np.max(self.T) + 2))[0]
            t = tdf / np.sum(tdf)
            x = np.arange(self.tmin, self.tmax + 1)
            y = (np.size(np.where(self.T == self.tmin + 6)[0]) / np.power(self.tmin + 6, -self.beta)) * \
                np.power(x, -self.beta)
            y = y / np.sum(tdf)

            ps = t[self.tmin:self.tmax]
            diffs = np.diff([ps, y[:-1]], axis = 0)
            above_idx = np.where(diffs < 0)[1]
            below_idx = np.where(diffs > 0)[1]
            prob_above = np.sum(ps[above_idx])
            prob_below = np.sum(ps[below_idx])

            kprob_t = prob_above / prob_below
            self.kprob_t = kprob_t
        except Exception:
            print('kprob_t not working')
            self.kprob_t = None

    def __gen_k3(self):
        pdf = np.histogram(self.burst, bins = np.arange(1, np.max(self.burst) + 2))[0]
        p = pdf / np.sum(pdf)
        x = np.arange(self.xmin, np.max(self.burst) + 1)
        y = (np.size(np.where(self.burst == self.xmin + 6)[0]) / np.power(self.xmin + 6, -self.alpha)) * \
            np.power(x, -self.alpha)
        y = y / np.sum(pdf)


        half = int(math.sqrt(self.xmin*np.max(self.burst)))
        xh = np.arange(x[half], np.max(self.burst) + 1)
        yh = (np.size(np.where(self.burst == self.xmin + 6)[0]) / np.power(self.xmin + 6, -self.alpha)) * \
            np.power(xh, -self.alpha)
        yh = yh / np.sum(pdf)


        ps = p[x[half]:]
        diffs = np.diff([ps, yh[:-1]], axis = 0)
        above_idx = np.where(diffs < 0)[1]
        below_idx = np.where(diffs > 0)[1]
        prob_above = np.sum(ps[above_idx])
        prob_below = np.sum(ps[below_idx])

        kprob_b = prob_above / prob_below

        self.k3b = kprob_b

        tdf = np.histogram(self.T, bins = np.arange(1, np.max(self.T) + 2))[0]
        p = tdf / np.sum(tdf)
        x = np.arange(self.tmin, np.max(self.T) + 1)
        y = (np.size(np.where(self.T == self.tmin + 6)[0]) / np.power(self.tmin + 6, -self.beta)) * \
            np.power(x, -self.beta)
        y = y / np.sum(tdf)


        half = int(math.sqrt(self.tmin*np.max(self.T)))
        xh = np.arange(x[half], np.max(self.T) + 1)
        yh = (np.size(np.where(self.T == self.tmin + 6)[0]) / np.power(self.tmin + 6, -self.beta)) * \
            np.power(xh, -self.beta)
        yh = yh / np.sum(tdf)


        ps = p[x[half]:]
        diffs = np.diff([ps, yh[:-1]], axis = 0)
        above_idx = np.where(diffs < 0)[1]
        below_idx = np.where(diffs > 0)[1]
        prob_above = np.sum(ps[above_idx])
        prob_below = np.sum(ps[below_idx])

        kprob_t = prob_above / prob_below

        self.k3t = kprob_t




