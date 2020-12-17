import numpy as np
from criticality_hlab import criticality as cr
from musclebeachtools_hlab import musclebeachtools as mbt
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

class Crit:
    """
    
    Class to look at criticality stuff

    tbd on init and what not
    
    """
    def __init__(self, spikewords, perc = 0.35, nfactor_bm = 0, nfactor_tm = 0, nfactor_bm_tail = 1, nfactor_tm_tail = 1, saveloc = '', pltname = '', plot = True):

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

        # all parameters set by run_crit
        self.burst = None
        self.T = None
        self.bm = None
        self.tm = None
        self.p_value_burst = None
        self.p_value_t = None
        self.dcc = None
        self.scaling_plot = None
        self.burst_cdf_plot = None
        self.t_cdf_plot = None
        self.pre = None
        self.fit = None
        self.xmin = None
        self.xmax = None
        self.tmin = None
        self.tmax = None
        self.alpha = None


    def __repr__(self):
        """
        mild descriptors of the object
        """
        return str(f'CRIT: {self.animal} -- {self.date} -- {self.time_frame}')

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

    def run_crit(self, flag = 2):
        """
        So this will do the whole calculation. I think i'm gonna break it up into methods even further for re-running purposes

        need to:
        - run get avalanches
        - run av_analyses
        """

        R = cr.get_avalanches(self.spikewords, perc = self.perc)
        burst = R['S']
        T = R['T']
        self.burst = burst
        self.T = T

        self.bm = int(np.max(burst) / 20)
        self.tm = int(np.max(T) / 20)

        crit_params = {
            'flag': flag,
            'pltname': self.pltname,
            'saveloc': self.saveloc,
            'plot': self.plot,
            'bm': self.bm,
            'tm': self.tm
        }

        Result = cr.AV_analysis(burst, T, crit_params, nfactor_bm = self.nfactor_bm, nfactor_tm = self.nfactor_tm, nfactor_bm_tail = self.nfactor_bm_tail,
                                nfactor_tm_tail = self.nfactor_tm_tail)

        if flag == 2:
            self.p_value_burst = Result['P_burst']
            self.p_value_t = Result['P_t']

        self.dcc = Result['df']

        if self.plot:
            self.scaling_plot = Result['scaling_relation_plot']
            if flag == 2:
                self.burst_cdf_plot = Result['burst_cdf']
                self.t_cdf_plot = Result['time_cdf']

        self.pre = Result['pre']
        self.fit = Result['fit']
        self.xmin = Result['xmin']
        self.xmax = Result['xmax']
        self.tmin = Result['tmin']
        self.tmax = Result['tmax']
        self.alpha = Result['alpha']
        self.beta = Result['beta']
        self.__gen_kappa()
        self.__gen_k2()
        self.__gen_kprob()

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


def get_cell_stats(cell):
    fr, xbins = cell.plotFR(binsz = 3600, start = False, end = False,
                            lplot = 0)

    isis = []
    bins = np.arange(cell.start_time, cell.end_time, 300)  # 5 min bins

    for i in range(len(bins)):
        if i == len(bins) - 1:
            end_bin = cell.end_time
        else:
            end_bin = bins[i + 1]
        spk_idxs = np.where(np.logical_and(cell.spike_time_sec > bins[i], cell.spike_time_sec < end_bin))
        isis.append(np.diff(cell.spike_time[spk_idxs]))

    means = np.array([np.mean(i) for i in isis])
    stds = np.array([np.std(i) for i in isis])
    cvs = stds / means
    binned_cvs = [cvs[i * 12:(i + 1) * 12] for i in np.arange(int(xbins[-1]))]
    cv = [np.mean(i) for i in binned_cvs]
    if len(fr) != len(cv):
        print('this isnt working, fix it')
    return fr, cv


def construct_fr_df(paths):
    bdays = {
        'caf01': dt(2019, 12, 24, 7, 30),
        'caf19': dt(2020, 1, 19, 7, 30),
        'caf22': dt(2020, 2, 17, 7, 30),
        'caf26': dt(2020, 2, 20, 7, 30),
        'caf34': dt(2020, 3, 18, 7, 30),
        'caf37': dt(2019, 8, 18, 7, 30),
        'caf40': dt(2020, 2, 20, 7, 30),
        'caf42': dt(2020, 2, 20, 7, 30),
        'caf48': dt(2020, 7, 20, 7, 30),
        'caf49': dt(2020, 7, 20, 7, 30),
        'caf50': dt(2020, 7, 20, 7, 30),
        'caf52': dt(2020, 4, 19, 7, 30),
        'caf54': dt(2020, 7, 11, 7, 30),
        'caf55': dt(2020, 7, 11, 7, 30),
        'caf58': dt(2020, 9, 23, 7, 30),
        'caf60': dt(2020, 9, 23, 7, 30),
        'caf61': dt(2019, 12, 11, 7, 30),
        'caf62': dt(2019, 11, 18, 7, 30),
        'eab52': dt(2019, 4, 19, 7, 30),
        'eab47': dt(2019, 2, 17, 7, 30),
        'eab': dt(2019, 2, 17, 7, 30),
        'eab50': dt(2019, 2, 15, 7, 30),
        'eab40': dt(2018, 12, 5, 7, 30)
    }
    seconds_in_day = 60 * 60 * 24
    with open('/media/HlabShare/clayton_sahara_work/criticality/cell_stats.csv', 'w', newline = '') as c:
        w = csv.DictWriter(c, fieldnames = ['animal', 'rstart_time', 'age_start', 'days_old', 'hours_old', 'cell_idx', 'quality', 'fr', 'cv', 'wf'])
        w.writeheader()

    done = []
    final = []
    for i, p in enumerate(paths):
        print(i)
        path_info = p[:p.find('_perc')]
        if path_info in done:
            print('already done')
        else:
            done.append(path_info)
            crit = None
            try:
                crit = np.load(p, allow_pickle = True)[0]
                birth = bdays[crit.animal]
                for cell in crit.cells:
                    start_time = crit.cells[0].rstart_time
                    start_time = dt.strptime(start_time, '%Y-%m-%d_%H-%M-%S')
                    age = start_time - birth

                    fr, cv = get_cell_stats(cell)

                    for i in range(len(fr)):
                        age_now = age + timedelta(hours = i)
                        days_old = age_now.total_seconds() / seconds_in_day
                        hours_old = age_now.total_seconds() / 3600

                        with open('/media/HlabShare/clayton_sahara_work/criticality/cell_stats.csv', 'a', newline = '') as c:
                            w = csv.DictWriter(c, fieldnames = ['animal', 'rstart_time', 'age_start', 'days_old', 'hours_old', 'cell_idx', 'quality', 'fr', 'cv', 'wf'])

                            d = {
                                'animal': crit.animal,
                                'rstart_time': start_time,
                                'age_start': age_now,
                                'days_old': days_old,
                                'hours_old': hours_old,
                                'cell_idx': cell.clust_idx,
                                'quality': cell.quality,
                                'fr': fr[i],
                                'cv': cv[i],
                                'wf': cell.waveform
                            }
                            w.writerow(d)
            except Exception:
                print('This object is final or weird. Go back and find these cells individually to add to the csv')
                final.append(p)
    np.save('/media/HlabShare/clayton_sahara_work/fr_csv_done.npy', done)
    np.save('/media/HlabShare/clayton_sahara_work/fr_csv_FINAL.npy', final)


def run_crit_from_start(obj, flag = 2, save = True):
    if obj.final:
        print('This crit object is final, there are no cells saved here. If youd like to rerun this block start from lilo_and_stitch')
        return
    total_time = __get_totaltime(obj.time_frame)
    num_bins = int(total_time / obj.hour_bins)
    bin_len = int((obj.hour_bins * 3600) / obj.ava_binsize)
    good_cells = [cell for cell in obj.cells if cell.quality in obj.qualities and cell.cell_type in obj.cell_types]
    spikewords = mbt.n_spiketimes_to_spikewords(good_cells, binsz = obj.ava_binsize, binarize = 1)
    idx = obj.block_num
    if idx == num_bins - 1:
        data = spikewords[:, (idx * bin_len):]
    else:
        data = spikewords[:, (idx * bin_len): ((idx + 1) * bin_len)]
    obj.spikewords = data
    param_str = __get_paramstr(obj.animal, obj.probe, obj.date, obj.time_frame, obj.hour_bins, obj.perc, obj.ava_binsize, obj.qualities, obj.cell_types, idx)
    obj.pltname = param_str
    obj.run_crit(flag = flag)
    print(f'BLOCK RESULTS: P_vals - {obj.p_value_burst}   {obj.p_value_t} \n DCC: {obj.dcc}')
    if save:
        to_save = np.array([obj])
        np.save(f'{obj.saveloc}Crit_{param_str}', to_save)
    return obj


def __get_totaltime(time_frame):
    start_time = int(time_frame[0:time_frame.find('_')])
    stop_time = int(time_frame[time_frame.find('_') + 1:])
    total_time = stop_time - start_time
    return total_time


def __get_paramstr(animal, probe, date, time_frame, hour_bins, perc, ava_binsize, quals, cells, idx):
    qual_str = '_'.join(map(str, quals))
    cell_str = '_'.join(cells)
    s = f'{animal}_{probe}_{date}_{time_frame}_{str(hour_bins)}hrs_perc{str(int(perc * 100))}_binsz{str(int(ava_binsize * 1000))}ms_q{qual_str}_cells{cell_str}_{idx}'
    return s


def generate_timeframes(start, end, blocksize):
    ts = np.arange(start, end + 1, blocksize)
    time_frames = [str(a) + "_" + str(b) for a, b in zip(ts[:-1], ts[1:])]
    return time_frames


def save_obj(crit):
    to_save = np.array([crit])
    np.save(f'{crit.saveloc}Crit_{crit.pltname}', to_save)


def signal_handler(signum, frame):
    print("timeout")
    raise Exception('timeout')


def get_info_from_path(path):
    animal_pattern = '((caf|eab)\d{2})'
    matches = re.findall(animal_pattern, path)
    animal = matches[0][0]

    date_pattern = '\d{4}'
    matches = re.findall(date_pattern, path)
    date = matches[0]

    time_frame_pattern = '/\d{1,}_\d{1,}'
    matches = re.findall(time_frame_pattern, path)
    time_frame = matches[0][1:]

    probe = path[path.find('probe'):path.find('probe') + 6]

    return animal, date, time_frame, probe


params = {
    'rerun': True,
    'flag': 2,  # 1 is DCC 2 is p_val and DCC
    'ava_binsz': 0.04,  # in seconds
    'hour_bins': 4,  # durration of block to look at
    'perc': 0.35,
    'nfactor_bm': 0,
    'nfactor_tm': 0,
    'nfactor_bm_tail': .8,  # upper bound to start exclude for burst
    'nfactor_tm_tail': .8,  # upper bound to start exclude for time 
    'cell_type': ['FS', 'RSU'],
    'plot': True
}

