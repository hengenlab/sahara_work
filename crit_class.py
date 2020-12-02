import numpy as np
import criticality as cr
import musclebeachtools as mbt
import re
import pandas as pd
import os 
import glob
import signal
import gc
from datetime import datetime as dt 
from datetime import timedelta
import csv

class Crit:
    """
    Class to look at criticality stuff

    tbd on init and what not
    """

    def __init__(self, spikewords, perc = 0.35, nfactor_bm = 0, nfactor_tm = 0, nfactor_bm_tail = 1, nfactor_tm_tail = 1, saveloc = '',pltname='', plot = True):

        #required parameters
        self.perc = perc
        self.spikewords = spikewords
        n,b = np.shape(spikewords)
        self.num_cells = n
        self.num_bins = b
        self.saveloc = saveloc
        self.pltname = pltname
        self.plot = plot
        self.nfactor_bm = nfactor_bm
        self.nfactor_tm = nfactor_tm
        self.nfactor_bm_tail = nfactor_bm_tail
        self.nfactor_tm_tail = nfactor_tm_tail

        #all parameters set by run_crit
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

        #all parameters that can be set by lilo and stitch
        self.time_frame = None
        self.block_num = None
        self.qualities = None
        self.cell_types = None
        self.hour_bins = None
        self.ava_binsize = None
        self.animal = None
        self.date = None
        self.final = False
        self.cells = []
        self.probe = None

    def __repr__(self):
        '''
        how to generate crit
        '''
        return str(f'CRIT: {self.animal} -- {self.date} -- {self.time_frame}')
    
    def __enter__(self):
        return self
    
    def __exit__(self, exception_type, exception_value, traceback):
        print('closing file')
        
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
    
    def finalize(self):
        self.cells = None
        self.final = True

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

        self.bm = int(np.max(burst)/20)
        self.tm = int(np.max(T)/20)


        crit_params = {
            'flag' : flag,
            'pltname': self.pltname,
            'saveloc' : self.saveloc,
            'plot' : self.plot,
            'bm' : self.bm,
            'tm' : self.tm
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
        self.gen_kappa()
        self.gen_k2()
        self.gen_kprob()

    def run_crit_from_start(self, flag = 2, save=False):
        if self.final:
            print('This crit object is final, there are no cells saved here. If youd like to rerun this block start from lilo_and_stitch')
            return
        total_time = __get_totaltime(self.time_frame)
        num_bins = int(total_time/self.hour_bins)
        bin_len = int((self.hour_bins * 3600) / self.ava_binsize)
        good_cells = [cell for cell in self.cells if self.quality in self.qualities and self.cell_type in self.cell_types]
        spikewords = mbt.n_spiketimes_to_spikewords(good_cells, binsz = self.ava_binsize, binarize = 1)
        idx = self.block_num
        if idx == num_bins - 1:
            data = spikewords[:, (idx * bin_len):]
        else:
            data = spikewords[:, (idx * bin_len): ((idx + 1) * bin_len)]
        self.spikewords = data
        param_str = __get_paramstr(self.animal, self.probe, self.date, self.time_frame, self.hour_bins, self.perc, self.ava_binsize, self.qualities, self.cell_types, idx)
        self.pltname = param_str
        self.run_crit(flag = flag)
        print(f'BLOCK RESULTS: P_vals - {self.p_value_burst}   {self.p_value_t} \n DCC: {self.dcc}')
        if save:
            to_save = np.array([obj])
            np.save(f'{obj.saveloc}Crit_{param_str}', to_save)

    def gen_beta(self):
        tm = int(np.max(self.T)/20)
        _, _, beta = \
            cr.EXCLUDE(self.T[self.T < np.power(np.max(self.T), self.nfactor_tm_tail)], tm,
                   nfactor=self.nfactor_tm)
        self.beta = beta
    
    def gen_k2(self, num=100):
        if self.burst is None:
            print("You must run_crit() before you can run this function")
            return
        if self.xmin is None:
            print('this block failed while running for some reason')
            return
        idx = np.where(np.logical_and(self.burst<=self.xmax, self.burst>=self.xmin)) 
        b = self.burst[idx]
        n = np.size(b)
        cdf = np.cumsum(np.histogram(b, np.arange(self.xmin, self.xmax+2))[0]/n)
        s = np.unique(b)
        A = 1/np.sum(np.power(s, -self.alpha))
        fit = np.cumsum(A*np.power(np.arange(self.xmin, self.xmax+1), -self.alpha)) 

        idxs = np.geomspace(1, np.size(cdf)-1, num=num, dtype=int)
        half = int(num/2)
        second_half = idxs[half:]
        diffs = [cdf[i]-fit[i] for i in second_half]
        mean_diff = np.mean(diffs)
        k2b = mean_diff

        self.k2b = k2b
        try:
            idx2 = np.where(np.logical_and(self.T<=self.tmax, self.T>=self.tmin)) 
            t = self.T[idx2]
            n = np.size(t)
            cdf = np.cumsum(np.histogram(t, np.arange(self.tmin, self.tmax+2))[0]/n)
            s = np.unique(t)
            A = 1/np.sum(np.power(s, -self.beta))
            fit = np.cumsum(A*np.power(np.arange(self.tmin, self.tmax+1), -self.beta)) 

            idxs = np.geomspace(1, np.size(cdf)-1, num=num, dtype=int)
            half = int(num/2)
            second_half = idxs[half:]
            diffs = [cdf[i]-fit[i] for i in second_half]
            mean_diff = np.mean(diffs)
            k2t = mean_diff

            self.k2t = k2t
        except AttributeError:
            self.k2t = None
        return self.k2b, self.k2t
            
    def gen_kappa(self, num = 100):
        if self.burst is None:
            print("You must run_crit() before you can run this function")
            return
        if self.xmin is None:
            print('block failed')
            return
        idx = np.where(np.logical_and(self.burst<=self.xmax, self.burst>=self.xmin)) 
        b = self.burst[idx]
        n = np.size(b)
        cdf = np.cumsum(np.histogram(b, np.arange(self.xmin, self.xmax+2))[0]/n)
        s = np.unique(b)
        A = 1/np.sum(np.power(s, -self.alpha))
        fit = np.cumsum(A*np.power(np.arange(self.xmin, self.xmax+1), -self.alpha)) 

        idxs = np.geomspace(1, np.size(cdf)-1, num=num, dtype=int)
        diffs = [fit[i]-cdf[i] for i in idxs]
        mean_diff = np.mean(diffs)
        kappa_burst = 1 + mean_diff

        self.kappa_burst = kappa_burst
        try:
            idx2 = np.where(np.logical_and(self.T<=self.tmax, self.T>=self.tmin)) 
            t = self.T[idx2]
            n = np.size(t)
            cdf = np.cumsum(np.histogram(t, np.arange(self.tmin, self.tmax+2))[0]/n)
            s = np.unique(t)
            A = 1/np.sum(np.power(s, -self.beta))
            fit = np.cumsum(A*np.power(np.arange(self.tmin, self.tmax+1), -self.beta)) 

            idxs = np.geomspace(1, np.size(cdf)-1, num=num, dtype=int)
            diffs = [fit[i]-cdf[i] for i in idxs]
            mean_diff = np.mean(diffs)
            kappa_t = 1 + mean_diff

            self.kappa_t = kappa_t
        except AttributeError:
            self.kappa_t = None
        return self.kappa_burst, self.kappa_t

    def gen_kprob(self):
        if self.burst is None:
            print("You must run_crit() before you can run this function")
            return
        if self.xmin is None:
            print('this block failed while running for some reason')
            return

        pdf = np.histogram(self.burst, bins = np.arange(1, np.max(self.burst) + 2))[0]
        p = pdf / np.sum(pdf)
        x = np.arange(self.xmin, self.xmax + 1)
        y = (np.size(np.where(self.burst == self.xmin + 6)[0]) / np.power(self.xmin + 6, -self.alpha)) *\
            np.power(x, -self.alpha)
        y = y / np.sum(pdf)

        ps=p[self.xmin:self.xmax]
        diffs = np.diff([ps,y[:-1]],axis=0)
        above_idx = np.where(diffs<0)[1]
        below_idx = np.where(diffs>0)[1]
        prob_above = np.sum(ps[above_idx])
        prob_below = np.sum(ps[below_idx])

        kprob_b = prob_above/prob_below
        self.kprob_b = kprob_b
        try:
            tdf = np.histogram(self.T, bins = np.arange(1, np.max(self.T) + 2))[0]
            t = tdf / np.sum(tdf)
            x = np.arange(self.tmin, self.tmax + 1)
            y = (np.size(np.where(self.T == self.tmin + 6)[0]) / np.power(self.tmin + 6, -self.beta)) *\
                np.power(x, -self.beta)
            y = y / np.sum(tdf)

            ps=t[self.tmin:self.tmax]
            diffs = np.diff([ps,y[:-1]],axis=0)
            above_idx = np.where(diffs<0)[1]
            below_idx = np.where(diffs>0)[1]
            prob_above = np.sum(ps[above_idx])
            prob_below = np.sum(ps[below_idx])

            kprob_t = prob_above/prob_below
            self.kprob_t = kprob_t
        except Exception:
            print('kprob_t not working')
            self.kprob_t = None

def get_cell_stats(cell):
    fr, xbins = cell.plotFR(binsz=3600, start=False, end=False,
               lplot=0)

    isis = []
    bins = np.arange(cell.start_time, cell.end_time, 300) # 5 min bins

    for i in range(len(bins)):
        if i == len(bins)-1:
            end_bin = cell.end_time
        else:
            end_bin = bins[i+1]
        spk_idxs = np.where(np.logical_and(cell.spike_time_sec > bins[i], cell.spike_time_sec < end_bin))
        isis.append(np.diff(cell.spike_time[spk_idxs]))

    means = np.array([np.mean(i) for i in isis])
    stds = np.array([np.std(i) for i in isis])
    cvs = stds/means
    binned_cvs = [cvs[i*12:(i+1)*12] for i in np.arange(int(xbins[-1]))]
    cv = [np.mean(i) for i in binned_cvs]
    if len(fr) != len(cv):
        print('this isnt working, fix it')
    return fr, cv

def construct_fr_df(paths):
    bdays = {
        'caf01':dt(2019, 12, 24, 7, 30),
        'caf19':dt(2020, 1, 19, 7, 30),
        'caf22':dt(2020, 2, 17, 7, 30),
        'caf26':dt(2020, 2, 20, 7, 30),
        'caf34':dt(2020, 3, 18, 7, 30),
        'caf37':dt(2019, 8, 18, 7, 30),
        'caf40':dt(2020, 2, 20, 7, 30),
        'caf42':dt(2020, 2, 20, 7, 30),
        'caf48':dt(2020, 7, 20, 7, 30),
        'caf49':dt(2020, 7, 20, 7, 30),
        'caf50':dt(2020, 7, 20, 7, 30),
        'caf52':dt(2020, 4, 19, 7, 30),
        'caf54':dt(2020, 7, 11, 7, 30),
        'caf55':dt(2020, 7, 11, 7, 30),
        'caf58':dt(2020, 9, 23, 7, 30),
        'caf60':dt(2020, 9, 23, 7, 30),
        'caf61':dt(2019, 12, 11, 7, 30),
        'caf62':dt(2019, 11, 18, 7, 30),
        'eab52':dt(2019, 4, 19, 7, 30),
        'eab47':dt(2019, 2, 17, 7, 30),
        'eab':dt(2019, 2, 17, 7, 30),
        'eab50':dt(2019, 2, 15, 7, 30),
        'eab40':dt(2018, 12, 5, 7, 30)
    }
    seconds_in_day = 60*60*24
    with open('/media/HlabShare/clayton_sahara_work/criticality/cell_stats.csv', 'w', newline='') as c:
        w =  csv.DictWriter(c, fieldnames=['animal', 'rstart_time', 'age_start', 'days_old', 'hours_old', 'cell_idx', 'quality', 'fr', 'cv', 'wf'])
        w.writeheader()
    
    done = []
    final = []
    for i,p in enumerate(paths):
        print(i)
        path_info = p[:p.find('_perc')]
        if path_info in done:
            print('already done')
        else:
            done.append(path_info)
            crit = None
            try:
                crit = np.load(p, allow_pickle=True)[0]
                birth = bdays[crit.animal]
                for cell in crit.cells:
                    start_time = crit.cells[0].rstart_time
                    start_time = dt.strptime(start_time, '%Y-%m-%d_%H-%M-%S')
                    age = start_time - birth

                    fr, cv = get_cell_stats(cell)

                    for i in range(len(fr)):

                        age_now = age + timedelta(hours=i)
                        days_old = age_now.total_seconds()/seconds_in_day
                        hours_old = age_now.total_seconds()/3600

                        with open('/media/HlabShare/clayton_sahara_work/criticality/cell_stats.csv', 'a', newline='') as c:
                            w =  csv.DictWriter(c, fieldnames=['animal', 'rstart_time', 'age_start', 'days_old', 'hours_old', 'cell_idx', 'quality', 'fr', 'cv', 'wf'])

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


def lil_helper_boi(p):
    err = False
    try:
        crit = np.load(p, allow_pickle=True)
        crit = crit[0]
    except Exception as er:
        print("won't load object", flush = True)
        err = True
        errors = [p, errors]
    try:
        to_append = [crit.animal, crit.probe, crit.date, crit.time_frame, crit.block_num, crit.p_value_burst, crit.p_value_t, crit.dcc, (crit.p_value_burst > 0.05 and crit.p_value_t > 0.05), crit.kappa_burst, crit.kappa_t, crit.k2b, crit.k2t, crit.kprob_b, crit.kprob_t]
    except Exception:
        try:
            to_append = [crit.animal, probe, crit.date, crit.time_frame, crit.block_num, crit.p_value_burst, crit.p_value_t, crit.dcc, (crit.p_value_burst > 0.05 and crit.p_value_t > 0.05), crit.kappa_burst, crit.kappa_t, crit.k2b, crit.k2t, crit.kprob_b, crit.kprob_t]
        except Exception as er:
            print(f"not going to work --- skipping this path {p}", flush=True)
            err = True
            errors = [p, er]
            return err, errors

    return err, to_append

def lil_helper_dude(crit):
    err = False
    
    birth = bdays[crit.animal]
    start_time = crit.cells[0].rstart_time
    start_time = dt.strptime(start_time, '%Y-%m-%d_%H-%M-%S')
    age = start_time - birth
    age = age + timedelta(hours = (crit.block_num*crit.hour_bins))
    geno = genos[crit.animal]
    info = [crit.animal, crit.probe, crit.date, crit.time_frame, crit.block_num, birth, start_time, age, geno, crit.p_value_burst, crit.p_value_t, crit.dcc, (crit.p_value_burst > 0.05 and crit.p_value_t > 0.05), crit.kappa_burst, crit.kappa_t, crit.k2b, crit.k2t, crit.kprob_b, crit.kprob_t]

    return info

def write_to_csv(data, cols, loc):
    d = dict(zip(cols, data))
    with open(loc, 'a', newline='') as c:
        w =  csv.DictWriter(c, fieldnames=cols)
        w.writerow(d)

def write_to_results_csv(crit, loc):
    cols = ['animal', 'probe', 'date', 'time_frame', 'block_num','bday','rstart_time', 'age', 'geno', 'p_val_b', 'p_val_t', 'dcc', 'passed', 'kappa_b', 'kappa_t', 'k2b', 'k2t', 'kprob_b', 'kprob_t']
    err, data = lil_helper_dude(crit)
    if err:
        print('this path failed, plz just fucking delete it and re-do this path ffs')
        return err, []
    write_to_csv(data, cols, loc)
    return err, data

def get_results(animal,probe='', paths = None, save=False, saveloc='', re_load = False):
    if paths is None:  
        paths = glob.glob(f'/media/HlabShare/clayton_sahara_work/criticality/{animal}*/*/{probe}*/Crit*')
    
    results = []
    print(f'Total # of paths: {len(paths)}')
    errs = []

    for i,p in enumerate(paths):
        good=True
        if i%5 == 0:
            print(f'#paths: {i}', flush = True)
        if i%100==0 and i!=0:
            gc.collect()

        if 'LOADED' not in p or re_load:

            err, to_append = lil_helper_boi(p)
            if err:
                print(to_append)
                errs.append(to_append)
            else:
                results.append(to_append)

            if 'LOADED' in p:
                base = p[:p.find('_LOADED')]
            else:
                base = p[:-4]

            n = base+'_LOADED.npy'
            os.rename(p, n)
            
    cols = ['animal', 'probe', 'date', 'time_frame', 'block_num', 'p_val_b', 'p_val_t', 'dcc', 'passed', 'kappa_b', 'kappa_t', 'k2b', 'k2t', 'kprob_b', 'kprob_t']
    df = pd.DataFrame(results, columns = cols)
    df_clean = df.drop_duplicates(subset=['animal','probe','date', 'time_frame', 'block_num'], keep = 'last')
    
    if not re_load:
        og_df = pd.read_pickle(f'{saveloc}/ALL_ANIMALS_all_results.pkl')
        df_clean = pd.concat([df_clean, og_df])
        df_clean = df_clean.drop_duplicates(subset=['animal','probe','date', 'time_frame', 'block_num'], keep = 'last')
        print(f'OG SIZE: {og_df.size()}', flush = True)
        print(f'NEW SIZE: {df_clean.size()}', flush = True)
    if save:
        if animal=='':
            df_clean.to_pickle(f'{saveloc}/ALL_ANIMALS_all_results.pkl')
        else:
            df_clean.to_pickle(f'{saveloc}/{animal}_all_results.pkl')

    return df_clean, errs

def run_crit_from_start(obj, flag = 2, save=True):
    if obj.final:
        print('This crit object is final, there are no cells saved here. If youd like to rerun this block start from lilo_and_stitch')
        return
    total_time = __get_totaltime(obj.time_frame)
    num_bins = int(total_time/obj.hour_bins)
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
    stop_time = int(time_frame[time_frame.find('_')+1:])
    total_time = stop_time-start_time
    return total_time

def __get_paramstr(animal,probe, date, time_frame, hour_bins, perc, ava_binsize, quals, cells, idx):
    qual_str = '_'.join(map(str, quals))
    cell_str = '_'.join(cells)
    s = f'{animal}_{probe}_{date}_{time_frame}_{str(hour_bins)}hrs_perc{str(int(perc*100))}_binsz{str(int(ava_binsize*1000))}ms_q{qual_str}_cells{cell_str}_{idx}'
    return s

def generate_timeframes(start, end, blocksize):
    ts = np.arange(start, end+1, blocksize)
    time_frames = [str(a)+"_"+str(b) for a, b in zip(ts[:-1], ts[1:])]
    return time_frames

def save_obj(crit):
    to_save = np.array([crit])
    np.save(f'{crit.saveloc}Crit_{crit.pltname}', to_save)

def finalize_all_paths(paths):
    for c in objs: 
        crit = np.load(c, allow_pickle=True)[0] 
        try: 
            if not crit.final: 
                crit.finalize() 
                print(f'{crit.time_frame} - finalized') 
                save_obj(crit) 
        except AttributeError: 
            print(f'{crit.time_frame} - old version no cells') 

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

    probe = path[path.find('probe'):path.find('probe')+6]

    return animal, date, time_frame, probe

params = {
    'rerun' : True,
    'flag': 2, # 1 is DCC 2 is p_val and DCC
    'ava_binsz': 0.04, # in seconds
    'hour_bins': 4,# durration of block to look at
    'perc': 0.35,
    'nfactor_bm':0, 
    'nfactor_tm':0,
    'nfactor_bm_tail':.8, # upper bound to start exclude for burst
    'nfactor_tm_tail': .8, # upper bound to start exclude for time 
    'cell_type': ['FS', 'RSU'], 
    'plot' : True
    }

def lilo_and_stitch_the_sequel(paths, params, save=True, delete=True):
    all_objs = []
    errors = []
    for p in paths:
        c = np.load(p, allow_pickle = True)
        from_start = False
        if c.cell_types != params['cell_type'] or c.qualities != params['quality'] or c.ava_binsize != params['ava_binsz'] or c.hour_bins != params['hour_bins']:
            from_start = True
            c.ava_binsize = params['ava_binsz']
            c.hour_bins = params['hour_bins']
            c.qualities = params['quality']
            c.cell_types = params['cell_type']
        c.perc = params['perc']
        c.nfactor_bm = params['nfactor_bm']
        c.nfactor_tm = params['nfactor_tm']
        c.nfactor_bm_tail = params['nfactor_bm_tail']
        c.nfactor_tm_tail = params['nfactor_tm_tail']

        signal.signal(signal.SIGALRM, signal_handler)
        signal.alarm(600)
        try:
            if from_start:
                c_new = run_crit_from_start(obj, flag = params['flag'], save=save)
                all_objs.append(c_new)
                
            else:
                c.run_crit()
                if save:
                    save_obj(c)
                all_objs.append(c)
            if delete:
                os.remove(p)
        except Exception:
            print('TIMEOUT or ERROR')
            errors.append(f'{c.animal} -- {c.date} -- {c.time_frame} -- {c.block_num} --- ERRORED')

    return all_objs, errors


def lilo_and_stitch(paths, params, rerun=False, save=True):
    all_objs = []
    errors = []
    for idx, path in enumerate(paths):
        basepath = path[:path.rfind('/')]
        if not os.path.exists(f'{basepath}/done.txt') or params['rerun']:
            
            print(f'\n\nWorking on ---- {path}', flush = True)
            animal, date, time_frame, probe = get_info_from_path(path)
            print(f'INFO: {animal} -- {date} -- {time_frame} -- {probe}')
            total_time = __get_totaltime(time_frame)
            saveloc = f'/media/HlabShare/clayton_sahara_work/criticality/{animal}/{date}/{probe}/'
            if not os.path.exists(saveloc):
                os.makedirs(saveloc)

            num_bins = int(total_time/params['hour_bins'])
            bin_len = int((params['hour_bins'] * 3600) / params['ava_binsz'])
            quals = [1,3]
            try:
                cells = np.load(path, allow_pickle = True)
                good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type']]
                if len(good_cells) < 10:
                    quals = [1,2,3]
                    good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type']]
                spikewords = mbt.n_spiketimes_to_spikewords(good_cells, binsz = params['ava_binsz'], binarize = 1)
            except Exception:
                print("Neuron File Won't Load")
                pass
            for idx in np.arange(0, num_bins):
                signal.signal(signal.SIGALRM, signal_handler)
                signal.alarm(600)
                noerr = True
                try:
                    print(f'Working on block {idx} --- hours {idx*params["hour_bins"]}-{(idx+1)*params["hour_bins"]}', flush = True)
                    if idx == num_bins - 1:
                        data = spikewords[:, (idx * bin_len):]
                    else:
                        data = spikewords[:, (idx * bin_len): ((idx + 1) * bin_len)]

                    param_str = __get_paramstr(animal,probe, date, time_frame, params['hour_bins'], params['perc'], params['ava_binsz'], quals, params['cell_type'], idx)
                    crit = Crit(data, perc = params['perc'], nfactor_bm = params['nfactor_bm'], nfactor_tm = params['nfactor_tm'],
                                nfactor_bm_tail = params['nfactor_bm_tail'], nfactor_tm_tail = params['nfactor_tm_tail'], saveloc = saveloc,
                                pltname=param_str, plot = params['plot'])

                    crit.run_crit(flag = params['flag'])
                    crit.time_frame = time_frame
                    crit.block_num = idx
                    crit.qualities = quals
                    crit.cell_types = params['cell_type']
                    crit.hour_bins = params['hour_bins']
                    crit.ava_binsize = params['ava_binsz']
                    crit.animal = animal
                    crit.date = date
                    crit.final = False
                    crit.cells = [cell for cell in cells if cell.quality < 4]
                    crit.probe = probe
                    
                except Exception:
                    print('TIMEOUT or ERROR', flush = True)
                    errors.append(f'{animal} -- {probe} -- {date} -- {time_frame} -- {idx} --- ERRORED')
                    noerr=False
                    signal.alarm(0)

                if rerun and noerr:
                    while crit.p_value_burst < 0.05 or crit.p_value_t < 0.05:
                        signal.signal(signal.SIGALRM, signal_handler)
                        signal.alarm(600)
                        print('\nRERUNNING BLOCK', flush = True)
                        if crit.nfactor_tm_tail < 0.75 or crit.nfactor_bm_tail < 0.75:
                            print('DONE RERUNNNING -- BLOCK WILL NOT PASS\n')
                            signal.alarm(0)
                            break
                        if crit.p_value_burst < 0.05:
                            crit.nfactor_bm_tail -= 0.05
                        if crit.p_value_t < 0.05:
                            crit.nfactor_tm_tail -= 0.05
                        try:
                            crit.run_crit(flag = params['flag'])
                        
                        except Exception:
                            print('TIMEOUT or ERROR', flush = True)
                            errors.append(f'{animal} -- {probe} -- {date} -- {time_frame} -- {idx} --- ERRORED')
                            signal.alarm(0)
                            noerr=False
                            break
                        signal.alarm(0)
                
                if noerr and save:
                    print(f'BLOCK RESULTS: P_vals - {crit.p_value_burst}   {crit.p_value_t} \n DCC: {crit.dcc}', flush = True)
                    to_save = np.array([crit])
                    np.save(f'{saveloc}Crit_{param_str}', to_save)
                    all_objs.append(crit)

            with open(f'{basepath}/done.txt', 'w+') as f:
                f.write('done')
        else:
            print(f'\n\n{path} -- ALREADY DONE')

    return all_objs, errors

