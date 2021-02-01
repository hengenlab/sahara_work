import glob
import sahara_work as s
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import seaborn as sns
import criticality as cr
from musclebeachtools_hlab import musclebeachtools as mbt
import csv
import os
import numpy as np
from datetime import datetime as dt
from datetime import timedelta
from sahara_work import Crit
from sahara_work.crit_hlab import Crit_hlab
import re
import pandas as pd
import os
import signal
import gc
from copy import deepcopy as cdc

def get_all_results(csvloc, loaded_file, re_load):
    """
    because python refuses to play nice with memory. This will load all crit info into
    a csv to be loaded into a pandas df at a later time

    it doesn't really work. Gonna eventually re-write this correctly but for now - a bandaid
    """
    paths = sorted(glob.glob(f'/media/HlabShare/clayton_sahara_work/criticality/*/*/*/Crit*'))
    print(f'Total # of paths: {len(paths)}')
    error_save = '/media/HlabShare/clayton_sahara_work/criticality/results_errors_TODEL.npy'
    if re_load:
        write_csv_header(csvloc)
        loaded = np.array([])
        errs = np.array([])
    else:
        loaded = np.load(loaded_file)
        errs = np.load(error_save, allow_pickle = True)
        print(f'Number of paths already loaded: {len(loaded)}')

    count = 0
    for i, p in enumerate(paths):
        if i % 5 == 0:
            print(f'#paths: {i}', flush = True)
            print(f'count: {count}', flush = True)

        if count == 100:
            print('saving progress', flush = True)
            np.save(loaded_file, loaded)
            np.save(error_save, errs)
            count = 0

        if (p not in loaded and p not in errs) or re_load:
            count += 1
            err, dat = s.write_to_results_csv_from_path(p, csvloc)
            if err:
                if len(dat) < 2:
                    dat = np.append(dat, p)
                errs = np.append(errs, dat)
            else:
                loaded = np.append(loaded, p)

            if 'LOADED' in p:
                base = p[:p.find('_LOADED')]
                os.rename(p, base + '.npy')

    print('saving final progress', flush = True)
    np.save(loaded_file, loaded)
    np.save(error_save, errs)

    return errs, loaded


def write_to_csv(data, cols, loc):
    d = dict(zip(cols, data))
    with open(loc, 'a', newline = '') as c:
        w = csv.DictWriter(c, fieldnames = cols)
        w.writerow(d)


def write_to_results_csv(crit, loc):
    cols = ['animal', 'probe', 'date', 'time_frame', 'block_num', 'scored', 'bday', 'rstart_time', 'age', 'geno', 'p_val_b', 'p_val_t', 'dcc', 'passed', 'kappa_b', 'kappa_t', 'k2b', 'k2t', 'kprob_b', 'kprob_t']
    err, data = s.lil_helper_boi(crit)
    if err:
        print('this path failed, plz just fucking delete it and re-do this path ffs')
        return err, data
    write_to_csv(data, cols, loc)
    return err, data


def write_csv_header(loc):
    cols = ['animal', 'probe', 'date', 'time_frame', 'block_num', 'scored', 'bday', 'rstart_time', 'age', 'geno', 'p_val_b', 'p_val_t', 'dcc', 'passed', 'kappa_b', 'kappa_t', 'k2b', 'k2t', 'kprob_b', 'kprob_t']
    with open(loc, 'w', newline = '') as c:
        w = csv.DictWriter(c, fieldnames = cols)
        w.writeheader()


def write_to_results_csv_from_path(p, loc):
    err = False
    crit = None
    try:
        crit = np.load(p, allow_pickle = True)
        crit = crit[0]
    except Exception as er:
        print("won't load object", flush = True)
        err = True
        errors = [p, er]
        return err, errors

    err, data = write_to_results_csv(crit, loc)
    return err, data


def lil_helper_boi(crit):
    err = False

    try:
        birth = s.get_birthday(crit.animal)
        start_time = crit.cells[0].rstart_time
        start_time = dt.strptime(start_time, '%Y-%m-%d_%H-%M-%S')
        age = start_time - birth
        age = age + timedelta(hours = int((crit.block_num * crit.hour_bins)))
        geno = s.get_genotype(crit.animal)
        info = [crit.animal, crit.probe, crit.date, crit.time_frame, crit.block_num, crit.scored_by, birth, start_time, age, geno, crit.p_value_burst, crit.p_value_t, crit.dcc, (crit.p_value_burst > 0.05 and crit.p_value_t > 0.05), crit.kappa_burst, crit.kappa_t, crit.k2b, crit.k2t, crit.kprob_b, crit.kprob_t]
    except Exception as e:
        print(f'error: {e}')
        err = True
        info = [e]
    return err, info

def get_paths(scorer = '', geno=None, animal = '', probe = ''):
    s = f'/media/HlabShare/clayton_sahara_work/clustering/{animal}*/*/*/{probe}*/co/'
    print(s)
    basepaths = [f for f in glob.glob(s)]
    print(f'total # of folders: {len(basepaths)}', flush = True)
    og = []
    errors = []
    for f in basepaths:
        neuron_files = glob.glob(f+'*neurons_group0*npy')
        scored_files = glob.glob(f+f'*scored_{scorer}*.npy')
        xgb_files = glob.glob(f+f'*neurons_group0.npy')
        animal, _, _, _ = get_info_from_path(f)
        thisgeno = get_genotype(animal)
        if geno is None or thisgeno in geno:
            if scorer == 'xgb':
                og = np.concatenate([og, xgb_files])
            elif len(scored_files) > 0:
                og = np.concatenate([og, scored_files])
            elif len(neuron_files) == 1 and (scorer=='xgb' or scorer == ''):
                og = np.concatenate([og, neuron_files])
            else:
                errors.append(f)
    return og

def get_birthday(animal):
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
        'caf66': dt(2020, 3, 22, 7, 30),
        'caf69': dt(2020, 1, 5, 7, 30),
        'caf71': dt(2019, 8, 30, 7, 30),
        'caf72': dt(2020, 11, 24, 7, 30),
        'caf73': dt(2020, 1, 5, 7, 30),
        'caf74': dt(2019, 12, 5, 7, 30),
        'caf75': dt(2020, 1, 5, 7, 30),
        'eab52': dt(2019, 4, 19, 7, 30),
        'eab47': dt(2019, 2, 17, 7, 30),
        'eab': dt(2019, 2, 17, 7, 30),
        'eab50': dt(2019, 2, 15, 7, 30),
        'eab40': dt(2018, 12, 5, 7, 30)
    }

    return bdays[animal]

def encode_animal(animal):
    animals = ['caf01', 'caf19', 'caf22', 'caf26', 'caf34', 'caf37', 'caf40', 'caf42', 
                'caf48', 'caf49', 'caf50', 'caf52', 'caf54', 'caf55', 'caf58', 'caf60', 
                'caf61', 'caf62', 'caf66', 'caf69', 'caf71', 'caf72', 'caf73', 'caf74', 
                'caf75','eab52', 'eab47', 'eab', 'eab50', 'eab40']
    nums = np.arange(len(animals))

    keys = dict(zip(animals, nums))

    return keys[animal]

def decode_animal(num):

    animals = ['caf01', 'caf19', 'caf22', 'caf26', 'caf34', 'caf37', 'caf40', 'caf42', 
                'caf48', 'caf49', 'caf50', 'caf52', 'caf54', 'caf55', 'caf58', 'caf60', 
                'caf61', 'caf62', 'caf66', 'caf69', 'caf71', 'caf72', 'caf73', 'caf74',
                 'caf75', 'eab52', 'eab47', 'eab', 'eab50', 'eab40']

    return animals[num]

def get_regions(animal):
    regions = {
        'caf01': ['CA1'],
        'caf19': ['CA1'],
        'caf22': ['V1','CA1'],
        'caf26': ['M1','CA1','S1'],
        'caf34': ['S1','M1','CA1','NAc'],
        'caf37': ['CA1'],
        'caf40': ['CA1'],
        'caf42': ['M1','ACaD','CA1','RSPv','V1'],
        'caf48': ['CA1'],
        'caf49': ['CA1'],
        'caf50': ['CA1'],
        'caf52': ['CA1'],
        'caf54': ['V1'],
        'caf55': ['V1'],
        'caf58': ['CA1'],
        'caf60': ['CA1'],
        'caf61': ['CA1'],
        'caf62': ['CA1'],
        'caf66': ['RSPv','V1','C_Pu','LGN','Sup_col','S1','CA1','M1'],
        'caf69': ['ACaD','RSPv','V1','CA1'],
        'caf71': ['V1','RSPv','CA1','ACaD'],
        'caf72': ['CA1'],
        'caf73': ['ACaD','CA1','RSPv','V1'],
        'caf74': ['ACaD','RSPv','CA1','V1'],
        'caf75': ['ACaD','CA1','RSPv','V1'],
        'eab52': ['CA1','V1'],
        'eab47': ['M1_M2','CA1','V2'],
        'eab': ['M1_M2','CA1','V2'],
        'eab50': ['C_Pu','C_Pu','M1_M2','CA1_DG','CA1_DG','S1','Sup_col','V1_V2'],
        'eab40': ['S1','CA1','M1','M2']
    }
    return regions[animal]


def get_genotype(animal):
    genos = {
        'caf01': 'e4',
        'caf19': 'te4',
        'caf22': 'te4',
        'caf26': 'wt',
        'caf34': 'wt',
        'caf37': 'te4',
        'caf40': 'wt',
        'caf42': 'wt',
        'caf48': 'te4',
        'caf49': 'te4',
        'caf50': 'e4',
        'caf52': 'te4',
        'caf54': 'myt1l',
        'caf55': 'myt1l',
        'caf58': 'e4',
        'caf60': 'te4',
        'caf61': 'e4',
        'caf62': 'te4',
        'caf66': 'wt',
        'caf69': 'wt',
        'caf71': 'app_ps1',
        'caf72': 'te4',
        'caf73': 'app_ps1',
        'caf74': 'app_ps1',
        'caf75': 'app_ps1',
        'eab52': 'te4',
        'eab47': 'te4',
        'eab': 'te4',
        'eab50': 'wt',
        'eab40': 'wt'
    }

    return genos[animal]


def load_crit(path):
    return np.load(path, allow_pickle = True)[0]

def update_object(old, save_new = False):
    new_obj = Crit_hlab(

        spikewords = old.spikewords, perc = old.perc, nfactor_bm = old.nfactor_bm,
        nfactor_tm = old.nfactor_tm, nfactor_bm_tail = old.nfactor_bm_tail, nfactor_tm_tail = old.nfactor_tm_tail,
        saveloc = old.saveloc, pltname = old.pltname, plot = old.plot, burst = old.burst, T = old.T, tm = old.tm,
        p_value_burst = old.p_value_burst, p_value_t = old.p_value_t, dcc = old.dcc, scaling_plot = old.scaling_plot,
        burst_cdf_plot = old.burst_cdf_plot, t_cdf_plot = old.t_cdf_plot, pre = old.pre, fit = old.fit,
        xmin = old.xmin, xmax = old.xmax, tmin = old.tmin, tmax = old.tmax, alpha = old.alpha, time_frame = old.time_frame,
        block_num = old.block_num, qualities = old.qualities, cell_types = old.cell_types, hour_bins = old.hour_bins,
        ava_binsize = old.ava_binsize, animal = old.animal, date = old.date, final = old.final, cells = old.cells,
        probe = old.probe, filename = old.filename, pathname = old.pathname, scored_by = old.scored_by
        )

    if save_new:
        np.save(old.filename, [new_obj])

    return new_obj

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

    date_pattern = '\d{8}'
    matches = re.findall(date_pattern, path)
    date = matches[0]
    date = date[0:4]+date[-2:]

    time_frame_pattern = '/\d{1,}_\d{1,}'
    matches = re.findall(time_frame_pattern, path)
    time_frame = matches[0][1:]

    probe = path[path.find('probe'):path.find('probe') + 6]

    return animal, date, time_frame, probe

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

def get_age_sec(start_time, birthday):
    start_time = dt.strptime(start_time, '%Y-%m-%d_%H-%M-%S')
    age = start_time - birthday
    seconds = age.total_seconds()
    return seconds


params = {
    'flag': 2,  # 1 is DCC 2 is p_val and DCC
    'ava_binsz': 0.04,  # in seconds
    'hour_bins': 4,  # durration of block to look at
    'perc': 0.35,
    'nfactor_bm': 0,
    'nfactor_tm': 0,
    'nfactor_bm_tail': .9,  # upper bound to start exclude for burst
    'nfactor_tm_tail': .9,  # upper bound to start exclude for time 
    'cell_type': ['FS', 'RSU'],
    'plot': True,
    'quals': None, 
    'saveloc': f'/media/HlabShare/clayton_sahara_work/criticality/{animal}/{date}/{probe}/'
}


def lilo_and_stitch(paths, params, rerun = False, save = True, overlap = False, verbose = True):
    all_objs = []
    errors = []
    for idx, path in enumerate(paths):
        basepath = path[:path.rfind('/')]
        
        print(f'\n\nWorking on ---- {path}', flush = True)
        animal, date, time_frame, probe = get_info_from_path(path)
        print(f'INFO: {animal} -- {date} -- {time_frame} -- {probe}')
        total_time = __get_totaltime(time_frame)
        saveloc = f'/media/HlabShare/clayton_sahara_work/criticality/{animal}/{date}/{probe}/'
        if not os.path.exists(saveloc):
            os.makedirs(saveloc)

        if path.find('scored') < 0:
            scorer = 'xgb'
        else:
            scorer = path[path.find('scored')+7:path.find('.npy')]

        num_bins = int(total_time / params['hour_bins'])
        bin_len = int((params['hour_bins'] * 3600) / params['ava_binsz'])


        quals = [1, 2, 3]
        fr_cutoff = 50
        try:
            cells = np.load(path, allow_pickle = True)
            good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]
            num_cells = len(good_cells)
    
            if len(good_cells) < 10:
                quals = [1, 2, 3]
                good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]

            elif len(good_cells) < 60:
                quals = [1, 2]
                good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]
            elif len(good_cells) >= 60:
                quals = [1]
                good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]

                # if len(good_cells) > 100:
                #     cell_idxs = np.random.choice(len(good_cells), 50, replace=False)
                #     good_cells = good_cells[cell_idxs]
            if overlap :
                start = 3600
            else:
                start = False
            spikewords = mbt.n_spiketimes_to_spikewords(good_cells, binsz = params['ava_binsz'], binarize = 1, start = start)
        except Exception as err:
            print("Neuron File Won't Load")
            print(err)
            errors.append([f'{animal} -- {probe} -- {date} -- {time_frame} -- ALL --- {scorer} --- ERRORED', path])
            continue
        for idx in np.arange(0, num_bins):
            signal.signal(signal.SIGALRM, signal_handler)
            signal.alarm(600)
            noerr = True
            try:
                print(f'Working on block {idx} --- hours {idx * params["hour_bins"]}-{(idx + 1) * params["hour_bins"]}', flush = True)
                if idx == num_bins - 1:
                    data = spikewords[:, (idx * bin_len):]
                else:
                    data = spikewords[:, (idx * bin_len): ((idx + 1) * bin_len)]

                param_str = __get_paramstr(animal, probe, date, time_frame, params['hour_bins'], params['perc'], params['ava_binsz'], quals, params['cell_type'], idx)
                crit = Crit_hlab(spikewords = data, perc = params['perc'], nfactor_bm = params['nfactor_bm'], nfactor_tm = params['nfactor_tm'],
                            nfactor_bm_tail = params['nfactor_bm_tail'], nfactor_tm_tail = params['nfactor_tm_tail'], saveloc = saveloc,
                            pltname = f'{param_str}_{scorer}', plot = params['plot'])

                crit.run_crit(flag = params['flag'], verbose = verbose)
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
                crit.scored_by = scorer
                crit.pathname = path
                crit.filename = f'{saveloc}Crit_{param_str}_{scorer}'

            except Exception as err:
                print('TIMEOUT or ERROR', flush = True)
                print(err)
                errors.append([f'{animal} -- {probe} -- {date} -- {time_frame} -- {idx} --- {scorer} --- ERRORED', path])
                noerr = False
                signal.alarm(0)

            if rerun and noerr:
                while crit.p_value_burst < 0.05 or crit.p_value_t < 0.05:
                    signal.signal(signal.SIGALRM, signal_handler)
                    signal.alarm(900)
                    print('\nRERUNNING BLOCK', flush = True)
                    if crit.nfactor_tm_tail < 0.75 or crit.nfactor_bm_tail < 0.75:
                        print('DONE RERUNNNING -- BLOCK WILL NOT PASS\n')
                        signal.alarm(0)
                        break
                    if crit.p_value_burst < 0.05:
                        crit.nfactor_bm_tail -= 0.05
                        #crit.bm += 5
                    if crit.p_value_t < 0.05:
                        crit.nfactor_tm_tail -= 0.05
                        #crit.tm += 5
                    try:
                        crit.run_crit(flag = params['flag'], verbose = verbose)

                    except Exception:
                        print('TIMEOUT or ERROR', flush = True)
                        errors.append([f'{animal} -- {probe} -- {date} -- {time_frame} -- {idx} --- {scorer} --- ERRORED', path])
                        signal.alarm(0)
                        noerr = False
                        break
                    signal.alarm(0)

            if noerr:
                print(f'BLOCK RESULTS: P_vals - {crit.p_value_burst}   {crit.p_value_t} \n DCC: {crit.dcc}', flush = True)
                if save:
                    to_save = np.array([crit])
                    np.save(crit.filename, to_save)
                all_objs.append(crit)

        with open(f'{basepath}/done.txt', 'w+') as f:
            f.write('done')

    return all_objs, errors



params = {
    'flag': 2,  # 1 is DCC 2 is p_val and DCC
    'ava_binsz': 0.04,  # in seconds
    'hour_bins': 4,  # durration of block to look at
    'perc': 0.35,
    'nfactor_bm': 0,
    'nfactor_tm': 0,
    'nfactor_bm_tail': .9,  # upper bound to start exclude for burst
    'nfactor_tm_tail': .9,  # upper bound to start exclude for time 
    'cell_type': ['FS', 'RSU'],
    'plot': True,
    'quals': None, 
    'saveloc': f'/scratch/sensley/crit_testing/results/',
    'animal': 'caf37',
    'probe': 'probe1',
    'time_frame':'0_12',
    'date':'01012021'
}
# on the chpc, lil annoying but what can you do
def lilo_and_stitch_on_blu_ray(paths, params, rerun = False, save = True, overlap = False, verbose = False):
    all_objs = []
    errors = []
    for idx, path in enumerate(paths):
        basepath = path[:path.rfind('/')]
        
        print(f'\n\nWorking on ---- {path}', flush = True)
        animal = params['animal']
        date = params['date']
        time_frame = params['time_frame']
        probe = params['probe']
        print(f'INFO: {animal} -- {date} -- {time_frame} -- {probe}')
        total_time = __get_totaltime(time_frame)
        saveloc = params['saveloc']
        if not os.path.exists(saveloc):
            os.makedirs(saveloc)

        if path.find('scored') < 0:
            scorer = 'xgb'
        else:
            scorer = path[path.find('scored')+7:path.find('.npy')]

        num_bins = int(total_time / params['hour_bins'])
        bin_len = int((params['hour_bins'] * 3600) / params['ava_binsz'])


        quals = [1, 2, 3]
        fr_cutoff = 50
        try:
            cells = np.load(path, allow_pickle = True)
            good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]
            num_cells = len(good_cells)
    
            if len(good_cells) < 10:
                quals = [1, 2, 3]
                good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]

            elif len(good_cells) < 60:
                quals = [1, 2]
                good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]
            elif len(good_cells) >= 60:
                quals = [1]
                good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]

                # if len(good_cells) > 100:
                #     cell_idxs = np.random.choice(len(good_cells), 50, replace=False)
                #     good_cells = good_cells[cell_idxs]
            if overlap :
                start = 3600
            else:
                start = False
            spikewords = mbt.n_spiketimes_to_spikewords(good_cells, binsz = params['ava_binsz'], binarize = 1, start = start)
        except Exception as err:
            print("Neuron File Won't Load")
            print(err)
            errors.append([f'{animal} -- {probe} -- {date} -- {time_frame} -- ALL --- {scorer} --- ERRORED', path])
            continue
        for idx in np.arange(0, num_bins):
            signal.signal(signal.SIGALRM, signal_handler)
            signal.alarm(600)
            noerr = True
            try:
                print(f'Working on block {idx} --- hours {idx * params["hour_bins"]}-{(idx + 1) * params["hour_bins"]}', flush = True)
                if idx == num_bins - 1:
                    data = spikewords[:, (idx * bin_len):]
                else:
                    data = spikewords[:, (idx * bin_len): ((idx + 1) * bin_len)]

                param_str = __get_paramstr(animal, probe, date, time_frame, params['hour_bins'], params['perc'], params['ava_binsz'], quals, params['cell_type'], idx)
                crit = Crit_hlab(spikewords = data, perc = params['perc'], nfactor_bm = params['nfactor_bm'], nfactor_tm = params['nfactor_tm'],
                            nfactor_bm_tail = params['nfactor_bm_tail'], nfactor_tm_tail = params['nfactor_tm_tail'], saveloc = saveloc,
                            pltname = f'{param_str}_{scorer}', plot = params['plot'])

                crit.run_crit(flag = params['flag'], verbose = verbose)
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
                crit.scored_by = scorer
                crit.pathname = path
                crit.filename = f'{saveloc}Crit_{param_str}_{scorer}'

            except Exception as err:
                print('TIMEOUT or ERROR', flush = True)
                print(err)
                errors.append([f'{animal} -- {probe} -- {date} -- {time_frame} -- {idx} --- {scorer} --- ERRORED', path])
                noerr = False
                signal.alarm(0)

            if rerun and noerr:
                while crit.p_value_burst < 0.05 or crit.p_value_t < 0.05:
                    signal.signal(signal.SIGALRM, signal_handler)
                    signal.alarm(900)
                    print('\nRERUNNING BLOCK', flush = True)
                    if crit.nfactor_tm_tail < 0.75 or crit.nfactor_bm_tail < 0.75:
                        print('DONE RERUNNNING -- BLOCK WILL NOT PASS\n')
                        signal.alarm(0)
                        break
                    if crit.p_value_burst < 0.05:
                        crit.nfactor_bm_tail -= 0.05
                        #crit.bm += 5
                    if crit.p_value_t < 0.05:
                        crit.nfactor_tm_tail -= 0.05
                        #crit.tm += 5
                    try:
                        crit.run_crit(flag = params['flag'], verbose = verbose)

                    except Exception:
                        print('TIMEOUT or ERROR', flush = True)
                        errors.append([f'{animal} -- {probe} -- {date} -- {time_frame} -- {idx} --- {scorer} --- ERRORED', path])
                        signal.alarm(0)
                        noerr = False
                        break
                    signal.alarm(0)

            if noerr:
                print(f'BLOCK RESULTS: P_vals - {crit.p_value_burst}   {crit.p_value_t} \n DCC: {crit.dcc}', flush = True)
                if save:
                    to_save = np.array([crit])
                    np.save(crit.filename, to_save)
                all_objs.append(crit)

        with open(f'{basepath}/done.txt', 'w+') as f:
            f.write('done')

    return all_objs, errors
