import glob
import sahara_work as s
import csv
import os
import numpy as np
from datetime import datetime as dt
from datetime import timedelta
from sahara_work import Crit
from sahara_work import *
from sahara_work.crit_hlab import Crit_hlab

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
        cols = ['animal', 'probe', 'date', 'time_frame', 'block_num', 'bday', 'rstart_time', 'age', 'geno', 'p_val_b', 'p_val_t', 'dcc', 'passed', 'kappa_b', 'kappa_t', 'k2b', 'k2t', 'kprob_b', 'kprob_t']
        with open(csvloc, 'w', newline = '') as c:
            w = csv.DictWriter(c, fieldnames = cols)
            w.writeheader()
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
        'eab52': dt(2019, 4, 19, 7, 30),
        'eab47': dt(2019, 2, 17, 7, 30),
        'eab': dt(2019, 2, 17, 7, 30),
        'eab50': dt(2019, 2, 15, 7, 30),
        'eab40': dt(2018, 12, 5, 7, 30)
    }

    return bdays[animal]


def get_genotype(animal):
    genos = {
        'caf01': 'te4',
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
        'caf61': '?',
        'caf62': 'te4',
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

def lilo_and_stitch(paths, params, rerun = False, save = True, overlap = False):
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

            scorer = path[path.find('scored')+7:path.find('.npy')]

            num_bins = int(total_time / params['hour_bins'])
            bin_len = int((params['hour_bins'] * 3600) / params['ava_binsz'])
            quals = [1, 2, 3]
            fr_cutoff = 10
            try:
                cells = np.load(path, allow_pickle = True)
                good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]

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
                if overlap:
                    start = 3600
                else:
                    start = False
                spikewords = mbt.n_spiketimes_to_spikewords(good_cells, binsz = params['ava_binsz'], binarize = 1, start = start)
            except Exception:
                print("Neuron File Won't Load")
                pass
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
                                pltname = param_str, plot = params['plot'])

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
                    crit.scored_by = scorer
                    crit.pathname = path
                    crit.filename = f'{saveloc}Crit_{param_str}_{scorer}'

                except Exception as err:
                    print('TIMEOUT or ERROR', flush = True)
                    print(err)
                    errors.append(f'{animal} -- {probe} -- {date} -- {time_frame} -- {idx} --- {scorer} --- ERRORED')
                    noerr = False
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
                            errors.append(f'{animal} -- {probe} -- {date} -- {time_frame} -- {idx} --- {scorer} --- ERRORED')
                            signal.alarm(0)
                            noerr = False
                            break
                        signal.alarm(0)

                if noerr and save:
                    print(f'BLOCK RESULTS: P_vals - {crit.p_value_burst}   {crit.p_value_t} \n DCC: {crit.dcc}', flush = True)
                    to_save = np.array([crit])
                    np.save(crit.filename, to_save)
                    all_objs.append(crit)

            with open(f'{basepath}/done.txt', 'w+') as f:
                f.write('done')
        else:
            print(f'\n\n{path} -- ALREADY DONE')

    return all_objs, errors
