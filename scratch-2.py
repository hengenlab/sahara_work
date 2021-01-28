import glob
import sahara_work as sw
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import seaborn as sns
from criticality_hlab import criticality as cr
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
    'quals': None
}
def __get_paramstr(animal, probe, date, time_frame, hour_bins, perc, ava_binsize, quals, cells, idx):
    qual_str = '_'.join(map(str, quals))
    cell_str = '_'.join(cells)
    s = f'{animal}_{probe}_{date}_{time_frame}_{str(hour_bins)}hrs_perc{str(int(perc * 100))}_binsz{str(int(ava_binsize * 1000))}ms_q{qual_str}_cells{cell_str}_{idx}'
    return s

def __get_totaltime(time_frame):
    start_time = int(time_frame[0:time_frame.find('_')])
    stop_time = int(time_frame[time_frame.find('_') + 1:])
    total_time = stop_time - start_time
    return total_time

def lilo_and_stitch(paths, params, rerun = True, save = False, overlap = False):
    all_objs = []
    errors = []
    for idx, path in enumerate(paths):
        basepath = path[:path.rfind('/')]
        
        print(f'\n\nWorking on ---- {path}', flush = True)
        animal, date, time_frame, probe = sw.get_info_from_path(path)
        print(f'INFO: {animal} -- {date} -- {time_frame} -- {probe}')
        total_time = sw.__get_totaltime(time_frame)
        saveloc = f'/media/HlabShare/clayton_sahara_work/criticality/param_testing/{animal}/{date}/{probe}/'
        if not os.path.exists(saveloc):
            os.makedirs(saveloc)

        if path.find('scored') < 0:
            scorer = 'xgb'
        else:
            scorer = path[path.find('scored')+7:path.find('.npy')]

        num_bins = int(total_time / params['hour_bins'])
        bin_len = int((params['hour_bins'] * 3600) / params['ava_binsz'])


        quals = params['quals']
        fr_cutoff = 50
        try:
            cells = np.load(path, allow_pickle = True)
            good_cells = [cell for cell in cells if cell.quality in quals and cell.cell_type in params['cell_type'] and cell.plotFR(binsz=cell.end_time, lplot=0, lonoff=0)[0][0] < fr_cutoff and cell.presence_ratio() > .99]
            num_cells = len(good_cells)

            print(f'WITH QUALS: {quals}, NCELLS: {num_cells}')
            spikewords = mbt.n_spiketimes_to_spikewords(good_cells, binsz = params['ava_binsz'], binarize = 1, start = start)
        except Exception as err:
            print("Neuron File Won't Load")
            print(err)
            pass
        for idx in np.arange(0, num_bins):
            signal.signal(signal.SIGALRM, sw.signal_handler)
            signal.alarm(600)
            noerr = True
            try:
                print(f'Working on block {idx} --- hours {idx * params["hour_bins"]}-{(idx + 1) * params["hour_bins"]}', flush = True)
                if idx == num_bins - 1:
                    data = spikewords[:, (idx * bin_len):]
                else:
                    data = spikewords[:, (idx * bin_len): ((idx + 1) * bin_len)]

                param_str = sw.__get_paramstr(animal, probe, date, time_frame, params['hour_bins'], params['perc'], params['ava_binsz'], quals, params['cell_type'], idx)
                crit = Crit_hlab(spikewords = data, perc = params['perc'], nfactor_bm = params['nfactor_bm'], nfactor_tm = params['nfactor_tm'],
                            nfactor_bm_tail = params['nfactor_bm_tail'], nfactor_tm_tail = params['nfactor_tm_tail'], saveloc = saveloc,
                            pltname = f'{param_str}_{scorer}', plot = params['plot'])

                crit.run_crit(flag = params['flag'], verbose = False)
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
                    signal.signal(signal.SIGALRM, sw.signal_handler)
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
                        crit.run_crit(flag = params['flag'])

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


saveloc = '/media/HlabShare/clayton_sahara_work/criticality/param_testing/'
subset = np.load('/media/HlabShare/clayton_sahara_work/criticality/paramtesting_subset.npy')

print('RUNNING WITH 1s')
params['quals'] = [1]
results1s, errors1s = lilo_and_stitch(subset, params)

np.save(saveloc + 'results1s.npy', results1s)
np.save(saveloc + 'errors1s.npy', errors1s)

print('RUNNING WITH 1s AND 2s')

params['quals'] = [1,2]
results12s, errors12s = lilo_and_stitch(subset, params)

np.save(saveloc + 'results12s.npy', results12s)
np.save(saveloc + 'errors12s.npy', errors12s)

print('RUNNING WITH 1s AND 2s AND 3s')

params['quals'] = [1,2,3]
results123s, errors123s = lilo_and_stitch(subset, params)

np.save(saveloc + 'results123s.npy', results123s)
np.save(saveloc + 'errors123s.npy', errors123s)
