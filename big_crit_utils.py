import numpy as np
import glob
import pandas as pd
from sahara_work import Crit
from sahara_work import lilo_and_stitch
from sahara_work import lilo_and_stitch_on_blu_ray
import sahara_work as sw
from datetime import datetime as dt
import signal
import sys
import os

def write_to_files(o, csvloc):
    err, appended = sw.write_to_results_csv(o, csvloc)
    if err:
        print('something weird happened, this should not have errored')
    else:
        new_path = o.pathname
        loaded = np.load('/media/HlabShare/clayton_sahara_work/criticality/loaded_paths_results.npy')
        loaded = np.append(loaded, new_path)
        np.save('/media/HlabShare/clayton_sahara_work/criticality/loaded_paths_results.npy', loaded)
    return appended

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
        'plot': True
    }

def run_testing_chpc(paths, params, jobnum=0, animal = '', probe = '', rerun = True, redo = False):
    status_file = f'/scratch/sensley/STATUS_{jobnum}_test.txt'
    csv_file = f'/scratch/sensley/results_{jobnum}.csv'

    all_objs, errors = lilo_and_stitch_on_blu_ray(paths, params, rerun = rerun, save = True, verbose=False)
    results = []
    for o in all_objs:
        appended = write_to_files(o, csv_file)
        results.append(appended)

    if len(all_objs) > 0:
        df = pd.DataFrame(results, columns = ['animal', 'probe', 'date', 'time_frame', 'block_num', 'scored', 'bday', 'rstart_time', 'age', 'geno', 'p_val_b', 'p_val_t', 'dcc', 'passed', 'kappa_b', 'kappa_t', 'k2b', 'k2t', 'kprob_b', 'kprob_t'])
        group = df.groupby(['animal', 'probe', 'date', 'scored'])
        strs = []
        for i, row in group:
            num_passed = row[row["passed"]].count()['passed']
            total_num = row.count()['passed']
            avg_dcc = row.mean()['dcc']
            animal = row['animal'].to_numpy()[0]
            date = row['date'].to_numpy()[0]
            probe = row['probe'].to_numpy()[0]
            scored = row['scored'].to_numpy()[0]
            s = f'{str(animal)} -- {probe} -- {date} -- {scored} -- passed {num_passed}/{total_num} -- avg dcc {avg_dcc}'
            strs.append(s)
    
    now = dt.now()
    with open(status_file, 'a+') as f:
        f.write(f'\n{now.strftime("%d/%m/%Y %H:%M:%S")} ------------ \n')
        f.write(f'{b} PATHS DONE - of this job\n')
        if len(all_objs) > 0: 
            for s in strs:
                f.write(f'{s}\n')
        if len(errors) > 0:
            f.write('\tERRORS:\n')
            for e in errors:
                f.write(f'\t{e[0]}\n')
    return 0
