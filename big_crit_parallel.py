import numpy as np
import glob
import pandas as pd
from sahara_work import Crit
from sahara_work import lilo_and_stitch
import sahara_work as sw
from datetime import datetime as dt
import signal
import sys
import os
import multiprocessing as mp
import time

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

def run_linear(paths, csvloc, redo, rerun):
    params = {
        'redo_paths': redo,
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

    all_objs, errors = lilo_and_stitch(paths, params, rerun = rerun, save = False)
    results = []
    for o in all_objs:
        appended = write_to_files(o, csvloc)
        results.append(appended)

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

    if len(all_objs) > 0:
        with open('/media/HlabShare/clayton_sahara_work/criticality/status_test.txt', 'a+') as f:
            f.write(f'\n{now.strftime("%d/%m/%Y %H:%M:%S")} ------------ \n')
            f.write(f'{b} PATHS DONE - of this job')
            for s in strs:
                f.write(f'{s}\n')
            f.write('\tERRORS:\n')
            for e in errors:
                f.write(f'\t{e}\n')


def run(paths, csvloc, b, redo = False, rerun = True):
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

    all_objs, errors = lilo_and_stitch(paths, params, rerun = rerun, save = True)
    results = []
    for o in all_objs:
        LOCK.acquire()
        appended = write_to_files(o, csvloc)
        LOCK.release()
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
    LOCK2.acquire()
    with open('/media/HlabShare/clayton_sahara_work/criticality/STATUS.txt', 'a+') as f:
        f.write(f'\n{now.strftime("%d/%m/%Y %H:%M:%S")} ------------ \n')
        f.write(f'{b} PATHS DONE - of this job\n')
        f.write(f'worker:\t{mp.current_process()}\n')
        if len(all_objs) > 0: 
            for s in strs:
                f.write(f'{s}\n')
        if len(errors) > 0:
            f.write('\tERRORS:\n')
            for e in errors:
                f.write(f'\t{e[0]}\n')
                errored = np.load('/media/HlabShare/clayton_sahara_work/criticality/errored_paths.npy')
                errored = np.append(errored, e[1])
                np.save('/media/HlabShare/clayton_sahara_work/criticality/errored_paths.npy', errored)

    LOCK2.release()
    return 0

def get_paths(animal, probe, redo, csvloc):

    s = f'/media/HlabShare/clayton_sahara_work/clustering/{animal}*/*/{probe}*/co/*scored_*.npy'
    print(s)
    og = [f for f in glob.glob(s)]
    print(f'total # of paths: {len(og)}', flush = True)
    if redo:
        paths = og
        sw.write_csv_header(csvloc)
        loaded = []
        np.save('/media/HlabShare/clayton_sahara_work/criticality/loaded_paths_results.npy', loaded)
        errors = []
        np.save('/media/HlabShare/clayton_sahara_work/criticality/errored_paths.npy', errors)
    else:
        loaded = np.load('/media/HlabShare/clayton_sahara_work/criticality/loaded_paths_results.npy')
        errored = np.load('/media/HlabShare/clayton_sahara_work/criticality/errored_paths.npy')
        paths = []
        for p in og:
            if p not in loaded and p not in errored:
                paths.append(p)
    paths = sorted(paths)

    print(f'Number of paths left to analyze: {len(paths)}', flush = True)

    return paths

def test_running(paths, csvloc, temppath, id):
    print(f'starting: {id}')
    results = []
    for p in paths:
        results.append(sw.get_info_from_path(p))

    time.sleep(2)
    for info in results:
        LOCK.acquire()
        sw.write_to_csv(info, ['animal', 'date', 'time_frame', 'probe'], csvloc)
        LOCK.release()

    LOCK2.acquire()
    with open(temppath, 'a+') as f:
        f.write(f'PROCESS: {id} ---\n')
        for info in results:
            f.write(f'\tinfo: {info}\n')
    LOCK2.release()

def init(l, l2):
    global LOCK
    LOCK = l
    global LOCK2
    LOCK2 = l2

def test_linear(paths, csvloc, temppath):
    results = []
    for p in paths:
        results.append(sw.get_info_from_path(p))

    time.sleep(2)
    for info in results:
        sw.write_to_csv(info, ['animal', 'date', 'time_frame', 'probe'], csvloc)

    with open(temppath, 'a+') as f:
        f.write(f'PROCESS: linear ---\n')
        for info in results:
            f.write(f'\tinfo: {info}\n')

if __name__ == '__main__':
    TEMP_PATH = '/media/HlabShare/clayton_sahara_work/criticality/STATUS.txt'
    l = len(sys.argv)
    animal = sys.argv[1]
    probe = sys.argv[2]
    rerun = True if sys.argv[3] == 'True' else False
    redo = True if sys.argv[4] == 'True' else False
    print(f'specifying animal -- {animal}')
    print(f'specifying probe -- {probe}')
    print(f'rerun -- {rerun}')
    print(f'redo paths -- {redo}')

    csvloc = '/media/HlabShare/clayton_sahara_work/criticality/all_results.csv'

    paths = get_paths(animal, probe, redo, csvloc)

    now = dt.now()
    with open(TEMP_PATH, 'a+') as f:
        f.write('\n\n\n------ JOB START -------- ')
        f.write(f'{now.strftime("%d/%m/%Y %H:%M:%S")} ------------ \n')
        f.write(f'{len(paths)} PATHS TO DO - of this job - lol fuck u right?\n')

    l = mp.Lock()
    l2 = mp.Lock()

    bins = np.arange(0, len(paths), 5)
    tic = time.time()
    results = []
    with mp.Pool(processes = 4, initializer=init, initargs=(l,l2)) as p:
        for i, b in enumerate(bins):
            if b == bins[-1]:
                ps = paths[b:]
            else:
                ps = paths[b:bins[i+1]]

            results.append(p.apply_async(run, args = (ps, csvloc, b, redo, rerun)))
        [r.get() for r in results]
    toc = time.time()
    print('parallel time: ', toc-tic)
