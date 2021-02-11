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
import shutil
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

def write_to_files_chpc(o, csvloc):
    err, appended = sw.write_to_results_csv(o, csvloc)
    if err:
        print('something weird happened, this should not have errored')
    else:
        new_path = o.pathname
        loaded = np.load('/scratch/khengen_lab/crit_sahara/loaded_paths_results.npy')
        loaded = np.append(loaded, new_path)
        np.save('/scratch/khengen_lab/crit_sahara/loaded_paths_results.npy', loaded)
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
        'plot': True,
        'quals': None, 
        'base_saveloc': f'/media/HlabShare/clayton_sahara_work/criticality/',
        'verbose':False,
        'timeout':5000,
        'none_fact':40, 
        'exclude':True, 
        'exclude_burst':50,
        'exclude_time':20,
        'exclude_diff_b':20,
        'exclude_diff_t':10,
        'save': True
    }

def run_testing_chpc(paths, params, JOBDIR, jobnum=0, jobname = '',animal = '', probe = '', rerun = True, redo = False):
    tic = time.time()
    status_file = f'{JOBDIR}/STATUS_{jobname}.txt'
    csv_file = f'{JOBDIR}/results_{jobname}.csv'

    all_objs, errors = sw.lilo_and_stitch(paths, params, rerun = rerun, save = params['save'], verbose=params['verbose'], timeout=params['timeout'])

    results = []
    for o in all_objs:
        appended = write_to_files_chpc(o, csv_file)
        results.append(appended)

    if len(all_objs) > 0:
        cols = sw.get_cols()
        df = pd.DataFrame(results, columns = cols)
        
        group = df.groupby(['animal', 'probe', 'date', 'scored'])
        strs = []
        for i, row in group:
            if params['flag'] == 1:
                num_passed = 0
            else:
                num_passed = row[row["passed"]].count()['passed']
            total_num = row.count()['dcc']
            avg_dcc = row.mean()['dcc']
            animal = row['animal'].to_numpy()[0]
            date = row['date'].to_numpy()[0]
            probe = row['probe'].to_numpy()[0]
            scored = row['scored'].to_numpy()[0]
            age = row['age'].astype(str).to_numpy()[0] # check this line
            s = f'{str(animal)} -- {probe} -- {date} -- {scored} -- {age}-- passed {num_passed}/{total_num} -- avg dcc {avg_dcc}'
            strs.append(s)
    toc = time.time()
    now = dt.now()
    with open(status_file, 'a+') as f:
        f.write(f'\n{now.strftime("%d/%m/%Y %H:%M:%S")} ------------ \n')
        f.write(f'{jobnum} PATHS DONE - of this job\n')
        f.write(f'{(toc-tic)/60/60} hours to complete these paths\n')
        if len(all_objs) > 0: 
            for s in strs:
                f.write(f'{s}\n')
        if len(errors) > 0:
            f.write('\tERRORS:\n')
            for e in errors:
                f.write(f'\t{e[0]}\n')
                errored = np.load('/scratch/khengen_lab/crit_sahara/errored_paths.npy')
                errored = np.append(errored, e[1])
                np.save('/scratch/khengen_lab/crit_sahara/errored_paths.npy', errored)
    return 0

def make_chpc_crit_jobs(paths_per_job, total_jobs=None):
    BASE = '/scratch/khengen_lab/crit_sahara/'
    print(f'base dir: ', BASE)
    all_paths = sorted(glob.glob('/scratch/khengen_lab/crit_sahara/DATA/media/HlabShare/clayton_sahara_work/clustering/*/*/*/*/co/*neurons_group0.npy'))
    print(f'total num paths: {len(all_paths)}', flush=True)
    all_animals = np.unique([sw.get_info_from_path(p)[0] for p in all_paths])
    print(f'total num animals: {len(all_animals)}', flush=True)
    pathcount = 0
    jobcount = 0
    for animal in all_animals:
        probe = sw.get_probe(animal, region = 'CA1')
        animal_paths = sorted(glob.glob(f'/scratch/khengen_lab/crit_sahara/DATA/media/HlabShare/clayton_sahara_work/clustering/{animal}*/*/*/{probe}/co/*neurons_group0.npy'))
        bins = np.arange(0, len(animal_paths), paths_per_job)
        for i, b in enumerate(bins):
            if total_jobs is not None and jobcount > total_jobs:
                print('Killing this, jobnum reached')
                return
            os.chdir(BASE)
            if i == len(bins)-1:
                these_paths = animal_paths[b:]
            else:
                these_paths = animal_paths[b:b+paths_per_job]
            newjobdir = os.path.join(BASE, 'JOBS', f'{animal}_job_{i}')
            print('newdir: ', newjobdir)
            if not os.path.exists(newjobdir):
                os.makedirs(newjobdir)
            shutil.copyfile(BASE+'qsub_criticality_chpc.sh', newjobdir+'/qsub_criticality_chpc.sh')
            shutil.copyfile(BASE+'criticality_script_test.py', newjobdir+'/criticality_script_test.py')
            
            os.chdir(newjobdir)
            with open('qsub_criticality_chpc.sh', 'r') as f:
                shellfile = f.read()
            shellfile = shellfile.replace('REPLACEJOBNAME', f'{animal}_job_{i}')
            shellfile = shellfile.replace('REPLACEBASE', newjobdir)
            shellfile = shellfile.replace('REPLACEOUT', newjobdir)
            shellfile = shellfile.replace('REPLACECOUNT', str(pathcount))

            with open('qsub_criticality_chpc.sh', 'w') as f:
                f.write(shellfile)

            with open('job_paths.txt', 'w') as pathfile:
                for p in these_paths:
                    pathfile.write(f'{p}\n')

            pathcount+=paths_per_job
            jobcount+=1

def run_linear(paths, params, jobnum, animal = '', probe = '', rerun = True, redo = False):
    paths = sw.get_paths(animal = animal, probe = probe)
    all_objs, errors = lilo_and_stitch(paths, params, rerun = rerun, save = True, verbose=False)
    results = []
    for o in all_objs:
        appended = write_to_files(o, csvloc)
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
    with open(f'/media/HlabShare/clayton_sahara_work/criticality/STATUS_{jobnum}.txt', 'a+') as f:
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

    return 0