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


def write_qsub_header(shfile):
    f = '''#!/bin/bash\n# Make sure ncpus in spikeinterface_currentall.py is same as ppn\n# Please change BASEDIR\nBASEDIR=/scratch/khengen_lab/crit_sahara/\
        \nOUTDIR=/scratch/khengen_lab/crit_sahara/\n# Get name for log file\nJOBID=`echo ${PBS_JOBID} | cut -c1-12`\
        \noutput_name=${PBS_JOBNAME}_${JOBID}.log\n# Load modules\nmodule purge all\nmodule load gcc-4.7.2\n\n'''
    with open(shfile, 'w') as t:
        t.write(f)



def write_to_files(o, csvloc):
    err, appended = sw.write_to_results_csv(o, csvloc)
    if err:
        print('something weird happened, this should not have errored')
    # else:
    #     new_path = o.pathname
    #     loaded = np.load('/media/HlabShare/clayton_sahara_work/criticality/loaded_paths_results.npy')
    #     loaded = np.append(loaded, new_path)
    #     np.save('/media/HlabShare/clayton_sahara_work/criticality/loaded_paths_results.npy', loaded)
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

def write_to_pkl_chpc(o, pkl_loc):
    err, appended = sw.write_to_results_pkl(o, pkl_loc)
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
    'bm':None,
    'tm':None,
    'nfactor_bm': 0,
    'nfactor_tm': 0,
    'nfactor_bm_tail': .9,  # upper bound to start exclude for burst
    'nfactor_tm_tail': .9,  # upper bound to start exclude for time 
    'cell_type': ['FS', 'RSU'],
    'quals':[1,2,3],
    'plot': True,
    'base_saveloc': f'/media/HlabShare/clayton_sahara_work/criticality/',
    'verbose':False,
    'timeout':5000,
    'none_fact':40, 
    'exclude':True, 
    'exclude_burst':50,
    'exclude_time':20,
    'exclude_diff_b':20,
    'exclude_diff_t':10,
    'fr_cutoff':50,
    'save':True,
    'start': None,
    'end': None,
    'shuffle':True,
    'subsample':False,
    'subsample_factor':None,
    'subsample_iter':None, 
    'subsample_replace':False
}

def run_testing_chpc(paths, params, JOBDIR, jobnum=0, jobname = '',animal = '', probe = '', rerun = True, redo = False):
    

    tic = time.time()
    basejobdir = JOBDIR[:JOBDIR.rfind('/')]

    errcols = ['animal', 'probe', 'date', 'time_frame', 'block_num', 'scored', 'file', 'error', 'now', 'when']
    errf = basejobdir+'/errored.pkl'


    status_file = f'{JOBDIR}/STATUS_{jobname}.txt'
    csv_file = f'{JOBDIR}/results_{jobname}.csv'
    pkl_file = f'{JOBDIR}/results_{jobname}.pkl'

    sw.write_csv_header(csv_file)

    all_objs, errors = sw.lilo_and_stitch(paths, params, save = params['save'], timeout=params['timeout'])

    results = []
    for o in all_objs:
        appended = write_to_pkl_chpc(o, pkl_file)
        appended2 = write_to_files_chpc(o, csv_file)
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
                num_passed = row[row["passed"]==True].count()['passed']
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
                f.write(f'\t{e[0]} --- {e[1]} --- {e[2]} --- {e[3]} --- {e[4]} --- {e[5]}: {e[-1]}\n')
                temp = pd.DataFrame([e], columns = errcols)
                errdf = pd.read_pickle(errf)
                errdf = errdf.append(temp)
                errdf.to_pickle(errf)
    return 0


def smol_2_big(animal):
    if len(animal)==5:
        a = animal[:3].upper() + '000' + animal[3:]
    else:
        a = animal[:3].upper() + '00' + animal[3:]
    return a

def get_all_paths(animal):
    all_paths = sorted(glob.glob(f'/scratch/khengen_lab/crit_sahara/DATA/media/HlabShare/Clustering_Data/{animal}*/*/*/*/co/*neurons_group0.npy'))
    all_paths = [p for p in all_paths if 'block' not in p]
    print(f'total num paths: {len(all_paths)}', flush=True)
    all_animals = np.unique([sw.get_info_from_path(p)[0] for p in all_paths])
    print(f'total num animals: {len(all_animals)}', flush=True)
    allpaths = []
    for animal in all_animals:
        probe = sw.get_probe(animal, region = 'CA1')
        geno = sw.get_genotype(animal)
        if geno == 'app_ps1' or (len(sw.get_regions(animal))>=4):
            a = smol_2_big(animal)
            animal_paths = sorted([p for p in all_paths if a in p])
            print(f'{animal}: {len(animal_paths)}')
            allpaths.append(animal_paths)
        elif probe != -1:
            a = smol_2_big(animal)
            animal_paths = sorted([p for p in all_paths if a in p and probe in p])
            print(f'{animal}: {len(animal_paths)}')
            allpaths.append(animal_paths)
    
    allpaths = np.concatenate(allpaths)
    return allpaths

def get_rand_subset(per_animal = 2):
    paths = []
    allpaths = get_all_paths('')
    all_animals = np.unique([sw.get_info_from_path(p)[0] for p in allpaths])
    for animal in all_animals:
        probe = sw.get_probe(animal, region = 'CA1')
        if probe != -1:
            a = smol_2_big(animal)
            animal_paths = np.sort([p for p in allpaths if a in p and probe in p])
            rand = np.random.randint(low=0, high = len(animal_paths), size=per_animal)
            ps = animal_paths[rand]
            paths.append(ps)
    paths = np.concatenate(paths)
    return paths
    

def make_chpc_crit_jobs(paths_per_job, jobname, total_jobs=None, paths = None, animal = '', resubmit = False):
    BASE = '/scratch/khengen_lab/crit_sahara/'
    print(f'base dir: ', BASE)

    if paths is None:
        paths = get_all_paths(animal)

    all_animals = np.unique([sw.get_info_from_path(p)[0] for p in paths])
    pathcount = 0
    jobcount = 0
    finalpaths = []
    for animal in all_animals:
        a = smol_2_big(animal)
        animal_paths = [p for p in paths if a in p]
        print(f'{animal}: {len(animal_paths)}')
        bins = np.arange(0, len(animal_paths), paths_per_job)
        for i, b in enumerate(bins):
            if total_jobs is not None and jobcount > total_jobs:
                print('Killing this, jobnum reached')
                break
            os.chdir(BASE)
            if i == len(bins)-1:
                these_paths = animal_paths[b:]
            else:
                these_paths = animal_paths[b:b+paths_per_job]
            if resubmit:
                newjobdir = os.path.join(BASE, 'JOBS', 'RERUN', jobname, f'{animal}_job_{i}')
            else:
                newjobdir = os.path.join(BASE, 'JOBS', jobname, f'{animal}_job_{i}')
            print('newdir: ', newjobdir)
            if not os.path.exists(newjobdir):
                os.makedirs(newjobdir)
            shutil.copyfile(BASE+f'qsub_criticality_chpc_{jobname}.sh', newjobdir+f'/qsub_criticality_chpc_{jobname}.sh')
            shutil.copyfile(BASE+f'criticality_script_{jobname}.py', newjobdir+f'/criticality_script_{jobname}.py')
            
            os.chdir(newjobdir)
            with open(f'qsub_criticality_chpc_{jobname}.sh', 'r') as f:
                shellfile = f.read()
            shellfile = shellfile.replace('REPLACEJOBNAME', f'{animal}_job_{i}_{jobname}')
            shellfile = shellfile.replace('REPLACEBASE', newjobdir)
            shellfile = shellfile.replace('REPLACEOUT', newjobdir)
            shellfile = shellfile.replace('REPLACECOUNT', str(pathcount))
            shellfile = shellfile.replace('SCRIPTNAME', f'criticality_script_{jobname}.py')

            with open(f'qsub_criticality_chpc_{jobname}.sh', 'w') as f:
                f.write(shellfile)

            with open('job_paths.txt', 'w') as pathfile:
                for p in these_paths:
                    pathfile.write(f'{p}\n')

            pathcount+=paths_per_job
            jobcount+=1
            finalpaths.append(newjobdir)
    os.chdir('/scratch/khengen_lab/crit_sahara')
    write_qsub_header('qsub_tosubmit.sh')
    with open('qsub_tosubmit.sh', 'a+') as sub:
        for f in finalpaths:
            qsub = glob.glob(f+'/qsub*')[0]
            sub.write(f'qsub {qsub}\n')
    print('qsub_tosubmit.sh written -- done')
    return finalpaths


def resubmit_jobs(efiles):

    ef = efiles[0]
    jobname = ef[ef.rfind('_')+1:ef.rfind('.e')]
    print(f'JOBNAME: {jobname} --- FIXING {len(efiles)} JOBS')
    edirs = [f'JOBS/{jobname}/{e[:e.find(jobname)-1]}' for e in efiles]

    paths_to_fix = []
    for d in edirs:
        with open(f'{d}/job_paths.txt', 'r') as f:
            for line in f:
                paths_to_fix.append(line.strip())
    print(f'TOTAL PATHS TO RERUN: {len(paths_to_fix)}')
    newdirs = make_chpc_crit_jobs(paths_per_job = 1, jobname = jobname, paths = paths_to_fix, resubmit = True)
    os.chdir('/scratch/khengen_lab/crit_sahara')
    write_qsub_header('qsub_tosubmit.sh')
    with open('qsub_tosubmit.sh', 'a+') as sub:
        for f in newdirs:
            qsub = glob.glob(f+'/qsub*')[0]
            sub.write(f'qsub {qsub}\n')
    print('qsub_tosubmit.sh written -- done')


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


def plot_dist(ax, burst, xmin, alpha, c, shuffled):
    pdf = np.histogram(burst, bins = np.arange(1, np.max(burst) + 2))[0]
    p = pdf / np.sum(pdf)
    p[p==0] = np.nan
    ax.plot(np.arange(1, np.max(burst) + 1), p, color = c, alpha = 0.75, linewidth=1)
    
    if shuffled is not None:
        pdfs = np.histogram(shuffled, bins = np.arange(1, np.max(shuffled) + 2))[0]
        ps = pdfs / np.sum(pdfs)
        ax.plot(np.arange(1, np.max(shuffled) + 1), ps, color = 'gray', alpha = 0.75, linewidth=1)
        
    x = np.arange(xmin, np.max(burst)**0.8)
    y = (np.size(np.where(burst == xmin + 6)[0]) / np.power(xmin + 6, -alpha)) *\
        np.power(x, -alpha)
    y = y / np.sum(pdf)
    ax.plot(x, y, color = 'red')
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1,10**4)

def scrub_dists(crit_objs):
    res = []
    for p in crit_objs:
        crit = saw.load_crit(p)
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[5, 5])
        plot_dist(ax, crit.burst, crit.xmin, crit.alpha, 'lightseagreen', None)
        fig.show()
        score = input('rating?: ')
        res.append([crit.animal, crit.date, crit.time_frame, crit.block_num, crit.probe, crit.burst, 'burst', score])

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[5, 5])
        plot_dist(ax, crit.T, crit.tmin, crit.beta, 'lightcoral', None)
        fig.show()
        score = input('rating?: ')
        res.append([crit.animal, crit.date, crit.time_frame, crit.block_num, crit.probe, crit.T, 'T', score])
