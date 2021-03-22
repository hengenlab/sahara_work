import sahara_work as s 
import glob
import sys
import time

def run(basedir, pathnum, jobname):
    paths = []
    with open('job_paths.txt', 'r') as f:
        for line in f:
            paths.append(line.strip())

    params = {
        'flag': 2,  # 1 is DCC 2 is p_val and DCC
        'ava_binsz': 0.02,  # in seconds
        'hour_bins': 4,  # durration of block to look at
        'perc': 0.35,
        'bm':50,
        'tm':20,
        'nfactor_bm': 0,
        'nfactor_tm': 0,
        'nfactor_bm_tail': .8,  # upper bound to start exclude for burst
        'nfactor_tm_tail': .8,  # upper bound to start exclude for time 
        'cell_type': ['FS', 'RSU'],
        'quals':[1,2],
        'plot': True,
        'base_saveloc': f'/scratch/khengen_lab/crit_sahara/RESULTS/----addjobname-----',
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
        'subsample_replace':False,
        'branching_ratio': True,
        'br_binsize': 0.004,
        'br_kmax': 500,
        'br_binlen': 5, # in minutes
        'br_numbins': 3 # begining middle end
    }

    tic = time.time()
    s.run_testing_chpc(paths, params, JOBDIR = basedir, jobnum = pathnum, jobname = jobname)
    toc = time.time()
    print(f'\nTOTAL JOB TIME: {(toc-tic)/60} min')

if __name__ == '__main__':
    basedir = sys.argv[1]
    pathnum = sys.argv[2]
    jobname = sys.argv[3]
    run(basedir, pathnum, jobname)

