import sahara_work as s 
import glob
import sys

def run(basedir, pathnum, jobname):
    paths = []
    with open('job_paths.txt', 'r') as f:
        for line in f:
            paths.append(line.strip())

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

    s.run_testing_chpc(paths, params, JOBDIR = basedir, jobnum = pathnum, jobname = jobname)

if __name__ == '__main__':
    basedir = sys.argv[1]
    pathnum = sys.argv[2]
    jobname = sys.argv[3]
    run(basedir, pathnum, jobname)

