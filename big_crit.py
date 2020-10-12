import numpy as np 
import glob
import pandas as pd
from sahara_work import Crit
from sahara_work import lilo_and_stitch
from sahara_work import lilo_and_stitch_the_sequel
from datetime import datetime as dt
import signal 
import sys

def rerun(p):
    #passing in end of glob string to rerun specific paths
    s = f'/media/HlabShare/clayton_sahara_work/criticality/{p}'
    paths = glob.glob(s)
    print(f'# of paths to analyze: {len(paths)}')

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
    'quality': [1,3], # all qualities would be [1,2,3]
    'cell_type': ['FS', 'RSU'], 
    'plot' : True
    }

    final_paths = []
    for p in paths:
        c = np.load(p, allow_pickle = True)
        if c.p_value_burst < 0.05 or c.p_value_t < 0.05:
            final_paths.append(p)
    
    all_objs, errors = lilo_and_stitch_the_sequel(final_paths, params)

    df = pd.DataFrame(results, columns = ['animal', 'date', 'time_frame', 'block_num', 'p_value_burst', 'p_value_t', 'dcc', 'passed'])

    group = df.groupby(['animal', 'date'])

    strs = []
    for i, row in group:
        num_passed = row[row["passed"]].count()['passed']
        total_num = row.count()['passed']
        avg_dcc = row.mean()['dcc']
        animal = row['animal'].to_numpy()[0]
        date = row['date'].to_numpy()[0]
        s = f'{str(animal)} -- {date} -- passed {num_passed}/{total_num} -- avg dcc {avg_dcc}'
        strs.append(s)

    now = dt.now()
    if len(all_objs) > 0 :
        with open('/media/HlabShare/clayton_sahara_work/criticality/STATUS.txt', 'a+') as f:
            f.write(f'{now.strftime("%d/%m/%Y %H:%M:%S")} ------------ \n')
            for s in strs:
                f.write(f'{s}\n')
            f.write('\tERRORS:\n')
            for e in errors:
                f.write(f'\t{e}\n')


def run(animal='', rerun=False):
    s = f'/media/HlabShare/clayton_sahara_work/clustering/{animal}*/*/probe*/co/*scored_clayton_spks_rm_new_mbt_caf.npy'
    print(s)
    paths = glob.glob(s)
    print(f'# of paths to analyze: {len(paths)}')

    params = {
        'rerun' : False,
        'flag': 2, # 1 is DCC 2 is p_val and DCC
        'ava_binsz': 0.04, # in seconds
        'hour_bins': 4,# durration of block to look at
        'perc': 0.35,
        'nfactor_bm':0, 
        'nfactor_tm':0,
        'nfactor_bm_tail':.8, # upper bound to start exclude for burst
        'nfactor_tm_tail': .8, # upper bound to start exclude for time
        'quality': [1,2,3], # all qualities would be [1,2,3]
        'cell_type': ['FS', 'RSU'], 
        'plot' : True
        }
    bins = np.arange(0,len(paths),10)

    for i,b in enumerate(bins):
        print(f"\n\n{b} ---- PATHS COMPLETE \n\n")
        if b == bins[-1]:
            ps = paths[b:]
        else:
            ps = paths[b: bins[i+1]]

        all_objs, errors = lilo_and_stitch(ps, params, rerun=rerun)

        results = []
        for o in all_objs:
            results.append([o.animal, o.probe, o.date, o.time_frame, o.block_num, o.p_value_burst, o.p_value_t, o.dcc, (o.p_value_burst > 0.05 and o.p_value_t > 0.05)])

        df = pd.DataFrame(results, columns = ['animal', 'probe','date', 'time_frame', 'block_num', 'p_value_burst', 'p_value_t', 'dcc', 'passed'])

        group = df.groupby(['animal','probe', 'date'])

        strs = []
        for i, row in group:
            num_passed = row[row["passed"]].count()['passed']
            total_num = row.count()['passed']
            avg_dcc = row.mean()['dcc']
            animal = row['animal'].to_numpy()[0]
            date = row['date'].to_numpy()[0]
            probe = row['probe'].to_numpy()[0]
            s = f'{str(animal)} -- {probe} -- {date} -- passed {num_passed}/{total_num} -- avg dcc {avg_dcc}'
            strs.append(s)

        now = dt.now()
        if len(all_objs) > 0 :
            with open('/media/HlabShare/clayton_sahara_work/criticality/STATUS.txt', 'a+') as f:
                f.write(f'{now.strftime("%d/%m/%Y %H:%M:%S")} ------------ \n')
                for s in strs:
                    f.write(f'{s}\n')
                f.write('\tERRORS:\n')
                for e in errors:
                    f.write(f'\t{e}\n')

if __name__ == "__main__":
    l = len(sys.argv)
    if l > 1:
        animal = sys.argv[1]
        rerun = True if sys.argv[2] == 'True' else False
        print(f'specifying animal -- {animal}')
        print(f'rerun -- {rerun}')
        run(sys.argv[1], sys.argv[2])
    else:
        run()


