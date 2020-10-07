import numpy as np 
import glob
import pandas as pd
from sahara_work import Crit
from sahara_work import lilo_and_stitch
from datetime import date as dt
import signal 

paths = glob.glob('/media/HlabShare/clayton_sahara_work/clustering/*/*/probe*/co/*scored_clayton_spks_rm_new_mbt_caf.npy')

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

for i in bins:
    print(f"\n\n{i} ---- PATHS COMPLETE \n\n")
    if i == bins[-1]:
        ps = paths[i:]
    else:
        ps = paths[i:(i+1)]

    all_objs, errors = lilo_and_stitch(ps, params)

    results = []
    for o in all_objs:
        results.append([o.animal, o.date, o.time_frame, o.block_num, o.p_value_burst, o.p_value_t, o.dcc, (o.p_value_burst > 0.05 and o.p_value_t > 0.05)])

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

    with open('/media/HlabShare/clayton_sahara_work/criticality/STATUS.txt', 'a+') as f:
        f.write(f'{dt.today()} ------------ \n')
        for s in strs:
            f.write(f'{s}\n')
        f.write('\tERRORS:\n')
        for e in errors:
            f.write(f'\t{e}\n')


