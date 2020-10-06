import numpy as np 
import glob
import pandas as pd
from sahara_work import Crit
from sahara_work import lilo_and_stitch
from datetime import date
import signal 

paths = ['/media/HlabShare/clayton_sahara_work/clustering/caf46_1001_1/80_88/probe1/co/H_2020-10-04_18-49-12_2020-10-05_02-44-12_neurons_group0_scored_clayton_spks_rm_new_mbt_caf.npy']

params = {
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
    'saveloc' : "/media/HlabShare/clayton_sahara_work/criticality/caf46/1001/",
    'plot' : True
    }

all_objs, errors = lilo_and_stitch(paths, params)

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

with open('STATUS.txt', 'a+') as f:
    f.write(f'{date.today()} ------------ \n')
    for s in strs:
        f.write(f'{s}\n')
    f.write('\tERRORS:\n')
    for e in errors:
        f.write(f'\t{e}\n')


