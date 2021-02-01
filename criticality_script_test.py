import sahara_work as s 

jobnum = 0
paths = glob.glob('/scratch/sensley/crit_testing/*')
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
    'saveloc': f'/scratch/sensley/crit_testing/results/',
    'animal': 'caf37',
    'probe': 'probe1',
    'time_frame':'0_12',
    'date':'01012021'
}

test_subset = paths[0:5]

s.run_testing_chpc(test_subset, params)

