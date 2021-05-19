import glob
import sahara_work as saw
import numpy as np
import pandas as pd

params = {
    'flag': 2,  # 1 is DCC 2 is p_val and DCC
    'ava_binsz': 0.02,  # in seconds
    'hour_bins': 2,  # durration of block to look at
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
    'base_saveloc': f'/media/HlabShare/AD_paper/clustering_test/results/',
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
    'branching_ratio': False,
    'br_binsize': 0.004,
    'br_kmax': 500,
    'br_binlen': 5, # in minutes
    'br_numbins': 3 # begining middle end - this isn't implemented right now, it's just in the middle
}

hours = [2, 4, 8, 12]
cols = np.concatenate([saw.get_cols(), ['clust_durr']])
pkl_loc = '/media/HlabShare/AD_paper/clustering_test/clust_dcc_results.pkl'
for h in hours[1:]:
    files = glob.glob(f'/media/HlabShare/AD_paper/clustering_test/*/*_{h}h_*/*/*/co/*neurons_group0.npy')
    print(f'\n\nStarting {h} hour tests ---- {len(files)} files to process')
    all_objs, errs = saw.lilo_and_stitch(files, params)
    for c in all_objs:
        err, info = saw.lil_helper_boi(c)
        info = np.concatenate([info, [h]])
        saw.write_to_pkl(info, cols, pkl_loc)