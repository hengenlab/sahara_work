from sahara_work import Crit
import sahara_work as saw
import numpy as np 
import pandas as pd
import signal
import glob
from copy import deepcopy as cdc
import pickle
import musclebeachtools as mbt

files = glob.glob('/media/HlabShare/Clustering_Data/EAB00047/*/*/*/co/*neurons_group0.npy') 
paths = [files[10]] 


params = { 
    'flag': 2,  # 1 is DCC 2 is p_val and DCC 
    'ava_binsz': None,  # in seconds 
    'hour_bins': 4,  # durration of block to look at 
    'perc': 0.35, 
    'bm':50, 
    'tm':20, 
    'nfactor_bm': 0, 
    'nfactor_tm': 0, 
    'nfactor_bm_tail': .9,  # upper bound to start exclude for burst 
    'nfactor_tm_tail': .9,  # upper bound to start exclude for time  
    'cell_type': ['FS', 'RSU'], 
    'quals':None, 
    'plot': True,  
    'base_saveloc': f'/media/HlabShare/clayton_sahara_work/criticality/testing/paramchecks/', 
    'verbose':False, 
    'timeout':False, 
    'none_fact':40,  
    'exclude':True,  
    'exclude_burst':50, 
    'exclude_time':20, 
    'exclude_diff_b':20, 
    'exclude_diff_t':10, 
    'fr_cutoff':50, 
    'save':False, 
    'start': 1, 
    'end': 2, 
    'shuffle': True,
    'subsample':False,
    'subsample_factor':None,
    'subsample_iter':None, 
    'subsample_replace':False 
} 

big = {}


    
binz = [0.004, 0.01, 0.02, 0.04]
quals = [[1], [1, 2], [1, 2, 3]] 
i = 0
for f in allpaths:   
    for b in binz: 
        print(f'binsize: {b}') 
        for q in quals: 
            print(f'quals: {q}') 
            R = {}  
            R['binsize'] = b
            R['quals'] = q
            params['ava_binsz']=b 
            params['quals']=q 
                
            results = saw.lilo_and_stitch([f], params, save = params['save'], verbose = params['verbose']) 
            if len(results[0]) > 0:
                crit = results[0][0]
                  
                R['burst'] = crit.burst 
                R['xmin'] = crit.xmin 
                R['xmax'] = crit.xmax 
                R['alpha'] = crit.alpha
                R['T'] = crit.T 
                R['tmin'] = crit.tmin 
                R['tmax'] = crit.tmax 
                R['beta'] = crit.beta
                R['dcc'] = crit.dcc
                R['ncells'] = crit.num_cells
                R['pvalb'] = crit.p_value_burst
                R['pvalt'] = crit.p_value_t
                R['burstS'] = crit.burstS
                R['TS'] = crit.TS
                R['animal'] = crit.animal
                R['failed'] = False
            else:
                R['failed'] = True
            
            big[str(i)] = R
            i+=1
f = open("R_all.pkl","wb") 
pickle.dump(big,f) 
f.close()   