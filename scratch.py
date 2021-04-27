from sahara_work import Crit
import sahara_work as saw
import numpy as np 
import pandas as pd
import signal
import glob
from copy import deepcopy as cdc
import pickle
import musclebeachtools as mbt

files = ['/media/HlabShare/Clustering_Data/CAF00022/caf22_05142020/24_32/probe2/co/H_2020-05-15_14-32-55_2020-05-15_22-27-55_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00026/caf26_08262020/0_8/probe2/co/H_2020-08-26_10-02-34_2020-08-26_17-57-34_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00034/caf34_08042020/0_8/probe3/co/H_2020-08-04_17-03-44_2020-08-05_00-58-44_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00037/caf37_08312020/16_24/probe1/co/H_2020-09-01_05-39-56_2020-09-01_13-34-56_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00040/caf40_09112020/0_8/probe1/co/H_2020-09-11_17-49-53_2020-09-12_01-44-54_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00042/caf42_09142020/72_84/probe3/co/H_2020-09-17_17-20-37_2020-09-18_05-15-37_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00048/caf48_11062020/0_12/probe1/co/H_2020-11-06_10-06-14_2020-11-06_22-01-15_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00049/caf49_01292021/8_12/probe1/co/H_2021-01-29_23-31-43_2021-01-30_03-26-43_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00050/caf50_01252021/64_68/probe1/co/H_2021-01-28_08-40-07_2021-01-28_09-40-07_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00052/caf52_11202020/24_36/probe1/co/H_2020-11-21_10-04-45_2020-11-21_21-59-45_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00060/caf60_12022020/8_16/probe1/co/H_2020-12-02_18-44-49_2020-12-03_02-34-49_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00061/caf61_12022020/16_24/probe1/co/H_2020-12-03_02-44-06_2020-12-03_10-34-06_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00062/caf62_11202020/12_24/probe1/co/H_2020-11-20_22-06-20_2020-11-21_09-56-20_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00069/caf69_12282020/36_48/probe4/co/H_2020-12-29_23-09-40_2020-12-30_10-59-40_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00072/02192021/96_100/probe1/co/H_2021-02-23_17-14-18_2021-02-23_21-09-18_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00077/02112021/0_12/probe1/co/H_2021-02-11_15-08-25_2021-02-12_03-03-26_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00078/02242021/108_112/probe1/co/H_2021-03-01_04-12-54_2021-03-01_08-07-54_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00079/02112021/40_44/probe1/co/H_2021-02-13_07-03-45_2021-02-13_10-58-45_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00080/02192021/72_76/probe1/co/H_2021-02-22_17-10-31_2021-02-22_21-00-31_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00081/02242021/12_24/probe3/co/H_2021-02-25_04-17-11_2021-02-25_16-12-11_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/CAF00082/02192021/60_72/probe1/co/H_2021-02-22_05-11-45_2021-02-22_17-01-45_neurons_group0.npy',
 '/media/HlabShare/Clustering_Data/EAB00047/eab47_07022019/64_72/probe2/co/H_2019-07-05_08-04-24_2019-07-05_13-59-28_neurons_group0.npy']

paths = files


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
for f in paths:   
    for b in binz: 
        print(f'binsize: {b}') 
        for q in quals: 
            print(f'quals: {q}') 
            R = {}  
            R['binsize'] = b
            R['quals'] = q
            params['ava_binsz']=b 
            params['quals']=q 
            animal, date, time_frame, probe = saw.get_info_from_path(f)
            results = saw.lilo_and_stitch([f], params, save = params['save']) 
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
                R['animal'] = animal
            
            big[str(i)] = R
            i+=1
f = open("R_all.pkl","wb") 
pickle.dump(big,f) 
f.close()   



# concat files
import glob
import pandas as pd
files = glob.glob('*.csv')
df = pd.read_csv(files[0])
for f in files[1:]:
    t = pd.read_csv(f)
    df = df.append(t)