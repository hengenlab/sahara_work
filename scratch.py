from sahara_work import Crit
import sahara_work as saw
import numpy as np 
import pandas as pd
import signal
import glob
import pickle

files = glob.glob('/media/HlabShare/Clustering_Data/CAF00037/*/*/*/co/*neurons_group0.npy') 
paths = [files[10]] 
 
params = { 
    'flag': 1,  # 1 is DCC 2 is p_val and DCC 
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
    'timeout':5000, 
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
    'shuffle': True 
} 



R = {
    '0.01' : { 
            '[1]' : { 
               
            }, 
            '[1, 2]' : {  
                    
            }, 
            '[1, 2, 3]' : {  
                    
            },   
    }, 
        '0.02' : { 
            '[1]' : { 
               
            }, 
            '[1, 2]' : {  
                    
            }, 
            '[1, 2, 3]' : {  
                    
            },   
    }, 
        '0.04' : {  
                '[1]' : {  
    
                },  
                '[1, 2]' : {   
                       
                },  
                '[1, 2, 3]' : {   
                        
                },    
            }  
} 
    
binz = [0.01, 0.02, 0.04]
quals = [[1], [1, 2], [1, 2, 3]]    
for b in binz: 
    print(f'binsize: {b}') 
    for q in quals: 
        print(f'quals: {q}') 
        
        params['ava_binsz']=b 
        params['quals']=q 
             
        results = saw.lilo_and_stitch(paths, params, save = params['save'], verbose = params['verbose']) 
        if len(results[0]) > 0:
            crit = results[0][0]     
            R[str(b)][str(q)]['burst'] = crit.burst 
            R[str(b)][str(q)]['xmin'] = crit.xmin 
            R[str(b)][str(q)]['xmax'] = crit.xmax 
            R[str(b)][str(q)]['alpha'] = crit.alpha
            R[str(b)][str(q)]['T'] = crit.T 
            R[str(b)][str(q)]['tmin'] = crit.tmin 
            R[str(b)][str(q)]['tmax'] = crit.tmax 
            R[str(b)][str(q)]['beta'] = crit.beta
            R[str(b)][str(q)]['dcc'] = crit.dcc
            R[str(b)][str(q)]['ncells'] = crit.num_cells
            R[str(b)][str(q)]['burstS'] = crit.burstS
            R[str(b)][str(q)]['TS'] = crit.TS
            R[str(b)][str(q)]['failed'] = False
        else:
            R[str(b)][str(q)]['failed'] = True
f = open("R_caf37.pkl","wb") 
pickle.dump(R,f) 
f.close()   