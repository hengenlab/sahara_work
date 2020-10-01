import numpy as np
import sahara_work as sw 

params = {
    'flag': 2, # 1 is DCC 2 is p_val and DCC
    'ava_binsz': 0.040, # in seconds
    'hour_bins': 4, # durration of block to look at
    'perc': 0.25,
    'nfactor_bm_tail':.8, # upper bound to start exclude for burst
    'nfactor_tm_tail':.8, # upper bound to start exclude for time
    'quality': [1,2,3],
    'cell_type': ['RSU', 'FS'], # all qualities would be [1,2,3]
    'animal' : 'CAF37',
    'saveloc' : "/media/HlabShare/clayton_sahara_work/criticality/caf37/0827/24_36/",
    'notes': 'no FS cells',
    'time_frame':'24_36', # if you're running multiple blocks and the paths are right - this should be None
    'plot' : True
    }

paths = ['/media/bs007s/caf/caf37/0827_12htest/24_36/probe1/co/scored_clayton.npy']

data = sw.lilo_and_stitch(paths, params)
