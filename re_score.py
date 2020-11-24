import numpy as np
import musclebeachtools_hlab.musclebeachtools as mbt 
import glob 
import collections
from matplotlib import pyplot as plt 
import os
import os.path as op
import time

s = f'/media/HlabShare/clayton_sahara_work/clustering/*/*/*/co/*scored_clayton_spks_rm_new_mbt_caf.npy'
print(s)
fs = [f for f in glob.glob(s)]
print(f'total # of paths: {len(fs)}', flush=True)

for i, p in enumerate(fs):
    new_name = p[0:p.find('neurons')] + 'neurons_rescored_for_xgb.npy'
    if op.exists(new_name):
        print("Already did this one")
    else:
        tic = time.time()
        print(p)
        print('Loading . . . . . . . . . . ')
        cells = np.load(p,allow_pickle=True)
        twos = len([cell for cell in cells if cell.quality==2])
        print(f'{twos} TWOS TO RESCORE')
        for cell in cells:
            if cell.quality == 2:
                fix_amp = cell.mean_amplitude*2.2
                cell.checkqual(fix_amp_ylim=fix_amp)
        twos_1 = len([cell for cell in cells if cell.quality==2])
        print(f'{twos} start --- {twos_1} fixed')
        np.save(new_name, cells)
        toc = time.time()
        print('########## It took you {} minutes to RESCORE this block ##########'.format((toc-tic)/60))
        print(f'\n\n Only {len(fs)-i} paths left lol')
