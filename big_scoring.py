import numpy as np
import musclebeachtools as mbt 
import glob 
import collections
from matplotlib import pyplot as plt 
import os
import os.path as op
import time


# scorer = 'xgb'
#scorer = 'clayton'
scorer = 'sahara'
autoqual = '/media/HlabShare/'
###     M:/caf\caf42_0914\60_72\probe5\co\H_2020-09-17_05-20-37_2020-09-17_17-15-37_neurons_group0.npy     ###
# fls = sorted(glob.glob('/media/*/caf/*/*/*/co/*_neurons_group0.npy'))
fls = sorted(glob.glob('/media/bs007s/caf/*/*/*/co/*_neurons_group0.npy'))
# fls = sorted(glob.glob('C:/Users/clayt/Desktop/clustering_files/*/*_neurons_group0.npy'))
# fls = sorted(glob.glob('M:/xyf/CAF40/914/*/*/*_neurons_group0.npy'))
i = len(fls)

fls2 = []
for fl in fls:
    if os.path.exists(str(fl[:-4]+'_scored_clayton.npy')) == True:
        pass
    elif os.path.exists(str(fl[:-4]+'_scored_xgb.npy')) == True:
        pass
    elif os.path.exists(str(fl[:-4]+'_scored_sahara.npy')) == True:
        pass
    else:
        fls2.append(fl)

ii = len(fls2)

print('###########You have scored {} outputs so far!###########'.format(i-ii))
print('###########You have {} files to cluster, have fun!###########'.format(ii))

print('scorer =',scorer)
if scorer == 'clayton':
    for fl2 in fls2:
        tic = time.time()
        print(fl2)
        print('Loading . . . . . . . . . . ')
        cells = np.load(fl2,allow_pickle=True)
        print('Boosting . . . . . . . . . . ')
        mbt.autoqual(cells,autoqual)
        o_ones = len([cell for cell in cells if cell.quality == 1])
        o_twos = len([cell for cell in cells if cell.quality == 2])
        o_threes = len([cell for cell in cells if cell.quality == 3])
        print('Boosted . . . . . . . . . . ')
        new_bad_list = []
        for cell in cells:
            if cell.quality in [1,2,3]:
                if cell.mean_amplitude > 900:
                    cell.quality = 4
                else:
                    fix_amp = cell.mean_amplitude*2.2
                    cell.checkqual(fix_amp_ylim=fix_amp)
                if cell.quality == 4:
                    new_bad_list.append(cell)
        if len(new_bad_list) > 0:
            print('Saving bad list . . . . . . . . . . ')
            np.save(str(fl2[:-4]+'_new_bad_list.npy'),cells)
        toc = time.time()
        print('########## It took you {} minutes to score this block ##########'.format((toc-tic)/60))
        print('delta ones = {} pre / {} post'.format(o_ones,len([cell for cell in cells if cell.quality in [1]])))
        print('delta twos = {} pre / {} post'.format(o_twos,len([cell for cell in cells if cell.quality in [2]])))
        print('delta threes = {} pre / {} post'.format(o_threes,len([cell for cell in cells if cell.quality in [3]])))
        print('{} good cells'.format(len([cell for cell in cells if cell.quality in [1,2,3]])))
        print('Saving . . . . . . . . . . ')
        np.save(str(fl2[:-4]+'_scored_clayton.npy'),cells)
        ii = ii - 1
        print('########## {} neuron outputs left ##########'.format(ii))

if scorer == 'sahara':
    for fl2 in fls2:
        tic = time.time()
        print(fl2)
        print('Loading . . . . . . . . . . ')
        cells = np.load(fl2,allow_pickle=True)
        print('Boosting . . . . . . . . . . ')
        mbt.autoqual(cells,autoqual)
        o_ones = len([cell for cell in cells if cell.quality == 1])
        o_twos = len([cell for cell in cells if cell.quality == 2])
        o_threes = len([cell for cell in cells if cell.quality == 3])
        print('Boosted . . . . . . . . . . ')
        new_bad_list = []
        for cell in cells:
            if cell.quality in [1,2,3]:
                if cell.mean_amplitude > 900:
                    cell.quality = 4
                else:
                    fix_amp = cell.mean_amplitude*2.2
                    cell.checkqual(fix_amp_ylim=fix_amp)
                if cell.quality == 4:
                    new_bad_list.append(cell)
        if len(new_bad_list) > 0:
            print('Saving bad list . . . . . . . . . . ')
            np.save(str(fl2[:-4]+'_new_bad_list.npy'),cells)
        toc = time.time()
        print('########## It took you {} minutes to score this block ##########'.format((toc-tic)/60))
        print('delta ones = {} pre / {} post'.format(o_ones,len([cell for cell in cells if cell.quality in [1]])))
        print('delta twos = {} pre / {} post'.format(o_twos,len([cell for cell in cells if cell.quality in [2]])))
        print('delta threes = {} pre / {} post'.format(o_threes,len([cell for cell in cells if cell.quality in [3]])))
        print('{} good cells'.format(len([cell for cell in cells if cell.quality in [1,2,3]])))
        print('Saving . . . . . . . . . . ')
        np.save(str(fl2[:-4]+'_scored_sahara.npy'),cells)
        ii = ii - 1
        print('########## {} neuron outputs left ##########'.format(ii))

elif scorer == 'xgb':
    for fl2 in fls2:
        tic = time.time()
        cells = np.load(fl2,allow_pickle=True)
        mbt.autoqual(cells,'Z:/models/xgboost_autoqual')
        o_ones = len([cell for cell in cells if cell.quality == 1])
        o_twos = len([cell for cell in cells if cell.quality == 2])
        o_threes = len([cell for cell in cells if cell.quality == 3])
        for cell in cells:
            if cell.quality in [1,2,3]:
                if cell.mean_amplitude < 45:
                    cell.quality = 4
                else:
                    fix_amp = cell.mean_amplitude*2.2
        toc = time.time()
        np.save(str(fl2[:-4]+'_scored_xgb.npy'),cells)
        ii = ii - 1