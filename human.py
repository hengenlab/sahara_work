from pickle import load
import numpy as np
import scipy
import scipy.io as sio
import h5py
import glob
from copy import deepcopy as cdc
import criticality as cr

def load_human(filename):
    f = h5py.File(filename)
    for k in f.keys():
        n = k
    data = f[n].value
    return data

files = glob.glob('MSC/msc*.mat')
mask = load_human('MSC/submask.mat')
networks = mask[0]
idxs = np.where(networks == 4)[0]

dfull = []
for f in files:
    data = load_human(f)
    d = data[:,idxs]
    dfull.append(d)
DATA = np.concatenate(dfull)
D = DATA.T  # SO! The current state of the data is that it's only the voxels from network 4 (who knows what network that is), and its VOXEL x TIME (so like spiketimes)

D_bin = []
for v in D:
    std = np.std(v) * 2
    temp = cdc(v) 
    temp[v >= std] = 1 
    temp[v < std] = 0 
    D_bin.append(temp)   

R = cr.get_avalanches(D_bin, perc = 0.25) 
burst = R['S'] 
T = R['T'] 

Result = cr.AV_analysis(burst, T, flag = 1, nfactor_bm_tail = 0.75, nfactor_tm_tail = 0.75, nfactor_bm = 50, bm = 200,  tm = 50,\
                        pltname = 'msc1_test', saveloc = '/media/HlabShare/AD_paper/human/', plot= True)

