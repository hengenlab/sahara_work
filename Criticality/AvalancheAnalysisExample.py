from musclebeachtools_hlab import musclebeachtools as mbt
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sahara_work import Criticality as cr

# Load data from mbt output; transform spiketimes into binary spikewords

goodCells = np.load('EAB50.npy',allow_pickle=True) 
ava_binsz = 0.045    # avalanches bin size 20ms  suggest to be 20-50ms
nrn_time = 16     # 12h recording block

# spks = mbt.getspikes(goodCells, 0, 3600*nrn_time)
# data_T = mbt.spiketimes_to_spikewords(spks,0,3600*nrn_time,ava_binsz,1)  # old version musclebeachtools
# spks = mbt.n_getspikes(goodCells, 0, 3600*nrn_time)   # new musclebeachtools and spk_interface output

# The last parameter: 1 means binarize the spike numbers into 1 or 0; 0 means calculate the actual number of spikes in each bin
data = mbt.n_spiketimes_to_spikewords(goodCells,ava_binsz,0,3600*nrn_time,1)   # new musclebeachtools and spk_interface output


 ############### get AVsize and AVduration ################ 
# the default  threshold is 25 percentile. But we could change it by tuning
# 'perc'. Make  sure 'perc' is in the range from 0 -- 1. When network
# activity is silent most time, should set 'perc' = 0, which means threshold
# is zero #########  

perc = 0.3   # threshold of neural avalanches   suggest to be 20-50%
r = cr.AV_analysis_BurstT(data, perc = perc)
x = r['S'] # x is AVsize
y = r['T'] # y is AVdura

################## Avalanche analysis including AVsize, AVduration
# distribution and scaling relation. burstM and tM are used to set the
# limits for lower boundary for AVsize and AVduration distributions.
# Result1 only returns exponents, lower/upper bounds limits and DCC value
# Result2 return pvalue for null hypothesis
# Result3 generate three figures for scaling relationship

burstM = 18  # suggest to be 8-20; adjust based on figures from Result3
tM = 4       # suggest to be 3-6; adjust based on figures from Result3

Result1 = cr.AV_analysis_ExponentErrorComments(x, y, burstM, tM)

Result2 = cr.AV_analysis_ExponentErrorComments(x, y, burstM, tM, flag = 2)

Result3, ax3 = cr.AV_analysis_ExponentErrorComments(x, y, burstM, tM, flag = 3)

