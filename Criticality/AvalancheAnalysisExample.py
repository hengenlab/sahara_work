from musclebeachtools_hlab import musclebeachtools as mbt
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import Criticality as cr

# Load data from mbt output; transform spiketimes into binary spikewords

goodCells = np.load('EAB50.npy',allow_pickle=True) 
ava_binsz = 0.02    # avalanches bin size 20ms  suggest to be 20-50ms
nrn_time = 12     # 12h recording block

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
burst = r['S'] #  AVsize
duration = r['T'] # AVdura

################## Avalanche analysis including AVsize, AVduration
# distribution and scaling relation. burstM and tM are used to set the
# limits for lower boundary for AVsize and AVduration distributions.
# Result1 only returns exponents, lower/upper bounds limits and DCC value
# Result2 return pvalue for null hypothesis
# Result3 generate three figures for scaling relationship

burstM = 12  # suggest to be 8-30; adjust based on figures from Result3
tM = 4       # suggest to be 3-12; adjust based on figures from Result3

Result1 = cr.AV_analysis_ExponentErrorComments(burst, duration, burstM, tM)

Result2 = cr.AV_analysis_ExponentErrorComments(burst, duration, burstM, tM, flag = 2)

Result3, ax3 = cr.AV_analysis_ExponentErrorComments(burst, duration, burstM, tM, flag = 3)

# Generate figures with shuffle data

# calcuate the size and duration of AVs from shuffled data 
burst_shuffle, duration_shuffle = cr.Genshuffle(goodCells,nrn_time,ava_binsz,perc, binary = 1, frame =0) 
# binary = 1, binarize the data; binary = 0, count the actual number of spikes in each bin
# frame = 0, generate random spiketimes (uniform distribution, keep the same FR for each neuron)
# frame = 1, do frame shuffling: for each neuron, swap all spikes in one time bin (30s window) with spikes at another randomly chosen time bin. 

cr.Plotshuffle(burstM,tM,burst,duration,burst_shuffle,duration_shuffle) 
# burst and duration is the size and duration of AVs from your experimental data, which you generated above.
