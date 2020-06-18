'''1. Load first file in block in ntk, bandpass it.
2. Load neuron list in mbt, all you need is len of neuron_list spike_times and quality
3. Do two for loops
loop 1. number of tetrodes  Plot tetrode 4 channels
        loop 2. Number of neuron and if neuron belong to that tetrode
                Plot spike time in sample above the data (copy the code to plot spikes as Keith wrote) based on quality'''



import musclebeachtools as mbt
import neuraltoolkit as ntk
import numpy as np
import os
import glob
import time
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as clr
import copy
plt.ion()
# datdir = '/media/bs005r/braingineers/2019-08-03/10uMGlutamate/t/final/'
# datdir = '/media/bs005r/braingineers/2019-08-03/10uMGlutamate/t_l5/final/'
# datdir = '/media/bs005r/braingineers/2019-08-03/10uMGlutamate/t_2/final/'
# datdir = '/media/bs005r/braingineers/2019-08-03/10uMGlutamate/t33/final/'
datdir = '/media/bs006r/CAF00022/CAF00022_2020-06-05_08-47-46/'
# rawdat = 'P_20190803-Gw17-DIV17-Glu10uM-002-190803-173035-190803-175810_190803_183526_-P_20190803-Gw17-DIV17-Glu10uM-002-190803-173035-190803-175810_190803_183526_times_300-timee_300_length_7500360_p_0_u_0_chg_1_spks_-_emg.bin'
# rawdat = 'P_20190803-Gw17-DIV17-Glu10uM-002-190803-173035-190803-175810_190803_183026_-P_20190803-Gw17-DIV17-Glu10uM-002-190803-173035-190803-175810_190803_183526_times_0-timee_300_length_15000720_p_0_u_0_chg_1_spks_-_emg.bin'
# base_name = 'P_20190803-Gw17-DIV17-Glu10uM-002-190803-173035-190803-175810_190803_183526_-P_20190803-Gw17-DIV17-Glu10uM-002-190803-173035-190803-175810_190803_183526_times_300-timee_300_length_7500360_p_0_u_0_chg_1_'
# base_name = 'P_20190803-Gw17-DIV17-Glu10uM-002-190803-173035-190803-175810_190803_183026_-P_20190803-Gw17-DIV17-Glu10uM-002-190803-173035-190803-175810_190803_183526_times_0-timee_300_length_15000720_p_0_u_0_chg_1_'
rawdat = 'Headstages_128_Channels_int16_2020-06-05_08-47-46.bin'
base_name = 'Headstages_128_Channels_int16_2020-06-05_08-47-46'
# Pull raw data into memory.
nchans = 128
nsec = 2
fs = 25000
f = open(datdir+'/'+rawdat, 'rb')
dr = np.fromfile(f, dtype=np.int16,  count=fs*nsec*nchans) # count is Fs x seconds x nchans
length = np.int64(np.size(dr)/nchans)
drr1 = np.reshape(dr, [nchans, length], order='F')
f.close()
drr = ntk.butter_bandpass(drr1, 250, 7500, 25000, 3)

spkclusts = np.load(datdir+'/'+base_name+'spike_clusters.npy')
spktimes = np.load(datdir+'/'+base_name+'spike_times.npy')
templates = np.load(datdir+'/'+base_name+'template_waveform.npy')
uclusts = np.load(datdir+'/'+base_name+'unique_clusters.npy')
maxchans = np.load(datdir+'/'+base_name+'max_channel.npy')
# quals = np.genfromtxt( datdir+'/../'+base_name+'cluster_group.tsv',delimiter = '\t',dtype = 'str',skip_header =1)
quals = np.genfromtxt( datdir+'/'+base_name+'cluster_group.tsv',delimiter = '\t',dtype = 'str',skip_header =1)
# thesechans = np.arange(17,26)
# thesechans = np.arange(0,16)
thesechans = np.arange(16,32)
# thesechans = np.array([1, 2, 5, 6, 14])
# thesechans = np.array([6, 14])
theseclusts = [i for i, x in enumerate(maxchans) if x in thesechans]
# Taken from SashaTrubetskoy
# dcolors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
#            '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
#            '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000', '#fff000']
current_palette = sns.color_palette()
# dcolors = sns.color_palette(palette="bright", n_colors=len(theseclusts))
dcolors = sns.color_palette(palette="deep", n_colors=len(theseclusts))
#ax_list = fig.axes
fig, ax = plt.subplots(ncols = 1, nrows = np.size(thesechans), figsize = [25,15], sharex=True, gridspec_kw={'hspace': 0})
for i,ch  in enumerate(thesechans):
    ax[i].plot(drr[ch,:], color = [0.4,0.4,0.4], linewidth = 0.5)
    sns.despine(left = True, bottom = True)
# for c,clust in enumerate(theseclusts):
#
#
#     temptimes = np.squeeze(spktimes[np.where(spkclusts == clust)[0]])
#     temptimes = temptimes[np.where(temptimes<=25000*2)]
#
#     for j,k in zip (fig.axes, thesechans ): # j is the axis, k is the channel (use for y data)
#
#         for t in temptimes:
#
#             if maxchans[clust][0] == k:
#                 alphaval = 1
#             else:
#                 alphaval = 0.3
#
#
#             if quals[clust][1] != 'noise': # don't plot noise clusters
#
#                 if quals[clust][1] == 'good':
#                     # solid line plotting for single units
#                     j.plot(np.arange(t-30,t+30), drr[k, int(t-30):int(t+30)], color = dcolors[c], alpha = alphaval)
#
#                 elif quals[clust][1] == 'mua':
#                     #scatter plot MUA
#                     sns.scatterplot(np.arange(t-30,t+30), drr[k, int(t-30):int(t+30)], color = dcolors[c],ax = j, alpha = alphaval)
#
for c,clust in enumerate(theseclusts):
    temptimes = np.squeeze(spktimes[np.where(spkclusts == clust)[0]])
    temptimes = temptimes[np.where(temptimes<=25000*nsec)]
    for j,k in zip (fig.axes, thesechans ): # j is the axis, k is the channel (use for y data)
        for t in temptimes:
            if maxchans[clust][0] == k:
                alphaval = 1
            else:
                alphaval = 0.3
            # if quals[c][1] != 'noise': # don't plot noise clusters
            if quals[c][1] == 'good':
                # solid line plotting for single units
                j.plot(np.arange(t-30,t+30), drr[k, int(t-30):int(t+30)], color = dcolors[c], alpha = alphaval)
            elif quals[c][1] == 'mua':
                #scatter plot MUA
                sns.scatterplot(np.arange(t-30,t+30), drr[k, int(t-30):int(t+30)], color = dcolors[c],ax = j, alpha = alphaval)
            # elif quals[c][1] == 'noise':
            #     #scatter plot MUA
            #     sns.scatterplot(np.arange(t-30,t+30), drr[k, int(t-30):int(t+30)], color = dcolors[c],ax = j, alpha = alphaval)
xlims = ax[0].get_xlim();
ax[0].set_xlim([0,nsec*fs])
plt.draw()
# labels = [item.get_text() for item in ax[-1].get_xticklabels()]
# print("labels ", labels)
# labels2 = [int(int(i)/25) for i in labels]
# print("labels2 ", labels2)
# ax[-1].set_xticklabels(labels2);
# ax[-1].set_xlabel('Time (msec)',fontsize = 16);
fig.text(0.09, 0.5, 'Voltage (uV)', va='center', rotation='vertical',fontsize = 16);
ax[0].set_title('Clusters {} on channels {}'.format(theseclusts, thesechans));
# make sure you're either in the directory you want to save into, or else add a savedir to the beginning of the filename
fn = ('clusters_ch{}_to_{}.pdf'.format(thesechans[0],thesechans[-1]))
plt.savefig(fn)