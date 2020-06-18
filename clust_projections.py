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
import matplotlib.backends.backend_pdf as mpdf


rawdat_dir = '/Volumes/carina/EAB00050/block3/'
save_loc = '/Volumes/HlabShare/clayton_sahara_work/clustering/'
animal_name = 'EAB50'

# rawdat_dir = '/Volumes/bs006r/CAF00022/CAF00022_2020-06-05_08-47-46/'
# save_loc = '/Volumes/HlabShare/clayton_sahara_work/clustering/caf22/'

chan_map = 'EAB50chmap_00'
probe_num = 0
num_chans = 512
#cells = ''
loc = 1
fs=25000
window_len = 25000

def clust_projections_loc(rawdat_dir, save_loc, chan_map, probe_num, num_chans, cells, loc, fs=25000, window_len=25000):
    num_probes = int(num_chans/64)

    chan_map_full = np.repeat(chan_map, num_probes)

    start_time = cells[0].rstart_time
    num_files_away = int(loc/300)
    all_files = sorted(glob.glob(rawdat_dir+'/*'))
    first_file = [f for f in all_files if start_time in f][0]
    index_first = all_files.index(first_file)
    this_file = all_files[index_first+num_files_away]


    ts = loc*fs 
    te = ts+window_len

    xbins = np.arange(ts,te)

    t, dat = ntk.load_raw_gain_chmap_1probe(this_file, num_chans,
                                chan_map_full, nprobes=int(num_chans/64),
                                lraw=1, ts=ts%(300*fs), te=te%(300*fs),
                                probenum=probe_num, probechans=64)


    clust_groups = np.reshape(np.arange(0,64), (16,4))

    with mpdf.PdfPages(os.path.join(save_loc + f"/clust_projections_{animal_name}_sec{loc}.pdf")) as pdf:
        for thesechans in clust_groups:

            figc, ax = plt.subplots(ncols=1, nrows=len(thesechans),figsize=[16, 10], sharex=True,gridspec_kw={'hspace': 0})
            thesecells = [cell for cell in cells if cell.peak_channel in thesechans]
            dcolors = sns.color_palette(palette="deep", n_colors=len(thesecells))

            for i, c in enumerate(thesechans):
                ch_dat = dat[c,:]

                ch_dat = ntk.butter_bandpass(ch_dat, 500, 10000, fs=fs)

                ax[i].plot(xbins, ch_dat, color=[0.4, 0.4, 0.4], linewidth=0.5)
                sns.despine(left=True, bottom=True)
                ax[i].set_ylabel(f"Channel: {c}")
                ylim = ax[i].get_ylim()
                if ylim[0]<-500:
                    ax[i].set_ylim(bottom=-500)
                if ylim[1]>100:
                    ax[i].set_ylim(top=100)
            
            for c, cell in enumerate(thesecells):
                peak_channel = cell.peak_channel
                color = dcolors[c]

                t0 = ts + 30
                t1 = te - 30
                

                temptimes=cell.spike_time
                idx = np.where(np.logical_and(temptimes<t1, temptimes>t0))
                tempspikes = temptimes[idx]
                spikes = tempspikes - ts

                for a, chan in zip(ax, thesechans):
                    
                    chan_dat = dat[chan, :]
                    chan_dat = ntk.butter_bandpass(chan_dat, 500, 10000, fs=fs)
                    if peak_channel == chan:
                        alpha = .75
                    else:
                        alpha=.5
                    
                    for si, s in enumerate(spikes):
                        a.plot(np.arange(s+ts-30, s+ts+30), chan_dat[s-30:s+30], color = color, alpha=alpha, linewidth=2, label=f'Unit {cell.clust_idx}' if  alpha == .75 and si == 0 else '')
            figc.legend()
            plt.tight_layout()
            ax[0].set_title(f'Clusters {[cell.clust_idx for cell in thesecells]} on channels {thesechans}')
            pdf.savefig(figc)
            plt.close(figc)
                





            
            

        

