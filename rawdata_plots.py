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
import sys
import scipy as sc


'''
    #load and set basic info
    cells = np.load("path/to/cells", allow_pickle=True)
    raw_dir = "path/to/directory/with/raw/files"
    num_chans= 128
    num_probes=2
    chan_map = np.repeat('EAB50chmap_00', num_probes)

    #need to find a cell to look at 
    good_cells = [cell for cell in cells if cell.quality in [1]]

    #pick one of those cells
    clust_num = good_cells[0].clust_idx

    #run the function
    def amp_spike_projection(cells, clust_num, num_chans, probe_num, raw_dir, chan_map, window_len=25000, fs=25000)
'''




def all_chans_projected_spikes(raw_file, nchans, hstype, loc_in_data, clusters, nsec_window=1, fs=25000):
    """
        plots all the raw data from a specific time point in the data with
        the cluster spikes highlighted 

        raw_file: path to raw data to plot
        loc_in_data: location in the data to plot (in sample points)
        clusters: cells to use to map spikes onto, either a path or the loaded cells
        nsec_window: number of seconds to plot, defaults to 1
        fs: sampling rate, 25000

        returns: the plot with all the channels 



        NOT DONE YET
    """

    # load the clusters
    if type(clusters)==str:
        if os.path.exists(clusters):
            cells = np.load(clusters, allow_pickle=True)
        else:
            print("--- That is not a valid path name for the clusters --- exiting")
            return
    else:
        cells = clusters
    

    #load raw data
    nps = int(nchans/64)

    # this isn't working, talk to kiran to figure out best way to load all channels
    tr, drr = ntk.load_raw_binary_gain_chmap(raw_path, number_of_channels = nchans, hstype=hstype, nprobes=nps)



def one_channel_projected_spikes(ax, chan_dat, all_cells, clust_num, chan_num, loc, window_len, fs=25000):
    '''
    plots one chan with spikes of that cell highlighted on it and the other cells on that channel plotted lighter

    chan_dat: channel data already loaded
    all_cells: all the clusters already loaded
    clust_num: main cluster in question
    chan_num: what channel this is
    loc: where in the data to start looking in sample points
    window_len: num sample points to look at
    fs: sample rate
    '''

    data = ntk.butter_bandpass(chan_dat, 500, 10000, fs=25000)
    st = loc
    if st<0:
        st = 0
    ax.plot(np.arange(st, loc+window_len), data, color='darkgrey')

    cells_on_chan = [cell for cell in cells if cell.peak_channel == chan_num]
    dcolors = sns.color_palette(palette="deep", n_colors=len(cells_on_chan))

    for i, cell in enumerate(cells_on_chan):
        spikes = cell.spike_time
        spikes = spikes[spikes<loc+window_len]
        spikes = spikes[spikes>loc]
        spikes = spikes - loc

        if cell.clust_idx == clust_num:
            alpha = 1
        else:
            alpha = 0.4
        
        for spike in spikes:
            ax.plot(np.arange((spike+loc)-50, (spike+loc)+50), data[spike-50 : spike+50], color=dcolors[i], alpha = alpha)


def find_5_min(start_string, raw_dir, sample_point, fs=25000):
    whole_dir = sorted(glob.glob(raw_dir+'/*'))
    first_file = [x for x in whole_dir if start_string in x]
    idx=whole_dir.index(first_file[0])

    sps_per_file = 300*fs
    num_files_away = int(sample_point/sps_per_file)
    print(num_files_away)

    st_base = start_string[:14]
    minu = int(start_string[14:16])
    new_min = minu+(num_files_away*5)
    new_time_str = st_base+str(new_min)
    print(new_time_str)

    new_file = whole_dir[idx+num_files_away]

    return new_file


def onclick(event):
    point[0]=int(event.xdata)

point = [0]
def amp_spike_projection(cells, clust_num, num_chans, probe_num, raw_dir, chan_map, window_len=25000, fs=25000):
    plt.ion()
    cell = [cell for cell in cells if cell.clust_idx == clust_num][0]

    # chan_map = ntk.find_channel_map(chan_map, num_chans)
    # peak_chan = chan_map[cell.peak_channel]
    peak_channel = cell.peak_channel

    start_string = cell.rstart_time

    fig = plt.figure(constrained_layout=True, figsize=(11, 5))
    gs = fig.add_gridspec(2, 3)

    #plot wf
    wf_ax = fig.add_subplot(gs[:-1, 0])
    amp_ax = fig.add_subplot(gs[0,1:])
    raw_ax = fig.add_subplot(gs[-1, :])

    line_ax = amp_ax.twinx()

    wf_sh = cell.waveform_tetrodes.shape[0]
    wf_ax.plot(np.linspace(0, (wf_sh * 1000.0) / cell.fs,
                                    wf_sh),
                        cell.waveform_tetrodes,
                        color='#6a88f7')
    wf_ax.plot(np.linspace(0, (wf_sh * 1000.0) / cell.fs,
                                    wf_sh),
                        cell.waveform,
                        'g*')
    
    #plot amp
    th_zs = 2
    z = np.abs(sc.stats.zscore(cell.spike_amplitude))
    amp_ax.plot((cell.spike_time / cell.fs),
                cell.spike_amplitude, 'bo',
                markersize=1.9, alpha=0.2)
    amp_ax.plot((cell.spike_time[np.where(z > th_zs)] / cell.fs),
                cell.spike_amplitude[np.where(z > th_zs)],
                'ro',
                markersize=1.9, alpha=0.2)
    amp_ax.set_xlabel('Time (s)')
    amp_ax.set_ylabel('Amplitudes', labelpad=-3)
    amp_ax.set_xlim(left=(cell.start_time),
                    right=(cell.end_time))
    amp_ax.set_ylim(bottom=0,
                    top=(min((np.mean(cell.spike_amplitude) +
                                (3*np.std(cell.spike_amplitude))),
                                np.max(cell.spike_amplitude))))
    amp_stats_str = '\n'.join((
        r'$Min: %d, Max: %d$' % (np.min(cell.spike_amplitude),
                                    np.max(cell.spike_amplitude), ),
        r'$Mean:%d, Med:%d, Std:%d$'
        % (np.mean(cell.spike_amplitude),
            np.median(cell.spike_amplitude),
            np.std(cell.spike_amplitude), )))
    props = dict(boxstyle='round', facecolor='wheat',
                    alpha=0.5)
    amp_ax.text(0.73, 0.27, amp_stats_str,
                transform=amp_ax.transAxes,
                fontsize=8,
                verticalalignment='top', bbox=props)
    
    done = False

    fig.canvas.mpl_connect('button_press_event', onclick)

    
    while(not done):
        plt.show()
        print("waiting")
        plt.waitforbuttonpress()
        #point = [int(input("What second do you want to look at?"))]
        print("samp point: ", point[0])
        
        file_to_load = find_5_min(start_string, raw_dir, point[0]*fs, fs=25000)

        ts = ((point[0]%300) * fs )
        if ts<0:
            ts=0
        
        te = ts + window_len
        if te>(300*fs):
            te=-1
        t, dat = ntk.load_raw_gain_chmap_1probe(file_to_load, num_chans,
                               chan_map, nprobes=int(num_chans/64),
                               lraw=1, ts=ts, te=te,
                               probenum=probe_num, probechans=64)

        ch_dat = dat[peak_channel,:]

        one_channel_projected_spikes(raw_ax, ch_dat, cells, clust_num, peak_channel, point[0]*fs, window_len, fs=25000)
        line_ax.plot([point[0], point[0]], [0, line_ax.get_ylim()[1]], color="black")
        plt.show()
        plt.pause(1)
        save = input("Would you like to save this figure? y / n")=='y'
        if save:
            fig.savefig(f"projected_spikes_w_amp_second-{point[0]}")
        done = input("Do you want to choose another segment to plot? y / n") == 'n'
        if not done:
            line_ax.clear()
            raw_ax.clear()
    





        


        







