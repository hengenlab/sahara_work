from sahara_work import Criticality as cr 
import musclebeachtools_hlab.musclebeachtools as mbt 
import neuraltoolkit as ntk
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.backends.backend_pdf as mpdf
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import seaborn as sns
from copy import deepcopy as cdc
import seaborn as sns
import csv
import os
import glob


def pull_crit_data(all_dicts, save_loc, animal, time_frame, paths=False):
    """
    helper method for crit_plots

    takes in an array of dictionaries and pulls all the dcc and p_value data, appends all that data into a nice block and saves it where you want
    """
    all_dicts_obj = []

    if paths:
        for path in all_dicts:
            print(f'loading dict: {path}')
            all_dicts_obj.append(np.load(path, allow_pickle=True).item())
    else:
        all_dicts_obj=all_dicts

    all_dccs=[]
    all_p_t=[]
    all_p_b=[]
    all_params=[]
    for i, data in enumerate(all_dicts_obj):
        print(f"working on block: {i}")
        all_dccs.append(np.array(data['all_dcc_values']))
        all_p_t.append(np.array(data['all_p_values_t']))
        all_p_b.append(np.array(data['all_p_values_burst']))
        all_params.append(data['parameters'])
    
    all_dccs = np.array(all_dccs).flatten()
    all_p_t = np.array(all_p_t).flatten()
    all_p_b = np.array(all_p_b).flatten()

    all_data = [all_dccs, all_p_b, all_p_t]

    np.save(save_loc+f"/dcc_pb_pt_{animal}_{time_frame}", all_data)
    np.save(save_loc+f'/all_params_{animal}_{time_frame}', all_params)
    return all_data, all_params

params={
    "animal": "caf22",
    "date": "0526",
    "time_range":"0-36_P2"
}

def crit_plots(dcc, p_b, p_t, labels, params, save=False):
    """
    makes pretty plots from arrays of all dcc and p value data. not from the dictionaries - make sure theres an extra 0 in the labels
    """
    plt.ion()
    fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1)

    bar_width=0.4
    size_x = np.arange(len(dcc))+1
    dur_x = size_x + bar_width
    color_b = np.where(p_b<0.05, "black", '#a6c875')
    color_t = np.where(p_t<0.05, "black", '#f1da7a')
    ax1.bar(size_x, p_b, bar_width, color = color_b, alpha = 0.7, label = 'DCC', zorder = 10)
    ax1.bar(dur_x, p_t, bar_width, color = color_t, alpha = 0.7, label = 'DCC', zorder = 10)
    ax1.set_xlim([0,len(dcc)+1]) 
    ax1.set_xticks(np.arange(np.size(size_x)+1)+bar_width/2)
    xlim = ax1.get_xlim()
    ax1.plot([xlim[0], xlim[1]], [0.05, 0.05], color = '#738595', linestyle = '--')
    ax1.set_ylabel('p value', fontsize = 20)
    ax1.set_xticklabels(labels, rotation=50)


    color_dcc = np.where(dcc>0.2, "lightcoral", '#464196')
    ax2.bar(np.arange(len(dcc))+1, dcc, color = color_dcc, alpha = 0.7, label = 'DCC', zorder = 10)
    ax2.set_ylim([0,1])
    ax2.set_xlim([0,len(dcc)+1]) 
    ax2.plot([0, len(dcc)+1], [0.2, 0.2], linestyle = '--', color = '#ff964f', zorder = 15)
    ax2.set_ylabel('DCC', fontsize = 20)
    ax2.set_xlabel("time-bin", fontsize=20)
    ax2.set_xticks(np.arange(np.size(size_x)+1))
    ax2.set_xticklabels(labels, rotation=50)

    ax1.set_title(f"{params['animal']} data for time range {params['time_range']} on {params['date']}")

    plt.tight_layout()
    plt.show()
    if(save):
        fig.savefig(f"criticality_figures_{params['animal']}_{params['date']}_{params['time_range']}")
    
    return fig, ax1, ax2


    # small_bin: original size of bin used to make the matrix, in seconds
    # hour_bins: how you want the matrix broken up, in hours 
        # ex: 4 - would break up a 12 hour matrix into 3 arrays of 4 hours


    small_bin=small_bin*1000
    total_time = FR_mat.shape[1]*(small_bin) # in ms
    num_cells = FR_mat.shape[0]
    new_time = hour_bins*3600*1000 # in ms
    num_arrays = int(total_time/new_time)
    new_binsz = int(new_time/(small_bin))


    try:
        reshaped_array = np.reshape(FR_mat, (num_arrays, num_cells, new_binsz))
    except:
        print("Your new binsize doesn't divide perfectly into the total time --- creating an array with the last block holding the remainder")
        reshaped_array=[]
        edges = np.arange(0, np.shape(FR_mat)[1], new_binsz)
        for e in edges:
            subset = FR_mat[:, e : e+new_binsz]
            reshaped_array.append(subset)
        final_bin_len=(reshaped_array[-1].shape[1]*small_bin)/1000/3600
        if final_bin_len < 1 :
            reshaped_array = reshaped_array[0:-1]
        print(f"your final bin is of length: {final_bin_len} hours")
    
    return reshaped_array


def large_plot(start_folder, end_folder, animal, save=False, save_loc=None):
    """
    this function is the main wrapper function for the large crit plot

    takes a start and an end folder, goes through every folder in-between
    and pulls out all the completed criticality data and plots it 

    any folder without crit data gets added as "..." in the plot until theres a completed folder 


    PARAMETERS:

    start_folder: first folder to pull crit data from 
    end_folder: last folder to pull crit data from 
    save: do you want to save the figure?
    save_loc: where to save it?
    animal: animal name

    """


def raster_avalanches(cells, av_binsize, perc, trange=10):
    '''
    returns a plot with the binarized spikes organized into avalanches

    cells: cells to look at (only good ones that would be used in crit)
    av_binsize: how they should be binarized
    perc: threshold for avalanche
    range: seconds to dispaly

    returns: plot and ax 
    '''
    FR_mat = mbt.n_spiketimes_to_spikewords(cells, binsz=av_binsize, binarize=1)
    end_bin = trange/av_binsize
    sm_fr = FR_mat[:,0:end_bin]
    n,m = np.shape(sm_fr)


    positions =[]
    for cell in sm_fr:
        positions.append(np.where(cell>0)[0])
    network = np.nansum(sm_fr, axis=0)
    sortN = np.sort(network)
    threshold = sortN[round(m*perc)]


    thresh_max = np.ma.masked_where(network<=threshold, network)

    zdata = cdc(network)
    zdata[~thresh_max.mask] = 1 #avalanches
    zdata[thresh_max.mask] = 0 #intervals

    
    edges = np.diff(zdata)
    ontimes = np.where(edges>0)[0]
    offtimes = np.where(edges<0)[0]

    if zdata[0]==1:
        ontimes = np.insert(ontimes, 0, 0)
    if zdata[-1]==1:
        offtimes =np.append(offtimes,len(edges))
   
    xys=[]
    widths=[]
    heights=[]
    for i, on in enumerate(ontimes):
        xys.append((on,0))
        widths.append(offtimes[i]-on)

    boxes = [Rectangle(xy, width, n) for xy, width in zip(xys, widths)]
    pc = PatchCollection(boxes,facecolor='moccasin', alpha=0.7, edgecolor='firebrick')

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(11,8))
    ax.eventplot(positions, orientation="horizontal")
    ax.set_ylabel("Cell")
    ax.set_xlabel("Bin")
    ax.set_title(f'Avalanches over {trange} seconds')
    
    ax.add_collection(pc)
        
    plt.show()

    return fig, ax

    