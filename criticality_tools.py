from sahara_work import Criticality as cr 
import musclebeachtools_hlab.musclebeachtools as mbt 
import neuraltoolkit as ntk
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import matplotlib.backends.backend_pdf as mpdf
import seaborn as sns
import csv
import os
import glob

params={
    'model_start_time' : 0,
    'model_stop_time' : 100,
    'ava_binsz' : 1, # in ms
    'perc' : .65
}
def run_model_data(model_data_path, params):

    model_data = np.load(model_data_path, allow_pickle=True)
    plt_name = model_data_path[model_data_path.rfind('/')+1:model_data_path.rfind('.npy')]

    model_data_seconds = model_data/1000
    data = cr.spiketimes_to_spikewords(model_data_seconds, startime=params['model_start_time'], stoptime=params['model_stop_time'], binsize=params['ava_binsz'], binarize=1)
    r = cr.AV_analysis_BurstT(data, perc = params['perc'])
    burst = r['S'] 
    duration = r['T'] 

    burstM = int(np.max(burst)/40)
    tM = int(np.max(duration)/20)
    if tM == 0:
        tM = int(np.percentile(duration,10))

    Result2 = cr.AV_analysis_ExponentErrorComments(burst, duration, burstM, tM, flag=2, pltname=plt_name)
    Result3 = cr.AV_analysis_ExponentErrorComments(burst, duration, burstM, tM, flag=3, pltname=plt_name)
    print(Result3[0]['df'][0])

    return Result2, Result3


def FR_plot(cells, binsz, rec_len, color=False):
    plt.ion()
    fig,ax = plt.subplots(nrows=1, ncols=1)
    fig.set_figheight(3)
    fig.set_figwidth(12)

    total_sps = rec_len*3600*25000 # total sample points. rec_len should be in hours
    bin_sps = binsz*25000 # total sample points in each bin, binsz should be in seconds

    bins = np.arange(0, (total_sps + bin_sps), bin_sps) # creates the bins for the firing rate plot. takes total recording length and breaks it up by binsize
    
    colors=['seagreen', 'gold', 'salmon']
    for cell in cells:
        spk_times=cells.spike_time
        fr=np.histogram(spk_times, bins=bins)
        fr=fr[0]/600
        a

def ratio_to_csv(alpha, beta, block, filename):

    if(os.path.exists(filename)):
        with open(filename, mode='a+') as ratio_file:
            writer = csv.writer(ratio_file, delimiter=',' , quotechar='"')
            a = str(alpha)
            b=str(beta)
            c=str(beta/alpha)
            writer.writerow([block, a,b,c])
    else:
        with open(filename, mode='w') as ratio_file:
            writer = csv.writer(ratio_file, delimiter=',' , quotechar='"')
            writer.writerow(['Block', 'Alpha', 'Beta', 'Beta/Alpha'])
            a = str(alpha)
            b=str(beta)
            c=str(beta/alpha)
            writer.writerow([block, a,b,c])

def pdf_output(neurons, rawdatdir, hstype, saveloc):
    base_time = neurons[0].rstart_time
    first_file = glob.glob(f'{rawdatdir}*_int16_{base_time}.bin')[0]
    num_channels = int(first_file[first_file.rfind("Headstages_")+11:first_file.rfind('_Channels')])
    raw = ntk.load_raw_binary_gain_chmap(first_file, num_channels, hstype)

    for cell in neurons:
        if(cell.rstart_time!=base_time):
            "loading new raw data file"
            new_file = glob.glob(f'{rawdatdir}*_int16_{cell.rstart_time}.bin')[0]
            raw = ntk.load_raw_binary_gain_chmap(first_file, num_channels, hstype)


def Genshuffle(neurons, nrn_time, ava_binsz, perc, binary = 1, frame = 0):

	spks = mbt.n_getspikes(neurons, 0, 3600*nrn_time)

	spks_shuffle = []

	if frame == 0:
	# generate random spiketimes (Uniform distribution, keep the same FR for each neuron)
		print('generate random spiketimes')
		for i in np.arange(0,len(spks)):
			spikes = spks[i]
			numspk = np.shape(spikes)[0]
			spikes_shuffle = np.random.random(size = numspk) * nrn_time * 3600
			spks_shuffle.append(spikes_shuffle)
	
	else:
	# frame shuffling; For each neuron, swap all spikes in one time bin (30s window) with spikes at another randomly chosen time bin.
		print('Frame Shuffling')

		endtime = nrn_time*3600
		for i in np.arange(0,len(spks)): 
			spikes = spks[i] 
			spikes_shuffle = frameshuffle(spikes,endtime) 
			spks_shuffle.append(spikes_shuffle)
	
	data_T = cr.spiketimes_to_spikewords(spks_shuffle,0,3600*nrn_time,ava_binsz*1000,binary) # 1 for binary
	data_shuffle = data_T.T

	return data_shuffle

def frameshuffle(spks,endtime):
    for j in np.arange(0,1000):
        frame1 = np.random.rand() * endtime
        frame2 = np.random.rand() * endtime
        idx1 = np.where(np.logical_and(spks <= frame1 + 30, spks >= frame1 ))[0]
        idx2 = np.where(np.logical_and(spks <= frame2 + 30, spks >= frame2 ))[0]
        spks[idx1] = spks[idx1] - frame1 + frame2
        spks[idx2] = spks[idx2] - frame2 + frame1
    spks = np.sort(spks)  
    return spks

def pull_crit_data(all_dicts, save_loc, animal, time_frame, paths=False):
    """
    takes in an array of dictionaries and pulls all the dcc and p_value data, appends all that data into a nice block and saves it where you want
    """
    all_dicts_obj=[]

    if paths:
        for path in all_dicts:
            print (f'loading dict: {path}')
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

# params={
#     "animal": "caf19",
#     "date": "0409",
#     "time_range":"108-168"
# }

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

def looped_crit(FR_mat, shuffled_FR_mat, params, plot_shuffled=True, plot=False):
    """
    FR_mat: binarized firing matrix binned by miliseconds and then reshaped into chunks
    binsz: how the binarized matrix is binned (in ms)
    quality: list of neuron qualities included, if you used qualities 1 and 2 then [1,2], only 1 would be [1]
    perc: percent of neurons needed for an avalance, 0-1 (.25 is 25%)
    burstM: usually 12
    tM: usually 4

    the master_dict that is returned contains the result of each criticality block: 
        master_dict['Result'] = {
            burst: 
            alpha:
            xmin:
            xmax:
            T:
            beta:
            tmin:
            tmax:
            pre:
            fit:
            df:
        }
    in addition it has a list of all the p_values:
        master_dict['all_p_values_burst/t] = [... , ... , ...]
    and dcc values:
        master_dict['all_dcc_values] = [... , ... , ...]
    and the parameters

    it also has all the figures produced by running criticality. 
        "dcc_ac_block0"

    """
    quality = params['quality']
    ava_binsz=params['ava_binsz']
    perc = params['perc']
    burstM=params['burstM']
    tM = params['tM']
    hour_bins=params['hour_bins']
    total_time=params['total_time']
    time_frame = params["time_frame"]
    animal = params["animal"]

    total_bins=len(FR_mat[0])
    total_seconds=total_bins*ava_binsz*1000 # converts to seconds from ms
    total_hours=total_seconds/3600

    num_bins = int(total_time/hour_bins)
    bin_len = int((hour_bins*3600)/ava_binsz)

    master_dict = {}
    all_p_values_burst = []
    all_p_values_t=[]
    all_dcc_values = []

    qual_str = '_'.join(map(str,quality))
    param_str = f'{time_frame}_{str(hour_bins)}hrs_perc{str(int(perc*100))}_binsz{str(int(ava_binsz*1000))}ms_bm{str(burstM)}_tm{str(tM)}_q{qual_str}'
    csv_filename = f'alpha_beta_ratios_{param_str}.csv'
    for idx, t_bin in enumerate(np.arange(0, num_bins)):
        print(f"working on block {idx+1} of {num_bins}")
        
        if idx == num_bins-1:
            # exec(f"data{idx+1}=FR_mat[:, (idx*bin_len):]")
            data=FR_mat[:, (idx*bin_len):]
            if plot_shuffled:
                data_shuffled = shuffled_FR_mat[:, (idx*bin_len):]
        else:
            # exec(f"data{idx+1}=FR_mat[:, (idx*bin_len) : ((idx+1)*bin_len)]")
            data=FR_mat[:, (idx*bin_len) : ((idx+1)*bin_len)]
            if plot_shuffled:
                data_shuffled = shuffled_FR_mat[:, (idx*bin_len) : ((idx+1)*bin_len)]

        Result = cr.AV_analysis_BurstT(data, perc=perc)
        if plot_shuffled:
            Result_shuffled = cr.AV_analysis_BurstT(data_shuffled, perc=perc)
            burst_shuffled=Result_shuffled['S']
            time_shuffled=Result_shuffled['T']
            
        else:
            burst_shuffled=None
            time_shuffled=None

        Result2, ax1, ax2 = cr.AV_analysis_ExponentErrorComments(Result["S"], Result["T"], burstM, tM, param_str+'_'+str(idx), flag = 2, EX_burst=1, EX_time=1)
        Result3, ax3 = cr.AV_analysis_ExponentErrorComments(Result["S"], Result["T"], burstM, tM, param_str+'_'+str(idx), flag = 3, burst_shuffled=burst_shuffled, T_shuffled=time_shuffled, plot_shuffled=plot_shuffled)
        
        ratio_to_csv(Result3['alpha'][0], Result3['beta'][0], f'{time_frame}_{idx+1}', csv_filename)
        master_dict["Result_block"+str(idx)] = Result2
        master_dict["dcc_ax_block"+str(idx)] = ax3
        master_dict["p_test_axs_block"+str(idx)] = (ax1, ax2)
        all_p_values_burst.append(Result2["P_burst"])
        all_p_values_t.append(Result2["P_t"])
        all_dcc_values.append(Result3["df"][0])
        print("P_VALUES: ", all_p_values_burst[-1], " -- ", all_p_values_t[-1])
        print("DCC: ", all_dcc_values[-1])

    master_dict["all_p_values_burst"]=all_p_values_burst
    master_dict["all_p_values_t"]=all_p_values_t
    master_dict["all_dcc_values"]=all_dcc_values
    master_dict["parameters"] = params

    if plot:
        fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1)
        # for row in ax:
        #     for col in row:
        #ax.bar(np.arange(np.shape(FR_mat)[0]), master_dict['all_dcc_values'], color = '#464196', alpha = 0.7, label = 'DCC', zorder = 10)
        bar_width=0.4
        size_x = np.arange(len(FR_mat))+1
        dur_x = size_x + bar_width
        ax1.bar(size_x, master_dict["all_p_values_burst"], bar_width, color = '#a6c875', alpha = 0.7, label = 'DCC', zorder = 10)
        ax1.bar(dur_x, master_dict["all_p_values_t"], bar_width, color = '#f1da7a', alpha = 0.7, label = 'DCC', zorder = 10)
        ax1.set_xlim([0,len(FR_mat)+1]) 

        ax1.set_xticks(np.arange(np.size(size_x)+1)+bar_width/2)
        xlim = ax1.get_xlim()
        ax1.plot([xlim[0], xlim[1]], [0.05, 0.05], color = '#738595', linestyle = '--')
        ax1.set_ylabel('p value', fontsize = 20)
        ax1.set_xlabel('Bin num', fontsize = 20)

        ax2.bar(np.arange(len(FR_mat))+1, master_dict["all_dcc_values"], color = '#464196', alpha = 0.7, label = 'DCC', zorder = 10)
        ax2.set_ylim([0,1])
        ax2.set_xlim([0,len(FR_mat)+1]) 
        ax2.plot([0, 16.5], [0.2, 0.2], linestyle = '--', color = '#ff964f', zorder = 15)
        ax2.set_ylabel('DCC', fontsize = 20)
        ax2.set_xlabel("hours", fontsize=20)

        fig.savefig("criticality_figures")
    return master_dict


# params = {
#     'ava_binsz': 0.045,
#     'hour_bins': 4,
#     'total_time':12,
#     'perc': 0.25,
#     'burstM': 10,
#     'tM': 4,
#     'quality': [1,2],
#     'time_frame': '0420',
#     'animal' : 'caf19',
#     'notes': 'local field test'
# }


def lilo_and_stitch(paths, params, overlap=0, plot=False, plot_shuffled=True):
    """
    paths: list of paths to neuron objects from clustering, full paths recomended
    params: dictionary with necessary parameters:
        params = {
            "ava_binsz": 0.045,
            "hour_bins": 4,
            "perc":.25,
            "burstM":10,
            "tM":4,
            "quality":[1,2]
        }
    overlap: does the continuous data have an overlap? if no: 0 if yes: number of hours in the overlap (1)
    plot: keep false for now, it's not done yet lol
    """
    all_data = []
    all_ps  = []
    all_dccs = []

    quality = params['quality']
    ava_binsz=params['ava_binsz']
    perc = params['perc']
    burstM=params['burstM']
    tM = params['tM']
    hour_bins=params['hour_bins']
    total_time=params['total_time']
    time_frame = params["time_frame"]
    animal = params["animal"]

    for idx, path in enumerate(paths):
        print("------- WORKING ON ", path, " --------")
        last_slash = path.rfind("/")
        base_path = path[0:last_slash+1]
        
        time_frame = base_path[0:base_path.find("/probe")]
        params['time_frame']=time_frame
        
        cells = np.load(path, allow_pickle=True)
        good_cells = [cell for cell in cells if cell.quality in quality]
        data = mbt.n_spiketimes_to_spikewords(good_cells, binsz=ava_binsz, binarize=1)

        if plot_shuffled:
            data_shuffled = Genshuffle(good_cells, total_time, ava_binsz, perc, binary = 1, frame = 0)
        else:
            data_shuffled=None
        # if overlap:
        #     ms_in_overlap = overlap * 3600 * 1000
        #     bins_in_overlap = ms_in_overlap/(ava_binsz*1000)
        #     fr_mat_trimmed = fr_mat[:, 0:-int(bins_in_overlap)]
        #     data = break_up_mat(fr_mat_trimmed, ava_binsz, hour_bins)
        # else:
        #     fr_mat_reshaped = break_up_mat(fr_mat, ava_binsz, hour_bins)
        #     data = fr_mat_reshaped
        
        master_dict = looped_crit(data, data_shuffled, params, plot_shuffled=plot_shuffled, plot=False)

        qual_str = '_'.join(map(str,quality))
        np.save(f'{base_path}{animal}_dict_{time_frame}_{str(hour_bins)}hrs_perc{str(int(perc*100))}_binsz{str(int(ava_binsz*1000))}ms_bm{str(burstM)}_tm{str(tM)}_q{qual_str}', master_dict)
        all_data.append(master_dict)
    
    if plot:
        # plot the giant dcc and p_value plot here
        print('gonna plot')

    
    return all_data












    

