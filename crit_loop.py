from sahara_work import Criticality_new as cr
import musclebeachtools_hlab.musclebeachtools as mbt
import neuraltoolkit as ntk
import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as mpdf
import seaborn as sns
import csv
import os
import glob


def get_timeframe(path):
    last_slash = path.rfind("/")
    base_path = path[0:last_slash+1]
    if base_path.find("/probe") < 0:
        print("The path is not complete. Please pass in paths that contain the time frame directory")
        return -1
    time_frame_1 = base_path[0:base_path.find("/probe")]
    time_frame = time_frame_1[time_frame_1.rfind('/')+1:]
    return time_frame, base_path


def get_totaltime(time_frame):
    start_time = int(time_frame[0:time_frame.find('_')])
    stop_time = int(time_frame[time_frame.find('_')+1:])
    total_time = stop_time-start_time
    return total_time


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

def looped_crit(FR_mat, shuffled_FR_mat, params,basepath, plot_shuffled=True):
    quality = params['quality']
    ava_binsz=params['ava_binsz']
    perc = params['perc']
    hour_bins=params['hour_bins']
    total_time=params['total_time']
    time_frame = params["time_frame"]
    animal = params["animal"]

    num_bins = int(total_time/hour_bins)
    bin_len = int((hour_bins*3600)/ava_binsz)

    master_dict = {}
    all_p_values_burst = []
    all_p_values_t=[]
    all_dcc_values = []

    qual_str = '_'.join(map(str,quality))
    param_str = f'{time_frame}_{str(hour_bins)}hrs_perc{str(int(perc*100))}_binsz{str(int(ava_binsz*1000))}ms_q{qual_str}'
    csv_filename = f'alpha_beta_ratios_{param_str}.csv'

    for idx, t_bin in enumerate(np.arange(0, num_bins)):
        print(f"working on block {idx+1} of {num_bins}")
        
        if idx == num_bins-1:
            data=FR_mat[:, (idx*bin_len):]
            if plot_shuffled:
                data_shuffled = shuffled_FR_mat[:, (idx*bin_len):]
        else:
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
            
        burstM = int(np.max(Result['S'])/20)
        tM = int(np.max(Result['T'])/20)

        Result = cr.AV_analysis_new(Result["S"], Result["T"], burstM, tM, param_str+'_'+str(idx),saveloc=basepath, flag = 2, burst_shuffled=burst_shuffled, T_shuffled=time_shuffled, plot_shuffled=plot_shuffled, plot=True)
        ratio_to_csv(Result['alpha'][0], Result['beta'][0], f'{time_frame}_{idx+1}', csv_filename)
        master_dict["Result_block"+str(idx)] = Result
        all_p_values_burst.append(Result["P_burst"])
        all_p_values_t.append(Result["P_t"])
        all_dcc_values.append(Result["df"][0])
        print("P_VALUES: ", all_p_values_burst[-1], " -- ", all_p_values_t[-1])
        print("DCC: ", all_dcc_values[-1])

    master_dict["all_p_values_burst"]=all_p_values_burst
    master_dict["all_p_values_t"]=all_p_values_t
    master_dict["all_dcc_values"]=all_dcc_values
    master_dict["parameters"] = params

    return master_dict
 

params = {
    'ava_binsz': 0.045,
    'hour_bins': 4,
    'perc': 0.25,
    'quality': [1,2],
    'animal' : 'EAB26',
    'notes': 'testing'
}

def lilo_and_stitch(paths, params, overlap=0, plot_shuffled=True):
    all_data = []
    all_ps  = []
    all_dccs = []

    quality = params['quality']
    ava_binsz=params['ava_binsz']
    perc = params['perc']
    hour_bins=params['hour_bins']
    animal = params["animal"]

    for idx, path in enumerate(paths):
        print("------- WORKING ON ", path, " --------")

        if get_timeframe(path) == -1:
            return 0
        params['time_frame'], basepath =get_timeframe(path)
        params['total_time']=get_totaltime(params['time_frame'])

        cells = np.load(path, allow_pickle=True)
        good_cells = [cell for cell in cells if cell.quality in quality]
        print(f"Number of cells: {len(good_cells)}")
        data = mbt.n_spiketimes_to_spikewords(good_cells, binsz=ava_binsz, binarize=1)

        if plot_shuffled:
            data_shuffled = Genshuffle(good_cells, params['total_time'], ava_binsz, perc, binary = 1, frame = 0)
        else:
            data_shuffled=None

        master_dict = looped_crit(data, data_shuffled, params,basepath, plot_shuffled=plot_shuffled)

        qual_str = '_'.join(map(str,quality))
        np.save(f'{basepath}{animal}_dict_{params["time_frame"]}_{str(hour_bins)}hrs_perc{str(int(perc*100))}_binsz{str(int(ava_binsz*1000))}ms_q{qual_str}', master_dict)
        all_data.append(master_dict)

    return all_data




    