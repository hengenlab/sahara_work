import criticality as cr
import musclebeachtools as mbt
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

def get_mean_isi(cells):
    isis = []
    for cell in cells:
        isi, _, _ = cell.isi_hist()
        isis.append(np.mean(isi))
    mean_isi = np.mean(np.asarray(isis))
    return mean_isi


def get_timeframe(path, time_frame):
    last_slash = path.rfind("/")
    base_path = path[0:last_slash+1]
    if time_frame is None:
        if base_path.find("/probe") < 0:
            print("The path is not complete. Please pass in paths that contain the time frame directory")
            return -1
        time_frame_1 = base_path[0:base_path.find("/probe")]
        time_frame = time_frame_1[time_frame_1.rfind('/')+1:]
        return time_frame, base_path
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


def looped_crit(FR_mat, shuffled_FR_mat, params, plot_shuffled=False):
    quality = params['quality']
    ava_binsz=params['ava_binsz']
    perc = params['perc']
    hour_bins=params['hour_bins']
    total_time=params['total_time']
    time_frame = params["time_frame"]
    animal = params["animal"]

    num_bins = int(total_time/hour_bins)
    bin_len = int((hour_bins*3600)/ava_binsz)

    if params['end_bin'] == -1:
        params['end_bin'] = num_bins

    master_dict = {}
    all_p_values_burst = []
    all_p_values_t=[]
    all_dcc_values = []

    qual_str = '_'.join(map(str,quality))
    cell_str = '_'.join(map(str, params['cell_type']))
    param_str = f'{time_frame}_{str(hour_bins)}hrs_perc{str(int(perc*100))}_binsz{str(int(ava_binsz*1000))}ms_q{qual_str}_cells{cell_str}'

    for idx in np.arange(params['start_bin'], params['end_bin']):
        print(f"working on block {idx+1} of {num_bins}")
        
        if idx == num_bins-1:
            data=FR_mat[:, (idx*bin_len):]
            if plot_shuffled:
                data_shuffled = shuffled_FR_mat[:, (idx*bin_len):]
        else:
            data=FR_mat[:, (idx*bin_len) : ((idx+1)*bin_len)]
            if plot_shuffled:
                data_shuffled = shuffled_FR_mat[:, (idx*bin_len) : ((idx+1)*bin_len)]

        Result = cr.get_avalanches(data, perc=perc)

        if plot_shuffled:
            Result_shuffled = cr.get_avalanches(data_shuffled, perc=perc)
            burst_shuffled=Result_shuffled['S']
            time_shuffled=Result_shuffled['T']
        else:
            burst_shuffled=None
            time_shuffled=None
        
        crit_params = {
            'flag': params['flag'],
            'pltname': param_str+'_'+str(idx),
            'saveloc': params['saveloc'],
            'burst_shuffled': burst_shuffled,
            'T_shuffled': time_shuffled,
            'plot_shuffled': False,
            'plot': params['plot'], 
            'nfactor_bm_tail':params['nfactor_bm_tail'], 
            'nfactor_tm_tail':params['nfactor_tm_tail']
        }


        crit_params['bm'] = int(np.max(Result['S'])/20)
        crit_params['tm'] = int(np.max(Result['T'])/20)

        Result = AV_analysis(Result["S"], Result["T"], crit_params, nfactor_bm = params['nfactor_bm'], nfactor_tm = params['nfactor_tm'], nfactor_bm_tail=crit_params['nfactor_bm_tail'], nfactor_tm_tail=crit_params['nfactor_tm_tail'])
        master_dict["Result_block"+str(idx)] = Result

        if params['flag'] == 2:
            all_p_values_burst.append(Result["P_burst"])
            all_p_values_t.append(Result["P_t"])
            print("P_VALUES: ", all_p_values_burst[-1], " -- ", all_p_values_t[-1])

        all_dcc_values.append(Result["df"])
        print("DCC: ", all_dcc_values[-1])

    master_dict["all_p_values_burst"]=all_p_values_burst
    master_dict["all_p_values_t"]=all_p_values_t
    master_dict["all_dcc_values"]=all_dcc_values
    master_dict["parameters"] = params

    return master_dict


# paths = [
#         '/media/HlabShare/clayton_sahara_work/clustering/eab47_0615_19/0_8/probe2/co/H_2019-06-15_20-16-59_2019-06-16_04-12-02_neurons_group0_scored_clayton.npy']

params = {
    'flag': 2, # 1 is DCC 2 is p_val and DCC
    'ava_binsz': 0.04, # in seconds
    'hour_bins': 4,# durration of block to look at
    'start_bin': 0, # start block to look at for criticality of the total blocks, 
                    #o if its a 12 hour block broken into 3x4hr bins and you only want 
                    #to look at the last block then the start block would be 2 (index at 0) 
                    #and the end block would be 3
    'end_bin': -1,
    'perc': 0.25,
    'nfactor_bm':0, 
    'nfactor_tm':0,
    'nfactor_bm_tail':.8, # upper bound to start exclude for burst
    'nfactor_tm_tail': .7, # upper bound to start exclude for time
    'quality': [1,2,3],
    'cell_type': ['RSU', 'FS'], # all qualities would be [1,2,3]
    'animal' : 'CAF19',
    'saveloc' : "/media/HlabShare/clayton_sahara_work/criticality/caf19/0420/",
    'notes': 'plz fucking work',
    'time_frame_list':['0_8', '8_16', '16_24', '24_32', '32_40', '40_48', '48_56',
                        '56_64', '64_72', '72_80', '80_88', '88_96', '96_104', '104_112',
                        '112_120', '120_128', '128_136', '136_144', '144_152', '152_160',
                        '160_168', '168_176'], # if you're running multiple blocks and the paths are right - this should be None
    'plot' : True
    }

params = {
    'flag': 2, # 1 is DCC 2 is p_val and DCC
    'ava_binsz': 0.04, # in seconds
    'hour_bins': 4,# durration of block to look at
    'start_bin': 0, # start block to look at for criticality of the total blocks, 
                    #o if its a 12 hour block broken into 3x4hr bins and you only want 
                    #to look at the last block then the start block would be 2 (index at 0) 
                    #and the end block would be 3
    'end_bin': -1,
    'perc': 0.35,
    'nfactor_bm':0, 
    'nfactor_tm':0,
    'nfactor_bm_tail':.8, # upper bound to start exclude for burst
    'nfactor_tm_tail': .8, # upper bound to start exclude for time
    'quality': [1,3], # all qualities would be [1,2,3]
    'cell_type': ['FS', 'RSU'], 
    'animal' : 'caf22',
    'saveloc' : "/media/HlabShare/clayton_sahara_work/criticality/caf22/0508/",
    'notes': 'lets try this again...',
    'time_frame_list':['0_8', '16_24', '24_32', '32_40', '40_48', '48_56',
       '56_64', '64_72', '72_80', '80_88', '88_96'], 
    'plot' : True
    }

# paths = ['/media/bs001s/caf/CAF19_0326/0_3/probe1/co/scored_clayton.npy', '/media/bs001s/caf/CAF19_0326/0_3/probe1/co/scored_clayton.npy']
def lilo_and_stitch(paths, params, plot_shuffled=False):
    all_data = []
    all_ps  = []
    all_dccs = []

    quality = params['quality']
    ava_binsz=params['ava_binsz']
    perc = params['perc']
    hour_bins=params['hour_bins']
    animal = params['animal']
    
    for idx, path in enumerate(paths):
        print("------- WORKING ON ", path, " --------")

                
        params['time_frame'] = params['time_frame_list'][idx]
        params['total_time']=get_totaltime(params['time_frame'])
        print(f'time frame: {params["time_frame"]}      ----      total time: {params["total_time"]}')

        cells = np.load(path, allow_pickle=True)
        good_cells = [cell for cell in cells if cell.quality in quality and cell.cell_type in params['cell_type']]
        print(f"Number of cells: {len(good_cells)}")
        data = mbt.n_spiketimes_to_spikewords(good_cells, binsz=ava_binsz, binarize=1)

        if plot_shuffled:
            data_shuffled = Genshuffle(good_cells, params['total_time'], ava_binsz, perc, binary = 1, frame = 0)
        else:
            data_shuffled=None
       
        master_dict = looped_crit(data, data_shuffled, params, plot_shuffled=plot_shuffled)
        qual_str = '_'.join(map(str,quality))
        np.save(f'{params["saveloc"]}{animal}_dict_{params["time_frame"]}_startbin{str(params["start_bin"])}_{str(hour_bins)}hrs_perc{str(int(perc*100))}_binsz{str(int(ava_binsz*1000))}ms_q{qual_str}', master_dict)
    
        all_data.append(master_dict)
    

    return all_data




    