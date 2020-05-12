from musclebeachtools_hlab import musclebeachtools as mbt
import numpy as np
import os
from copy import deepcopy as cdc
import neuraltoolkit as ntk
import glob
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as colors
import matplotlib.cm as cmx
import Criticality as cr

def spiketimes_to_spikewords(spiketimes,startime,stoptime,binsize,binarize): 
    ### ARGUMENTS
    #spiketimes - list of spiketime arrays where each element of list is for a different neuron (e.g. the output of n_getspikes())
    #startime,stoptime - bounds of time to convert to spikewords (in seconds)
    #binsize - size of bin in milliseconds
    #binarize - '1' to convert to 0's and 1's, '0' to keep as histogram counts
    ### RETURNS
    #array of spikewords with each column as a cell and each rows as time in bins 
    
    #get sec_time to bin conversion factor
    #startime in bins
	startime_ms = startime * 1000
	stoptime_ms = stoptime * 1000
	binrange = np.arange(start = startime_ms,stop = stoptime_ms+1, step = binsize)
	n_cells = len(spiketimes)

	spikewords_array = np.zeros([n_cells,binrange.shape[0]-1])
	for i in range(n_cells):
		spiketimes_cell = np.asarray(spiketimes)[i] * 1000. #spiketimes in seconds * 1000 msec/sec
		counts, bins = np.histogram(spiketimes_cell,bins = binrange)
		if binarize == 1:
	        #binarize the counts
			counts[counts>0] = 1
        # print(counts.astype(np.int))
		spikewords_array[i,:] = counts
	return(spikewords_array.astype(np.int8))

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
	print('generate shuffled spikewords')

	perc = perc     

	r_shuffle = cr.AV_analysis_BurstT(data_shuffle, perc = perc)
	burst_shuffle = r_shuffle['S'] 
	duration_shuffle = r_shuffle['T'] 

	return burst_shuffle,duration_shuffle

def Plotshuffle(bm, tm, burst, duration, burst_shuffle, duration_shuffle):

	Result = cr.AV_analysis_ExponentErrorComments(burst, duration, bm, tm)

	burst_min = Result['xmin']
	burst_max = Result['xmax']
	duration_min = Result['tmin']
	duration_max =  Result['tmax']
	alpha = Result['alpha']
	beta = Result['beta']

	pdf = np.histogram(burst, bins = np.arange(1, np.max(burst)+2))[0]
	pdf_shuffle = np.histogram(burst_shuffle, bins = np.arange(np.min(burst_shuffle), np.max(burst_shuffle)+2))[0]
	tdf = np.histogram(duration,bins = np.arange(1, np.max(duration)+2))[0]
	tdf_shuffle = np.histogram(duration_shuffle,bins = np.arange(np.min(duration_shuffle), np.max(duration_shuffle)+2))[0]

	fig1, ax1 = plt.subplots(nrows = 1, ncols = 2, figsize = [12, 8])
	ax1[0].plot(np.arange(1,np.max(burst)+1),pdf/np.sum(pdf), marker = 'o', markersize = 3, linestyle = 'None', color = '#2138ab', alpha = 0.75)
	ax1[0].set_yscale('log')
	ax1[0].set_xscale('log')

	x = np.arange(burst_min, burst_max+1)
	y = (np.size(np.where(burst == burst_min+6)[0])/np.power(burst_min+6, -alpha))*np.power(x, -alpha) # calculate the fitted probability
	y = y/np.sum(pdf)

	x_shuffle = np.arange(np.min(burst_shuffle),np.max(burst_shuffle)+1)
	y_shuffle = pdf_shuffle/np.sum(pdf_shuffle)   # shuffled probability

	ax1[0].plot(x,y, color = '#c5c9c7',label = "data",linestyle = '--',linewidth = 2)
	ax1[0].plot(x_shuffle,y_shuffle,color = '#c5c9c7',linestyle = '-', label = "shuffle")
	ax1[0].set_xlabel('AVsize')
	ax1[0].set_ylabel('PDF(D)')
	ax1[0].set_title('AVsize PDF, ' + str(np.round(alpha[0], 3)))

	# plot duration distribution
	ax1[1].plot(np.arange(1,np.max(duration)+1),tdf/np.sum(tdf), marker = 'o', markersize = 3, linestyle = 'None', color = '#7bfdc7', alpha = 0.75)
	ax1[1].set_yscale('log')
	ax1[1].set_xscale('log')

	x = np.arange(duration_min, duration_max+1)
	y = (np.size(np.where(duration == duration_min+4)[0])/np.power(duration_min+4, -beta))*np.power(x, -beta)
	y = y/np.sum(tdf)

	x_shuffle = np.arange(np.min(duration_shuffle),np.max(duration_shuffle)+1)
	y_shuffle = tdf_shuffle/np.sum(tdf_shuffle)

	ax1[1].plot(x,y, color = '#c5c9c7',label = "data",linestyle = '--',linewidth = 2)
	ax1[1].plot(x_shuffle,y_shuffle,color = '#c5c9c7',linestyle = '-', label = "shuffle")
	ax1[1].set_xlabel('AVduration')
	ax1[1].set_ylabel('PDF(T)')
	ax1[1].set_title('AVduration PDF, ' + str(np.round(beta[0], 3)))
	
	# plt.legend()
	plt.show()

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

def check_SW_comp(nrn, state, mat_time, SW_binsz):
	''' This function finds the percentage of a certain time bin that is spent in a desired state.
		INPUTS:
		nrn: any neuron object from the desired animal
		state: desired state (1 = Wake, 2 = NREM, 3 = REM)
		mat_time: total amount of time that the avalanche matrix spans
		SW_binsz: size of bins ou want for assessing state concentration
		'''
	perc = np.zeros(int(mat_time/SW_binsz))
	for i, p in enumerate(np.arange(0, mat_time, SW_binsz)):
		time_idx = np.where(np.logical_and(nrn.behavior[0] > p, nrn.behavior[0] < p + SW_binsz))[0]
		if np.size(time_idx) == 0:
			s_idx = np.where(nrn.behavior[0] < p)[0]
			sleep_states = nrn.behavior[:, s_idx[-1]]
			desired_state = np.where(sleep_states == state)[0]
			perc[i] = np.size(desired_state)/np.size(sleep_states)
		else:
			sleep_states = nrn.behavior[:, time_idx]
			ts = sleep_states[0]
			ss = sleep_states[1]
			desired_state = np.where(ss == state)[0]
			if nrn.behavior[1,time_idx[-1]+1] == state:
				desired_state = np.append(desired_state, desired_state[-1]+1)
			ts = np.insert(ts, 0, p)
			ts = np.append(ts, p + SW_binsz)
			t_len = [int(ts[d]+1-ts[d]) for d in desired_state]
			total_len = np.sum(t_len)
			perc[i] = total_len/SW_binsz
	return perc

def bin_FR_mat(datadir, scrubbed, file_startclust = [], binsz = 40, qual = 1, rawdatadir = [], multi_probe = False, probeNumber = 1, start_block = 0, end_block = 1):
	''' INPUTS:
		datadir: t folder for desired clustering job
		binsz: desired binsize (ms)
		qual: what quality to include (if more than one enter as list)
		rawdatadir: data with sleep info if applicable
		multi_probe: True or False
		probeNumber: if multiprobe, which one
		start_block/end_block: if choosing a particular block from clustering
		OUTPUT:
		matrix with columns as different timepoints and rows as different neurons
		'''
	os.chdir(datadir)
	fileList = mbt.makeFileList(datadir, file_startclust, rawdatadir, multi_probe, start_block, end_block, probeNumber)
	if scrubbed:
		qual_list = np.load('scrubbed_quality_' + str(start_block) + '.npy')
	else:
		qual_list = fileList[5][0]

	if np.size(qual) > 1:
		cell_list = []
		for q in qual:
			cell_list.append(np.where(qual_list == q)[0])
		cell_list = np.concatenate(cell_list)
	else:
		cell_list = np.where(qual_list == qual)[0]
	cell_list = np.sort(cell_list)

	cells = []
	for a,i in enumerate(cell_list):
		print(a, ' ', i)
		cells.append(mbt.neuron(datadir, datatype='npy', clust_idx=i, file_list = fileList, rawdatadir = rawdatadir, start_block = start_block, end_block = end_block, probenumber = probeNumber))

	#make a function that spikes up waking spike times and sleeping spike times
	spks = mbt.getspikes(cells, 0, 3600*24) # gets spiketimes
	nrn_time = int(cells[0].time[-1]/3600)
	data_T = mbt.spiketimes_to_spikewords(spks,0,3600*nrn_time,binsz,1) #binarizes spikes
	data = data_T.T #transposes matrix so each column is a time bin and each row is a cell
	np.save('cr_mat.npy',data)

	return data, cells

def getspikes(neuron_list, starttime, stoptime):
	"""
	returns an array of all the spike times cells, each row is a new cell

	originally in spikeword_tools but moved here for the new MBT
	"""
	n_cells = np.shape(neuron_list)[0]
	spiketimes_allcells = list()
	for i in np.arange(n_cells):
		print('Getting spiketimes for cell ' + str(i))
		#get all spiketimes for cell
		spiketimes = neuron_list[i].spike_time
		spiketimes = spiketimes/25000
		spiketimes = spiketimes[(spiketimes>starttime)&(spiketimes<stoptime)]
		spiketimes_allcells.append(spiketimes)
	return (spiketimes_allcells)

def bin_FR_mat_MS(neurons, bin_size):
	""" 
	returns a matrix of spike times by cell based on bin size 

	neurons: list of mbt neurons. this list needs to be the final list of neurons you want to include
				there is no neuron currating done in this function. I would suggest running a scrubbing
				function prior to this and creating a numpy object with the cells you are including.
	bin_size: the size of the bin you're analyzing 
	"""

	spks = getspikes(neurons, 0, 3600*24) #gets spiketimes
	print("---got spiketimes---")
	nrn_time = int(neurons[0].spike_time[-1]/3600)
	print("---converted to neuron time---")
	data_T = spiketimes_to_spikewords(spks, 0, 3600*nrn_time, bin_size, 1) #binarizes spikes
	print("---binarized---")
	data = data_T.This

	return data


def align_motion(rawdatadir, motion_dir, file_startclust):
	os.chdir(rawdatadir)
	file_startdir =  np.sort(glob.glob('Head*.bin'))[0]
	t1 = ntk.load_raw_binary_gain_chmap(file_startdir, 512, 'PCB_tetrode', nprobes=8, t_only=1)
	t2 = ntk.load_raw_binary_gain_chmap(file_startclust, 512, 'PCB_tetrode', nprobes=8, t_only=1)
	ts_start = file_startclust.find('.bin')-19
	ts_end = file_startclust.find('.bin')-3
	# TS = file_startclust[ts_start:ts_end]
	# TS = TS.replace('_', 'T')
	# TS = TS.replace('-15', '*')
	# TS = TS.replace('-', '', 2)
	os.chdir(motion_dir)
	move_list = []
	time_list = []
	move = np.sort(glob.glob('*_tmove.npy'))
	for i in move:
		move_list.append(np.load(i)[0])
		time_list.append(np.load(i)[1])
	move_arr = np.concatenate(move_list)
	time_arr = np.concatenate(time_list)
	time_arrnsec = time_arr*3600*1e9
	offset = (t2-t1)
	start_idx = np.where(time_arrnsec > offset)[0][0]

	aligned_arr = move_arr[start_idx:]
	time = time_arr[start_idx:]
	aligned_time = time - time[0]
	tnsec = aligned_time*3600

	return tnsec, aligned_arr

def chk_FR(cells, binsz = 1800):
	fig1, ax1 = plt.subplots(ncols = 1, nrows = 1, figsize= [10, 10])
	plt.yticks(fontsize = 15)
	sns.despine()
	plt.xticks(fontsize = 15)
	sns.despine()
	plt.ion()

	cm = plt.get_cmap('coolwarm')
	values = np.arange(0, np.size(cells))
	cNorm = colors.Normalize( vmin = 0, vmax = values[-1])
	scalarMap = cmx.ScalarMappable( norm = cNorm, cmap = cm)

	for idx, cell in enumerate(cells):

		edges   = np.arange(0,max(cell.time),binsz)
		bins    = np.histogram(cell.time,edges)
		hzcount = bins[0]
		hzcount = hzcount/binsz
		hzcount[hzcount==0] = 'NaN'
		xbins   = bins[1]
		xbins   = xbins/binsz
		colorVal = scalarMap.to_rgba(idx)
		ax1.plot(xbins[:-1]/(3600/binsz),hzcount, color = colorVal )
	ylim = ax1.get_ylim()[1]
	ax1.set_xlim([0, xbins[-2]/(3600/binsz)])
	ax1.set_ylim([-5, ylim])
	ax1.set_ylabel('Firing Rate (Hz)', fontsize = 20)
	ax1.set_xlabel('Time (Hours)', fontsize = 20)







