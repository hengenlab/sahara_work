import musclebeachtools_hlab.musclebeachtools as mbt
import numpy as np
import neuraltoolkit as ntk
import matplotlib.pyplot as plt
import glob

cells = np.load('/Volumes/Scratch/caf/test/clusters_scored/single_ref_clusters_scored.npy', allow_pickle=True)
cell = cells[19]

raw_dat = glob.glob(f'/Volumes/bs004r/D1/CAF00001_2020-06-22_12-28-05/*{cell.rstart_time}.bin')
t, load_dat = ntk.load_raw_binary_gain_chmap_range(raw_dat[0], 64, ['EAB50chmap_00'], te=(25000*30))
chdat = load_dat[cell.peak_channel, :]

data = ntk.butter_bandpass(chdat, 500, 10000, fs=25000)

large_amps = np.where(cell.spike_amplitude > 200)[0]
small_amps = np.where(cell.spike_amplitude < 100)[0]

lg_spk_dat = data[cell.spike_time[large_amps[0]]-50 : cell.spike_time[large_amps[0]]+50]
sm_spk_dat = data[cell.spike_time[small_amps[20]]-50 : cell.spike_time[small_amps[20]]+50]
sm_color = sm_spk_dat[45:56]
lg_color = lg_spk_dat[45:56]

fig = plt.figure(figsize=(15,6))
gs = fig.add_gridspec(2,3)

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15,2))

ax[1].plot(sm_spk_dat, color='grey')
ax[1].plot(np.arange(45,56),sm_color, color='green')

ax[2].plot(lg_spk_dat, color='grey')
ax[2].plot(np.arange(45,56), lg_color, color='r')

#ax[0].plot(cell.wf_b, alpha=0.1, color=[0.5, 0.5, 0.5])

ax[0].plot(cell.waveform_tetrodes)
ax[0].set_title('Cluster 12 Waveform')
ax[1].set_title("'Small' spike (<100mV)")
ax[2].set_title("'Large' spike (>200mV)")
