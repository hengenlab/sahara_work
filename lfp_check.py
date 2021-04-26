import scipy.signal as signal
import matplotlib.pyplot as plt
import numpy as np
from sys import platform
import matplotlib.backends.backend_pdf as mpdf
import time

def lfp_check(files, fs, num_checks, plt_dir = '', animal=''):
    if platform == "darwin":
        plt.switch_backend('TkAgg')
    else:
        plt.switch_backend('Agg')
    indx = np.random.randint(low=0, high=len(files), size=num_checks)
    count = 0
    num_plots = int(np.ceil((num_checks/4)))
    for i in range(num_plots):
        fig, ax = plt.subplots(nrows = 4, ncols = 1, figsize=(16,10))

        for row in ax:
            if count < num_checks:
                fl = files[indx[count]]
                eeg = np.load(fl)
                downdatlfp = np.mean(eeg[:,0:3600*fs], axis = 0)   # generate spectrogram from the 1st hour
                f, t_spec, x_spec = signal.spectrogram(downdatlfp, fs=fs, window='hanning', nperseg=1000, noverlap=1000-1, mode='psd')
                #f, t_spec, x_spec = signal.spectrogram(downdatlfp, fs=fs, window='hanning', nfft=400, detrend=False, noverlap=200, mode='psd') # lower quality but faster
                fmax = 64
                fmin = 1
                x_mesh, y_mesh = np.meshgrid(t_spec, f[(f<fmax) & (f>fmin)])
                row.pcolormesh(x_mesh, y_mesh, np.log10(x_spec[(f<fmax) & (f>fmin)]), cmap='jet')
                del(x_mesh)
                del(y_mesh)
                row.set_ylim(1,64)
                row.set_xlim(0,3600)
                row.set_yscale('log')
                row.set_title(f'{animal}--{fl}')
                count+=1
        plt.tight_layout()
        fig.savefig(f'{plt_dir}/LFP_check_{animal}_{i}.png')