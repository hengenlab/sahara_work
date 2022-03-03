from pickle import load
import numpy as np
import scipy
import scipy.io as sio
import h5py
import glob
from copy import deepcopy as cdc
import criticality as cr
import matplotlib.pyplot as plt


# copy this function before using code
def plot_shuff(fig, burst, T):
    # burst PDF
    pdf = np.histogram(burst, bins = np.arange(1, np.max(burst) + 2))[0]
    p = pdf / np.sum(pdf)
    #xlimburst = fig.axes[0].get_xlim()
    fig.axes[0].plot(np.arange(1, np.max(burst) + 1), p, color = 'darkgrey', alpha = 0.75)
    #fig.axes[0].set_xlim(xlimburst)
    # ax1[0].set_yscale('log')
    # ax1[0].set_xscale('log')

    
    # ax1[0].set_xlabel('AVsize')
    # ax1[0].set_ylabel('PDF(S)')
    # ax1[0].set_title('AVsize PDF, ' + str(np.round(alpha, 3)))

    # time pdf
    tdf = np.histogram(T, bins = np.arange(1, np.max(T) + 2))[0]
    t = tdf / np.sum(tdf)
    fig.axes[1].plot(np.arange(1, np.max(T) + 1), t, color = 'darkgrey', alpha = 0.75)
    
    # ax1[1].set_yscale('log')
    # ax1[1].set_xscale('log')
    # sns.despine()

    # ax1[1].set_xlabel('AVduration')
    # ax1[1].set_ylabel('PDF(D)')
    # ax1[1].set_title('AVdura PDF, ' + str(np.round(beta, 3)))


    # figure out how to plot shuffled data

    # plt.tight_layout()
    # plt.legend()
    # plt.savefig(saveloc + "/" + pltname + 'scaling_relations' + '.svg', format='svg')

    return fig

# loading function to load the matlab data
def load_human(filename):
    f = h5py.File(filename)
    for k in f.keys():
        n = k
    data = np.array(f[n])
    return data

subsample = True
smooth = False
binsize = 5
nsub = 100
perc = 0.25

def check():
    files = glob.glob('/media/HlabShare/AD_paper/human/DIAN/MN/*.mat')
    all_lens = []
    for f in files:
        dat = load_human(f)
        len = dat.shape[0]
        all_lens.append(len)
    

def DIAN(file):
    subsample = True
    smooth = False
    binsize = 5
    nsub = 100
    n_its = 1000
    perc = 0.25


    mask = load_human('/media/HlabShare/AD_paper/human/MSC/submask.mat')
    networks = mask[0]
    idxs = np.where(networks == 0)[0]
    data = load_human(file)
    d = data[:,idxs]
    D = d.T

    R = {'S' : [], 'T' : []}
    Rs = {'S' : [], 'T' : []}
    for n in range(n_its):
        if subsample:
            voxs = np.random.randint(low=0, high = D.shape[0]-1, size=nsub)
            D = D[voxs]

        D_bin = []
        D_bin_shuff = []
        for v in D:
            std = np.std(v) * 1.5
            temp = cdc(v) 

            if smooth:
                idxs = np.where(v>=std)[0]
                binrange = np.arange(0, (len(v) + binsize), binsize)
                counts, bins = np.histogram(idxs, bins=binrange) 
                counts[counts > 0] = 1
                D_bin.append(counts)

                shuff_idxs = np.random.randint(low=0, high=len(temp), size=int(len(idxs)))
                countsS, binsS = np.histogram(shuff_idxs, bins=binrange) 
                countsS[countsS > 0] = 1
                D_bin_shuff.append(countsS) 

            else:
                temp[v >= std] = 1 
                temp[v < std] = 0 
                D_bin.append(temp)

                numspks = sum(temp)
                shuff_idxs = np.random.randint(low=0, high=len(temp), size=int(numspks))
                temp_shuff = np.zeros(shape=len(temp))
                temp_shuff[shuff_idxs] = 1
                D_bin_shuff.append(temp_shuff) 

        r = cr.get_avalanches(D_bin, perc = perc) 
        burst = r['S'] 
        T = r['T'] 

        R['S'].append(burst)
        R['T'].append(T)

        rs = cr.get_avalanches(D_bin_shuff, perc = perc) 
        bursts = rs['S'] 
        Ts = rs['T'] 
        Rs['S'].append(bursts)
        Rs['T'].append(Ts)

    R['S'] = np.concatenate(R['S'])     
    R['T'] = np.concatenate(R['T'])
    Rs['S'] = np.concatenate(Rs['S'])
    Rs['T'] = np.concatenate(Rs['T']) 

    Result = cr.AV_analysis(R['S'], R['T'], flag = 1, nfactor_bm_tail = 1, nfactor_tm_tail = 1, nfactor_bm = 0, bm = 15,  tm = 20,\
                            pltname = 'dian_test', saveloc = '/media/HlabShare/AD_paper/human/', plot = True)


    fig = plot_shuff(Result['scaling_relation_plot'], Rs['S'], Rs['T'])
    fig.savefig(f'DIAN_network_{n}.jpeg')
    

def MSC():
    files = glob.glob('/media/HlabShare/AD_paper/human/MSC/msc*.mat')
    mask = load_human('/media/HlabShare/AD_paper/human/MSC/submask.mat')
    networks = mask[0]

    for n in np.unique(networks):
        idxs = np.where(networks == n)[0]

        dfull = []
        for f in files:
            data = load_human(f)
            d = data[:,idxs]
            dfull.append(d)
        DATA = np.concatenate(dfull)
        D = DATA.T  # SO! The current state of the data is that it's only the voxels from network 4 (who knows what network that is), and its VOXEL x TIME (so like spiketimes)

        if subsample:
            voxs = np.random.randint(low=0, high = D.shape[0]-1, size=nsub)
            D = D[voxs]

        D_bin = []
        D_bin_shuff = []
        for v in D:
            std = np.std(v) * 2.5
            temp = cdc(v) 

            if smooth:
                idxs = np.where(v>=std)[0]
                binrange = np.arange(0, (len(v) + binsize), binsize)
                counts, bins = np.histogram(idxs, bins=binrange) 
                counts[counts > 0] = 1
                D_bin.append(counts)

                shuff_idxs = np.random.randint(low=0, high=len(temp), size=int(len(idxs)))
                countsS, binsS = np.histogram(shuff_idxs, bins=binrange) 
                countsS[countsS > 0] = 1
                D_bin_shuff.append(countsS) 

            else:
                temp[v >= std] = 1 
                temp[v < std] = 0 
                D_bin.append(temp)

                numspks = sum(temp)
                shuff_idxs = np.random.randint(low=0, high=len(temp), size=int(numspks))
                temp_shuff = np.zeros(shape=len(temp))
                temp_shuff[shuff_idxs] = 1
                D_bin_shuff.append(temp_shuff) 

        
        try:
            R = cr.get_avalanches(D_bin, perc = perc) 
            burst = R['S'] 
            T = R['T'] 

            Rs = cr.get_avalanches(D_bin_shuff, perc = perc) 
            bursts = Rs['S'] 
            Ts = Rs['T'] 
            
            Result = cr.AV_analysis(burst, T, flag = 1, nfactor_bm_tail = 0.7, nfactor_tm_tail = 0.8, nfactor_bm = 5, bm = 15,  tm = 20,\
                                    pltname = 'msc_test', saveloc = '/media/HlabShare/AD_paper/human/', plot= True)


            fig = plot_shuff(Result['scaling_relation_plot'], bursts, Ts)
            fig.savefig(f'msc_network_{n}.jpeg')
        except Exception as err:
            print(f'Error on network {n} --- {err}')
