import numpy as np
import scipy
import scipy.io as sio
from sahara_work import Criticality_final as cr

DCC = np.zeros((19,4))
perc = 0.0
burstM = 10
tM = 5
# m=1
# n=1
# name = "sub_eig" + str(m) + "_num" + str(n) + ".mat"
# pltname = "sub_eig" + str(m) + "_num" + str(n)
# datamat = sio.loadmat(name)
# datamat = datamat['Data']
# data = scipy.sparse.csr_matrix.toarray(datamat)
# r = cr.AV_analysis_BurstT(data, perc = perc)
# x = r['S']  
# y = r['T'] 
# Result3 = cr.AV_analysis_new(x, y, burstM, tM, pltname, saveloc='/media/bs001s/caf/model_stuff/', plot=True) 

# DCC[m-1,n-1] = Result3['df']
params = {
    'flag': 1,
    'bm': None,
    'tm': None,
    'pltname': "testing_ncc",
    'saveloc': "/media/bs001s/caf/model_stuff/figures/output_figs/",
    'burst_shuffled': None,
    'T_shuffled': None,
    'plot_shuffled': False,
    'plot': True
}

for m in np.arange(1,19):
    print(f'EIG {m}')
    for n in np.arange(1,5):
        
        name = "10000_cells/super_10000_long_sparse_eig" + str(m) + "_num" + str(n) + ".mat"
        pltname = "super_10000_long_sparse_eig" + str(m) + "_num" + str(n)
        params['pltname'] = pltname
        datamat = sio.loadmat(name)
        datamat = datamat['Data']
        data = scipy.sparse.csr_matrix.toarray(datamat)
        # rand_nums = np.random.randint(low=0, high=2000, size=1000)
        # d = data[rand_nums, :]
        r = cr.get_avalanches(data, perc = perc)
        x = r['S']  
        y = r['T'] 

        Result3  = AV_analysis(x, y, params)
        DCC[m-1,n-1] = Result3['df']
        del Result3
        del data
        # except:
        #     print(f'no avalanches for eig {m}')
        #     DCC[m-1,n-1] = np.nan
np.save('super_10000_long_spaese.npy',DCC)

# single file
name = '2000_cells/super_2000_long_eig1_num2.mat'
datamat = sio.loadmat(name)
datamat = datamat['Data']
data = scipy.sparse.csr_matrix.toarray(datamat)
r = cr.get_avalanches(data, perc = perc) 
x = r['S']  
y = r['T'] 

Result3  = AV_analysis(x, y, params)




DCC = np.zeros(len(np.arange(10,2000,10)))
name = "2000_cells/sub_2000_long_eig1_num1.mat"
datamat = sio.loadmat(name)
datamat = datamat['Data']
data = scipy.sparse.csr_matrix.toarray(datamat)

for i, n in enumerate(np.arange(10,2000, 10)):
    print(i)
    pltname = "num_subsample_subtest_"+str(n)+"_"
    params['pltname'] = pltname
    rand_nums = np.random.randint(low=0, high=2000, size=n)
    d = data[rand_nums, :]
    r = cr.get_avalanches(d, perc = perc) 
    x = r['S']  
    y = r['T'] 
    try:
        Result3  = AV_analysis(x, y, params)
        DCC[i] = Result3['df']
    except ValueError:
        print('no avs')
        DCC[i] = np.NaN
    

# name = "matlab_files/super_subsampled.mat"
# datamat = sio.loadmat(name)
# datamat = datamat['Data']
# data = scipy.sparse.csc_matrix.toarray(datamat)

# for m in np.arange(1,20):
#     print(f'EIG {m}')
#     for n in np.arange(1,5):
#         pltname = "super_subsample300_perc0_eig" + str(m) + "_num" + str(n)
#         d = data[m+n, :]
#         r = cr.AV_analysis_BurstT(d, perc = perc)
#         x = r['S']  
#         y = r['T'] 
#         Result3  = AV_analysis_new(x, y, burstM, tM, pltname, flag = 1, saveloc='/media/bs001s/caf/model_stuff/figures/output_figs/', plot=True) 
        
#         DCC[m-1,n-1] = Result3['df']
# np.save('super_subsampled_spikes.npy',DCC)



# <CV> plot
mean_cvs_sub = np.zeros(19)
for m in np.arange(1,20):
        print(f"eig {m}")
        name = "5000_cells/sub_5000_eig" + str(m) + "_num1" + ".mat"
        datamat = sio.loadmat(name)
        datamat = datamat['Data']
        data = scipy.sparse.csr_matrix.toarray(datamat)

        spks = [np.where(x==1)[0] for x in data]
        isi = [np.diff(x) for x in spks]
        cvs = [scipy.stats.variation(x) for x in isi]

        mean_cvs_sub[m-1] = np.mean(cvs)



        

        
        
np.save('super_2000_long.npy',DCC)