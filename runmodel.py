import numpy as np
import scipy
import scipy.io as sio
import h5py
import criticality_hlab.criticality as cr
import time


# perc = 0.0
# burstM = 10
# tM = 10
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


def get_data_from_sparse(name, idxs):
    f = h5py.File(name)
    datamat = f['Data']
    d = scipy.sparse.csc_matrix((datamat['data'], datamat['ir'], datamat['jc']))
    if idxs is not None:
        d = d[idxs,:]
    data = np.asarray(d.todense())
    return data

def get_fucking_burried(f, datname, row, col):
    name = h5py.h5r.get_name(f[datname][row][col], f.id)
    dat = f[name]
    d = scipy.sparse.csc_matrix((dat['data'], dat['ir'], dat['jc']))
    return d.toarray()

def get_data_normal(name, idxs):
    datamat = sio.loadmat(name)
    datamat = datamat['Data']
    data = scipy.sparse.csr_matrix(datamat)
    if idxs is not None:
        data = data[idxs,:]
    data = data.toarray()
    return data

def get_txt_data(name):
    data = []
    with open(name, 'r') as f:
        for line in f:
            data.append(int(float(line.strip())))

    data = np.asarray(data)
    return data


params = {
    'flag': 1,
    'pltname': 'testing',
    'saveloc': "/media/HlabShare/clayton_sahara_work/model_shit/figures/output_figs",
    'plot': True,
    'perc':0.35,
    'bm': 20,
    'tm': 5,
    'none_fact': 20,
    'exclude':False,
    'verbose':False,
    'nfactor_tm':0,
    'nfactor_bm':0,
    'nfactor_tm_tail':1,
    'nfactor_bm_tail':1,
    'exclude_burst':50,
    'exclude_time':20,
    'exclude_diff_b':20,
    'exclude_diff_t':10
}

# Result = cr.AV_analysis(x, y, params, nfactor_bm = params['nfactor_bm'], nfactor_tm = params['nfactor_tm'], 
#                         nfactor_bm_tail = params['nfactor_bm_tail'], nfactor_tm_tail = params['nfactor_tm_tail'], 
#                         none_fact = params['none_fact'], verbose = params['verbose'], exclude = params['exclude'], 
#                         exclude_burst = params['exclude_burst'], exclude_time = params['exclude_time'], 
#                         exclude_diff_b = params['exclude_diff_b'], exclude_diff_t = params['exclude_diff_t'])

DCC = np.zeros((19,4))

fulldata = sio.loadmat('Sub_Super/SubSummary.mat')
fulldata = fulldata['Data_sub']
for m in np.arange(1, 20):
    eig = 20-m
    print(f'EIG {eig}')
    for n in np.arange(1, 5):
        pltname = f'sub_eig{m}_{n}'
        params['pltname'] = pltname
        
        d = fulldata[n-1, eig-1]
        d = scipy.sparse.csr_matrix(d)
        data = d.toarray()
        r = cr.get_avalanches(data, perc = params['perc'])
        x = r['S']
        y = r['T']
        worked = True
        try:
            Result3 = cr.AV_analysis(x, y, params, nfactor_bm = params['nfactor_bm'], nfactor_tm = params['nfactor_tm'], 
                            nfactor_bm_tail = params['nfactor_bm_tail'], nfactor_tm_tail = params['nfactor_tm_tail'], 
                            none_fact = params['none_fact'], verbose = params['verbose'], exclude = params['exclude'], 
                            exclude_burst = params['exclude_burst'], exclude_time = params['exclude_time'], 
                            exclude_diff_b = params['exclude_diff_b'], exclude_diff_t = params['exclude_diff_t'])
        except Exception as err:
            print(err)
            worked = False
        if worked:
            DCC[m - 1, n - 1] = Result3['df']
            del Result3
        else:
            DCC[m - 1, n - 1] = np.nan

        del data

params = {
    'flag': 1,
    'pltname': 'testing',
    'saveloc': "/media/HlabShare/clayton_sahara_work/model_shit/figures/output_figs",
    'plot': True,
    'perc':0.35,
    'none_fact': 20,
    'exclude':False,
    'verbose':False,
    'bm':20,
    'tm':5,
    'nfactor_tm':0,
    'nfactor_bm':0,
    'nfactor_tm_tail':1,
    'nfactor_bm_tail':1,
    'exclude_burst':50,
    'exclude_time':20,
    'exclude_diff_b':20,
    'exclude_diff_t':10
}

DCC = np.zeros((19,4))

fulldata = h5py.File('Sub_Super/SuperSummary.mat')
for m in np.arange(1, 20):
    print(f'EIG {m}')
    for n in np.arange(1, 5):
        pltname = f'super_eig{m}_{n}'
        params['pltname'] = pltname
        data = get_fucking_burried(fulldata, 'Data_super', m-1, n-1)
        r = cr.get_avalanches(data, perc = params['perc'])
        x = r['S']
        y = r['T']
        worked = True
        try:
            Result3 = cr.AV_analysis(x, y, params, nfactor_bm = params['nfactor_bm'], nfactor_tm = params['nfactor_tm'], 
                            nfactor_bm_tail = params['nfactor_bm_tail'], nfactor_tm_tail = params['nfactor_tm_tail'], 
                            none_fact = params['none_fact'], verbose = params['verbose'], exclude = params['exclude'], 
                            exclude_burst = params['exclude_burst'], exclude_time = params['exclude_time'], 
                            exclude_diff_b = params['exclude_diff_b'], exclude_diff_t = params['exclude_diff_t'])
        except Exception as err:
            print(err)
            worked = False
        if worked:
            DCC[m - 1, n - 1] = Result3['df']
            del Result3
        else:
            DCC[m - 1, n - 1] = np.nan

        del data

#np.save('sub_10000_5mil_sparse_300subsampled.npy', DCC)


for m in np.arange(19, 20):
    print(f'EIG {m}')
    for n in np.arange(1, 5):
        name = "2000_cells/super_2000_eig" + str(m) + "_num" + str(n) + ".mat"
        pltname = "super_2000_long_eig" + str(m) + "_num" + str(n)
        params['pltname'] = pltname
        datamat = sio.loadmat(name)
        datamat = datamat['Data']
        data = scipy.sparse.csr_matrix.toarray(datamat)
        rand_nums = np.random.randint(low = 0, high = 2000, size = 1000)
        d = data[rand_nums, :]
        r = cr.get_avalanches(data, perc = perc)
        x = r['S']
        y = r['T']

        Result3 = AV_analysis(x, y, params)
        DCC[m - 1, n - 1] = Result3['df']
        del Result3
        del data
        # except:
        #     print(f'no avalanches for eig {m}')
        #     DCC[m-1,n-1] = np.nan
np.save('super_2000_long.npy', DCC)

# run mauricio's data

name = 'subsampled_n_100_g_1.485.txt'
data = get_txt_data(name)
r = cr.get_avalanches(data, perc = perc, ncells=100)
x = r['S']
y = r['T']

params['bm'] = int(np.max(x)/20)
params['tm'] = int(np.max(y)/20)
Result3 = AV_analysis(x, y, params, nfactor_bm_tail=0.8, nfactor_tm_tail=0.8)








# single file
name = '500_lambda11_long.mat'
datamat = sio.loadmat(name)
datamat = datamat['Data']
data = scipy.sparse.csr_matrix.toarray(datamat)
r = cr.get_avalanches(data, perc = perc) 
x = r['S']  
y = r['T'] 

Result3  = AV_analysis(x, y, params)




DCC = np.zeros(len(np.arange(100,10000,100)))
name = "10000_cells/sub_10000_5mil_sparse_eig3_num1.mat"


for i, n in enumerate(np.arange(10,10000, 100)):
    print(i)
    pltname = "num_subsample_subtest_"+str(n)+"_"
    params['pltname'] = pltname
    rand_nums = np.random.randint(low=0, high=10000, size=n)
    data = get_data_from_sparse(name, rand_nums)
    r = cr.get_avalanches(data, perc = perc)
    x = r['S']  
    y = r['T'] 
    try:
        Result3  = AV_analysis(x, y, params)
        DCC[i] = Result3['df']
    except ValueError:
        print('no avs')
        DCC[i] = np.nan
    except RuntimeError:
        print("optimize not successful - no avs")
        DCC[i] = np.nan

np.save("subsample_test_subcritical_10000_5mil_60.npy", DCC)
    

name = "matlab_files/super_subsampled.mat"
datamat = sio.loadmat(name)
datamat = datamat['Data']
data = scipy.sparse.csc_matrix.toarray(datamat)

for m in np.arange(1,20):
    print(f'EIG {m}')
    for n in np.arange(1,5):
        pltname = "super_subsample300_perc0_eig" + str(m) + "_num" + str(n)
        d = data[m+n, :]
        r = cr.AV_analysis_BurstT(d, perc = perc)
        x = r['S']
        y = r['T']
        Result3  = AV_analysis_new(x, y, burstM, tM, pltname, flag = 1, saveloc='/media/bs001s/caf/model_stuff/figures/output_figs/', plot=True)
        
        DCC[m-1,n-1] = Result3['df']
np.save('super_subsampled_spikes.npy',DCC)



# <CV> plot
mean_cvs_super = np.zeros(19)
for m in np.arange(1,20):
        print(f"eig {m}")
        tic = time.time()
        name = "10000_cells/super_10000_long_sparse_eig" + str(m) + "_num1" + ".mat"
        try:
            data = get_data_from_sparse(name)
        except OSError:
            data = get_data_normal(name)
        data = data.astype(np.int8)
        print('loaded')
        bool_spikes = (data == 1)
        spks = [np.where(x)[0] for x in bool_spikes]
        isi = [np.diff(x) for x in spks]
        cvs = [scipy.stats.variation(x) for x in isi]

        mean_cvs_super[m-1] = np.mean(cvs)
        toc = time.time()
        print(f'T: {tic-toc}')

#overnight code to run
mean_cvs_sub = np.zeros(19)
for m in np.arange(1,20):
        print(f"eig {m}")
        tic = time.time()
        name = "10000_cells/sub_10000_long_sparse_eig" + str(m) + "_num1" + ".mat"
        try:
            data = get_data_from_sparse(name)
        except OSError:
            data = get_data_normal(name)
        data = data.astype(np.int8)
        print('loaded')
        bool_spikes = (data == 1)
        spks = [np.where(x)[0] for x in bool_spikes]
        isi = [np.diff(x) for x in spks]
        cvs = [scipy.stats.variation(x) for x in isi]

        mean_cvs_sub[m-1] = np.mean(cvs)
        toc = time.time()
        print(f'T: {tic-toc}')
np.save("10000_cells_cvs_sub", mean_cvs_sub)




        

        
        
