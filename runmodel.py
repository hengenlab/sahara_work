import numpy as np
import scipy
import scipy.io as sio
from sahara_work import Criticality_new as cr
DCC = np.zeros((19,4))
perc = 0.3       
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


for m in np.arange(1,2):
    print(f'EIG {m}')
    for n in np.arange(1,5):
        
        name = "super_eig" + str(m) + "_num" + str(n) + ".mat"
        pltname = "super_eig_tmax" + str(m) + "_num" + str(n)
        datamat = sio.loadmat(name)
        datamat = datamat['Data']
        data = scipy.sparse.csr_matrix.toarray(datamat)
        r = cr.AV_analysis_BurstT(data, perc = perc)
        x = r['S']  
        y = r['T'] 
        Result3  = AV_analysis_new(x, y, burstM, tM, pltname, flag = 1, saveloc='/media/bs001s/caf/model_stuff/', plot=True) 
        
        DCC[m-1,n-1] = Result3['df']
np.save('super_full_tmax.npy',DCC)