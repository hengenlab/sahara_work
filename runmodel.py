import numpy as np
import scipy
import scipy.io as sio
from sahara_work import Criticality_new as cr
DCC = np.zeros((19,4))
perc = 0.3       
burstM = 10
tM = 5
for m in np.arange(1,20):
    for n in np.arange(1,5):
        
        name = "sub_eig" + str(m) + "_num" + str(n) + ".mat"
        pltname = "sub_eig" + str(m) + "_num" + str(n)
        datamat = sio.loadmat(name)
        datamat = datamat['Data']
        data = scipy.sparse.csr_matrix.toarray(datamat)
        r = cr.AV_analysis_BurstT(data, perc = perc)
        x = r['S']  
        y = r['T'] 
        Result3, ax3 = AV_analysis_new(x, y, burstM, tM, pltname, saveloc='', plot=True) 
        
        DCC[m-1,n-1] = Result3['df']
np.save('subdcc.npy',DCC)