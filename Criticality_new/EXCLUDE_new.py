import numpy as np
from sahara_work import Criticality_new as cr
from copy import deepcopy as cdc
import powerlaw as plaw 
# def EXCLUDE(burst, num = 1):
#     KS = 1
#     dKS=1
#     xmin=1

#     while KS > np.min([num/np.sqrt(np.size(burst[burst>xmin])), 0.1]) and dKS > 0.0005:
#         fit = plaw.Fit(burst) # use package to determine the powerlaw fit and all the statistics
#         alpha = fit.alpha # shape of powerlaw
#         xmin = fit.xmin # systematically determined min based on minimizing the ks value

#         cdf = fit.cdf() # cumulative distribution function for our data
#         fit_cdf = fit.power_law.cdf() # cumulative distribution function for the powerlaw fit
        
#         KS_old = KS # save og KS value
#         KS = fit.power_law.KS() # find KS value from this fit
#         dKS = np.abs(KS_old-KS) # difference between the two

#         if dKS > 0.0005:
            
#             burst = burst[burst < np.max(burst)] # drops the max by one each time because its < not <=
#             burstMax = np.max(burst) # save the max value
#         else:
#             print("quitting")

#     burstMin = xmin # save the min value

#     return burstMax, burstMin, alpha

def EXCLUDE(burst, setmin, num=1):
    KS = 1
    dKS = 1
    xmin = 1



    while KS > np.min([num/np.sqrt(np.size(burst[burst>xmin])), 0.1]) and dKS > 0.0005:
        alpha, xmin, ks, Loglike = cr.tplfit(burst, setmin)
        alpha = alpha[0]
        N = np.size(burst)
        xmax = np.max(burst)
        k=0
        z = burst[burst>=xmin]
        n= np.size(z)
        cdf = np.cumsum(np.histogram(z, bins = np.arange(xmin,xmax+2))[0]/n)

        idx = np.where(np.logical_and(xmin<=burst, burst<=xmax))[0]
        s = np.unique(burst[idx]) 
        smin = np.min(s)
        A = 1/np.sum(np.power(s, -alpha))
        fit = np.cumsum(A*np.power(np.arange(xmin, xmax+1), -alpha))
        KS_old = KS
        KS = np.max(np.abs(cdf - fit))
        dKS = np.abs(KS_old-KS)
        burst = burst[burst < np.max(burst)]
        burstMax = np.max(burst)

    burstMin = xmin
    return burstMax, burstMin, alpha
