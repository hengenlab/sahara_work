import numpy as np
from sahara_work import Criticality_new as cr
from copy import deepcopy as cdc


def EXCLUDE(burst, setmin, num=1):
    KS = 1
    dKS = 1
    xmin = 1

    while KS > np.min([num/np.sqrt(np.size(burst[burst>xmin])), 0.1]) and dKS > 0.0005:
        alpha, xmin, ks, Loglike = cr.tplfit(burst, setmin)
        alpha = alpha[0]
        xmax = np.max(burst)

        z = burst[burst>=xmin]
        n = np.size(z)
        cdf = np.cumsum(np.histogram(z, bins = np.arange(xmin,xmax+2))[0]/n)

        idx = np.where(np.logical_and(xmin<=burst, burst<=xmax))[0]
        s = np.unique(burst[idx]) 
        A = 1/np.sum(np.power(s, -alpha))
        fit = np.cumsum(A*np.power(np.arange(xmin, xmax+1), -alpha))

        KS_old = cdc(KS)

        KS = np.max(np.abs(cdf - fit))
        dKS = np.abs(KS_old-KS)
        burst = burst[burst < np.max(burst)]
        burstMax = np.max(burst)

    burstMin = xmin
    return burstMax, burstMin, alpha
