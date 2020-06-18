import numpy as np
from copy import deepcopy as cdc
import scipy.optimize
from scipy.stats import powerlaw
from sahara_work import Criticality_new as cr
def tplfit(burst,limit):
    """
        this code runs a loop to determine a good minima for your data to fit the most ideal powerlaw,
        it runs from 1 to your boundary and then chooses the smallest ks value, the minima at that point
        is your xmin

        only going to use for determining the minima, not anything else 
    """
    KS = []
    alpha = []
    Loglike = []
    alpha = []
    xmax = np.max(burst)

    for xmin in np.arange(1, limit+1):
        idx = np.where(np.logical_and(burst >= xmin, burst <= xmax))[0]
        n = np.size(burst[idx])
        s=np.unique(burst[idx])
        smin=np.min(s)
        smax=np.max(s)

        LL = lambda x: x*np.sum(np.log(burst[idx])) - n*np.log(1/np.sum(np.power(s,-x)))
        a = scipy.optimize.fmin(func=LL, x0=2.3, disp = False)
        fval = LL(a)
        Loglike.append(-fval)

        cdf = np.cumsum(np.histogram(burst[idx],bins = np.arange(xmin,xmax+2))[0]/n)
        A = 1/(np.sum(np.power(s, -a)))
        alpha.append(a)
        fit = np.cumsum(A*np.power(np.arange(xmin,xmax+1), -a))
        KS.append(np.max(np.abs(cdf-fit)))

    xmin = int(np.where(KS==np.min(KS))[0])
    alpha = alpha[xmin]
    Loglike = Loglike[xmin] 
    ks = np.min(KS)
    xmin=xmin+1

    return alpha, xmin, ks, Loglike