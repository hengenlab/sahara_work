import numpy as np
from copy import deepcopy as cdc
import scipy.optimize
from sahara_work import Criticality as cr
def tplfit(burst,limit):
	KS = []
	alpha = []
	Loglike = []
	alpha = []

	for x0 in np.arange(1, limit+1):

		X = cdc(burst)
		xmax = np.max(burst)
		idx = np.where(np.logical_and(x0<=X, X<=xmax))[0]
		n = np.size(burst[idx])
		s = np.unique(X[idx]) 
		smin = np.min(s) # minimum avalanche
		smax = np.max(s) # maximum avalanche

		LL = lambda x: x*np.sum(np.log(burst[idx])) - n*np.log(1/np.sum(np.power(s,-x)))
		a = scipy.optimize.fmin(func=LL, x0=2.3, disp = False) # start search from 2.3
		fval = LL(a) #value of the minimum likelihood
		Loglike.append(-fval)

		N = np.size(burst)
		k = 0

		z = cdc(burst)
		z = z[z>=x0]
		n = np.size(z)
		cdf = np.cumsum(np.histogram(z,bins = np.arange(x0,xmax+2))[0]/n) # cdf
		A = 1/(np.sum(np.power(s, -a)))
		alpha.append(a)
		fit = np.cumsum(A*np.power(np.arange(x0,xmax+1), -a))
		KS.append(np.max(np.abs(cdf-fit)))

	xmin = int(np.where(KS==np.min(KS))[0])
	alpha = alpha[xmin]
	Loglike = Loglike[xmin];
	ks = np.min(KS)
	xmin = xmin + 1

	return alpha, xmin, ks, Loglike
