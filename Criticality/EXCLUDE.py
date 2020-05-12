import numpy as np
from sahara_work import Criticality as cr
from copy import deepcopy as cdc
def EXCLUDE(burst, setmin = 10, num = 1, flag = False):
# [burstMax, burstMin, alpha] = EXCLUDE(burst, varargin) is a function that
# determine both the lower and upper boundaries with a small KS.
# varargin return different inputs for setmin, num and flag.
                                                                                                                                            
# burst could be avalanche sizes or durations, which will be applied to the
# significant test and compared with power law distribution. 

# Copyright @Zhengyu Ma 2016
# Translated by Lizzie Tilden 5/29/2019	
# Edited by Yifan Xu 11/18/19
	KS = 1
	dKS = 1
	xmin = 1

	while KS > np.min([num/np.sqrt(np.size(burst[burst>xmin])), 0.1]) and dKS > 0.0005:
		alpha, xmin, ks, Loglike = cr.tplfit(burst, setmin)
		# print(xmin)
		# print(alpha)
		alpha = alpha[0]
		#xmin  = xmin[0]
		N = np.size(burst)
		xmax = np.max(burst)
		k = 0
		z = cdc(burst)
		z = z[z>=xmin]            
		n = np.size(z)
		cdf = np.cumsum(np.histogram(z, bins = np.arange(xmin,xmax+2))[0]/n)
		
		idx = np.where(np.logical_and(xmin<=burst, burst<=xmax))[0]
		s = np.unique(burst[idx])    
		smin = np.min(s) # minimum avalanche
		A = 1/np.sum(np.power(s, -alpha)) # constant factor for perfect power law
		fit = np.cumsum(A*np.power(np.arange(xmin, xmax+1), -alpha)) # CDF for perfect power law
		KS_old = KS # assign previous KS as KS_old
		KS = np.max(np.abs(cdf - fit)) # calculate the new KS
		dKS = np.abs(KS_old-KS) # difference between the current KS and the KS in the previous step

		burst = burst[burst < np.max(burst)] # lower the upper boundary by one since here we used '<' but not '<='
		burstMax = np.max(burst)
		# print(KS)

	burstMin = xmin
	# print('This is new')
	return burstMax, burstMin, alpha
