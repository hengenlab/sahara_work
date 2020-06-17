import numpy as np
from copy import deepcopy as cdc
import numpy.ma as ma
def AV_analysis_BurstT(data, perc = 0.25):
# data is a matrix with each row as a neuron and each column is a
# time bin. Or data could be a vector (network activity)

# Result is a structure. Result.S is avalanche sizes, and Result.T
# is avalanche durations.

# The threshold for network activity is 25 percentile by default,
# but could be set up manually. When network activity is dominated 
# by silent period, perc could be zero. Otherwise, could try from roughly 20% to 50%.
# Threshold method is based on Poil, Simon-Shlomo, et al. "Critical-state
# dynamics of avalanches and oscillations jointly emerge from balanced
# excitation/inhibition in neuronal networks." Journal of Neuroscience
# 32.29 (2012): 9817-9823.

# Copyright @ Zhengyu Ma 2017
# Tranlsated by Lizzie Tilden:

# get get neurons# and frame#
	n, m = np.shape(data)

	############# get network activity #######################
	if n == 1:
		network = cdc(data)
	else:
		network = np.nansum(data, axis = 0)

	################ define threshold ##################

	if perc > 0:
		sortN      = np.sort(network)
		threshold2 = sortN[round(m*perc)]
	else:
		threshold2 = 0

	############

	zdata  = cdc(network)
	z2data = cdc(zdata)

	thresh_max = ma.masked_where(zdata<=threshold2, zdata)

	zdata[~thresh_max.mask] = 1 #avalanches
	zdata[thresh_max.mask] = 0 #intervals

	Z1 = np.where(~thresh_max.mask)[0] #avalanches
	Z  = np.where(thresh_max.mask)[0] #intervals

	z1data = cdc(zdata)

	z1data = np.delete(zdata, Z[np.where(np.diff(Z) == 1)[0]]) #use a single 0 to separate avalanches (a series of 1s)
	z0data = np.delete(zdata, Z1[np.where(np.diff(Z1) == 1)[0]]) #use 1 to separate intervals (some study focused on interval distributions)
	avalanches = np.delete(z2data, Z[np.where(np.diff(Z) == 1)[0]])
	avalanches[z1data == 0] = 0 #  use a single 0 to separate network activities in each avalanche

	J  = np.where(z1data == 0)[0]
	J1 = np.where(z0data == 1)[0]

	data2 = cdc(data)
	data2 = np.delete(data2, Z[np.where(np.diff(Z)==1)[0]], axis = 1)
	


	#################### Find Spike and AV sizes ########################
	burst = []
	for i in np.arange(0,np.size(J)-1):
		#print(J[i+1]-J[i]-2)
		fired = np.sum(avalanches[J[i]+1:J[i+1]]) - threshold2*(J[i+1]-J[i]-1)
		burst.append(fired)

	######### Duration Distribution ########################################
	T = np.diff(J)-1 # AVduration
	T1 = np.diff(J1)
	T = T[T>0] # Duration should be positive

	#################### Get final result ######################################
	Result = {}  
	Result['S'] = np.asarray(burst)
	Result['T'] = T

	return Result










