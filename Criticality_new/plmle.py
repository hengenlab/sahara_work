# PLMLE - maximum likelihood estimation for power law distributed data
# Estimates slope of power law distributed data using the method of maximum
# likelihood.

# Syntax: [tau, xmin, xmax, L] = plmle(x, varargin)

# Input:
#   x (vector double) - random data to be fitted to the power law 
#     distribution p(x) ~ x^(-tau) for x >= xmin and x <= xmax. 

# Variable Inputs:
#   (..., 'xmin', xmin) - specifies the lower truncation of distribution
#     for the fit (scalar double) (default: min(x))
#   (..., 'xmax', xmax) - specifies the upper truncation of distribution
#     for the fit (scalar double) (default: max(x))
#   (..., 'tauRange', tauRange) - sets the range of taus to test for the
#     MLE fit (vector double) (default: [1, 5])
#   (..., 'precision', precision) - sets the decimal precision for the MLE
#     search (scalar double (power of ten)) (default: 10^-3)

# Outputs:
#   tau (scalar double) - slope of power law region
#   xmin (scalar double) - lower truncation of distribution
#   xmax (scalar double) - upper truncation of distribution
#   L (vector double) - log-likelihood that we wish to maximize for the MLE

# Example:
#   x = pldist(10^4);
#     % generates perfectly non-truncated power-law distributed data with
#     % slope of 1.5
#   tau = plmle(x)
#     % computes tau by MLE for x
#   tau = plmle(x, 'precision', 10^-4)
#     % computes tau to 4 decimal places
#   x = pldist(10^4, 'upper', 50, 1.5);
#     % generates perfectly power-law distributed data for x <= 50
#   tau = plmle(x, 'xmax', 50)
#     % computes tau for truncated region

# Other m-files required: none
# Subfunctions: none
# MAT-files required: none

# See also: BENTPL, PLPLOT, PLPARAMS

# Author: Najja Marshall and Nick Timme
# Email: njm2149@columbia.edu and nicholas.m.timme@gmail.com
# May 2014; Last revision: 8-Apr-2016

# ==============================================================================
# Copyright (c) 2014, The Trustees of Indiana University
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#   1. Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer.

#   2. Redistributions in binary form must reproduce the above copyright notice,
#      this list of conditions and the following disclaimer in the documentation
#      and/or other materials provided with the distribution.

#   3. Neither the name of Indiana University nor the names of its contributors
#      may be used to endorse or promote products derived from this software
#      without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# Tranlsated by Lizzie Tilden 6/03/2019
import sys
import numpy as np
import numpy.matlib
from sahara_work import Criticality as cr
def plmle(x, tauRange = [1,5], precision = 1e-3, xmin = False, xmax = False):
	
	if xmin == False:
		xmin = np.min(x)
	if xmax == False:
		xmax = np.max(x)
	
	# Error check the precision
	if np.log10(precision) != np.round(np.log10(precision)):
		sys.exit('The precision must be a power of ten.')

	### Process data
	# Reshape
	x = np.reshape(x, [np.size(x), 1])

	# Determine data type
	if  np.count_nonzero(np.abs(x-np.round(x)) > 3*np.spacing(1)) >0:
		dataType = 'CONT'
	else:
		dataType = 'INTS'
		x = np.round(x)
		# Truncate
	x_idx = np.where(np.logical_and(x >= xmin, x <= xmax))[0]
	z = x[x_idx]
	unqZ = np.unique(z)
	nZ = np.size(z)
	nUnqZ = np.size(unqZ)
	allZ = np.arange(xmin,xmax+1)
	nallZ = np.size(allZ)

	## MLE calculation
	r = xmin/xmax

	nIterations = -np.log10(precision)   

	for iIteration in np.arange(1, nIterations+1):
		spacing = 1*(10**-(iIteration)) 
		if iIteration == 1:
			taus = np.arange(tauRange[0], tauRange[1]+spacing, spacing)
		else:
			if tauIdx == 1:
				taus = np.arange(taus[0], taus[-1]+spacing, spacing)
			elif tauIdx == np.size(taus):
				taus = np.arange(taus[-2], taus[-1]+spacing, spacing)
			else:
				taus = np.arange(taus[tauIdx-1], taus[tauIdx+1]+spacing, spacing)
		nTaus = np.size(taus)

		if dataType == 'INTS':
			# replicate arrays to equal size
			allZMat = np.matlib.repmat(np.reshape(allZ,[nallZ,1]),1, nTaus)
			tauMat = np.matlib.repmat(taus,nallZ,1)

			# compute the log-likelihood function
			L = -np.log(np.sum(np.power(allZMat,(-tauMat)), 0))-(taus/nZ) * np.sum(np.log(z))

		elif dataType == 'CONT':
			# compute the log-likelihood function (method established via
			# Deluca and Corral 2013)
			L = np.log((taus-1)/(1-np.power(r,(taus - 1))))-taus*(1/nZ)*np.sum(np.log(z))-(1-taus)*np.log(xmin);

			if 1 in taus:
				L[taus == 1] = -np.log(np.log(1/r))-(1/nZ)*np.sum(np.log(z))

		tauIdx = np.where(L == np.max(L))[0]
	tau = taus[tauIdx]

	return tau, xmin, xmax, L












