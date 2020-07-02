'''
%% PLPARAMS - automated computation of power law parameters using MLE
% Higher level macro and automated algorithm that computes the exponent for
% power-law distributed data and its uncertainty, the p-value and ks
% statistic for the fit, and searches for optimal support [xmin,xmax]. The
% function utilizes a "smart" greedy search algorithm to settle on a
% support pair. Prior to initiating the greedy search, all support pairs
% for which xmin and xmax differ by 1 are removed; the remainder are then
% sorted by log(xmax/xmin).
%
% Syntax: [tau, xmin, xmax, sigmaTau, p, pCrit, ks] = plparams(x, varargin)
%
% Input:
%   x (vector double) - random data that we would like to fit to the power
%     law distribution p(x) ~ x^(-tau) for x >= xmin and x <= xmax.
%
% Variable Inputs:
%   (..., 'samples', nSamples) - the number of sample distributions to draw.
%     Sets the resolution of the p-value (scalar double) (default: 500)
%   (..., 'threshold', pCrit) - for computational efficiency, a critical
%     p-value can be used to halt the computation if the likelihood of a
%     successful trial (given by the binomial distribution, see likelihood
%     below) drops below a pre-determined likelihood value. If pCrit is set
%     to 1, the computation will execute in full. Note: this only affects 
%     the greedy search process; final p-value is computed in full (scalar 
%     double) (default: .2)
%   (..., 'likelihood', likelihood) - likelihood threshold for binomial
%     process (scalar double) (default: 10^(-3))
%
% Outputs:
%   tau (scalar double) - slope of power law region
%   xmin (scalar double) - lower truncation of distribution
%   xmax (scalar double) - upper truncation of distribution
%   sigmaTau (scalar double) - error of the tau original fit estimated 
%     using the samples drawn from the model fit. If p is small, this error
%     is not valid.
%   p (scalar double) - proportion of sample distributions with KS
%     statistics larger than the KS statistic between the empirical pdf and
%     the model distribution. p is bounded by 0 and 1. A p-value of 1 
%     indicates that the KS statistic between the empirical pdf and the 
%     model was smaller than the KS statistic between the sample 
%     distributions and the model. Conversely, a p-value of 0 indicates 
%     that KS(simulated) > KS(empirical) for every sample.
%   pCrit (scalar double) - critical p-value used to truncate computation
%   ks (scalar double) - Kolmogorov-Smirnov statistic between the 
%     empirical pdf and the model
%
% Example:
%   x = bentpl(10000);
%     % generates random data distributed between two disjoint power law
%       regions with slopes 1.5 and 2.5 and first upper truncation of 50
%   [tau, xmin, xmax, sigma, p, pCrit, ks] = plparams(x);
%     % compute all power law parameters
%   [tau,~,~,~,p] = plparams(x, 'threshold', 1);
%     % compute power law exponent and full p-value
%   [~,~,~,~,p] = plparams(x, 'threshold', 1, 'samples', 1000);
%     % compute full p-value using 1000 simulated sets (increases
%       resolution of p-value)
%
% Other m-files required: PLMLE, PVCALC, MYMNRND
% Subfunctions: PLMLE, PVCALC
% MAT-files required: none
%
% See also: BENTPL, PLMLE, PVCALC, MYMNRND

% Author: Najja Marshall and Nick Timme
% Email: njm2149@columbia.edu and nicholas.m.timme@gmail.com
% May 2014; Last revision: 13-Apr-2016
%==============================================================================
% Copyright (c) 2014, The Trustees of Indiana University
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
%   1. Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
% 
%   2. Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
% 
%   3. Neither the name of Indiana University nor the names of its contributors
%      may be used to endorse or promote products derived from this software
%      without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


---------------------------------
adapted from NCC MATLAB toolbox - translated 6/26/2020 by Sahara Ensley

'''

import numpy as np
from itertools import combinations
from sahara_work import Criticiality_new as cr
def plparams(x, nSamples=500, pCrit=.05, likelihood=10e-3):
    #ensure data is oriented vertically
    nX = np.size(x)
    x_reshaped = np.reshape(x, [nX, 1]) # this line seems pointless as it's immediatly reverted when we take the unique array, but probably less pointless in matlab

    # finds all the combinations of avalanches, these will be the "min max" combos to test
    unqX = np.unique(x_reshaped)
    support = np.array(list(combinations(unqX, 2)))

    # removes all adjacent pairs of combinations
    for unique_elm in unqX:
        idx = np.where(support[:,0]==unique_elm)[0]
        if len(idx) > 0:
            first = idx[0]
            support = np.delete(support, first, axis=0)
    
    # compute 1/r = xmax/xmin for all support pairs
    rInv = support[:,1]/support[:,0]

    #get the amount of data covered by each support pair (nData)
    nSupport = np.size(support, 0)

    # computer and normalize natural log of 1/r
    lnRInv = np.log(rInv) / np.log(np.max(rInv))

    # rank support pairs by normalized ln(1/r)
    rank = np.power(lnRInv, 2)

    # sort support pairs by rank in descending order
    idxs = np.argsort(rank, axis=0)[::-1]
    support = support[idxs,:]

    ### greedy search for optimal support

    # try support pairs in ranked order until p = pCrit
    sweepFlag = True 
    iSupport = 0
    
    while sweepFlag and iSupport <= nSupport:
        if not iSupport % 10:
            print(f'{iSupport} pairs tried')
        current_pair = support[iSupport]
        xmin = current_pair[0]
        xmax = current_pair[1]

        # MLE of truncated distribution
        tau, _ , _ , _ = cr.plmle(x, xmin=xmin, xmax=xmax)

        # p-value for MLE
        p,_,_ = pvcalc(x, tau, xmin=xmin, xmax=xmax, samples=nSamples, threshold=pCrit, likelihood=likelihood)

        # halt search if p-value reaches critical p-value
        if p>= pCrit:
            sweepFlag=False
        else:
            iSupport = iSupport + 1
    
    # final statistics and full p-value
    p, ks, sigmaTau = pvcalc(x, tau, xmin=xmin, xmax=xmax, samples=nSamples, threshold=1)

    return tau, xmin, xmax, sigmaTau, p, pCrit, ks

