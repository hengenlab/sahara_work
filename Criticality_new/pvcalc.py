'''
%% PVCALC - compute p-value for the pdf fit
% Calculates the p-value by Monte Carlo. The method uses a raw pdf with a 
% model that assumes a power law distribution with exponent tau for the 
% section of the data bounded by support [xmin,xmax]. The function 
% generates many sample distributions using a model that is composed of a 
% power-law within [xmin,xmax]. p is the proportion of sample distributions 
% with KS statistics larger than the KS statistic between the original pdf 
% and the model distribution. For computational efficiency, the function 
% continually updates the likelihood of successful results using the 
% binomial distribution and halts for statistically unlikely results.
%
% Syntax: [p, ks, sigmaTau] = pvcalc(x, tau, varargin)
%
% Inputs:
%   x (vector double) - empirical data, assumed to be power-law distributed
%   tau (scalar double) - the exponent of the power-law distribution
%
% Variable Inputs:
%   (..., 'xmin', xmin) - sets lower truncation of distribution (scalar
%     double) [default: min(x)]
%   (..., 'xmax', xmax) - sets upper truncation of distribution (scalar
%     double) [default: max(x)]
%   (..., 'samples', nSamples) - the number of sample distributions to draw.
%     Sets the resolution of the p-value (scalar double) (default: 500)
%   (..., 'threshold', pCrit) - for computational efficiency, a critical
%     p-value can be used to halt the computation if the likelihood of a
%     successful trial (given by the binomial distribution) drops below a
%     pre-determined likelihood value. If pCrit is set to 1, the 
%     computation will execute in full (scalar double) (default: 1)
%   (..., 'likelihood', likelihood) - likelihood threshold for binomial
%     process (scalar double) (default: 10^(-3))
%
% Outputs:
%   p (scalar double) - proportion of sample distributions with KS
%     statistics larger than the KS statistic between the empirical pdf and
%     the model distribution. p is bounded by 0 and 1. A p-value of 1 
%     indicates that the KS statistic between the empirical pdf and the 
%     model was smaller than the KS statistic between the sample 
%     distributions and the model. Conversely, a p-value of 0 indicates 
%     that KS(simulated) > KS(empirical) for every sample.
%   ks (scalar double) - Kolmogorov-Smirnov statistic between the 
%     empirical pdf and the model
%   sigmaTau (scalar double) - error of the tau original fit estimated 
%     using the samples drawn from the model fit. If p is small, this error
%     is not valid. 
%
% Example:
%   x = pldist(10^4);
%     % generate non-truncated, power-law distributed data
%   tau = plmle(x);
%   [p, ks, sigma] = pvcalc(x, tau)
%     % compute (truncated) p-value, KS statistic, and error on power law 
%     % exponent
%   p = pvcalc(x, tau, 'threshold', 1)
%     % compute full p-value
%   p = pvcalc(x, tau, 'likelihood', 10^-5)
%     % compute p-value with threshold for binomial process set to 10^-5
%
% Other m-files required: MYMNRND
% Subfunctions: MYMNRND
% MAT-files required: none
%
% See also: PLDIST, PLMLE, MYMNRND

% Author: Najja Marshall and Nick Timme
% Email: njm2149@columbia.edu and nicholas.m.timme@gmail.com
% November 2013; Last revision: 1-Jun-2015

% Version Information
%
%   1.0: 11/4/13 - Creation of the original program. (Nick Timme)
%
%   2.0: 11/8/13 - Modification of the ks statistic calculation algorithm.
%   Also, the pdf is cut to only include the region between xmin and xmax.
%   (Nick Timme)
%
%   3.0: 6/1/15 - Addition of exponential modification functionality. (Nick
%   Timme)
%

%==============================================================================
% Copyright (c) 2013, The Trustees of Indiana University
% All rights reserved.
% 
% Authors: Nick Timme (nicholas.m.timme@gmail.com)
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



--------------------------
translated from NCC MATLAB toolbox to python - Sahara Ensley 6/26/2020
'''
import numpy as np
import scipy.stats
def pvcalc(x, tau, xmin=-1, xmax=-1, nSamples=500, pCrit = .05, likelihood = 10e-3)
    if xmin < 0 :
        xmin = np.min(x)
    if xmax < 0:
        xmax = np.max(x)
    
    # ensure data is vertical 
    x = np.reshape(x, [np.size(x), 1])

    diff = list(set(x) - set(np.floor(x)))

    if len(diff) == 0: # integer data

        # calculate pdf of data using unique bins
        pdfX = np.histogram(x, bins=np.arange(xmin, xmax+1))[0]
        nSupportEvents = sum(pdfX)
        pdfX = pdfX / nSupportEvents

        # calculate pdf of fit
        pdfFit = np.arange(xmin, xmax+1)**-tau
        pdfFit = pdfFit / np.sum(pdfFit)

        # caclulate cdfs of empirical data and spliced data
        cdfX = 1 - np.cumsum(pdfX)
        cdfFit = 1 - np.cumsum(pdfFit)

        # ensure cdfs are oriented vertically 
        cdfX = np.reshape(cdfX, [np.size(cdfX), 1])
        cdfFit = np.reshape(cdfFit, [np.size(cdfFit), 1])
    
    else: # continious data
        idxs = np.where(np.logical_and(x>xmin, x<xmax))[0]
        x = x[idxs]

        #sort values of data
        sortedx = np.sort(x)

        #calculate the data cdf
        cdfX = np.arange(1, np.size(x)) / np.size(x)

        # calculate the y value of the fit cdf
        cdfFit = ( (1 - (sortedx ** (1-tau))) - (1 - xmin**(1-tau)) ) / ((1 - xmax**(1-tau)) - (1 - xmin**(1-tau)))

        # ensure cdfs are oriented vertically 
        cdfX = np.reshape(cdfX, [np.size(cdfX), 1])
        cdfFit = np.reshape(cdfFit, [np.size(cdfFit), 1])
    
    # calculate Kolmogorov-Smirnov statistic for empirical data

    empiricalKS = np.max(np.abs(cdfX - cdfFit))

    # save result
    ks = {}
    ks['empiricalKS'] = empiricalKS
    ks['samples_KS'] = np.zeros((nSamples, 1))

    #### calculate the p-value
    
    successCounts = np.zeros(1, nSamples)

    # computation reduces to binomial process
    # carry out conditional on the likelihood of the number of succesess given the number of trials

    nSuccesses = 0
    thisLikelihood = 1
    binomialFlag = True 
    criticalThreshold = nSamples*pCrit

    iSample = 0
    sampleTau = np.zeros(nSamples, 1)

    while iSample <= nSamples and thisLikelihood > likelihood and binomialFlag:
        if len(diff) == 0: # integer data

            # generate sample dara using fit of the real data
            xSampleNo = cr.mymnrnd(nSupportEvents, pdfFit)
            xSample = rldecode(xSampleNo, np.arange(xmin, xmax+1))

            # computer sample pdf using unique bins
            pdfsample = np.histogram(xSample, bins = np.arange(xmin, xmax+1))[0] / nSupportEvents

            #fit sample to pre-determined support range (only compute if we'll use it)
            if pCrit == 1
                thisTau = cr.plmle(xSample, xmin=xmin, xmax=xmax)
                sampleTau[iSample] = thisTau
            
            # compute cdf of sample
            cdfSample = 1 - np.cumsum(pdfSample)

            # ensure cdfs are oriented vertically
            cdfSample = np.reshape(cdfSample, [np.size(cdfSample), 1])

        else: # continuous 
            print(" i'm pretty sure this should never happen and for the love of god i don't want to translate another \
                    matlab script so if you get this message sucks for you -- go back to the matlab and translate this section\
                    -- Sahara")

        # calculate the KS statistic for simulated data
        sampleKS = np.max(np.abs(cdfSample - cdfFit))

        # record sample ks
        ks['samples_ks'][iSample] = sampleKS

        # record a success if the empirical KS is bounded above by the sample KS
        if empiricalKS <= sampleKS:
            successCounts[iSample] = 1
            nSuccesses+=1
        
        # stop computation if number of successes reaches threshold
        if nSuccesses == criticalThreshold:
            binomialFlag = False 
        
        # update likelihood if critical threshold less than 1
        if pCrit not 1:
            thisLikelihood = 1 - scipy.stats.binom.cdf(criticalThreshold - nSuccesses -1, nSamples-iSample, pCrit)
        
        iSample+=1
    
    p = np.sum(successCounts) / nSamples
    
    if iSample == (nSamples+1):
        sigmaTau = np.std(sampleTau)
    else:
        sigmaTau = None

    return p, ks, sigmaTau