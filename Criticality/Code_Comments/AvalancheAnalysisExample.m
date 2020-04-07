clear all; close all; clc
% Load example data. Getting from a deprived hemisphere during baseline
% slot. Have already binned to binary dataset.
load('Data.mat'); Data = NewData;


% ############### get AVsize and AVduration ################ 
%the default  threshold is 25 percentile. But we could change it by tuning
%'perc'. Make  sure 'perc' is in the range from 0 -- 1. When network
%activity is silent most time, should set 'perc' = 0, which means threshold
%is zero #########  
r = AV_analysis_BurstT(NewData,'perc',0.25);
x = r.S'; % x is AVsize
y = r.T; % y is AVdura


% ################## Avalanche analysis including AVsize, AVduration
% distribution and scaling relation. burstM and tM are used to set the
% limits for lower boundary for AVsize and AVduration distributions.
% Result1 only returns exponents and lower/upper bounds limits
% Result2 return pvalue for null hypothesis in addition Result2.P
% Result3 generate a figure including the distributions
burstM = 10; tM = 2;
Result1 = AV_analysis_ExponentErrorComments(x, y, burstM, tM);
Result2 = AV_analysis_ExponentErrorComments(x, y, burstM, tM, 'flag', 2);
display(['Pvalue for size distribution is : ', num2str(Result2.P.burst)])
display(['Pvalue for duration distribution is : ', num2str(Result2.P.t)])
display('For this example, pvalues for both size and duration distributions are larger than 0.05')
display('Null hypothesis could not be rejected. Dataset follows power law distribution')

Result3 = AV_analysis_ExponentErrorComments(x, y, burstM, tM, 'flag', 3);



