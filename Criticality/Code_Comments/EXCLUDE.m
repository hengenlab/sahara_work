function [burstMax, burstMin, alpha] = EXCLUDE(burst, varargin)
% [burstMax, burstMin, alpha] = EXCLUDE(burst, varargin) is a function that
% determine both the lower and upper boundaries with a small KS.
% varargin return different inputs for setmin, num and flag.
%                                                                                                                                             
% burst could be avalanche sizes or durations, which will be applied to the
% significant test and compared with power law distribution. 
%
% Copyright @Zhengyu Ma 2016
iVarArg = 1;
KS = 1; setmin = 10;
num = 1;
while iVarArg <= length(varargin)
    argOkay = true;
    switch varargin{iVarArg},
        case 'setmin',      setmin = varargin{iVarArg+1}; iVarArg = iVarArg + 1
        case 'num',       num = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'flag',   flag = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(AV_analysis) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end

dKS = 1; xmin = 1;
while KS > min(num/sqrt(length(burst(burst>xmin))),0.1) & dKS > 0.0005
    setmin
    [alpha xmin] = tplfit(burst, setmin) % get exponent and lower boundary
    alpha = alpha(1);
    xmin = xmin(1);
    xmax = max(burst);
    N = length(burst);
    k = 0;
    z = burst;
    z = z(z>=xmin(1));            n    = length(z);
    cdf = cumsum(hist(z,xmin(1):xmax)./n);    
    s = unique(burst(xmin<=burst & burst<=xmax));
    smin = min(s);      %smax = max(s);  s = smin:smax;
    A = 1/ sum(s.^-alpha ); % constant factor for perfect power law
    fit = cumsum ( A*(xmin:xmax).^-alpha); % CDF for perfect power law
    KS_old  = KS; % assign previous KS as KS_old
    KS = max(abs(cdf - fit)); % calculate the new KS
    dKS = abs(KS_old - KS); % difference between the current KS and the KS in the previous step

    burst = burst(burst < max(burst)); % lower the upper boundary by one since here we used '<' but not '<='
    burstMax = max(burst); 
end
burstMin = xmin;

