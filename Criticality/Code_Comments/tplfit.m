function [alpha xmin ks Loglike] = tplfit(burst,limit, x_start, vargarin)

% function fit distribution of burst to power law, and return the
% parameters which give the best power fit (with minimum KS compared with
% perfect power law
%
% alpha: exponent for power law fitting
% xmin: truncated point
% ks: value for Kolmogorov–Smirnov test
%
% parameters:
% burst: could be either avalanche sizes or durations (a vector)
% limit set up the range from which we could get the xmin
% x_start fix the xmin as x_start

KS = []; alpha = [];  Loglike = [];

switch nargin
    case 2
        range = 1:limit;
    case 3
        range = x_start;
end


for x0 = range
    
    % [alpha xmin]=plfit(burst);
    X = burst;
    xmax = max(burst);
    n = length(burst(x0<=burst & burst<=xmax)) ; % number of avalanches
    s = unique(X(x0<=X & X<=xmax)); 
    smin = min(s); % minimum avalanche
    smax = max(s); % maximum avalanche
    
    LL = @(x) x*sum( log( burst(x0<=burst & burst<=xmax) ) ) - n*log( 1/sum(s.^-x ) ) ; % Maximum likelihood
    [a,fval] = fminsearch(LL , 2.3 ); % start search from 2.3    
    Loglike = [Loglike , -fval]; % value of the minimum likelihood
    
    
    % X = burst;
    %
    % for i = 1:1000
    %
    %     X = Shuffle(X);
    %     I = randi(length(X),[1 length(X)]);
    %     AV = X(I);
    %     [sth xmin]=plfit(AV);
    % %     xmax = max(AV);
    %     n = length(AV(xmin<=AV & AV<=xmax)) ;
    %     s = unique(AV(xmin<=AV & AV<=xmax));
    %     smin = min(s);      smax = max(s);
    %     LL = @(x) x*sum( log( AV(xmin<=AV & AV<=xmax) ) ) - n*log( 1/sum(s.^-x ) ) ;
    % %     options = optimset('PlotFcns',@optimplotfval);
    %     [a2,fval] = fminsearch(LL , 2.8);
    %     alpha = [alpha , a2];
    %     clear a2 I
    % end
    %
    % SD = std(alpha);   
    %% ########################## Hypothesis Testing ###############################
    % alpha = a;
    % [sth xmin] = plfit(burst);
    xmax = max(burst);    
    N = length(burst);
     k = 0;    
    z = burst;
    z = z(z>=x0);            n    = length(z);
    cdf = cumsum(hist(z,x0:xmax)./n); % cdf
    
    s = unique(X(x0<=X & X<=xmax));
     smin = min(s);      smax = max(s);
    A = 1/ sum(s.^-a ); % constant for perfect power law
    alpha = [alpha , a];
    fit{x0} = cumsum (A*(x0:xmax).^-a ); % perfect power law PDF
    KS = [KS , max(abs(cdf - fit{x0})) ] ; % calculate Kolmogorov–Smirnov

    
end

switch nargin
    case 2
        xmin = find(KS==min(KS)) ;  alpha = alpha( xmin );  Loglike = Loglike(xmin);
    case 3
        xmin = x_start;
end


ks = min(KS);

 





