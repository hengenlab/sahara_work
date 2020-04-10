function [P_value ks]= pvaluenew(burst)
% Null hypothesis test. When P_value is very large, we could not reject
% the Null hypothesis that the distribution of burst follows power law.
% Usually, we use 0.05 as criteria.

[alpha xmin] = tplfit(burst,40);%max(5,min(burst)))
xmax = max(burst);
N = length(burst);
k = 0;
z = burst;
z = z(z>=xmin);            n    = length(z);
cdf = cumsum(hist(z,xmin:xmax)./n);

s = unique(burst(xmin<=burst & burst<=xmax));
smin = min(s);      %smax = max(s);  s = smin:smax;
A = 1/ sum(s.^-alpha );
fit = cumsum ( A*(xmin:xmax).^-alpha );
KS = max(abs(cdf - fit))
ks = [];        j = 1;
h = figure; semilogx(xmin:xmax,fit,'*',xmin:xmax,cdf,'+r'); hold on; xlabel('Size','Fontsize',16);   ylabel('Prob(size < S)','Fontsize',16);
title('CDF,field view # ','Fontsize',17);
A = lognfit(z); mu = A(1); sig = A(2);
% plot( ( erfc( (log(xmin)-mu)/sig/sqrt(2) ) - erfc( (log(s+1)-mu)/sig/sqrt(2) ) ) ...
%     ...
%                                                  / ( erfc( (log(xmin)-mu)/sig/sqrt(2) ) - erfc( (log(xmax)-mu)/sig/sqrt(2) ) ),'g' );
%
% plot( 1 - erfc( (log(s)-mu)/sqrt(2)/sig) / erfc( (log(xmin)-mu)/sqrt(2)/sig) ,'m');

legend('Power law CDF','Experimental CDF','Location','Southeast'); legend boxoff;
% text(18, 0.6, ['\alpha',' = ',num2str(alpha) , '\pm', num2str(SD)] ,'Fontsize',17 );
% pubgraph(h,14,2,'w')

cdfplot(z); grid off

Niter = 1000;

while j< Niter
    j
    %% ########################## Inverse Method  #############################################
    % % xmax = max(burst); N = length(burst(burst>xmin));
    % xmin=1;
    % generate a data set following power law with the same exponent.
    N = 20*length(burst(burst>=xmin)); % generate N data points
    syn_data = floor( (xmin-1/2)*(1-rand(1,N)).^(1/(1-alpha)) + 1/2 ) ;
    syn_data = floor(heaviside(xmax-syn_data)) .* syn_data;
    syn_data(syn_data==0)=[];
    syn_data = syn_data(1:length(burst(burst>=xmin))); % choose the same size as burst
    %% ########################## Accept-Reject Method  #############################################
    % N = 2*length(burst(burst>=xmin));
    % beta = alpha - 1;
    % umax = xmin^-beta;      u = umax * rand(1,N);       r = floor(u.^-(beta^-1) );
    % syn_data = r .* floor(heaviside( (P(r,xmin,alpha) .* Q(xmin,xmin,beta))./(P(xmin,xmin,alpha) .* Q(r,xmin,beta)   ) - rand(1,N) ));
    % syn_data = floor(heaviside(xmax-syn_data)) .* syn_data;
    % % syn_data(syn_data==0)=[];
    % Ind = find(syn_data > 0);
    % syn_data = syn_data(Ind(1:length(burst(burst>=xmin)))) ;
    
    
    %% ############################################################################################
    X = syn_data;
    n = length(syn_data(xmin<=syn_data & syn_data<=xmax)) ; 
    a = tplfit(X,xmin); % calculate exponent for surrogated data    
    if abs(a-alpha)<=0.1 & a >1.0   % Make sure we use the right surrogated data
        X = X(X>=xmin);            n    = length(X);
        cdf = cumsum(hist(X,xmin:xmax)./n);
        s = unique(X(xmin<=X & X<=xmax));
        smin = min(s);      smax = max(s);
        A = 1/ sum(s.^-alpha );
        fit = cumsum ( A*(xmin:xmax).^-alpha );
        %  A = 1/ sum(s.^-a );
        %  fit = cumsum ( A*(xmin:xmax).^-a );
        %  ks = max(abs(cdf - fit));
        
        % ######## ks is for surrogated and perfect power law with -alpha
        % as exponent #######################
        ks = [ks , max(abs(cdf - fit)) ];
        j = j + 1
    end
end
ks
P_value = sum(sign( ks(ks>=KS) ) )/Niter;

display(['P_value = ',num2str(P_value),])
display(['KS = ',num2str(KS),])

AA=[P_value,KS];
end

