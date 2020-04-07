function Result  = AV_analysis_ExponentErrorComments(burst, T, bm, tm, varargin)
% This function calculate exponents for PDF(AVsize), PDF(AVdura), and
% scaling relation. When flag == 1 (default),
iVarArg = 1;
flag = 1; % default value
while iVarArg <= length(varargin)
    argOkay = true;
    switch varargin{iVarArg},
        case 'flag',   flag = varargin{iVarArg+1}; iVarArg = iVarArg + 1; % switch flag
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(AV_analysis) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end


L = max(burst);
s = 1:round(L);     s = s';
[pdf junk] = histc(full(burst),s);

% calculate exponent and other parameters for AVsize.
[burstMax, burstMin, alpha] = EXCLUDE(burst(burst < max(burst)^0.8),'setmin',bm) % get burstMin and burstMax
[alpha ,xmin] = plmle(burst( burst > burstMin  & burst > burstMin),'limit',burstMin);
xmin = burstMin;
Result.burst = burst;
Result.alpha = alpha;
Result.xmin = xmin;
Result.xmax = burstMax;

if flag == 2
    % ##### calculate pvalue for null hypothesis test #####
    Result.P.burst  = pvaluenew(burst(burst <burstMax & burst > burstMin));
    
elseif flag == 3
    %  plot distribution together with fitted distribution
    % ############### Plot PDF ####################
    subplot(1,3,1)
    pdf = hist(burst,1:max(burst));
    loglog(pdf/sum(pdf),'o', 'markerfacecolor', [.8 .2 .2], 'color', [.8 .2 .2]);
    box off; hold on;
    
    % ############### Plot fitted PDF #################
    x = burstMin:burstMax;
    y = length(find(burst == xmin))/(xmin^-alpha)*x.^-alpha;
    y = y/sum(pdf);
    loglog(x,y, 'k');
    xlabel('AVsize'); ylabel('PDF(S)')
    title(['AVsize PDF, ', num2str(alpha)])
end

%% ######### Duration Distribution ########################################

r = 1:max(T);
[tdf junk] = histc(T, r);
% calculate exponent and other parameters for AVdura.
[tMax, tMin, beta] = EXCLUDE(T,'setmin',tm)
[beta tmin] = plmle(T(T < tMax & T > tMin),'limit',tMin);
%     [beta tmin] = plmle(T,'limit',tMin);
tmin = tMin;
tMax = max(T);
Result.T = T;
Result.beta = beta;
Result.tmin = tMin;
Result.tmax = tMax;

if flag == 2
    % ##### calculate pvalue for null hypothesis test #####
    Result.P.t  = pvaluenew(T(T < tMax & T > tMin));
elseif flag == 3
    %  plot distribution together with fitted distribution
    % ############### Plot PDF ####################
    subplot(1,3,2)
    tdf = hist(T,1:max(T));
    loglog(tdf/sum(tdf),'o', 'markerfacecolor', [.8 .2 .2], 'color', [.8 .2 .2]);
    box off; hold on;
    
    % ############### Plot fitted PDF #################
    x = tMin:tMax;
    y = length(find(T == (tmin+2)))/((tmin+2)^-beta)*x.^-beta; y = y/sum(tdf);y
    loglog(x,y, 'k');
        xlabel('AVduration'); ylabel('PDF(D)')
    title(['AVdura PDF, ', num2str(beta)])
end

%% ################## scaling relation #####################
TT = 1:max(T); Sm = []; clear Smshape
% ########### Calculate average size for each duration #########
for i=1:length(TT)
    Sm=[Sm mean(burst(find(T==TT(i))))];
end

% ################## get fit and pre exponent ##################
Loc=find(Sm>0);TT=TT(Loc);Sm=Sm(Loc);
fit_sigma = polyfit(log(TT(intersect(find(TT>tMin), find(TT<tMin+30)))), log(Sm(intersect(find(TT>tMin), find(TT<tMin+30)))), 1);
sigma = (beta - 1)/(alpha - 1);
Result.pre = sigma;
Result.fit = fit_sigma;
Result.df = abs(sigma - fit_sigma(1));
Result.TT = TT;
Result.Sm = Sm;

if flag == 3
    % ############# Plot scaling relation/fitted/predicted ###########
    subplot(1,3,3)
    loglog(TT, (TT).^sigma/(10^sigma)*s(10)); hold on; box off
    loglog(TT, (TT).^fit_sigma(1)/(10^fit_sigma(1))*s(10),'b')
    loglog(TT, Sm,'o', 'markerfacecolor', [.8 .2 .2], 'color', [.8 .2 .2]);
    legend('pre', 'fit'); legend boxoff;
    xlabel('Duration'); ylabel('<S>')
    title(['Difference = ', num2str(Result.df)])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.8]);
    
end
