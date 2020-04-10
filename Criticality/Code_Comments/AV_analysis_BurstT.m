function Result =AV_analysis_BurstT(data, varargin)
% data is a matrix with each row as a neuron and each column is a
% time bin. Or data could be a vector (network activity)
%
% Result is a structure. Result.S is avalanche sizes, and Result.T
% is avalanche durations.
%
% The threshold for network activity is 25 percentile by default,
% but could be set up manually. When network activity is dominated
% by silent period, perc could be zero. Otherwise, could try from roughly 20% to 50%.
% Threshold method is based on Poil, Simon-Shlomo, et al. "Critical-state
% dynamics of avalanches and oscillations jointly emerge from balanced
% excitation/inhibition in neuronal networks." Journal of Neuroscience
% 32.29 (2012): 9817-9823.
%
% Copyright @ Zhengyu Ma 2017
iVarArg = 1; perc = 0.25;
while iVarArg <= length(varargin)
    argOkay = true;
    switch varargin{iVarArg},
        case 'perc',   perc = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(AV_analysis) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end

[n m]=size(data); %get neurons# and frame#
% data = sign(data); % choose to use this line when the inputs are not
% binary values.

% ############# get network activity #######################
if size(data,1) == 1
    network = data;
else
    network=nansum(data);  % define the network
end

if network(1) == 0
    network(1) = 1;
end

% ################ define threshold ##################
if perc > 0
    sortN = sort(network);
    threshold2 = sortN(round(m*perc));
else
    threshold2  = 0;
end

% ############
m = length(network);
zdata = network;
z2data=zdata;

zdata=sign(zdata-threshold2);
% if min(zdata) == 0
%     zdata(find(zdata>0.25))=1;
% else % in this situation, min(network) = threshold2 = 0
zdata(find(zdata <1)) = 0;
% end
z1data = zdata;
Z = find(zdata==0); % intervals
Z1 = find(zdata~=0); % avalanches
zdata(Z(diff(Z)==1))=[]; % use a single 0 to separate avalanches (a series of 1s)
z1data(Z1(diff(Z1)==1))=[]; % use 1 to separate intervals (some study focused on interval distributions)
z2data(Z(diff(Z)==1))=[]; % use a single 0 to separate network activities in each avalanche
J = find(zdata==0);
J1 = find(z1data == 1);
data2 = data; data2(:,Z(diff(Z)==1))=[];

%% #################### Find Spike and AV sizes ########################
% figure(chunk)
burst = [];
clear burstshape
for i = 1:(length(J)-1)
    % #### Use this line if AVsize = # of Spikes. ##########
    fired = sum(z2data(J(i)+1: J(i+1)-1))-threshold2*(J(i+1)-J(i)-2);
    % fired = sum(z2data(J(i)+1: J(i+1)-1));
    
    % #### Use this line if AVsize = # of participants (neurons or units)
    %  fired= sum(sign(nansum(data2(:,J(i)+1:J(i+1)-1),2)));
    
    burst = [burst ; fired]; % stack burst
    %  burstshape{i} = z2data(J(i)+1: J(i+1)-1); % get burst shape for  shape collapse
end

%% ######### Duration Distribution ########################################
T=diff(J)+1; % AVduration
T1 = diff(J1);
T = T(T>0); % Duration should be positive

%% #################### Get final result ######################################
Result.S = burst;
Result.T = T;
