% Extract data to KH__CirCadianSignbhvControl and KH__CirCadianSignbhvDeprived data,
% for control hemisphere, data from all chunks except chunk 6 are applied.
% for deprived hemisphere, only chunk 1-5, although actually both
% hemisphere share the same behavior states
% 4 behavior, 4 hour windows (3 windows in each chunk)
close all; clear all; clc
addpath('F:\Research\Project with Keith\Codes\anova_rm')
cd('F:\Research\Project with Keith\Codes')
set(0, 'defaultaxesfontsize', 15)
set(0, 'defaulttextfontsize', 12)
set(0, 'defaultaxeslinewidth', 2)
Sub = [42 45 47	49 50 63 36 40 67  73 75 ]; % for the first six and last 3 animals, we have both the contiguous sorted spikes and the chunk data,
%for KH36 and KH40, we have only contiguous sorted data.
% ############# basic information #################
% Chk shows how many half days for each subject, colorbhv defines what color
% each behavior uses. bhvr include 4 behavior states.
Chk = [18 18 18 17 18 21 18 18 21  21 21];
colorbhv = [0.9058    0.2785    0.9706; ...
    0.1270    0.5469    0.9572; ...
    0               0               0            ; ...
    0.9134    0.9575    0.4854; ...
    0.7524    0.9649    0.8003];
bhvr = {'REM', 'NREM', '' ,'Active', 'Quiet'} % four behaviors, the middle (3rd) one is empty
days = 1:9; % define which days are included in the analysis
day_length = length(days);
bhv = 4
for sub_loc =[1:6 8]%7:numel(Sub) % choose animals we want to include
    clear Result % clear Result to avoid influence from last subject
    subject = Sub(sub_loc);
    behavior = load(['F:\Research\Project with Keith\Data\MLS\KH',num2str(subject),'_STATETIMES_FINAL']);
    % ####### first column of behavior is behavior labels, second column is
    % the  begining of the behavior #########
    behavior = behavior.statetimes; 
    % ##### calculate duration for each behavior slot. the 3rd column is
    % the end of each behavior, and the 4th are the durations #######
    behavior(:,3) =[behavior(2:end,2);734000]; % end of each behavior
    behavior(:,4) = behavior(:,3)-behavior(:,2);    
    behavior = behavior(2:end-1,:); % eliminate the first and last behavior slots. 
    i = 0
    clear perc12
    for time = 1:4*3600:9*24*3600
        i = i+1;
        % ############# find all locations within 4 hour window ###########        
        loc =  intersect(find(behavior(:,2) > time), find(behavior(:,3)< time+4*3600));
        behavior_window = behavior(loc,:);
        if bhv == 4
            % #### Active sometimes include 3 as well. #####
            % #### 12 in loc12(perc12) represent single behavior.
            % #### 45 in loc45 represent all behaviors
            loc12 = [find(behavior_window(:,1) == bhv);find(behavior_window(:,1) == 3)]; 
        else
            loc12 = [find(behavior_window(:,1) == bhv)];
        end
        loc45 = [find(behavior_window(:,1) == 4);find(behavior_window(:,1) == 3);...
            find(behavior_window(:,1) == 5); find(behavior_window(:,1) == 2); find(behavior_window(:,1) == 1)];
        if isempty(behavior_window) == 0
            perc12(i) =nansum(behavior_window(loc12,4))/nansum(behavior_window(loc45,4));
            %perc12(i) = nansum(behavior_window(loc12,4))/(behavior_window(end,3)-behavior_window(1,2));%/(nansum(behavior_window(loc12,4))+nansum(behavior_window(loc45,4)));
        else
            perc12(i) = nan
        end
    end
    perc = reshape(perc12(1:54), 6, 9)*100; % change ratio to percentage
    Perc_Bhv_All{bhv}{sub_loc} = perc;
    perc(perc == 0) = nan;
    perc(4:6, 3) = nan; perc_sub(:,sub_loc) = nanmean(perc(:,days),2); %first 3 days or all days
    %     perc_sub(:,sub_loc) = nanmean(perc(:,[1:2 4:8]),2); % all days
    
    % ########## plot percentage bars for each individual ########
    subplot(2,4,sub_loc-heaviside(sub_loc-7))
    c = bar(2:4:12, perc_sub(1:3,sub_loc)); hold on % day (1st 12 hours)
    c.FaceColor = colorbhv(bhv,:);
    c = bar(14:4:24, perc_sub(4:6,sub_loc)); hold on % night (2nd 12 hours)
    c.FaceColor = [1 1 1];
    c.EdgeColor = colorbhv(bhv,:);
    errorbar(2:4:24, perc_sub(:, sub_loc), nanstd(perc')/sqrt(day_length),'k.'); % sqrt(days we used)
    axis([0 25 0 ceil(max(perc_sub(:))/10)*10+5]); box off
    title([bhvr{bhv}, ': Animal KH ', num2str(Sub(sub_loc))]); xlabel('Time (hour)'); ylabel('Asleep percentage')
    set(gca,'xtick',[2:4:24])
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.75, .6]*1.5);

%% ########## plot percentage bars for all animal ########
figure()
c = bar(2:4:12, mean(perc_sub(1:3,[1:6 8]),2)); hold on
c.FaceColor = colorbhv(bhv,:);
c = bar(14:4:24, mean(perc_sub(4:6,[1:6 8]),2)); hold on
c.FaceColor = [1 1 1];
c.EdgeColor = colorbhv(bhv,:);
errorbar(2:4:24, mean(perc_sub(:,[1:6 8]),2), std((perc_sub(:,[1:6 8]))')/sqrt(7),'k.'); % sqrt(# of subjects)
box off
xlim([0 25])
title([bhvr{bhv}, ': average percentage']); xlabel('Time (hour)'); ylabel('Asleep percentage')
set(gca,'xtick',[2:4:24])

% ############### add Anova test ######################
[p,~,stats]=anova2((perc_sub(:,[8 1:6]))') % p(1) is windows, p(2) is subjects
k = multcompare(stats);close
k2  = anova_2((perc_sub(:,[8 1:6]))')
loc_sig = find(k2(:,3)<0.05);
for loc_sig_num = 1:length(loc_sig);    groups{loc_sig_num} = [k2(loc_sig(loc_sig_num),[1 2])*4-2]; end
% close
if isempty(loc_sig_num) == 0
    H=sigstar(groups,k2(loc_sig,end));
end
box off
% saveas(gcf,[bhvr{bhv},' perc_alldays.pdf'])
%%  ############################## Save data in an excel table  #############################
Mean_ind = table(perc_sub(:,8),perc_sub(:,1),perc_sub(:,2),perc_sub(:,3),perc_sub(:,4),perc_sub(:,5),perc_sub(:,6),...
    'VariableNames',{'KH40','KH42', 'KH45', 'KH47', 'KH49', 'KH50', 'KH63'},'RowNames',{'h1to4','h5to8','h9to12','h13to16','h17to20','h21to24'})
Mean_Sem = table(mean(perc_sub(:,[1:6 8]),2),(std((perc_sub(:,[1:6 8]))'))'/sqrt(7),...
    'VariableNames',{'mean','sem'},'RowNames',{'h1to4','h5to8','h9to12','h13to16','h17to20','h21to24'})
Anova = table(k2(:,1),k2(:,2),k2(:,3),'VariableNames',{'Period1','Period2','Pvalue'})
Days ='First3days ' %'AllDays '% 
writetable(Mean_ind,'Sleep_percentage_Mean_individual.xls','Sheet',[Days, bhvr{bhv}],'WriteRowNames',true);
writetable(Mean_Sem,'Sleep_percentage_Mean_Sem.xls','Sheet',[Days,bhvr{bhv}],'WriteRowNames',true);
writetable(Anova,'Sleep_percentage_Anova.xls','Sheet',[Days ,bhvr{bhv}],'WriteRowNames',true);
