% Population Visual Response
% Code below creates population summary for Significnat Visually Reponsive Neurons
% Written by Seth Konig  written Seth Konig 9/1/16, updated 1/16/2017
% Code does the following
% 1) Summarizes visual responses to images for short and long windows
% 2) Determines if neurons may be sequentially organized in time
% 3) Determines whether place cells are also visually responsive
% 4) Determines if visually responsive neurons are also modulated by novel/repeat images
% 5) Tracks AP location, unit counts, and which monkey (not currently used)
% 6) Copies relevant figures to summary directory

%Code rechecked by SDK on 1/16/2017


clar %clear,clc

twin1 = 200;% how much time to take before an event
twin2 = 1000;%how much time to look at after reward and ITI start
twin3 = 500;%minimum fixation on cross hair
twin4 = 2200; %for long window on image on,should max at 2100 but cortex could have some lag
twin5 = 500;%aligned to appearance of crosshair
twin6 = 700;%maximum response duration
smval2 =60;

task = 'cvtnew';
min_blks = 1; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency

%---Misc. Parameters---%
monkey_all_unit_count = zeros(2,2);%row 1 place row 2 non-place, column by monkey
all_multi_unit_count = zeros(1,2); %all_multi_unit_countunits
all_unit_names = {}; %place cell unit names
all_monkeys = []; %1s and 2s
AP_location = []; %AP location of recorded place cell

%---Firing Rate Curves---%
cross_on_firing_rates = [];%appearance of cross
cross_fixation_firing_rates = [];%fixation on across
dot_on_firing_rates  = [];%dot appearance
all_dot_on_firing_rates = [];%not significant units
dot_color_change_firing_rates = [];
response_firing_rate = [];%when monkey releases bar
reward_firing_rate = [];%reward period
ITI_firing_rate = [];%firing rate aligned to ITI start

monkeys = {'Vivian','Tobii'};
figure_dir = {};
for monk =2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir{1} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir{2} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    end
    
    for sess = 1:length(session_data)
        %read in task file data
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        %load task file data
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials',...
            'stability_attribute','cfg','hdr','data');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        if num_units == 0
            continue
        end
        
        num_trials = length(cfg.trl); %number of trials
        valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
        if all(isnan(valid_trials(1,:)))
            continue
        end
        
        disp(task_file(1:8))
        if exist([data_dir task_file(1:8) '-cvtnew-time_locked_results.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-cvtnew-time_locked_results.mat']);
        else
            if num_units ~= 0
                error('Where is this file')
            end
            continue
        end
        
        
        for unit = 1:num_units
            if ~isempty(time_lock_firing{unit,1})
                
                %---Misc. Parameters---%
                monkey_all_unit_count(1,monk) = monkey_all_unit_count(1,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                
                
                %---for crosshair on---%
                if (epoch_data.rate_prctile(unit,1) > 95) && (epoch_data.temporalstability_prctile(unit,1) > 95)
                    firing_rate = time_lock_firing{unit,1};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate(:,1:twin1));
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    cross_on_firing_rates = [cross_on_firing_rates; firing_rate];
                end
                
                %---for fixation on crosshair---%
                if (epoch_data.rate_prctile(unit,2) > 95) && (epoch_data.temporalstability_prctile(unit,2) > 95)
                    firing_rate = time_lock_firing{unit,2};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate(:,1:twin1));
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    cross_fixation_firing_rates = [cross_fixation_firing_rates; firing_rate];
                end
                
                %---for dot on---%
                firing_rate = time_lock_firing{unit,3};
                num_not_nans = sum(~isnan(firing_rate));
                indeces = find(num_not_nans > 33);
                indeces(indeces > twin4) = [];
                firing_rate = firing_rate(:,indeces);
                
                [firing_rate,~]= nandens3(firing_rate,smval2,Fs);%going to smooth slightl lesss than before
                firing_rate = firing_rate-nanmean(firing_rate(:,1:twin1));
                if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                    firing_rate = firing_rate/max(abs(firing_rate));
                else%normalize to min %could have some neurons that only show supression
                    firing_rate = firing_rate/min(firing_rate);
                end
                if length(firing_rate) <twin4
                    firing_rate = [firing_rate NaN(1,twin4-length(firing_rate))];
                end
                
                if (epoch_data.rate_prctile(unit,3) > 95) && (epoch_data.temporalstability_prctile(unit,3) > 95)
                    dot_on_firing_rates = [dot_on_firing_rates; firing_rate];
                else
                    all_dot_on_firing_rates =[all_dot_on_firing_rates; firing_rate];
                end
                
                %---for dot color change---%
                if (epoch_data.rate_prctile(unit,4) > 95) && (epoch_data.temporalstability_prctile(unit,4) > 95)
                    firing_rate = time_lock_firing{unit,4};
                    num_not_nans = sum(~isnan(firing_rate));
                    indeces = find(num_not_nans > 33);
                    indeces(indeces > twin1+twin6) = [];
                    firing_rate = firing_rate(:,indeces);
                    [firing_rate,~]= nandens3(firing_rate,smval,Fs);%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate(:,1:twin1));
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    
                    if length(firing_rate) <twin1+twin6
                        firing_rate = [firing_rate NaN(1,twin1+twin6-length(firing_rate))];
                    end
                    
                    dot_color_change_firing_rates = [dot_color_change_firing_rates; firing_rate];
                end
                
                %---for monkey response---%
                if (epoch_data.rate_prctile(unit,5) > 95) && (epoch_data.temporalstability_prctile(unit,5) > 95)
                    firing_rate = time_lock_firing{unit,5};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate); %whole period since not sure what to look for yet
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    response_firing_rate = [response_firing_rate; firing_rate];
                end
                
                %---for reward period---%
                if (epoch_data.rate_prctile(unit,6) > 95) && (epoch_data.temporalstability_prctile(unit,6) > 95)
                    firing_rate = time_lock_firing{unit,6};
                    [firing_rate,~]= nandens(firing_rate,smval2,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate(:,1:twin1)); %whole period since not sure what to look for yet
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    reward_firing_rate = [reward_firing_rate; firing_rate];
                end
                
                %---for ITI period---%
                if (epoch_data.rate_prctile(unit,7) > 95) && (epoch_data.temporalstability_prctile(unit,7) > 95)
                    firing_rate = time_lock_firing{unit,7};
                    [firing_rate,~]= nandens(firing_rate,smval2,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate(:,1:twin1)); %whole period since not sure what to look for yet
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    ITI_firing_rate = [ITI_firing_rate; firing_rate];
                end
            end
        end
    end
end
%% Population Plots

%---Cross on---%
t13 = -twin1:twin3-1;
figure
subplot(1,2,1)
plot(t13,nanmean(cross_on_firing_rates))
xlabel('Time From Cross On (ms)')
ylabel('Normalized Firing Rate')
title('Population Average')
xlim([-twin1 twin4-twin1])
box off

[m,i] = max(cross_on_firing_rates,[],2);
[mm,ii] = sort(i);
subplot(1,2,2)
imagesc([-twin1:twin3-1],[1:size(cross_on_firing_rates,1)],cross_on_firing_rates(ii,:))
colormap('jet')
vals = cross_on_firing_rates(:,1:twin1);
caxis([-std(vals(:)) 1]);
title('Time Cell Plot')
ylabel('Cell #')
xlabel('Time From Cross on (ms)')
subtitle('Cross on Aligned Activity')

%---Fixation on Cross---%
figure
subplot(1,2,1)
plot(t13,nanmean(cross_fixation_firing_rates))
xlabel('Time From Fixation (ms)')
ylabel('Normalized Firing Rate')
title('Population Average')
xlim([-twin1 twin3])
box off

[m,i] = max(cross_fixation_firing_rates,[],2);
[mm,ii] = sort(i);
subplot(1,2,2)
imagesc([-twin1:twin4-1],[1:size(cross_fixation_firing_rates,1)],cross_fixation_firing_rates(ii,:))
colormap('jet')
vals = cross_fixation_firing_rates(:,1:twin1);
caxis([-std(vals(:)) 1]);
title('Time Cell Plot')
ylabel('Cell #')
xlabel('Time From Fixation (ms)')
subtitle('Fixation on Cross Algined Activity')
xlim([-twin1 twin4-twin1])

%---Dot on---%
t14 = -twin1:twin4-twin1-1;

num_not_nans = sum(~isnan(dot_on_firing_rates));
indeces = find(num_not_nans > 10);

figure
subplot(2,2,1)
plot(indeces-twin1,nanmean(dot_on_firing_rates(:,indeces)))
xlabel('Time From Dot On (ms)')
ylabel('Normalized Firing Rate')
title('Significant Population Average')
xlim([-twin1 twin4-twin1])
box off

[m,i] = max(dot_on_firing_rates,[],2);
[mm,ii] = sort(i);
subplot(2,2,2)
imagesc([-twin1:twin4-twin1-1],[1:size(dot_on_firing_rates,1)],dot_on_firing_rates(ii,:))
colormap('jet')
vals = dot_on_firing_rates(:,1:twin1);
caxis([-std(vals(:)) 1]);
title('Significant Time Cell Plot')
ylabel('Cell #')
xlabel('Significant Time From Dot on (ms)')
subtitle('Dot on Algined Activity')


subplot(2,2,3)
plot(indeces-twin1,nanmean(all_dot_on_firing_rates(:,indeces)))
xlabel('Time From Dot On (ms)')
ylabel('Normalized Firing Rate')
title('Non-Significant Population Average')
xlim([-twin1 twin4-twin1])
box off

[m,i] = max(all_dot_on_firing_rates,[],2);
[mm,ii] = sort(i);
subplot(2,2,4)
imagesc([-twin1:twin4-twin1-1],[1:size(all_dot_on_firing_rates,1)],all_dot_on_firing_rates(ii,:))
colormap('jet')
vals = all_dot_on_firing_rates(:,1:twin1);
caxis([-std(vals(:)) 1]);
title('Non-significant Time Cell Plot')
ylabel('Cell #')
xlabel('Significant Time From Dot on (ms)')
subtitle('Dot on Algined Activity')

%%

%---Dot Color Change---%
t16 = -twin1:twin6-1;
figure
subplot(1,2,1)
plot(t16,nanmean(dot_color_change_firing_rates))
xlabel('Time From Fixation (ms)')
ylabel('Normalized Firing Rate')
title('Population Average')
xlim([-twin1 twin3])
box off

[m,i] = max(dot_color_change_firing_rates,[],2);
[mm,ii] = sort(i);
subplot(1,2,2)
imagesc([-twin1:twin6-1],[1:size(dot_color_change_firing_rates,1)],dot_color_change_firing_rates(ii,:))
colormap('jet')
vals = dot_color_change_firing_rates(:,1:twin1);
caxis([-std(vals(:)) 1]);
title('Time Cell Plot')
ylabel('Cell #')
xlabel('Time From Fixation (ms)')
subtitle('Fixation on Cross Algined Activity')
%%

%---Bar Response---%
t51 = -twin5:twin1-1;
figure
subplot(1,2,1)
plot(t51,nanmean(response_firing_rate))
xlabel('Time From Bar Release (ms)')
ylabel('Normalized Firing Rate')
title('Population Average')
xlim([-twin5 twin1])
box off

[m,i] = max(response_firing_rate,[],2);
[mm,ii] = sort(i);
subplot(1,2,2)
imagesc([-twin5:twin1-1],[1:size(response_firing_rate,1)],response_firing_rate(ii,:))
colormap('jet')
vals = response_firing_rate(:,1:twin1);
caxis([-std(vals(:)) 1]);
title('Time Cell Plot')
ylabel('Cell #')
xlabel('Time From Bar Release (ms)')
subtitle('Bar Release Algined Activity')
%%

%---Reward---%
t12 = -twin1:twin2-1;
figure
subplot(1,2,1)
plot(t12,nanmean(reward_firing_rate))
xlabel('Time From Reward Start (ms)')
ylabel('Normalized Firing Rate')
title('Population Average')
xlim([-twin1 twin1+twin2])
box off

[m,i] = max(reward_firing_rate,[],2);
[mm,ii] = sort(i);
subplot(1,2,2)
imagesc([-twin1:twin1+twin2-1],[1:size(reward_firing_rate,1)],reward_firing_rate(ii,:))
colormap('jet')
vals = reward_firing_rate(:,1:twin1);
caxis([-std(vals(:)) 1]);
title('Time Cell Plot')
ylabel('Cell #')
xlabel('Time From Reward Start (ms)')
subtitle('Reward Algined Activity')
%%
%---ITI---%
t12 = -twin1:twin2-1;
figure
subplot(1,2,1)
plot(t12,nanmean(ITI_firing_rate))
xlabel('Time From ITI Start (ms)')
ylabel('Normalized Firing Rate')
title('Population Average')
xlim([-twin1 twin2])
box off

[m,i] = max(ITI_firing_rate,[],2);
[mm,ii] = sort(i);
subplot(1,2,2)
imagesc([-twin1:twin2-1],[1:size(ITI_firing_rate,1)],ITI_firing_rate(ii,:))
colormap('jet')
vals = ITI_firing_rate(:,1:twin1);
caxis([-std(vals(:)) 1]);
title('Time Cell Plot')
ylabel('Cell #')
xlabel('Time From ITI Start (ms)')
subtitle('ITI Algined Activity')