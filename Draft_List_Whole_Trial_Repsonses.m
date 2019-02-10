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
%where to store spatial analysis figure copies
summary_directory = 'C:\Users\seth.koenig\Desktop\Significant Units\Visual Response\';
if ~isdir(summary_directory)
    mkdir(summary_directory)
end

task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size

spatialness = [];
all_rates = [];
count = 0;
monkeys = {'Vivian','Tobii'};
figure_dir = {};
for monk =1:2
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
        if exist([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat']) %visual response analysis
            load([data_dir task_file(1:8) '-ListSQ-Visual_Response_Memory_results']) %memory visual response analysis
            load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'spatial_info') %spatial analysis
        else
            if num_units ~= 0
                error('Where is this file')
            end
            continue
        end
        
        
        for unit = 1:num_units
            if ~isempty(time_lock_firing{unit,1})
                
                
                
                %---for fixation on cross hair---%
                if ((epoch_data.rate_prctile(unit,1) > 95) && (epoch_data.temporalstability_prctile(unit,1) > 95)) ...
                        || ((epoch_data.rate_prctile(unit,2) > 95) && (epoch_data.temporalstability_prctile(unit,2) > 95)) ...
                        || ((epoch_data.rate_prctile(unit,3) > 95) && (epoch_data.temporalstability_prctile(unit,3) > 95)) ...
                        || ((epoch_data.rate_prctile(unit,4) > 95) && (epoch_data.temporalstability_prctile(unit,4) > 95))...
                        || ((epoch_data.rate_prctile(unit,5) > 95) && (epoch_data.temporalstability_prctile(unit,5) > 95))...
                        
                    %---Unit Test Significance and Firing Rate Curves---%
                    if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                        spatialness = [spatialness 1]; %place cell
                    else
                        spatialness = [spatialness 0]; %non place cell
                    end
                    
                    firing_rate = time_lock_firing{unit,2}(:,1:twin1+twin3); %for now only look at first 500 ms after fixation
                    
                    
                    %---for long image on 5 second period---%
                    firing_rate3 = time_lock_firing{unit,5};
                    if size(firing_rate3,2) > 5200;
                        firing_rate3 = firing_rate3(:,1:5200);
                    end
                    
                    firing_rate4 = time_lock_firing{unit,4};
                    
                    try
                    combined_rates = [firing_rate firing_rate3 firing_rate4];
                    [combined_rates,~]= nandens(combined_rates,100,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    
                    combined_rates = combined_rates-nanmean(combined_rates);
                    if max(combined_rates) > abs(min(combined_rates)) %normalize to max
                        combined_rates = combined_rates/max(abs(combined_rates));
                    else%normalize to min %could have some neurons that only show supression
                        combined_rates = combined_rates/min(combined_rates);
                    end
                    all_rates = [all_rates; combined_rates];
                    count = count + 1
                    catch
                        disp('Skipped something')
                        continue
                    end
                end
            end
        end
    end
end
%%
vals = all_rates(:);
vals(vals > 0) = [];

[mx,mxi] = max(all_rates');
[~,smx] = sort(mxi);

figure
imagesc(all_rates(smx,:))
xlabel('Trial Time (ms)')
ylabel('Neuron #')
caxis([-std(vals) 1])

figure
plot(mean(all_rates))
%%
[U,S,V] = pca(all_rates);
T = kmeans(U,8);
hold all
for i = 1:max(T)
    if sum(T == i) > 1
        plot(mean(all_rates(T == i,:)))
    else
        plot(all_rates(T == i,:));
    end
end
hold off
title('Clustered Curves')

subtitle('Population Averages')