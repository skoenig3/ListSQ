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
figure_dir =  {};
for monk = 1:2
    monkey = monkeys{monk};
    
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
        
        for session = 1:length(session_data)
            [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
                get_task_data(session_data{session},task);
            if isempty(task_file)
                continue
            end
            
            
            load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials','cfg','hdr','data','fixationstats');
            
            %get unit data
            [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
                multiunit,unit_confidence,sorting_quality);
            
            if exist([data_dir task_file(1:8) '-Fixation_locked_Sequence_results.mat'],'file') %want to remove later
                load([data_dir task_file(1:8) '-Fixation_locked_Sequence_results.mat']);
                load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'spatial_info')
            else
                continue
            end
            
            disp(['Importing ' task_file(1:8)])
            
            for unit = 1:num_units
                if ~isempty(fixation_locked_firing{unit})
                    if fixation_info.rate_prctile(unit) > 95
                        col = 1;%significant modulation
                        sig_count(1) = sig_count(1)+1;
                        %                 elseif fixation_info2.rate_prctile(unit) > 95
                        %                     col = 1;%significant modulation
                        %                     sig_count(2) = sig_count(2)+1;
                    else
                        col = 2;%not saccade modulated
                        sig_count(3) = sig_count(3)+1;
                        continue
                    end
                    
                    if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                        spatialness = [spatialness 1]; %place cell
                    else
                        spatialness = [spatialness 0]; %non place cell
                    end
                    
                    
                    
                    fixation_firing = fixation_locked_firing{unit};
                    fixation_firing(nansum(fixation_firing,2) == 0,:) = [];%remove fixations without spikes
                    
                    %%
                    all = [fixation_locked_firing{unit}(item_nums{unit} == 1,:) ... 
                    fixation_locked_firing{unit}(item_nums{unit} == 2,:) ...
                    fixation_locked_firing{unit}(item_nums{unit} == 3,:) ...
                    fixation_locked_firing{unit}(item_nums{unit} == 4,:)];
                    
                    [firing_rate,~]= nandens(fixation_firing,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    %                 firing_rate = firing_rate-nanmean(firing_rate);
                    %                 fr = firing_rate;
                    %                 if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                    %                     firing_rate = firing_rate/max(firing_rate);
                    %                 else%normalize to min %could have some neurons that only show supression
                    %                     firing_rate = firing_rate/min(firing_rate);
                    %                 end
                    
                    
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