% Code below creates population summary for Significnat Saccade Direction Cells
% Written by Seth Konig
% Code does the following
% 1) Summarizes MRLs (mean resultant vector length) for place and non-place cells
% 2) Summarize circular non-uniformity p-value (biased by fixation count and firing rate)
% 3) Copies relevant figures for place cells to summary directory

%Code rechecked by SDK on 1/5/2017

clar %clear,clc

%where to store spatial analysis figure copies
summary_directory = 'C:\Users\seth.koenig\Desktop\Significant Units\Saccade Direction\';
if ~isdir(summary_directory)
    mkdir(summary_directory)
end

task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size

bin_deg = 1; %number of degrees per bin
bin_deg2 = 45; %initial degree bins to determine time of peak direction modualtion
smval_deg = 18; %9 degrees std
smval = 30; %smoothing in time



%---Values All Fixations out2out and in2in---%
all_mrls = []; %all observed MRLs (mean resultant vector length) ignoring out2in and in2out
all_mrl_pctiles = []; %observed MRLs shuffled percentile ignoring out2in and in2out

%---Values All Fixations out2out only---%
all_mlrs_out = []; %observed MRLs for out2out fixations only
all_mrl_out_pctiles = []; %observed MRLs shuffled percentile ignoring out2in and in2out

%---Other values---%
spatialness = []; %1 for place cell, 0 for non place cell
all_unit_names = {};
all_monkeys = []; %1s and 2s for monkeys

prefered_firing_rate_curves = [];
anti_prefered_firing_rate_curves = [];
ratio_prefered = [];
all_windows  = [];

figure_dir = {};
all_dirs = []; %saccade directions for significant units
for monkey = 2:-1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir{1} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML[
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir{2} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for sess =1:length(session_data)
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
        if exist([data_dir task_file(1:8) '-Saccade_Direction_Analysis.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-Saccade_Direction_Analysis.mat'])
            load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'spatial_info')
            load([data_dir  task_file(1:8) '-Place_Cell_Analysis.mat']);
        else
            if num_units ~= 0
                error('Where is this file')
            end
            continue
        end
        
        num_units = size(unit_stats,2);
        for unit = 1:num_units
            if ~isnan(mrls.all_fixations(unit)) %if unit was processed
                
                %---Values All Fixations out2out and in2in---%
                all_mrls = [all_mrls mrls.all_fixations(unit)]; %all observed MRLs (mean resultant vector length) ignoring out2in and in2out
                all_mrl_pctiles = [all_mrl_pctiles mrls.all_fixations_shuffled_prctile(unit)];%observed MRLs shuffled percentile ignoring out2in and in2out
                
                all_dirs = [all_dirs saccade_direction{unit}];
                
                %---Values All Fixations out2out only---%
                all_mlrs_out = [all_mlrs_out mrls.out2out(unit)]; %observed MRLs for out2out fixations only
                all_mrl_out_pctiles = [all_mrl_out_pctiles mrls.out2out_shuffled_prctile(unit)]; %observed MRLs shuffled percentile ignoring out2in and in2out
                
                %---Other values---%
                all_unit_names = [all_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}];
                all_monkeys = [all_monkeys monkey];
                if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                    spatialness = [spatialness 1]; %place cell
                else
                    spatialness = [spatialness 0]; %non place cell
                end
                
                if  mrls.all_fixations_shuffled_prctile(unit) > 95
                    
                    %---Get Firing Rate Curves for Prefered vs AntiPreferredDirections---%
                    %unfortunately not enough was saved
                    degrees = [0:bin_deg:360]-180;
                    egrees = degrees(2:end);
                    degrees = degrees*pi/180;
                    % degrees = [degrees degrees(1)];
                    
                    fix_aligned = list_fixation_locked_firing{unit}; %fixation aligned firing rate
                    sac_dirs = saccade_direction{unit}; %saccade directions organized the same way as fix_algined
                    fr = nandens(fix_aligned,smval,'gauss',1000,'nanflt'); %firing rate curve aligned to fixations
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%---Method 2: Determine time window of interest based on Time of Peak Firing Rate Modulation---%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %prefered since we currently assume most units show a preference for
                    %certain stimulus features at a particular time relative to eye
                    %movements. Method simply finds the maximum firing rate deviation from
                    %the mean when we bin the firing rate curves into arbitrary directions.
                    %Note: the exact bins are not important for this method just to remove
                    %the anti-prefered direction since this seems to often cause inhbition.
                    degrees2 = [0:bin_deg2:360]-180;
                    curves = NaN(length(degrees2),twin1+twin2);
                    for bin = 2:length(degrees2)
                        these_directions = (sac_dirs < degrees2(bin) & sac_dirs >= degrees2(bin-1));
                        curves(bin,:) = nandens(fix_aligned(these_directions,:),30,'gauss',Fs);
                    end
                    mean_of_curves = nanmean(fr);
                    crude_direction_curves = curves; %save for later
                    curves = abs(curves-mean_of_curves); %negatvie and positive work
                    %don't include more than 100 ms before fixation and 200 ms after, have
                    %yet to see a good example with prominent tuning outside this anyway
                    time_of_max_modulation = find(curves(:,100:400) == max(max(curves(:,100:400))));
                    [~,time_of_max_modulation] = ind2sub(size(curves(:,100:400)),time_of_max_modulation); %got to convert back into time index
                    time_of_max_modulation = time_of_max_modulation+100; %since ingored the first 100 ms before
                    window = time_of_max_modulation-((window_width/2)-1):time_of_max_modulation+window_width/2; %take window around peak
                    
                    if (spatial_info.shuffled_rate_prctile(unit) > 95) || (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                        fix_aligned = fix_aligned(in_out{unit} == 2 | in_out{unit} == 4,:); %fixations in2in or out2out
                        sac_dirs = sac_dirs(in_out{unit} == 2 | in_out{unit} == 4); %directions for fixations in2in or out2out
                    end
                    
                    all_windows = [all_windows {window}];
                    
                    bin_deg = 12;
                    [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,fix_aligned,window,sac_dirs);
                    binned_firing_rate_curves{1,unit} = mean_binned_firing_rate; %binned firing rates
                    [prefered_dirs,anti_prefered_dirs,smoothed_direction_curve] = ...
                        select_prefred_indeces(binned_firing_rate_curves{1,unit},degrees,sac_dirs,2);
                    
                    antipref = nandens3(fix_aligned(anti_prefered_dirs,:),smval,1000);
                    pref = nandens3(fix_aligned(prefered_dirs,:),smval,1000);
                    a = antipref;
                    a(a < 0.1) = mean(a);
                    a(a == 0) = 1;
                    ratio_prefered = [ratio_prefered; pref./a];
                    
%                     if max(pref./a) > 10
%                        disp('now') 
%                     end
                    
                    if any(isnan(ratio_prefered(end,:)))
                        disp('now')
                    elseif sum(isinf(ratio_prefered(:))) > 0
                        disp('now')
                    end
                    
                    pref = pref-mean(pref);
                    pref = pref/max(abs(pref));
                    antipref = antipref-mean(antipref);
                    antipref = antipref/max(abs(antipref));
                    
                    prefered_firing_rate_curves = [prefered_firing_rate_curves; pref];
                    anti_prefered_firing_rate_curves = [anti_prefered_firing_rate_curves; antipref];
                end
                
                %%
            end
        end
    end
end
%%
clc
disp([num2str(nansum(spatialness)) ' place cells'])
disp([num2str(sum(all_mrl_pctiles > 95)) ' directionally modulated cells for "all fixations"'])
disp([num2str(sum(all_mrl_out_pctiles > 95)) ' directionally modulated cells for OUT2OUT fixations'])
disp([num2str(sum(all_mrl_out_pctiles > 95 & all_mrl_pctiles > 95)) ' directionally modulated cells for OUT2OUT fixations & "all fixations"'])
disp('--------------------------------------------------------------')
disp([num2str(sum(all_mrl_pctiles > 95 & spatialness == 1)) ' directionally modulated place cells for "all fixations"'])
disp([num2str(sum(all_mrl_out_pctiles > 95 & spatialness == 1)) ' directionally modulated place cells for OUT2OUT fixations'])
%%
%%---Copy Relevant Figures to Summary Directory---%
for unit = 1:length(all_unit_names)
    if all_mrl_pctiles(unit) > 95
        sub_dir1 = 'Saccade Direction\';
        name1 = [all_unit_names{unit} '_Saccade_Direction_Analysis.png'];
        if spatialness(unit) == 1 %place cell
            copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name1],...
                [summary_directory 'Place\' name1])
        elseif spatialness(unit) == 0 %non place cell
            copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name1],...
                [summary_directory 'Non Place\' name1])
        end
    end
end

%% Distribution of Saccade Directions
smval_deg =18; %9 degrees std

binned_dirs = zeros(1,360);
degrees = [-180:180];
for b = 1:length(degrees)-1
    binned_dirs(b) = sum(all_dirs < degrees(b+1) & all_dirs >= degrees(b));
end
degrees = degrees*pi/180;

bf = [binned_dirs(end-(3*smval_deg):end) binned_dirs binned_dirs(1:3*smval_deg)];%so don't get edege artifacts
binned_dirs = nandens(bf,smval_deg,'gauss',1); %smooth to display firing rate by saccade direction
binned_dirs = binned_dirs(3*smval_deg+2:end-3*smval_deg);%remove buffers
polarplot(degrees,[binned_dirs(end) binned_dirs])

%%
tm = -twin1:twin2-1;
figure
subplot(2,2,1)
hold on
plot(tm,nanmean(anti_prefered_firing_rate_curves),'r')
plot(tm,nanmean(prefered_firing_rate_curves),'g')
plot([tm(1) tm(end)],[0 0],'k')
plot([0 0],[-1 1],'k--')
hold off
xlabel('Time From Fixation Start (ms)')
ylabel('Normalized Firing Rate')
title('Population Average Mormalzed Firing Rates')
legend('Anti-prefered','Prefered')

subplot(2,2,2)
[~,mxi] = max(prefered_firing_rate_curves');
[~,si] = sort(mxi);
imagesc(tm,1:size(prefered_firing_rate_curves,1),prefered_firing_rate_curves(si,:))
xlabel('Time From Fixation Start (ms)')
ylabel('Neuron #')
colormap('jet')
title('Time Cell Plot')
%%
subplot(2,2,3)
plot(tm,nandens3(ratio_prefered,20,1))
xlabel('Time From Fixation Start (ms)')
ylabel('Ratio')
title('Average Ratio of Prefered/Anti-Prefered')
box off

%%
window_conc = zeros(1,twin1+twin2);
for w = 1:length(all_windows);
    window_conc(all_windows{w}) = window_conc(all_windows{w})+1;
end

subplot(2,2,4)
plot(tm,100*window_conc/length(all_windows))
xlabel('Time From Fixation Start (ms)')
ylabel('% of Units')
title('Distribution of 100 ms Window with Greatest Direction Modulation')
box off
