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

figure_dir = {};
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
                
                %---Values All Fixations out2out only---%
                all_mlrs_out = [all_mlrs_out mrls.out2out(unit)]; %observed MRLs for out2out fixations only
                all_mrl_out_pctiles = [all_mrl_out_pctiles mrls.out2out_shuffled_prctile(unit)]; %observed MRLs shuffled percentile ignoring out2in and in2out
                
                %---Other values---%
                all_unit_names = [all_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}];
                all_monkeys = [all_monkeys monkey];
                if strcmpi([task_file(1:8) '_' unit_stats{1,unit}],'PW141008_sig002c'); %was found to be unreliable spatial unit so don't use
                    spatialness = [spatialness NaN];
                elseif (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                    spatialness = [spatialness 1]; %place cell
                else
                    spatialness = [spatialness 0]; %non place cell
                end
            end
        end
    end
end
%%
clc
disp([num2str(nansum(spatialness)) ' place cells'])
disp([num2str(sum(all_mrl_pctiles > 95)) ' direcitonally modulated cells for "all fixations"'])
disp([num2str(sum(all_mrl_out_pctiles > 95)) ' direcitonally modulated cells for OUT2OUT fixations'])
disp([num2str(sum(all_mrl_out_pctiles > 95 & all_mrl_pctiles > 95)) ' direcitonally modulated cells for OUT2OUT fixations & "all fixations"'])
disp('--------------------------------------------------------------')
disp([num2str(sum(all_mrl_pctiles > 95 & spatialness == 1)) ' direcitonally modulated place cells for "all fixations"'])
disp([num2str(sum(all_mrl_out_pctiles > 95 & spatialness == 1)) ' direcitonally modulated place cells for OUT2OUT fixations'])
%%
%%---Copy Relevant Figures to Summary Directory---%
for unit = 1:length(all_unit_names)
    if all_mrl_pctiles(unit) > 95
        sub_dir1 = 'Saccade Direction\';
        name1 = [all_unit_names{unit} '.png'];
        if spatialness(unit) == 1 %place cell
            copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name1],...
                [summary_directory 'Place\' name1])
        elseif spatialness(unit) == 0 %non place cell
            copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name1],...
                [summary_directory 'Non Place\' name1])
        end
    end
end