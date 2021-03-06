% Code below creates population summary for Significnat Place Cells
% Written by Seth Konig
% Code does the following
% 1) Summarizes spatial correlations for place cell and non-place cells
% 2) Calculates fixation aligned firing rate curves for place cells
% 3) Calculates eye coverage and place field coverage for place cells
% 4) Tracks AP location, unit counts, and which monkey (not currently used)
% 5) Contextual differences between list and sequence task
% 6) Copies relevant figures for place cells to summary directory

%Code rechecked by SDK on 1/5/2017
%Small bug fixed on 2/6/18 on designating inside and outside sequence
%trials SDK

clar %clear,clc

%where to store spatial analysis figure copies
summary_directory = 'C:\Users\seth.koenig\Desktop\Significant Units\Spatial Analysis\';
if ~isdir(summary_directory)
    mkdir(summary_directory)
end

task = 'ListSQ';
min_blks = 2; %only anaedlyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size

%---Misc. Parameters (Mostly Place Cells)---%
monkey_all_unit_count = zeros(2,2);%row 1 place row 2 non-place, column by monkey
all_multi_unit_count = zeros(1,2); %all_multi_unit_countunits
all_place_cell_unit_names = {}; %place cell unit names
all_place_cell_monkeys = []; %1s and 2s
place_cell_AP_location = []; %AP location of recorded place cell
place_cell_subregion = []; %e.g. CA3 or CA1
non_place_cell_subregion = [];

n_in = []; %number of fixations ?->in
n_out2in = [];%number of fixaitons out->in

%---Spatial Correlation Values for all cells---%
all_place_cell_spatial_corrs = []; %place cell spatial correlations
all_place_cell_spatial_skagg_percents= []; %place cell skagg information score percentiles
all_non_place_cell_spatial_corrs = [];%non-place cell spatial correlations
all_non_place_cell_skagg_percents = [];%non-place cell skagg information score percentiles
non_place_skagg_names = [];

%---place field properties---%
place_field_area = []; %area of place field
place_coverage = {}; %place field locations
coverage = {}; %eye data coverage for place cells

%---Place Cell Firing Rate Curves for List---%
all_in_rates = [];%fixations out-> in
all_in_minus_out_rates = [];%out-> in minus out-> out
all_out_rates = []; %fixaiton out-> out normalized by max of out->out
all_fixation_rates = [];%all fixaions normalized by max of all fixations
all_list_peak_times = [];%time of peak firing rate of place cells for out-> in fixations
sig_p_list = []; %whether fixation firing rates were reliably different for in vs out for list
all_peak_list = []; %peak firing rate during list
all_list_peak_time = [];%peak firing time list

%---Firing Rate Curve Properties for Sequence Task---%
sig_p_seq = []; %whether firing rates were significantly different for in vs out for seq
all_peak_seq = []; %peak firing rate during seq
all_seq_peak_times = [];%peak firing time seq
all_seq_firing_curves = [];

%---Firing Rates Between Sequence and List Tasks---%
all_context_gain = []; %change in firing rate (list_fr-seq_fr)/list_fr
all_context_gain2 = [];%gain list_fr/seq_fr


monkeys = {'Vivian','Tobii'};
figure_dir = {};
for monk =2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
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
        
        %load spatial analysis data
        disp(task_file(1:8))
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
            'spatial_info','spike_times','eyepos','binsize','filter_width')
        
        %load Place Cell Fixation Analysis data
        load([data_dir task_file(1:8) '-Place_Cell_Analysis.mat'])
        if numshuffs < 1000
            error('Should have 1000 shuffles')%make sure file has right number of shuffles
        end
        if smval ~= 30
            error('Smoothing Value (2xStd) does not match expectations!')
        end
        
        for unit = 1:num_units
            if ~isempty(spatial_info.shuffled_info_rate{unit})&& length(spatial_info.shuffled_info_rate{unit}) < 1000
                error('Should have 1000 shuffles') %make sure file has right number of shuffles
            elseif isempty(spatial_info.shuffled_info_rate{unit})
               continue  
            end
            
            if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                    && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                
                
                %---Misc Parameters---%
                monkey_all_unit_count(1,monk) = monkey_all_unit_count(1,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                all_place_cell_unit_names = [all_place_cell_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}]; %unit name
                all_place_cell_monkeys = [all_place_cell_monkeys monk]; %1s and 2s for monkey
                place_cell_AP_location = [place_cell_AP_location chamber_zero(1)+ session_data{sess}.location(1)]; %AP location of recorded place cell
                
                subregion = session_data{sess}.subregion;
                vals = textscan(subregion,'%s','Delimiter',',');
                if length(vals{1}) == 1
                    place_cell_subregion = [place_cell_subregion vals{1}];
                else
                    chan = str2double(unit_stats{1,unit}(6));
                    place_cell_subregion = [place_cell_subregion vals{1}(chan)];
                end
                
                %---Spatial Correlation Values for Place cells---%
                all_place_cell_spatial_corrs = [all_place_cell_spatial_corrs spatial_info.spatialstability_halves(unit)];
                all_place_cell_spatial_skagg_percents = [all_place_cell_spatial_skagg_percents spatial_info.shuffled_rate_prctile(unit)];
                

                %---place field properties---%
                H = define_spatial_filter(filter_width); %spatial filter
                trial_data{1} = eyepos{unit}; %eye data
                trial_data{2} = spike_times{unit}; %spike times
                [~,timemaps] = get_firing_rate_map(trial_data,imageX,imageY,binsize,H,Fs,'all'); %get firing rate map
                coverage = [coverage {timemaps}];  %eye data coverage for place cells
                place_coverage = [place_coverage all_place_field_matrix(unit)]; %place field location
                place_field_area = [place_field_area area(unit)]; %size of place field relative to image size
                
                
                %---Place Cell Firing Rate Curves for List---%
                n_in = [n_in sum(in_out{unit} == 1 | in_out{unit} == 2)]; %number of fixations ?->in
                n_out2in = [n_out2in sum(in_out{unit} == 1)];%number of fixaitons out->in
                
                            
                %firing rate out-> in
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 1,:); %get spike trains
                in_curve = nandens(firing_rate,smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
                in_curve2 = in_curve;
                in_curve = in_curve-nanmean(in_curve(1:twin1)); %remove base line
                if ~isnan(stats_across_tasks(1,unit))%peak detected
                    in_curve = in_curve/in_curve(stats_across_tasks(1,unit));%divide by peak firing rate
                else %no peak so normalize by max
                    in_curve = in_curve/max(in_curve);
                end
                all_in_rates = [all_in_rates; in_curve];%fixations out-> in
                all_peak_list = [all_peak_list stats_across_tasks(2,unit)]; %peak firing rate during list
                all_list_peak_times = [all_list_peak_times stats_across_tasks(1,unit)];%time of peak firing rate of place cells for out-> in fixations
                
                %firing rate out-> out
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 4,:); %get spike trains
                out_curve = nandens(firing_rate,smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
                out_curve2 = out_curve;
                out_curve = out_curve-nanmean(out_curve(1:twin1)); %remove base line
                out_curve = out_curve/nanmax(out_curve); %not sure what peak would be so normalize by max
                all_out_rates = [all_out_rates; out_curve];%fixation out->out
                
                in_minus_out = in_curve2-out_curve2;
                %in_minus_out = in_minus_out-mean(in_minus_out(1:twin1));
                in_minus_out = in_minus_out/max(in_minus_out);
                all_in_minus_out_rates = [all_in_minus_out_rates; in_minus_out];
                
                %firing rate for all fixations
                firing_rate = list_fixation_locked_firing{unit};%get spike trains
                firing_rate = nandens(firing_rate,smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
                firing_rate = firing_rate-nanmean(firing_rate(1:twin1)); %remove base line
                firing_rate = firing_rate/nanmax(firing_rate);%not sure what peak would be so normalize by max
                all_fixation_rates = [all_fixation_rates; firing_rate]; %all fixations
                
                
                %find when firing rate for in field fixations is higher than out of field
                %for fixation out->in vs out->out
                if sum(list_sig_times{1,unit}) > 0
                    pass = 1;
                else
                    pass = 0;
                end
                
                %find when firing rate for in field fixations is higher than out of field
                %for fixation ?->in vs ?->out
                if sum(list_sig_times{2,unit}) > 0
                    pass2 = 1;
                else
                    pass2 = 0;
                end
                
                if pass == 1 && pass2 == 1 %both show sig difference
                    sig_p_list = [sig_p_list 3];
                elseif pass == 1 %just out->in vs out-> out shows difference
                    sig_p_list = [sig_p_list 1];
                elseif pass2 == 1 %just ?->in vs ?-> out shows difference
                    sig_p_list = [sig_p_list 2];
                else %no difference shown :(
                    sig_p_list = [sig_p_list 0];
                end
                
                %---Firing Rate Curve Properties for Sequence Task---%
                if ~isnan(stats_across_tasks(4,unit))
                    all_peak_seq = [all_peak_seq stats_across_tasks(4,unit)]; %peak firing rate during seq
                end
                all_seq_peak_times = [all_seq_peak_times stats_across_tasks(3,unit)];%peak firing time seq

                %---Firing Rates Between Sequence and List Tasks---%
                %change in firing rate (list_fr-seq_fr)/list_fr
                all_context_gain = [all_context_gain (stats_across_tasks(2,unit)-stats_across_tasks(4,unit))/stats_across_tasks(2,unit)];
                all_context_gain2 = [all_context_gain2 stats_across_tasks(2,unit)/stats_across_tasks(4,unit)];%gain list_fr/seq_fr
                
                %indexed weirdly so have to do it this way
                all_which_sequence = laundry(all_which_sequence);
                fixation_firing = [];
                for c = 1:4
                    for seq = 1:2
                        fixation_firing = [fixation_firing; sequence_fixation_locked_firing{c,unit}(all_which_sequence{unit}(trial_nums{c,unit}) == seq,:)];
                    end
                end
                
                if ~isempty(in_out_sequence{unit}) %at least 1 item in field and 1 item out of field
                    if  isempty(sequence_sig_times{unit}) %no significant difference
                        sig_p_seq = [sig_p_seq 0];
                        all_seq_firing_curves = [all_seq_firing_curves; NaN(1,twin1+twin2)];
                    else
                        seq_in_out = in_out_sequence{unit};
                        in_curve = nandens(fixation_firing(seq_in_out == 1,:),smval,'gauss',Fs,'nanflt');%smoothed firing rates for in field fixations
                        out_curve = nandens(fixation_firing(seq_in_out == 0,:),smval,'gauss',Fs,'nanflt');%smoothed firing rates for out of field fixations
                        pos_ind = (in_curve-out_curve) > 0; %find expected in field > out of field
                        pos_ind(sequence_sig_times{unit} == 0) = 0; %remove times that are not significant
                        pos_ind = find(pos_ind);
                        neg_ind =  (in_curve-out_curve) < 0 ; %find unexpected in field < out of field
                        neg_ind(sequence_sig_times{unit} == 0) = 0; %remove times that are not significant
                        neg_ind = find(neg_ind); 
                        
                        if isnan(stats_across_tasks(4,unit))
                            all_peak_seq = [all_peak_seq max(in_curve)]; %peak firing rate during seq
                        end
                        
                        in_curve = in_curve-nanmean(in_curve(1:twin1));
                        in_curve = in_curve/max(in_curve);
                        all_seq_firing_curves = [all_seq_firing_curves; in_curve];
                        
                        %find contiguous positive significant regions that
                        %at least 2std in length
                        rmv = [];
                        gaps_pos = findgaps(pos_ind); %find breaks
                        total_pos = 0;%total postive time
                        if ~isempty(gaps_pos)
                            for g = 1:size(gaps_pos,1)
                                gp = gaps_pos(g,:);
                                gp(gp == 0) = [];
                                if length(gp) < 50;%1.5*smval %3  standard deviations
                                    rmv = [rmv g];
                                else
                                    total_pos = total_pos+length(gp);
                                end
                            end
                        end
                        if ~isempty(rmv)
                            gaps_pos(rmv,:) = []; %remove time points that are too short
                        end
                        
                        %find contiguous positive significant regions that
                        %at least 2std in length
                        rmv = [];
                        gaps_neg = findgaps(neg_ind);  %find breaks
                        total_neg = 0;%total negative time
                        if ~isempty(gaps_neg)
                            for g = 1:size(gaps_neg,1)
                                gp = gaps_neg(g,:);
                                gp(gp == 0) = [];
                                if length(gp) < 50;%1.5*smval %3  standard deviations
                                    rmv = [rmv g]; %remove time points that are too short
                                else
                                    total_neg = total_neg+length(gp);
                                end
                            end
                        end
                        if ~isempty(rmv)
                            gaps_neg(rmv,:) = [];
                        end
                        
                        if total_pos > 0 && total_neg == 0 %if only expected result
                            sig_p_seq = [sig_p_seq 1];
                        elseif total_pos == 0 && total_neg == 0 %no result
                            sig_p_seq = [sig_p_seq 0];
                        elseif total_pos > 0 && total_neg > 0 %need to probe further since mixed result
                            
                            %find if negative result is before the start of
                            %the fixation since this could be from ITI or
                            %from another fixation (which could be
                            %influenced by this one). Both cases have been
                            %visually confirmed!
                            rmv = [];
                            for g = 1:size(gaps_neg,1)
                                gp = gaps_neg(g,:);
                                gp(gp == 0) = [];
                                if any(gp < twin1)  %contigous negative started before fixation
                                    rmv = [rmv g];
                                end
                            end
                            gaps_neg(rmv,:) = [];
                            total_neg = sum(sum(gaps_neg > 0));
                            
                            if isempty(gaps_neg) %so only negative time is before fixation so call postive effect
                                sig_p_seq = [sig_p_seq 1];
                            elseif total_pos > 2*total_neg %much more positive than neg so call positive effect
                                sig_p_seq = [sig_p_seq 1];
                                disp('much more positive')
                            elseif 2*total_pos < total_neg %much more negative than positive so call negative effect
                                sig_p_seq = [sig_p_seq -1];
                                disp('much more negative')
                            else %not sure so remove
                                sig_p_seq = [sig_p_seq NaN];
                                disp('unsure so removed')
                            end
                        elseif total_neg > 0 %if only negative result
                            sig_p_seq = [sig_p_seq -1];
                        else
                            error('What else could happen')
                        end
                    end
                else
                    sig_p_seq = [sig_p_seq NaN]; %no shapes in and out of field so keep parallel structure
                    all_peak_seq = [all_peak_seq NaN];
                    all_seq_firing_curves = [all_seq_firing_curves; NaN(1,twin1+twin2)];
                end
                
                if~isnan(sig_p_seq(end))
                    if isnan(all_context_gain2(end)) %can happen if no definable peaks so set as ratio of maximum firing rates
                        all_context_gain2(end) =  all_peak_list(end)/all_peak_seq(end);
                    elseif all_context_gain2(end) ~= all_peak_list(end)/all_peak_seq(end);
                        error('Contextual gain 2 should be list peak/sequence peak')
                    end
                    
                    if isnan(all_context_gain(end)) %can happen if no definable peaks so set as ratio of maximum firing rates
                        all_context_gain(end) = (all_peak_list(end)-all_peak_seq(end))/all_peak_list(end);
                    elseif  all_context_gain(end) ~= (all_peak_list(end)-all_peak_seq(end))/all_peak_list(end)
                        error('Context gain should be the change in firing rate')
                    end
                end
                
                
            elseif ~isnan(spatial_info.shuffled_rate_prctile(unit))
                
                %---Misc Parameters---%
                monkey_all_unit_count(2,monk) = monkey_all_unit_count(2,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                
                %---Spatial Correlation Values for Non-Place cells---%
                all_non_place_cell_spatial_corrs = [all_non_place_cell_spatial_corrs spatial_info.spatialstability_halves(unit)];
                all_non_place_cell_skagg_percents = [all_non_place_cell_skagg_percents spatial_info.shuffled_rate_prctile(unit)];
                non_place_skagg_names = [non_place_skagg_names  {[task_file(1:end-11) '_' unit_stats{1,unit}]}]; %unit name
                
                subregion = session_data{sess}.subregion;
                vals = textscan(subregion,'%s','Delimiter',',');
                if length(vals{1}) == 1
                    non_place_cell_subregion = [non_place_cell_subregion vals{1}];
                else
                    chan = str2double(unit_stats{1,unit}(6));
                    non_place_cell_subregion = [non_place_cell_subregion vals{1}(chan)];
                end
            end
        end
    end
end
%% Copy Place Cell Figures to Summary Direction
clc

%---First Display Summary Results---%
disp(['Found ' num2str(length(all_list_peak_times)) ' place cells']);%number of place cells
disp(['Found ' num2str(sum(sig_p_list == 0)) ' place cells do not show significant increase in firing rate in field...removing'])
nans = find(sig_p_list == 0);
for i = 1:length(nans)
    disp(['Removing ' all_place_cell_unit_names{nans(i)}])
end
disp(['Found ' num2str(sum(isnan(all_list_peak_times))) ' place cells do not show a peak in time...removing'])
nans = find(isnan(all_list_peak_times));
for i = 1:length(nans)
    disp(['Removing ' all_place_cell_unit_names{nans(i)}])
end
%%
%---Second Copy Relevant Figures to Summary Directory---%
% for unit = 1:length(all_place_cell_unit_names)
%     sub_dir1 = 'Place Cells Fixation Analysis\Best Place Cells Fixation Analysis\';
%     sub_dir2 = 'Spatial Analysis\';
%     
%     name1 = [all_place_cell_unit_names{unit} '_place_cell_fixation_analysis.png'];
%     name2 = [all_place_cell_unit_names{unit} '_place_cell_Sequence_InSideOutside.png'];
%     name3 = [all_place_cell_unit_names{unit} '_List_spatial_analysis.png'];
%     
%     if sig_p_list(unit) == 0 %not reliable
%         copyfile([figure_dir{all_place_cell_monkeys(unit)} sub_dir1 name1],...
%             [summary_directory 'Not reliable\' name1])
%         if ~isnan(sig_p_seq(unit))
%             copyfile([figure_dir{all_place_cell_monkeys(unit)} sub_dir1 name2],...
%                 [summary_directory 'Not reliable\' name2])
%         end
%         copyfile([figure_dir{all_place_cell_monkeys(unit)} sub_dir2 name3],...
%             [summary_directory 'Not reliable\' name3])
%     elseif isnan(all_list_peak_times(unit))%no peak in firing rate detected
%         copyfile([figure_dir{all_place_cell_monkeys(unit)} sub_dir1 name1],...
%             [summary_directory 'No peak\' name1])
%         if ~isnan(sig_p_seq(unit))
%             copyfile([figure_dir{all_place_cell_monkeys(unit)} sub_dir1 name2],...
%                 [summary_directory 'No peak\' name2])
%         end
%         copyfile([figure_dir{all_place_cell_monkeys(unit)} sub_dir2 name3],...
%             [summary_directory 'No peak\' name3])
%     else %passes all criterion
%         copyfile([figure_dir{all_place_cell_monkeys(unit)} sub_dir1 name1],...
%             [summary_directory name1])
%         if ~isnan(sig_p_seq(unit))
%             copyfile([figure_dir{all_place_cell_monkeys(unit)} sub_dir1 name2],...
%                 [summary_directory name2])
%         end
%         copyfile([figure_dir{all_place_cell_monkeys(unit)} sub_dir2 name3],...
%             [summary_directory name3])
%     end
% end
%% Remove Un-Reliable Neurons and Neurons wthout a definitve peak in firing rate

%remove neurons without significant difference in firnig rate between in
%and out of field
all_in_rates(sig_p_list == 0,:) = [];
all_fixation_rates(sig_p_list == 0,:) = [];
all_out_rates(sig_p_list == 0,:) = [];
all_place_cell_unit_names(sig_p_list == 0) = [];
coverage(sig_p_list == 0) = [];
place_coverage(sig_p_list == 0) = [];
place_field_area(sig_p_list == 0) = [];
all_context_gain(sig_p_list == 0) = [];
all_context_gain2(sig_p_list == 0) = [];
all_peak_list(sig_p_list == 0) = [];
all_peak_seq(sig_p_list == 0) = [];
all_seq_peak_times(sig_p_list == 0) = [];
sig_p_seq(sig_p_list == 0) = [];
all_list_peak_times(sig_p_list == 0) = [];
all_in_minus_out_rates(sig_p_list == 0,:) = [];
place_cell_AP_location(sig_p_list == 0) = [];
place_cell_subregion(sig_p_list == 0) = [];
all_seq_firing_curves(sig_p_list == 0,:) =[];
sig_p_list(sig_p_list == 0) = [];

%remove neurons without a definitive peak
all_in_rates(isnan(all_list_peak_times),:) = [];
all_fixation_rates(isnan(all_list_peak_times),:) = [];
all_out_rates(isnan(all_list_peak_times),:) = [];
all_place_cell_unit_names(isnan(all_list_peak_times)) = [];
coverage(isnan(all_list_peak_times)) = [];
place_coverage(isnan(all_list_peak_times)) = [];
place_field_area(isnan(all_list_peak_times)) = [];
all_context_gain(isnan(all_list_peak_times)) = [];
all_context_gain2(isnan(all_list_peak_times)) = [];
all_peak_seq(isnan(all_list_peak_times)) = [];
all_seq_peak_times(isnan(all_list_peak_times)) = [];
sig_p_seq(isnan(all_list_peak_times)) = [];
sig_p_list(isnan(all_list_peak_times)) = [];
all_peak_list(isnan(all_list_peak_times)) = [];
all_in_minus_out_rates(isnan(all_list_peak_times),:) = [];
place_cell_AP_location(isnan(all_list_peak_times)) = [];
place_cell_subregion(isnan(all_list_peak_times)) = [];
all_seq_firing_curves(isnan(all_list_peak_times),:) = [];
all_list_peak_times(isnan(all_list_peak_times)) = [];
%% Plot Psuedo-Population Firing Rate Curves
figure

%---Firing Rate Curves for out->in fixations---%
subplot(2,2,1)
vals = all_in_rates(:,1:twin1); %"baseline" out of field firing rate
[~,place_order] = sort(all_list_peak_times); %sort order by peak firing time
imagesc([-twin1:twin2-1],[1:size(all_in_rates,1)],all_in_rates(place_order,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_in_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('out -> in')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar

%---Firing Rate Curves for out->out fixations---%
%sorted same order as above
subplot(2,2,2)
vals = all_out_rates(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_out_rates,1)],all_out_rates(place_order,:)) %sorted same order as above
colormap('jet')
hold on
plot([0 0],[1 size(all_out_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('out -> out')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
colorbar

%---Firing Rate Curves for all fixations---%
%sorted same order as above
subplot(2,2,3)
vals = all_fixation_rates(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_fixation_rates,1)],all_fixation_rates(place_order,:)) %sorted same order as above
colormap('jet')
hold on
plot([0 0],[1 size(all_fixation_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('All Fixations')
caxis([-std(vals(:)) 1])  %set minumun to standard deviation of baseline since some neurons are greatly inhibited
colorbar

%---Firing Rate Curves for all fixations---%
%sorted on own
subplot(2,2,4)
[~,i] = max(all_fixation_rates,[],2); %sort by time of maximum firing rate
[~,all_order] = sort(i);
imagesc([-twin1:twin2-1],[1:size(all_fixation_rates,1)],all_fixation_rates(all_order,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_fixation_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('All Fixations Sorted Seperately')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
colorbar

%% Population Average in Field Firing Rate Curve
[m,i] = max(nanmean(all_in_rates));
figure
plot([-twin1:twin2-1],nanmean(all_in_rates))
hold on
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--');
plot([-44 -44],[yl(1) yl(2)],'k--')
plot([-twin1 twin2],[0 0],'k')
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Normalized Firing Rate')
title(['Average Place Cell Fixation Aligned Firing Rate, peak @ ' num2str(i-twin1) ' ms'])
xlim([-twin1 twin2])
box off
axis square

%% Distribution of Peak Times and Calculate FWHM

FWHM = NaN(1,size(all_in_rates,1));
for n = 1:size(all_in_rates,1)
    cv = all_in_rates(n,:);
    peak = all_list_peak_times(n);
    tims = find(cv < 0.5*cv(peak));
    diffs = tims-peak;
    neg_diffs = diffs(diffs < 0);
    pos_diffs = diffs(diffs > 0);
    neg = min(abs(neg_diffs));
    pos = min(pos_diffs);
    
    if isempty(pos)
        pos = twin2;
        FWHM(n) = NaN;
    else
        FWHM(n) = neg+pos;
    end
    
end
%%

figure
subplot(1,2,1)
histogram(all_list_peak_times-twin1,22)
xlabel('Peak Delay from Fixation Start (ms)')
ylabel('Neruon Count')
hold on
plot([median(all_list_peak_times-twin1) median(all_list_peak_times-twin1)],[0 10],'k--')
title(['Median Peak Delay From Fixation Start = ' num2str(round(median(all_list_peak_times-twin1))) ' ms'])
xlim([-twin1 twin2])
box off
axis square

subplot(1,2,2)
histogram(FWHM,25)
xlabel('FWHM @ peak')
ylabel('Count')
title(['Median FWHM = ' num2str(nanmedian(FWHM),3) ' ms'])
xlim([25 275])
box off
axis square


%% Place Cells: Eye Coverage, Field Coverage, Field Sizes
all_coverage = zeros(size(coverage{1}));
coverage_count =  zeros(size(coverage{1}));
for c = 1:length(coverage)
    cv = coverage{c};
    coverage_count(~isnan(cv)) =  coverage_count(~isnan(cv))+1;
    cv(isnan(cv)) = 0;
    all_coverage = all_coverage+cv;
end
all_coverage = all_coverage./coverage_count;
all_coverage = all_coverage/sum(sum(all_coverage));

all_place_coverage  = zeros(imageY,imageX);
coverage_count =  zeros(imageY,imageX);
for c = 1:length(place_coverage);
    cv = place_coverage{c};
    cov = coverage{c};
    cov = imresize(cov,[imageY,imageX],'method','nearest');%upsample to images size in pixel coordinates
    coverage_count(~isnan(cov)) =  coverage_count(~isnan(cov))+1;
    %      all_place_coverage2 = all_place_coverage2+cv;
    cv(isnan(cov)) = 0;
    all_place_coverage = all_place_coverage+cv;
end
all_place_coverage  = all_place_coverage./coverage_count;
all_place_coverage = all_place_coverage/sum(sum(all_place_coverage));
%%
figure
subplot(2,2,1)
imagesc(all_coverage)
axis off
colormap('jet')
title('Eye Data Coverage')
axis equal
colorbar

subplot(2,2,2)
imagesc(all_place_coverage)
axis off
colormap('jet')
title('View Cell Field Coverage')
p99 = prctile(all_place_coverage(:),99);
clims = caxis;
caxis([clims(1) p99]);
axis equal
colorbar

subplot(2,2,3)
histogram(place_field_area,15)
yl = ylim;
hold on
plot([median(place_field_area) median(place_field_area)],[0 yl(2)],'r--')
hold off
xlabel('Field Area (% of Screen Space)')
ylabel('Count')
title(['Median Field area = ' num2str(median(place_field_area),2)])

%% Plot Distribution of Spatial Correlations for Place and Non-Place Cells
figure
histogram(all_place_cell_spatial_corrs,'facecolor','k','facealpha',1,'edgecolor','none','Normalization','probability')
hold on
histogram(all_non_place_cell_spatial_corrs,'facecolor','g','facealpha',1,'edgecolor','none','Normalization','probability');
yl = ylim;
plot([median(all_place_cell_spatial_corrs) median(all_place_cell_spatial_corrs)],[0 yl(2)],'k--')
plot([nanmedian(all_non_place_cell_spatial_corrs) nanmedian(all_non_place_cell_spatial_corrs)],[0 yl(2)],'m--')
hold off
xlabel('Spatial Correlation First vs Second Half')
ylabel('Proportion')
legend('Place Cells','Non-Place Cells')
title(['Median_{place cell} \rho_{1/2} = ' num2str(median(all_place_cell_spatial_corrs),2) ...
    ', Median_{non-place cell} \rho_{1/2} = ' num2str(nanmedian(all_non_place_cell_spatial_corrs),2)])
box off
axis square
[~,p] = ttest2(all_place_cell_spatial_corrs,all_non_place_cell_spatial_corrs);
%% Population Contextual Gain

disp([num2str(sum(~isnan(sig_p_seq))) ' View Cells with 1+ item in field & 1+ item out of field'])
disp([num2str(sum(sig_p_seq == 1)) ' (' num2str(100*sum(sig_p_seq == 1)/sum(~isnan(sig_p_seq)),2) '%) View Cells show expected effect'])
disp([num2str(sum(sig_p_seq == 0)) ' (' num2str(100*sum(sig_p_seq == 0)/sum(~isnan(sig_p_seq)),2) '%) View Cells show no effect'])
disp([num2str(sum(sig_p_seq == -1)) ' (' num2str(100*sum(sig_p_seq == -1)/sum(~isnan(sig_p_seq)),2) '%) View Cells show opposite effect'])

%%
figure
subplot(2,2,1)
hist(100*all_context_gain,24)
xlim([-300 100])
xlabel('Task Preference (% Change)')
ylabel('Count')
box off
axis square
title(['Median = ' num2str(nanmedian(100*all_context_gain),2) ' %'])

acg2 = all_context_gain2;
acg2(acg2 > 4) = 5;
acg2(acg2 < 1) = -1./acg2(acg2 <1);

subplot(2,2,2)
histogram(acg2,25)
xlabel('Task Gain (Relative)')
ylabel('Count')
set(gca,'Xtick',[-4 -3 -2 -1 0 1 2 3 4 5])
box off
axis square
title(['|Median| = ' num2str(nanmedian(abs(acg2)),2)])
xlim([-4 5])


%%
%---Firing Rate Curves for out->in fixations---%
vals = all_in_minus_out_rates(:,1:twin1); %"baseline" out of field firing rate
vals = vals(:);
vals(vals > 0) = [];

figure
[~,place_order] = sort(all_list_peak_times); %sort order by peak firing time
imagesc([-twin1:twin2-1],[1:size(all_in_minus_out_rates,1)],all_in_minus_out_rates(place_order,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_in_minus_out_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('out -> in')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar

figure
[~,mxi] = max(all_in_minus_out_rates');
[~,place_order] = sort(mxi); %sort order by peak firing time
imagesc([-twin1:twin2-1],[1:size(all_in_minus_out_rates,1)],all_in_minus_out_rates(place_order,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_in_minus_out_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('out -> in')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar

names = all_place_cell_unit_names(place_order);
%%
for n = 1:length(names)
    if strcmpi(names{n},'TO160515_3_sig004d')
        disp(n)
    end
end
    
%% AP axis vs Latency of Response
avg_delay = zeros(1,30);
for b = 0.5:0.5:15
    if sum(place_cell_AP_location == b) >= 1
    avg_delay(2*b) = mean(all_list_peak_times(place_cell_AP_location == b));
    else
         avg_delay(2*b) = NaN;
    end
end
    
figure
plot(place_cell_AP_location,all_list_peak_times-twin1,'k.')
hold on
plot(0.5:0.5:15,avg_delay-twin1,'*r')
hold off
%% Subregion vs Latency of Response
place_region_id = [];
non_place_region_id = [];

for p = 1:length(place_cell_subregion)
    if strcmpi(place_cell_subregion{p},'DG')
        place_region_id(p) = 0;
    elseif strcmpi(place_cell_subregion{p},'CA1')
        place_region_id(p) = 1;
    elseif strcmpi(place_cell_subregion{p},'CA2')
        place_region_id(p) = 2;
    elseif strcmpi(place_cell_subregion{p},'CA3')
        place_region_id(p) = 3;
    elseif strcmpi(place_cell_subregion{p},'CA4')
        place_region_id(p) = 4;
    elseif strcmpi(place_cell_subregion{p},'ProSub') || strcmpi(place_cell_subregion{p},'ProS')
        place_region_id(p) = 5;
    else
        disp('what')
        disp(place_cell_subregion{p})
    end
end
for p = 1:length(non_place_cell_subregion)
    if strcmpi(non_place_cell_subregion{p},'DG')
        non_place_region_id(p) = 0;
    elseif strcmpi(non_place_cell_subregion{p},'CA1')
        non_place_region_id(p) = 1;
    elseif strcmpi(non_place_cell_subregion{p},'CA2')
        non_place_region_id(p) = 2;
    elseif strcmpi(non_place_cell_subregion{p},'CA3')
        non_place_region_id(p) = 3;
    elseif strcmpi(non_place_cell_subregion{p},'CA4')
        non_place_region_id(p) = 4;
    elseif strcmpi(non_place_cell_subregion{p},'ProSub') || strcmpi(non_place_cell_subregion{p},'ProS')
        non_place_region_id(p) = 5;
    else
        disp('what')
        disp(non_place_cell_subregion{p})
    end
end
%%
%combine DG & CA4
place_region_id(place_region_id == 4) = 0;

%combine CA2 & CA3
place_region_id(place_region_id == 2) = 3;

%combine ProSub and CA1
place_region_id(place_region_id == 5) = 1;
%%
%% Plot Relationship between Peak times in List and Peak Times in Sequence
x = all_list_peak_times(sig_p_seq == 1)-twin1;
y = all_seq_peak_times(sig_p_seq == 1)-twin1;

% fr_too_slow = find((all_peak_list(sig_p_seq == 1) < 2) | (all_peak_seq(sig_p_seq == 1) < 2));
% x(fr_too_slow) = NaN;
% y(fr_too_slow) = NaN;

x(isnan(y)) = [];
y(isnan(y)) = [];
[r,p] = corrcoef(x,y);
% [r,p] = corr(x,y);

figure
plot(x,y,'.k')
xlabel('List Delay to Peak from Fixation Start (ms)')
ylabel('Sequence Delay to Peak from Fixation Start (ms)')
title('Significant Response in Sequence and List')

%%
%---Plot visually now---%
[~,place_order] = sort(all_list_peak_times); %sort order by peak firing time

list_ordered = all_in_rates(place_order,:);
seq_ordered = all_seq_firing_curves(place_order,:);
sig_p_ordered = sig_p_seq(place_order);

list_ordered = list_ordered(sig_p_ordered == 1,:);
seq_ordered = seq_ordered(sig_p_ordered == 1,:);

ordered_peaks = all_seq_peak_times(place_order);
ordered_peaks = ordered_peaks(sig_p_ordered == 1);
nanseq = find(isnan(ordered_peaks));
list_ordered(nanseq,:) = [];
seq_ordered(nanseq,:) = [];

listvals = list_ordered(:,1:twin1); %"baseline" out of field firing rate
seqvals = seq_ordered(:,1:twin1); %"baseline" out of field firing rate

figure
subplot(2,2,1)
imagesc([-twin1:twin2-1],[1:size(list_ordered,1)],list_ordered)
colormap('jet')
hold on
plot([0 0],[1 size(list_ordered,1)],'w--');
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Neuron #')
title('Images sorted by Image Peak Times')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar
 axis square

subplot(2,2,2)
imagesc([-twin1:twin2-1],[1:size(seq_ordered,1)],seq_ordered)
colormap('jet')
hold on
plot([0 0],[1 size(seq_ordered,1)],'w--');
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Neuron #')
title('Sequences sorted by Image Peak Times')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar
 axis square
 
[~,place_order] = sort(all_seq_peak_times); %sort order by peak firing time


list_ordered = all_in_rates(place_order,:);
seq_ordered = all_seq_firing_curves(place_order,:);
sig_p_ordered = sig_p_seq(place_order);

list_ordered = list_ordered(sig_p_ordered == 1,:);
seq_ordered = seq_ordered(sig_p_ordered == 1,:);

ordered_peaks = all_seq_peak_times(place_order);
ordered_peaks = ordered_peaks(sig_p_ordered == 1);
nanseq = find(isnan(ordered_peaks));
list_ordered(nanseq,:) = [];
seq_ordered(nanseq,:) = [];

listvals = list_ordered(:,1:twin1); %"baseline" out of field firing rate
seqvals = seq_ordered(:,1:twin1); %"baseline" out of field firing rate


subplot(2,2,3)
imagesc([-twin1:twin2-1],[1:size(list_ordered,1)],list_ordered)
colormap('jet')
hold on
plot([0 0],[1 size(list_ordered,1)],'w--');
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Neuron #')
title('Images sorted by Sequence Peak Times')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar
axis square

subplot(2,2,4)
imagesc([-twin1:twin2-1],[1:size(seq_ordered,1)],seq_ordered)
colormap('jet')
hold on
plot([0 0],[1 size(seq_ordered,1)],'w--');
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Neuron #')
title('Sequences sorted by Sequence Peak Times')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar
axis square

subtitle(['Correlation between peak response: r = ' num2str(r(2),2) ', p = ' num2str(p(2),2)])
%%
