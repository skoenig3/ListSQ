%get recording locations for view cells
clar
monkeys = {'Vivian','Tobii'};

task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000;
imageX = 800;
imageY = 600;

%---Spatial Correlation Values for all cells---%
all_spatial_corrs = []; %place cell spatial correlations
all_not_spatial_corrs = [];%non-place cell spatial correlations
all_peak_times = [];
sig_p_list = [];

%---Other Parameters (Mostly Place Cells)---%
monkey_count = zeros(2,2);%row 1 place row 2 non-place
all_unit_names = {}; %place cell unit names
all_monkeys = []; %1s and 2s
AP_location = []; %AP location of recorded place cell

%---Novel Repeat Difference---%
all_place_cell_short = [];
all_place_cell_long = [];
not_place_cell_short = [];
not_place_cell_long = [];
early_short = [];
early_long = [];

for monk = 2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
        
        
%         listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
%         listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    end
    
    for sess = 1:length(session_data)
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials',...
            'stability_attribute','cfg','hdr','data');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        
        if num_units == 0
            continue
        end
        
        num_trials = length(cfg.trl);
        valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
        if all(isnan(valid_trials(1,:)))
            continue
        end
        
        disp(task_file(1:8))
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
            'spatial_info','spike_times','eyepos','binsize','filter_width')
        
        load([data_dir task_file(1:8) '-Place_Cell_Analysis.mat'],...
            'twin1','twin2','list_95_curve','list_fixation_locked_firing','smval',...
            'in_out','all_place_field_matrix','numshuffs','contextual_gain',...
            'sequence_fixation_locked_firing','sequence_95_curve','in_out_sequence')
        
        if numshuffs < 1000
            error('Should have 1000 shuffles')%make sure file has right number of shuffles
        end
        
        
        load([data_dir task_file(1:8) '-ListSQ-Visual_Response_Memory_results.mat']);
        if numshuffs < 5000
            error('Should have 5000 shuffles for memory')%make sure file has right number of shuffles
        end
        
        for unit = 1:num_units
            if ~isempty(spatial_info.shuffled_info_rate{unit})&& length(spatial_info.shuffled_info_rate{unit}) < 1000
                error('Should have 1000 shuffles') %make sure file has right number of shuffles
            end
            
            if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                    && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                %                     && (spatial_info.spatialstability_even_odd_prctile(unit) > 95) % ... %spatial consistency
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%---Firing Rate Curves During List Image Viewing---%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %---Misc Parameters---%
                monkey_count(1,monk) = monkey_count(1,monk)+1;
                all_unit_names = [all_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}];
                all_monkeys = [all_monkeys monk]; %1s and 2s
                AP_location = [AP_location chamber_zero(1)+ session_data{sess}.location(1)]; %AP location of recorded place cell
                all_spatial_corrs = [all_spatial_corrs spatial_info.spatialstability_halves(unit)];
                
                
                
                %1) first fixation in: out-> in
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 1,:);
                in_curve = nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 4,:);
                out_curve = nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                
                sig_ind = find((in_curve-out_curve) > list_95_curve{1,unit});
                gaps = findgaps(sig_ind);
                pass = 0;
                if ~isempty(gaps)
                    for g = 1:size(gaps,1)
                        gp = gaps(g,:);
                        gp(gp == 0) = [];
                        if length(gp) > 30
                            pass = 1;
                        end
                    end
                end
                
                %find peaks
                in_curve = in_curve-nanmean(in_curve(1:twin1));
                [PKS,LOCS]= findpeaks(in_curve,'MinPeakWidth',40);
                if ~isempty(LOCS)
                    %remove peaks less than 1/2 the max
                    LOCS(PKS < 0.66*max(in_curve)) = [];
                    PKS(PKS < 0.66*max(in_curve)) = [];
                end
                if isempty(LOCS)
                    LOCS = NaN;
                end
                if isnan(LOCS)
                    all_peak_times = [all_peak_times NaN];
                else
                    %take the first peak
                    all_peak_times = [all_peak_times LOCS(1)];
                end
                
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 1 | in_out{unit} == 2,:);
                in_curve = nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 3 | in_out{unit} == 4,:);
                out_curve = nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                sig_ind2 = find((in_curve-out_curve) > list_95_curve{2,unit});
                pass2 = 0;
                gaps = findgaps(sig_ind2);
                if ~isempty(gaps)
                    for g = 1:size(gaps,1)
                        gp = gaps(g,:);
                        gp(gp == 0) = [];
                        if length(gp) > 30
                            pass2 = 1;
                        end
                    end
                end
                
                if pass2 == 1;
                    sig_p_list = [sig_p_list 1];
                elseif pass == 1;
                    sig_p_list = [sig_p_list 2];
                else
                    sig_p_list = [sig_p_list 0];
                end
                
                %Novel/Repeat Differences Short
                gaps = findgaps(find(sig_short{unit}));
                rmv = [];
                too_early = 0;
                if ~isempty(gaps)
                    for g = 1:size(gaps,1)
                        gp = gaps(g,:);
                        gp(gp == 0) = [];
                        if length(gp) < 60
                            rmv = [rmv g];
                        elseif max(gp) < 200
                            rmv = [rmv g];
                            too_early = 1;
                        end
                    end
                end
                gaps(rmv,:) = [];
                if ~isempty(gaps)
                    all_place_cell_short = [all_place_cell_short 1];
                else
                    all_place_cell_short = [all_place_cell_short 0];
                end
                
                if too_early
                    early_short = [early_short 1];
                else
                    early_short = [early_short 0];
                end
                
                %Novel/Repeat Differences long
                gaps = findgaps(find(sig_long{unit}));
                rmv = [];
                too_early = 0;
                if ~isempty(gaps)
                    for g = 1:size(gaps,1)
                        gp = gaps(g,:);
                        gp(gp == 0) = [];
                        if length(gp) < 300
                            rmv = [rmv g];
                        elseif max(gp) < 200
                            rmv = [rmv g];
                            too_early = 1;
                        end
                    end
                end
                gaps(rmv,:) = [];
                if ~isempty(gaps)
                    all_place_cell_long = [all_place_cell_long 1];
                else
                    all_place_cell_long = [all_place_cell_long 0];
                end
                
                if too_early
                    early_long = [early_long 1];
                else
                    early_long = [early_long 0];
                end
                
            elseif ~isnan(spatial_info.shuffled_rate_prctile(unit))
                monkey_count(2,monk) = monkey_count(2,monk)+1;
                all_not_spatial_corrs = [all_not_spatial_corrs spatial_info.spatialstability_halves(unit)];
                
                
                %Novel/Repeat Differences Short
                gaps = findgaps(find(sig_short{unit}));
                rmv = [];
                too_early = 0;
                if ~isempty(gaps)
                    for g = 1:size(gaps,1)
                        gp = gaps(g,:);
                        gp(gp == 0) = [];
                        if length(gp) < 30
                            rmv = [rmv g];
                        elseif max(gp) < 200
                            rmv = [rmv g];
                            too_early = 1;
                        end
                    end
                end
                gaps(rmv,:) = [];
                if ~isempty(gaps)
                    not_place_cell_short = [not_place_cell_short 1];
                else
                    not_place_cell_short = [not_place_cell_short 0];
                end
                
                if too_early
                    early_short = [early_short 1];
                else
                    early_short = [early_short 0];
                end
                
                %Novel/Repeat Differences long
                gaps = findgaps(find(sig_long{unit}));
                rmv = [];
                too_early = 0;
                if ~isempty(gaps)
                    for g = 1:size(gaps,1)
                        gp = gaps(g,:);
                        gp(gp == 0) = [];
                        if length(gp) < 60
                            rmv = [rmv g];
                        elseif max(gp) < 200
                            rmv = [rmv g];
                            too_early = 1;
                        end
                    end
                end
                gaps(rmv,:) = [];
                if ~isempty(gaps)
                    not_place_cell_long = [not_place_cell_long 1];
                else
                    not_place_cell_long = [not_place_cell_long 0];
                end
                
                if too_early
                    early_long = [early_long 1];
                else
                    early_long = [early_long 0];
                end
                
            end
        end
    end
end
%%
clc
disp(['Found ' num2str(length(all_peak_times)) ' place cells'])
disp(['Found ' num2str(sum(isnan(all_peak_times))) ' place cells do not show a peak in time...removing'])
nans = find(isnan(all_peak_times));
for i = 1:length(nans)
    disp(['Removing ' all_unit_names{nans(i)}])
end
disp(['Found ' num2str(sum(sig_p_list == 0)) ' place cells do not show significant increase in firing rate in field...removing'])
nans = find(sig_p_list == 0);
for i = 1:length(nans)
    disp(['Removing ' all_unit_names{nans(i)}])
end
%%
% %remove neurons without a definitve peak in time
all_in_rates(isnan(all_peak_times),:) = [];
all_fixation_rates(isnan(all_peak_times),:) = [];
all_out_rates(isnan(all_peak_times),:) = [];
all_unit_names(isnan(all_peak_times)) = [];
coverage(isnan(all_peak_times)) = [];
place_coverage(isnan(all_peak_times)) = [];
all_areas(isnan(all_peak_times)) = [];
all_context_gain(isnan(all_peak_times)) = [];
all_context_gain2(isnan(all_peak_times)) = [];
all_peak_seq(isnan(all_peak_times)) = [];
all_seq_peak_time(isnan(all_peak_times)) = [];
sig_p_seq(isnan(all_peak_times)) = [];
all_list_peak_time(isnan(all_peak_times)) = [];
sig_p_list(isnan(all_peak_times)) = [];
all_peak_list(isnan(all_peak_times)) = [];
all_peak_times(isnan(all_peak_times)) = [];


%
% %remove neurons without significant difference in firnig rate between in
% %and out of field
all_in_rates(sig_p_list == 0,:) = [];
all_fixation_rates(sig_p_list == 0,:) = [];
all_out_rates(sig_p_list == 0,:) = [];
all_unit_names(sig_p_list == 0) = [];
coverage(sig_p_list == 0) = [];
place_coverage(sig_p_list == 0) = [];
all_areas(sig_p_list == 0) = [];
all_peak_times(sig_p_list == 0) = [];
all_context_gain(sig_p_list == 0) = [];
all_context_gain2(sig_p_list == 0) = [];
all_peak_list(sig_p_list == 0) = [];
all_peak_seq(sig_p_list == 0) = [];
all_seq_peak_time(sig_p_list == 0) = [];
sig_p_seq(sig_p_list == 0) = [];
all_list_peak_time(sig_p_list == 0) = [];
sig_p_list(sig_p_list == 0) = [];