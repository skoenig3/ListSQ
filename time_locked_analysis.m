function time_locked_analysis(data_dir,preprocessed_data_file,figure_dir,task)
% written by Seth Konig August, 2014
% Function analyses spike times locked to events that occur on the monitor
% at the time dicated by cortex event codes. Analysis does not analyze eye
% movements directly but when cortex says the eye had entered the fixation
% window. Updates to come on this.
%
% Inputs:
%   1) data_dir: direction where preprocessed data is located
%   2) preprocessed_data_file: preprocessed ListSQ data file with LFPs
%   divided by channel and trial.
%   3) figure_dir: location of where to put save figures
%   4) task: what task was this data come from with i.e. 'cvtnew','ListSQ'
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-time_locked_results'

twin = 500;% how much time to take before and after an event.
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
numshuffs = 100; %recommend this is between 100 & 1000, for bootstrapping to
% dtermine significance

%CVTNEW variable trial lenght
%time is + 300 ms
short = [700 1133]; %short duration trials
mediumm = [1134 1567];
long = [1568 2500];%cap should be 2000 but cortex can have some lag

switch task
    case 'cvtnew'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import & reformat data so that spikes are locked to events---%%%
        
        load([data_dir preprocessed_data_file],'data','cfg','multiunit','num_units');
        num_trials = length(cfg.trl);
        
        disp('Aligning spike times to trial events')
        
        %get important task specific information
        [event_names,event_codes,event_durs,event_t0] = define_events(task,preprocessed_data_file(3:8));
        
        %preallocate space and parallel structure of cfg
        time_lock_firing = make_time_lock_cell(event_durs,num_units,num_trials,twin);
        trial_duration = NaN(1,num_trials);
        for t = 1:num_trials
            if any(cfg.trl(t).allval == event_codes(8)); %rewarded
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(1));
                trial_duration(t) = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(6))-cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(5));
                event_times = get_cvtnew_event_times(cfg.trl(t),event_codes,event_durs,trial_start,twin);
                for unit = 1:num_units
                    time_lock_firing = put_in_time_matrix2(data,unit,t,event_times,time_lock_firing,twin,event_names);
                end
            end
        end
        %remove excess NaNs associated with error trials
        time_lock_firing = laundry(time_lock_firing);
        trial_duration = laundry(trial_duration);
        
        if size(time_lock_firing,2) < long(2)
            long(2) = size(time_lock_firing,2);
        end
        
        
        %commented out SDK 1/8/16
        %         %remove long times for dot (event 4) since highly variable
        %         for n = 1:size(time_lock_firing,2)
        %             total_trials = size(time_lock_firing{4,n},1);%total number of trials
        %             trials_per_time =  sum(~isnan(time_lock_firing{4}));%number of trials per sample
        %             rmv_time_points = find(trials_per_time < 0.25*total_trials);
        %             time_lock_firing{4,n}(:,rmv_time_points) = NaN;
        %             time_lock_firing{4,n} = laundry(time_lock_firing{4,n})
        %         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Perform Temporal analysis and signifance testing---%%%
        disp('Determining if spikes are locked to trial events')
        Fs = data(1).fsample; %should be 1000
        
        smval = 60; %temporal 1/2 width of gaussian smoothing filter
        
        % Collect spikes per epoch to determine if neurons fire more within
        % a given epoch ignore the 1st twin as this is from the previous event
        epoch_data.firing_rates = cell(size(time_lock_firing,1),num_units);
        epoch_data.dur = NaN(size(time_lock_firing,1),num_units);
        epoch_data.num_trials = NaN(size(time_lock_firing,1),num_units);
        epoch_data.trial_duration = trial_duration; %length of the dot on
        
        %info per spike is equivalent to dividing information rate by average firing rate
        temporal_info.rate = NaN(size(time_lock_firing,1),num_units); %the observed information rate in bits/sec
        temporal_info.shuffled_info_rate = cell(size(time_lock_firing,1),num_units); %bootstrapped information rate in bits/sec expected by chance from spike train
        temporal_info.shuffled_95_percentile = NaN(size(time_lock_firing,1),num_units);
        temporal_info.shuffled_90_percentile = NaN(size(time_lock_firing,1),num_units);
        info_type = 'temporal';
        
        % perform mutual information theory analysis to determine
        % if neurons encode info about a particular time within an epoch
        for event = 1:size(time_lock_firing,1);
            for unit = 1:num_units
                [temporal_info.rate(event,unit),temporal_info.shuffled_info_rate{event,unit}]...
                    = estimated_mutual_information(time_lock_firing{event,unit},numshuffs,info_type,smval,Fs);
                epoch_data.firing_rates{event,unit} = nansum(time_lock_firing{event,unit}(:,twin:end),2);
                epoch_data.dur(event,unit) = (size(time_lock_firing{event,unit},2)-twin)/1000; %duration in seconds
                epoch_data.num_trials(event,unit) = sum(~isnan(epoch_data.firing_rates{event,unit}));
            end
        end
        
        %estimate the significance threshold has the 95-percentile of the
        %bootstrapped data
        for event = 1:size(time_lock_firing,1);
            for unit = 1:num_units
                if ~isempty(temporal_info.shuffled_info_rate{event,unit})
                    temporal_info.shuffled_95_percentile(event,unit) =...
                        prctile(temporal_info.shuffled_info_rate{event,unit},95);
                    temporal_info.shuffled_90_percentile(event,unit) =...
                        prctile(temporal_info.shuffled_info_rate{event,unit},90);
                end
            end
        end
        
        % store event descriptions, codes, durs, and t0 in epoch structure
        epoch_data.event_names = event_names;
        epoch_data.event_codes = event_codes;
        epoch_data.event_durs = event_durs;
        epoch_data.event_t0 = event_t0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Plot and Save Figures of Results---%%%
        unit_names.name = cfg.channel(1:num_units);
        unit_names.multiunit = multiunit;
        task_data.task_type = task;
        time_locked_plots(figure_dir,preprocessed_data_file,time_lock_firing,epoch_data,temporal_info,task_data,unit_names,smval)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Finally save all the data---%%%
        save([data_dir preprocessed_data_file(1:8) '-cvtnew-time_locked_results.mat'],...
            'time_lock_firing','smval','temporal_info','epoch_data');
        disp(['Time Locked Data Analyis for ' preprocessed_data_file(1:8) ' saved']);
        
    case 'Sequence'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import & reformat data so that spikes are locked to events---%%%
        
        load([data_dir preprocessed_data_file],'data','cfg','multiunit','num_units','item_set');
        num_trials = length(cfg.trl);
        
        disp('Aligning spike times to trial events')
        
        %get important task specific information
        [itmlist,sequence_items,sequence_locations] = read_Sequence_itm_and_cnd_files(item_set);
        [event_names,event_codes,event_durs,event_t0] = define_events(task,preprocessed_data_file(3:7));
        
        %preallocate space and parallel structure of cfg
        time_lock_firing = make_time_lock_cell(event_durs,num_units,num_trials,twin);
        which_sequence = NaN(1,length(cfg.trl));
        for t = 1:num_trials
            if any(cfg.trl(t).allval == event_codes(end)); %rewarded
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(1));
                event_times = get_Sequence_event_times(cfg.trl(t),event_codes,event_durs,trial_start);
                for unit = 1:num_units
                    time_lock_firing = put_in_time_matrix2(data,unit,t,event_times,time_lock_firing,twin,event_names);
                end
                which_sequence(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
            end
        end
        %remove excess NaNs associated with error trials
        time_lock_firing = laundry(time_lock_firing);
        which_sequence = laundry(which_sequence);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Perform Temporal analysis and signifance testing---%%%
        disp('Determining if spikes are locked to trial events')
        Fs = data(1).fsample; %should be 1000
        
        smval = 60; %temporal 1/2 width of gaussian smoothing filter
        
        % Collect spikes per epoch to determine if neurons fire more within
        % a given epoch ignore the 1st twin as this is from the previous event
        epoch_data.firing_rates = cell(size(time_lock_firing,1),num_units);
        epoch_data.dur = NaN(size(time_lock_firing,1),num_units);
        epoch_data.num_trials = NaN(size(time_lock_firing,1),num_units);
        
        %info per spike is equivalent to dividing information rate by average firing rate
        temporal_info.rate = NaN(size(time_lock_firing,1),num_units); %the observed information rate in bits/sec
        temporal_info.shuffled_info_rate = cell(size(time_lock_firing,1),num_units); %bootstrapped information rate in bits/sec expected by chance from spike train
        temporal_info.shuffled_95_percentile = NaN(size(time_lock_firing,1),num_units);
        temporal_info.shuffled_90_percentile = NaN(size(time_lock_firing,1),num_units);
        info_type = 'temporal';
        
        % perform mutual information theory analysis to determine
        % if neurons encode info about a particular time within an epoch
        for event = 1:size(time_lock_firing,1);
            for unit = 1:num_units
                [temporal_info.rate(event,unit),temporal_info.shuffled_info_rate{event,unit}]...
                    = estimated_mutual_information(time_lock_firing{event,unit},numshuffs,info_type,smval,Fs);
                epoch_data.firing_rates{event,unit} = nansum(time_lock_firing{event,unit}(:,twin:end),2);
                epoch_data.dur(event,unit) = (size(time_lock_firing{event,unit},2)-twin)/1000; %duration in seconds
                epoch_data.num_trials(event,unit) = sum(~isnan(epoch_data.firing_rates{event,unit}));
            end
        end
        
        %estimate the significance threshold has the 95-percentile of the
        %bootstrapped data
        for event = 1:size(time_lock_firing,1);
            for unit = 1:num_units
                if ~isempty(temporal_info.shuffled_info_rate{event,unit})
                    temporal_info.shuffled_95_percentile(event,unit) =...
                        prctile(temporal_info.shuffled_info_rate{event,unit},95);
                    temporal_info.shuffled_90_percentile(event,unit) =...
                        prctile(temporal_info.shuffled_info_rate{event,unit},90);
                end
            end
        end
        
        % store event descriptions, codes, durs, and t0 in epoch structure
        epoch_data.event_names = event_names;
        epoch_data.event_codes = event_codes;
        epoch_data.event_durs = event_durs;
        epoch_data.event_t0 = event_t0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Plot and Save Figures of Results---%%%
        unit_names.name = cfg.channel(1:num_units);
        unit_names.multiunit = multiunit;
        task_data.task_type = task;
        task_data.which_sequence = which_sequence;
        time_locked_plots(figure_dir,preprocessed_data_file,time_lock_firing,epoch_data,temporal_info,task_data,unit_names,smval)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Finally save all the data---%%%
        save([data_dir preprocessed_data_file(1:8) '-Sequence-time_locked_results.mat'],...
            'time_lock_firing','smval','temporal_info','epoch_data','which_sequence');
        disp(['Time Locked Data Analyis for ' preprocessed_data_file(1:8) ' saved']);
        
    case 'ListSQ'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import & reformat data so that spikes are locked to events---%%%
        
        load([data_dir preprocessed_data_file],'data','cfg','multiunit','item_set','num_units');
        num_trials = length(cfg.trl);
        
        disp('Aligning spike times to trial events')
        
        %get important task specific information
        [itmlist,sequence_items] = read_ListSQ_itm_and_cnd_files(item_set);
        [event_names,event_codes,event_durs,event_t0] = define_events(task,preprocessed_data_file(3:7));
        [which_img,novel_vs_repeat] = get_image_numbers(cfg,itmlist,sequence_items,event_codes(16));
        
        %preallocate space and parallel structure of cfg
        time_lock_firing = make_time_lock_cell(event_durs,num_units,num_trials,twin);
        trial_type = NaN(1,length(cfg.trl)); %sequences = 1, images = 2
        which_sequence = NaN(1,length(cfg.trl));
        broke_row = zeros(1,num_units); %since unknown number of break fixation errors :(
        for t = 1:num_trials
            if any(cfg.trl(t).allval == event_codes(2)); %in which image was displayed or 1st item in sequence was displayed
                if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(end) && isempty(find(cfg.trl(t).allval == event_codes(13)))% sequence trials and only want rewarded ones
                    continue % go to the next trial
                end
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(1));
                if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(end)
                    event_times = get_sequence_event_times(cfg.trl(t),event_codes,event_durs,trial_start);
                else
                    event_times = get_List_event_times(cfg.trl(t),event_codes,event_durs,trial_start);
                end
                for unit = 1:num_units
                    [time_lock_firing,broke_row] = put_in_time_matrix(data,unit,t,event_times,time_lock_firing,twin,broke_row,event_names);
                end
            end
            if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(end) %sequence trial
                trial_type(t) = 1;
                which_sequence(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
            else %image trial
                trial_type(t) =2;
            end
        end
        %remove excess NaNs associated with error trials
        time_lock_firing = laundry(time_lock_firing);
        which_sequence = laundry(which_sequence);
        trial_type = laundry(trial_type);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Perform Temporal analysis and signifance testing---%%%
        disp('Determining if spikes are locked to trial events')
        Fs = data(1).fsample; %should be 1000
        
        smval = 60; %temporal 1/2 width of gaussian smoothing filter
        smval_novrep = 250; %temporal 1/2 width of gaussian smoothing filter for image trials
        
        %each sequence own condition + novel + repeat images + 2 for
        %combined across images and combined across sequences
        max_conditions = length(sequence_items)+2+2;
        
        % Collect spikes per epoch to determine if neurons fire more within
        % a given epoch ignore the 1st twin as this is from the previous event
        epoch_data.firing_rates = cell(size(time_lock_firing,1),num_units,max_conditions);
        epoch_data.dur = NaN(size(time_lock_firing,1),num_units,max_conditions);
        epoch_data.num_trials = NaN(size(time_lock_firing,1),num_units,max_conditions);
        
        %info per spike is equivalent to dividing information rate by average firing rate
        temporal_info.rate = NaN(size(time_lock_firing,1),num_units,max_conditions); %the observed information rate in bits/sec
        temporal_info.shuffled_info_rate = cell(size(time_lock_firing,1),num_units,max_conditions); %bootstrapped information rate in bits/sec expected by chance from spike train
        temporal_info.shuffled_95_percentile = NaN(size(time_lock_firing,1),num_units,max_conditions);
        temporal_info.shuffled_90_percentile = NaN(size(time_lock_firing,1),num_units,max_conditions);
        info_type = 'temporal';
        
        % perform mutual information theory analysis to determine
        % if neurons encode info about a particular time within an epoch
        for event = 1:size(time_lock_firing,1);
            for unit = 1:num_units
                if ~isempty(time_lock_firing{event,unit})
                    if event == 1 || event == 18 || event == 19 %break fixations not really broken into novel vs repeat
                        condition = 1; %so i don't have to change the format of the data
                        temp_time_lock = time_lock_firing{event,unit};
                        [temporal_info.rate(event,unit,condition),temporal_info.shuffled_info_rate{event,unit,condition}]...
                            = estimated_mutual_information(temp_time_lock,numshuffs,info_type,smval_novrep,Fs);
                        
                        epoch_data.firing_rates{event,unit,condition} = nansum(temp_time_lock(:,twin:end),2);
                        epoch_data.dur(event,unit,condition) = (size(temp_time_lock,2)-twin)/1000; %duration in seconds
                        epoch_data.num_trials(event,unit,condition) = sum(~isnan(epoch_data.firing_rates{event,unit,condition}));
                    else
                        if event <= 13 %sequence trials
                            for condition = 1:length(sequence_items)+1;
                                if condition <= length(sequence_items) %sequence trials by sequence
                                    temp_time_lock = time_lock_firing{event,unit}(which_sequence  == condition,:);
                                elseif condition == length(sequence_items)+1 %all sequence trials combined
                                    temp_time_lock = time_lock_firing{event,unit};
                                end
                                [temporal_info.rate(event,unit,condition),temporal_info.shuffled_info_rate{event,unit,condition}]...
                                    = estimated_mutual_information(temp_time_lock,numshuffs,info_type,smval_novrep,Fs);
                                
                                epoch_data.firing_rates{event,unit,condition} = nansum(temp_time_lock(:,twin:end),2);
                                epoch_data.dur(event,unit,condition) = (size(temp_time_lock,2)-twin)/1000; %duration in seconds
                                epoch_data.num_trials(event,unit,condition) = sum(~isnan(epoch_data.firing_rates{event,unit,condition}));
                            end
                        else %image trials
                            for condition = length(sequence_items)+2:max_conditions
                                if condition == length(sequence_items)+2 %novel image presentations
                                    temp_time_lock = time_lock_firing{event,unit}(novel_vs_repeat == 1,:);
                                elseif condition == length(sequence_items)+3 %repeat image presentations
                                    temp_time_lock = time_lock_firing{event,unit}(novel_vs_repeat == 2,:);
                                else %combined image trials
                                    temp_time_lock = time_lock_firing{event,unit};
                                end
                                if event == 16 || event == 17 %for image trials since they have long events use different smoothing
                                    [temporal_info.rate(event,unit,condition),temporal_info.shuffled_info_rate{event,unit,condition}]...
                                        = estimated_mutual_information(temp_time_lock,numshuffs,info_type,smval_novrep,Fs);
                                else
                                    [temporal_info.rate(event,unit,condition),temporal_info.shuffled_info_rate{event,unit,condition}]...
                                        = estimated_mutual_information(temp_time_lock,numshuffs,info_type,smval_novrep,Fs);
                                end
                                epoch_data.firing_rates{event,unit,condition} = nansum(temp_time_lock(:,twin:end),2);
                                epoch_data.dur(event,unit,condition) = (size(temp_time_lock,2)-twin)/1000; %duration in seconds
                                epoch_data.num_trials(event,unit,condition) = sum(~isnan(epoch_data.firing_rates{event,unit,condition}));
                            end
                        end
                    end
                end
            end
        end
        
        %estimate the significance threshold has the 95-percentile of the
        %bootstrapped data
        for event = 1:size(time_lock_firing,1);
            for unit = 1:num_units
                for condition = 1:max_conditions
                    if ~isempty(temporal_info.shuffled_info_rate{event,unit,condition})
                        temporal_info.shuffled_95_percentile(event,unit,condition) =...
                            prctile(temporal_info.shuffled_info_rate{event,unit,condition},95);
                        temporal_info.shuffled_90_percentile(event,unit,condition) =...
                            prctile(temporal_info.shuffled_info_rate{event,unit,condition},90);
                    end
                end
            end
        end
        
        %above analysis shows temporal encoding in epochs so going to
        %analyze whether firing rate is different in defferent epochs
        %across sequences and across novel vs repeat images
        epoch_data.pvals = NaN(size(time_lock_firing,1),num_units);
        for unit = 1:num_units
            for event = 2:size(time_lock_firing,1) %event 1 is ITI
                if event <= 13 %for sequence trials
                    if ~isempty(epoch_data.firing_rates{event,unit,1})
                        [~,epoch_data.pvals(event,unit)] = ...
                            kstest2(epoch_data.firing_rates{event,unit,1},epoch_data.firing_rates{event,unit,2});
                    end
                else
                    if ~isempty(epoch_data.firing_rates{event,unit,end-2}) && ~isempty(epoch_data.firing_rates{event,unit,end-1})
                        [~,epoch_data.pvals(event,unit)] = ...
                            kstest2(epoch_data.firing_rates{event,unit,end-2},epoch_data.firing_rates{event,unit,end-1});
                    end
                end
            end
        end
        
        
        % store event descriptions, codes, durs, and t0 in epoch structure
        epoch_data.event_names = event_names;
        epoch_data.event_codes = event_codes;
        epoch_data.event_durs = event_durs;
        epoch_data.event_t0 = event_t0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Plot and Save Figures of Results---%%%
        unit_names.name = cfg.channel(1:num_units);
        unit_names.multiunit = multiunit;
        task_data.task_type = 'ListSQ_Sequence';
        task_data.which_sequence = which_sequence;
        task_data.novel_vs_repeat = novel_vs_repeat;
        time_locked_plots(figure_dir,preprocessed_data_file,time_lock_firing,epoch_data,temporal_info,task_data,unit_names,smval)
        task_data.task_type = 'ListSQ_List';
        
        time_locked_plots(figure_dir,preprocessed_data_file,time_lock_firing,epoch_data,temporal_info,task_data,unit_names,[smval,smval_novrep])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Finally save all the data---%%%
        save([data_dir preprocessed_data_file(1:8) '-ListSQ-time_locked_results.mat'],'time_lock_firing',...
            'which_sequence','trial_type','smval','smval_novrep','temporal_info','epoch_data');
        disp(['Time Locked Data Analyis for ' preprocessed_data_file(1:8) ' saved']);
end
end

function [event_names,event_codes,event_durs,event_t0] = define_events(task,recording_date)
% define events, their names, and their durations.
switch task
    case 'cvtnew'
        % epoch 1 ITI period
        % epoch 2 fixation spot turns on
        % epoch 3 fixation on fixation spot occurs
        % epoch 4 dot turns on
        % epoch 5 dot changes color
        % epoch 6 responds to color change
        % epoch 7 reward period
        
        % There is one more event code than there are blocks in
        % time_locked_firing since ITI period is defined by 2 events (15 & 16)
        event_names = {{'ITIstart','ITIend'},'fixspot on','fixation occurs',...
            'Dot on','Dot Colorchange','Response','reward'};
        
        event_codes = [15 16 35 8 25 27 4 3];
        event_durs = [1000 zeros(1,2) 2000 zeros(1,2) 1000]; %event durations
        event_t0 = 500*ones(1,7); %in the matrix what index is equivalent to the time the event occured
    case 'Sequence'
        % epoch 1 ITI period
        % epoch 2 item 1 in sequence is displayed
        % epoch 3 item 1 in sequence is fixated
        % epoch 4 item 1 turned off/III (inter item interval) 1
        % epoch 5 item 2 in sequence is displayed
        % epoch 6 item 2 in sequence is fixated
        % epoch 7 item 2 turned off/III (inter item interval) 2
        % epoch 8 item 3 in sequence is displayed
        % epoch 9 item 3 in sequence is fixated
        % epoch 10 item 3 turned off/III (inter item interval) 3
        % epoch 11 item 4 in sequence is displayed
        % epoch 12 item 4 in sequence is fixated
        % epoch 13 reward period during sequence trials
        event_names = {{'ITIstart','ITIend'},'item1 on','fix item1','III 1',....
            'item2 on','fix item2','III 2','item3 on','fix item3',...
            'III 3','item4 on','fix item4','reward'};
        event_codes = [15 23 8 24 25 8 26 27 8 28 29 8 3 ];
        event_durs = [1000 zeros(1,11) 1000]; %event durations
        event_t0 = [500*ones(1,13)]; %in the matrix what index is equivalent to the time the event occured
        
    case 'ListSQ'
        % epoch 1 ITI period
        % epoch 2 item 1 in sequence is displayed
        % epoch 3 item 1 in sequence is fixated
        % epoch 4 item 1 turned off/III (inter item interval) 1
        % epoch 5 item 2 in sequence is displayed
        % epoch 6 item 2 in sequence is fixated
        % epoch 7 item 2 turned off/III (inter item interval) 2
        % epoch 8 item 3 in sequence is displayed
        % epoch 9 item 3 in sequence is fixated
        % epoch 10 item 3 turned off/III (inter item interval) 3
        % epoch 11 item 4 in sequence is displayed
        % epoch 12 item 4 in sequence is fixated
        % epoch 13 reward period during sequence trials
        % epoch 14 fixation cross is displayed
        % epoch 15 fixation on crosshair to initialize image trials
        % epoch 16 image viewing period start
        % epoch 17 image viewing period end
        % epoch 18 broke fixation while viewing image
        % epoch 19 started viewing image again after breaking fixation
        
        % img_on and img_off are equivalent to test0 on and test0 off i.e.
        % item 1 on and item 1 off. inter_item_interval[1-3] is equivalent to
        % Item[1-3]_off. Item4_off should roughly cooccur with reward.
        % There is one more event code than there are blocks in
        % time_locked_firing since ITI period is defined by 2 events (15 & 16)
        event_names = {{'ITIstart','ITIend'},'item1 on','fix item1','III 1',....
            'item2 on','fix item2','III 2','item3 on','fix item3',...
            'III 3','item4 on','fix item4','reward','fixspot on',...
            'fixation','img on','img off','break_fixation','refixated'};
        event_codes = [15 23 8 24 25 8 26 27 8 28 29 8 3 35 8 23 24 209 8];
        
        %image duration changed part way through from 7 secs to 5 secs of
        %cumulative looking times
        if recording_date < 140805
            imgdur = 7000;
        else
            imgdur = 5000;
        end
        event_durs = [1000 zeros(1,11) 1000 zeros(1,2) imgdur imgdur zeros(1,2)]; %event durations
        event_t0 = [500*ones(1,16) imgdur+500 500*ones(1,2)]; %in the matrix what index is equivalent to the time the event occured
end
end

function time_lock_firing = make_time_lock_cell(event_durs,num_units,num_trials,twin)
% preallocate space for time locked data
time_lock_firing = cell(length(event_durs),num_units); % unit by epoch
for i = 1:size(time_lock_firing,1);%by event
    for ii = 1:size(time_lock_firing,2) %by unit
        time_lock_firing{i,ii} = NaN(num_trials,event_durs(i)+2*twin);
    end
end
end

function event_times = get_sequence_event_times(cfgtrl,event_codes,event_durs,trial_start)
% function determines event times in a given trial
event_times = cell(1,size(event_codes-1,2));
for events = 1:13
    if events == 1 %this is the ITI period
        event_times{events} = [cfgtrl.alltim(cfgtrl.allval == event_codes(events))...
            cfgtrl.alltim(cfgtrl.allval == event_codes(events)+1)]-trial_start;
        event_times{events}(event_times{events} > event_durs(events)) = event_durs(events);
    elseif events == 3 || events == 6 || events == 9 || events == 12 %% fixations on items in sequence
        temp_times = cfgtrl.alltim(cfgtrl.allval == event_codes(events))-trial_start;
        previous_event_time  = event_times{events-1};
        if isempty(previous_event_time) %mean no event occured, this is an image trial
            continue
        end
        temp_times(temp_times <  previous_event_time ) = [];
        temp_times = temp_times(1);
        event_times{events} = temp_times;
    elseif events == 13
        temp_times = cfgtrl.alltim(cfgtrl.allval == event_codes(events))-trial_start;
        event_times{events} = [temp_times(1) temp_times(end)];
        break % all other events are for image trials
    else
        event_times{events} = cfgtrl.alltim(cfgtrl.allval == event_codes(events))-trial_start;
    end
end
end

function event_times = get_List_event_times(cfgtrl,event_codes,event_durs,trial_start)
% function determines event times in a given trial
event_times = cell(1,size(event_codes,2));
for events = 14:size(event_codes-1,2)
    if events == 15 %fixation on crosshair
        temp_times = cfgtrl.alltim(cfgtrl.allval == event_codes(events))-trial_start;
        event_times{events} = temp_times(1);
    elseif events == 16 %locked to image on
        temp_times = cfgtrl.alltim(cfgtrl.allval == event_codes(events))-trial_start;
        event_times{events} = [temp_times temp_times + event_durs(events)];
    elseif events == 17 %locked to image off
        temp_times = cfgtrl.alltim(cfgtrl.allval == event_codes(events))-trial_start;
        if ~isempty(temp_times)
            event_times{events} = [temp_times - event_durs(events) temp_times];
        end
    elseif events == 18 || events == 19
        temp_times = cfgtrl.alltim(cfgtrl.allval == event_codes(events))-trial_start;
        temp_times(temp_times > event_durs(16)*1.5+event_times{16}(1)) = []; %don't really care if they're really not paying attention to the image
        temp_times(temp_times < event_times{16}(1)+50) = []; %don't include breaks or fixations before image is turned on
        event_times{events} = temp_times;
    else
        event_times{events} = cfgtrl.alltim(cfgtrl.allval == event_codes(events))-trial_start;
    end
end

if ~isempty(event_times{18}) %so there are some break fixations and refixations
    if length(event_times{18})-1 == length(event_times{19}) %occurs if the refixation is late in the trial
        event_times{18}(end) = [];
    elseif length(event_times{19})-1 == length(event_times{18})
        event_times{19}(end) = [];
    elseif length(event_times{18})-1 > length(event_times{19})
        disp('error break fixations not matching refixations')
    end
    duration_outside = event_times{19}-event_times{18};
    %only take trials in which the monkey looked outside for long
    %durations setting to 500 since this is clearly more than 1 fixation
    event_times{18}(duration_outside < 500) = [];
    event_times{19}(duration_outside < 500) = [];
end
end

function event_times = get_cvtnew_event_times(cfgtrl,event_codes,event_durs,trial_start,twin)
% function determines event times in a given trial
event_times = cell(1,size(event_codes,2)-1);
for events = 1:length(event_times)
    if events == 1 %this is the ITI period
        event_times{events} = [cfgtrl.alltim(cfgtrl.allval == event_codes(events))...
            cfgtrl.alltim(cfgtrl.allval == event_codes(events)+1)]-trial_start;
        event_times{events}(event_times{events} > event_durs(events)) = event_durs(events);
    elseif events == length(event_times)%reward period
        temp_times = cfgtrl.alltim(cfgtrl.allval == event_codes(events+1))-trial_start;
        event_times{events} = [temp_times(1) temp_times(end)];
    elseif event_durs(events)== 2000% event 4 dot on for cvtnew_task
        temp_times = cfgtrl.alltim(cfgtrl.allval == event_codes(events+1))-trial_start;
        temp_times(2) = cfgtrl.alltim(cfgtrl.allval == event_codes(events+2))-trial_start-twin;%when color change occurs
        event_times{events} = temp_times;
    else
        event_times{events} = cfgtrl.alltim(cfgtrl.allval == event_codes(events+1))-trial_start;
    end
end
end

function event_times = get_Sequence_event_times(cfgtrl,event_codes,event_durs,trial_start)
% function determines event times in a given trial
event_times = cell(1,size(event_codes,2));
for events = 1:length(event_times)
    if events == 1 %this is the ITI period
        event_times{events} = [cfgtrl.alltim(cfgtrl.allval == event_codes(events))...
            cfgtrl.alltim(cfgtrl.allval == event_codes(events)+1)]-trial_start;
        event_times{events}(event_times{events} > event_durs(events)) = event_durs(events);
    elseif events == length(event_times)%reward period
        temp_times = cfgtrl.alltim(cfgtrl.allval == event_codes(events))-trial_start;
        event_times{events} = [temp_times(1) temp_times(end)];
    else
        temp_times = cfgtrl.alltim(cfgtrl.allval == event_codes(events))-trial_start;
        event_times{events} = temp_times(1); %in case there are multiple events like refixating a item
    end
end
end


function [time_lock_firing,broke_row] = put_in_time_matrix(data,unit,t,event_times,time_lock_firing,twin,broke_row,event_names)
% function finds spikes that are time locked to events
spikes = find(data(unit).values{t});
for event = 1:length(event_times)-2;
    if ~isempty(event_times{event})
        if length(event_times{event}) > 1 %so events that have a time period i.e. img display, ITI, and rewad
            if event == 1
                if t == 1 %only can look forward in time for spikes
                    event_spikes = spikes(spikes <= event_times{event}(2)+twin)+twin;
                    timevec = [1 event_times{event}(2)+twin]+twin;
                    tempvec = [NaN(1,twin) zeros(1,timevec(2)-timevec(1)+1)];
                    tempvec(event_spikes) = 1;
                    tempvec(1:twin) = [];
                else %look back in time to end of previous trial
                    prespikes = find(data(unit).values{t-1});
                    prespikes = prespikes(prespikes > length(data(unit).values{t-1})-twin) - length(data(unit).values{t-1});
                    event_spikes = [prespikes spikes(spikes <= event_times{event}(2)+twin)]+twin;
                    timevec = [1-twin event_times{event}(2)+twin]+twin;
                    tempvec =  zeros(1,timevec(2)-timevec(1)+1);
                    tempvec(event_spikes) = 1;
                end
                time_lock_firing{event,unit}(t,timevec(1):timevec(2)) = tempvec;
            elseif strcmpi(event_names{event},'reward') %whihc event is reward changes by task
                if t == length(data(unit).values);
                    if length(data(unit).values{t}) > event_times{event}(2)+twin;
                        rewardstop = event_times{event}(2)+twin;
                    else
                        rewardstop = length(data(unit).values{t});
                    end
                    event_spikes = spikes(spikes > event_times{event}(1)-twin & spikes <= rewardstop)-event_times{event}(1)+twin;
                    timevec = [1-twin rewardstop-event_times{event}(1)]+twin;
                else %look to the next trial
                    if length(data(unit).values{t}) > event_times{event}(2)+twin;
                        rewardstop = event_times{event}(2)+twin;
                        event_spikes = spikes(spikes > event_times{event}(1)-twin & spikes <= rewardstop)-event_times{event}(1)+twin;
                        timevec = [1-twin rewardstop-event_times{event}(1)]+twin;
                    else
                        rewardstop = twin-(length(data(unit).values{t})-event_times{event}(2));
                        postspikes = find(data(unit).values{t+1});
                        postspikes = postspikes(postspikes <= rewardstop);
                        event_spikes = [spikes(spikes > event_times{event}(1)-twin)-event_times{event}(1)...
                            postspikes+(length(data(unit).values{t})-event_times{event}(2))+(event_times{event}(2)-event_times{event}(1))]+twin;
                        timevec = [1-twin event_times{event}(2)-event_times{event}(1)+twin]+twin;
                    end
                end
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                time_lock_firing{event,unit}(t,timevec(1):timevec(2)) = tempvec;
            else %image trials
                event_spikes = spikes(spikes > event_times{event}(1)-twin & spikes <= event_times{event}(2)+twin)-event_times{event}(1)+twin;
                timevec = [event_times{event}(1)-twin+1 event_times{event}(2)+twin]-event_times{event}(1)+twin;
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                time_lock_firing{event,unit}(t,timevec(1):timevec(2)) = tempvec;
            end
        else %events that occur at a single instance
            event_spikes = spikes(spikes > event_times{event}-twin & spikes <= event_times{event}+twin)-event_times{event}+twin;
            timevec = [event_times{event}-twin+1 event_times{event}+twin]-event_times{event}+twin;
            tempvec = zeros(1,timevec(2)-timevec(1)+1);
            tempvec(event_spikes) = 1;
            time_lock_firing{event,unit}(t,timevec(1):timevec(2)) = tempvec;
        end
    end
end
if ~isempty(event_times{18}) %if there are break fixations during the image trials
    for b = 1:length(event_times{18})
        broke_row(unit) = broke_row(unit)+1;
        
        event_spikes = spikes(spikes > event_times{18}(b)-twin & spikes <= event_times{18}(b)+twin)-event_times{18}(b)+twin;
        timevec = [1-twin twin]+twin;
        tempvec = zeros(1,timevec(2)-timevec(1)+1);
        tempvec(event_spikes) = 1;
        time_lock_firing{18,unit}(broke_row(unit),timevec(1):timevec(2)) = tempvec;
        
        event_spikes = spikes(spikes > event_times{19}(b)-twin & spikes <= event_times{19}(b)+twin)- event_times{19}(b)+twin;
        timevec = [1-twin twin]+twin;
        tempvec = zeros(1,timevec(2)-timevec(1)+1);
        tempvec(event_spikes) = 1;
        time_lock_firing{19,unit}(broke_row(unit),timevec(1):timevec(2)) = tempvec;
    end
end
end

function time_lock_firing = put_in_time_matrix2(data,unit,t,event_times,time_lock_firing,twin,event_names)
spikes = find(data(unit).values{t});
for event = 1:length(event_times);
    if ~isempty(event_times{event})
        if length(event_times{event}) > 1 %so events that have a time period i.e. img display, ITI, and rewad
            if event == 1 %ITI period
                if t == 1 %only can look forward in time for spikes since 1st ITI
                    event_spikes = spikes(spikes <= event_times{event}(2)+twin)+twin;
                    timevec = [1 event_times{event}(2)+twin]+twin;
                    tempvec = [NaN(1,twin) zeros(1,timevec(2)-timevec(1)+1)];
                    tempvec(event_spikes) = 1;
                    tempvec(1:twin) = [];
                else %look back in time to end of previous trial
                    prespikes = find(data(unit).values{t-1});
                    prespikes = prespikes(prespikes > length(data(unit).values{t-1})-twin) - length(data(unit).values{t-1});
                    event_spikes = [prespikes spikes(spikes <= event_times{event}(2)+twin)]+twin;
                    timevec = [1-twin event_times{event}(2)+twin]+twin;
                    tempvec =  zeros(1,timevec(2)-timevec(1)+1);
                    tempvec(event_spikes) = 1;
                end
                time_lock_firing{event,unit}(t,timevec(1):timevec(2)) = tempvec;
            elseif strcmpi(event_names{event},'reward') %which event is reward changes by task
                if t == length(data(unit).values);
                    if length(data(unit).values{t}) > event_times{event}(2)+twin;
                        rewardstop = event_times{event}(2)+twin;
                    else
                        rewardstop = length(data(unit).values{t});
                    end
                    event_spikes = spikes(spikes > event_times{event}(1)-twin & spikes <= rewardstop)-event_times{event}(1)+twin;
                    timevec = [1-twin rewardstop-event_times{event}(1)]+twin;
                else %look to the next trial
                    if length(data(unit).values{t}) > event_times{event}(2)+twin;
                        rewardstop = event_times{event}(2)+twin;
                        event_spikes = spikes(spikes > event_times{event}(1)-twin & spikes <= rewardstop)-event_times{event}(1)+twin;
                        timevec = [1-twin rewardstop-event_times{event}(1)]+twin;
                    else
                        rewardstop = twin-(length(data(unit).values{t})-event_times{event}(2));
                        postspikes = find(data(unit).values{t+1});
                        postspikes = postspikes(postspikes <= rewardstop);
                        event_spikes = [spikes(spikes > event_times{event}(1)-twin)-event_times{event}(1)...
                            postspikes+(length(data(unit).values{t})-event_times{event}(2))+(event_times{event}(2)-event_times{event}(1))]+twin;
                        timevec = [1-twin event_times{event}(2)-event_times{event}(1)+twin]+twin;
                    end
                end
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                time_lock_firing{event,unit}(t,timevec(1):timevec(2)) = tempvec;
            else %dot on
                event_spikes = spikes(spikes > event_times{event}(1)-twin & spikes <= event_times{event}(2)+twin)-event_times{event}(1)+twin;
                timevec = [event_times{event}(1)-twin+1 event_times{event}(2)+twin]-event_times{event}(1)+twin;
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                time_lock_firing{event,unit}(t,timevec(1):timevec(2)) = tempvec;
            end
        else %events that occur at a single instance
            event_spikes = spikes(spikes > event_times{event}-twin & spikes <= event_times{event}+twin)-event_times{event}+twin;
            timevec = [event_times{event}-twin+1 event_times{event}+twin]-event_times{event}+twin;
            tempvec = zeros(1,timevec(2)-timevec(1)+1);
            tempvec(event_spikes) = 1;
            time_lock_firing{event,unit}(t,timevec(1):timevec(2)) = tempvec;
        end
    end
end
end
