function Sequence_Saccade_Analysis(data_dir,preprocessed_data_file,figure_dir,task,predicted_rt)
% written by Seth Konig September, 2014
% updated by Seth Konig April 7, 2015 to include estimating temporal
% information locked to saccades and fixations on each item.
% Significantly modified original version (Find in OLd and Misc Code) to
% include newer main function for anlalyzing all sequence trial data. SDK
% 3/3/15. SDK Removed LFP from analysis 1/6/2016
%
% Function analyses spike times correlated with eye movements in the
% sequence task.
%
% Inputs:
%   1) data_dir: direction where preprocessed data is located
%   2) preprocessed_data_file: preprocessed ListSQ data file with LFPs
%   divided by channel and trial.
%   3) figure_dir: location of where to put save figures
%   4) task: what task was this data come from with i.e. 'ListSQ' or 'Sequence'
%   5) predicted_rt: maximum "reaction time" for what constitutes as predictive, everything else is reactive

%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-sequence_saccade_locked_results'


twin = 500;% how much time to take before and after saccade.
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
fixwin = 5;%size of fixation window on each crosshair
event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
trial_start_code = 15;
smval = 120;%gaussian 1/2 width for smoothing
predicted_thresh = 10;% percent of saccades that must be < 150 ms to item to constitue as a learned-predicted item
% for Vivian 5% was the false positive rate for 150 ms so wanted at least double the False positive rate
numshuffs = 10; %recommend this is between 100 & 1000, for bootstrapping to
% dtermine significance
info_type = 'temporal';
Fs = 1000;

switch task
    case 'ListSQ'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import data & get successful trials---%%%
        
        load([data_dir preprocessed_data_file]);
        
        %get important task specific information
        [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
        overlap = find((sequence_locations{1}(1,:) ==  sequence_locations{2}(1,:)) & ...
            (sequence_locations{1}(2,:) ==  sequence_locations{2}(2,:)));
        
        %preallocate space and parallel structure of cfg
        successful_sequence_trials = NaN(1,length(cfg.trl));
        which_sequence = NaN(1,length(cfg.trl));
        for t = 1:length(cfg.trl);
            if sum(cfg.trl(t).allval == 3) >= 6; %in which sequence trials were rewarded
                which_sequence(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
                successful_sequence_trials(t) = t;
            end
        end
        successful_sequence_trials = laundry(successful_sequence_trials);
        which_sequence = laundry(which_sequence);
        
        num_trials = length(successful_sequence_trials);
        event_times = NaN(length(successful_sequence_trials),8);
        for t = 1:num_trials
            trial_start = cfg.trl(successful_sequence_trials(t)).alltim(cfg.trl(successful_sequence_trials(t)).allval == trial_start_code);
            for event = 1:length(event_codes);
                event_times(t,event) = cfg.trl(successful_sequence_trials(t)).alltim(cfg.trl(successful_sequence_trials(t)).allval == event_codes(event))-trial_start;
            end
        end
        
    case 'Sequence'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import data & get successful trials---%%%
        
        load([data_dir preprocessed_data_file]);
        
        %get important task specific information
        [itmlist,sequence_items,sequence_locations] = read_Sequence_itm_and_cnd_files(item_set);
        overlap = find((sequence_locations{1}(1,:) ==  sequence_locations{2}(1,:)) & ...
            (sequence_locations{1}(2,:) ==  sequence_locations{2}(2,:)));
        
        %preallocate space and parallel structure of cfg
        successful_sequence_trials = NaN(1,length(cfg.trl));
        which_sequence = NaN(1,length(cfg.trl));
        for t = 1:length(cfg.trl);
            if sum(cfg.trl(t).allval == 3) >= 6; %in which sequence trials were rewarded
                which_sequence(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
                successful_sequence_trials(t) = t;
            end
        end
        successful_sequence_trials = laundry(successful_sequence_trials);
        which_sequence = laundry(which_sequence);
        
        num_trials = length(successful_sequence_trials);
        event_times = NaN(length(successful_sequence_trials),8);
        for t = 1:num_trials
            trial_start = cfg.trl(successful_sequence_trials(t)).alltim(cfg.trl(successful_sequence_trials(t)).allval == trial_start_code);
            for event = 1:length(event_codes);
                event_times(t,event) = cfg.trl(successful_sequence_trials(t)).alltim(cfg.trl(successful_sequence_trials(t)).allval == event_codes(event))-trial_start;
            end
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---process eye data locked to trial events---%%%
fixationstats = fixationstats(successful_sequence_trials);
cfg.trl = cfg.trl(successful_sequence_trials);

saccade_start_time = NaN(length(fixationstats),4);%when did saccade to item start
fixation_start_time = NaN(length(fixationstats),4);%when did fixation on item start
reaction_time = NaN(length(fixationstats),4); %how long after a stimulus disappears do they saccade
time_to_fixation = NaN(length(fixationstats),4); %how fast did they get to it the first time around
fixation_accuracy = NaN(length(fixationstats),4); %how far off
fixation_duration = NaN(length(fixationstats),4); %fixation duration
extrafixations = NaN(length(fixationstats),4); %how many noncross hair fixations do they mak
for trial = 1:num_trials
    locs = sequence_locations{which_sequence(trial)};
    
    %convert to DVA for this analysis
    locs(1,:) = (locs(1,:)-400)/24;
    locs(2,:) = (locs(2,:)-300)/24;
    fixationstats{trial}.fixations(1,:) = (fixationstats{trial}.fixations(1,:)-400)/24;
    fixationstats{trial}.fixations(2,:) = (fixationstats{trial}.fixations(2,:)-300)/24;
    
    event_codes = cfg.trl(trial).allval;
    event_codes(event_codes == 100)= 0;
    event_codes(1) = 100;%eye data starts for recording right away
    event_times = cfg.trl(trial).alltim;
    
    trialdata = analyze_sequence_trial(fixationstats{trial},locs,fixwin,...
        event_codes,event_times);
    
    time_to_fixation(trial,:) = trialdata.t2f;
    fixation_duration(trial,:) = trialdata.fixation_duration;
    reaction_time(trial,:)  = trialdata.time_to_leave;
    extrafixations(trial,:) = trialdata.extrafixations;
    fixation_accuracy(trial,:) =  trialdata.accuracy;
    
    fixation_numbers = trialdata.fixationnums; %fixation number for each item
    fixationtimes = fixationstats{trial}.fixationtimes;
    saccadetimes = fixationstats{trial}.saccadetimes;
    if saccadetimes(1,1) < fixationtimes(1,1); %then started with a saccade
        sacind = 0;
    else%started with a fixation
        sacind = -1;
    end
    for item = 1:4
        if ~isnan(fixation_numbers(item))
            fixation_start_time(trial,item) = fixationtimes(1,fixation_numbers(item));
            if fixation_numbers(item)+sacind  >= 1
                try
                    saccade_start_time(trial,item) = saccadetimes(1,fixation_numbers(item)+sacind);
                catch
                    saccade_start_time(trial,item) = NaN;%should only occur if saccade was off screen or a blink occured otherwise should work
                    disp('Saccade Start Time not found')
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Firing Rate Locked to Eye Data---%%%
saccade_locked_firing = cell(4,num_units);
fixation_locked_firing = cell(4,num_units);
for c = 1:4
    for unit = 1:num_units
        saccade_locked_firing{c,unit} = NaN(num_trials,twin*2);
        fixation_locked_firing{c,unit} = NaN(num_trials,twin*2);
    end
end
for trial = 1:num_trials
    for unit = 1:num_units;
        spikes = find(data(unit).values{successful_sequence_trials(trial)});
        for c = 1:4;
            fixt = fixation_start_time(trial,c);
            sact = saccade_start_time(trial,c);
            if ~isnan(fixt)
                fix_spikes = spikes(spikes > fixt-twin & spikes <= fixt+twin)-fixt+twin;
                temp = zeros(1,twin*2);
                temp(fix_spikes) = 1;
                fixation_locked_firing{c,unit}(trial,:) = temp;
            end
            if ~isnan(sact)
                sac_spikes = spikes(spikes > sact-twin & spikes <= sact+twin)-sact+twin;
                temp = zeros(1,twin*2);
                temp(sac_spikes) = 1;
                saccade_locked_firing{c,unit}(trial,:) = temp;
            end
        end
    end
end

fixation_info = NaN(4,num_units);
fixation_info_90 = NaN(4,num_units);
fixation_info_95 = NaN(4,num_units);
fixation_shuffled_info = cell(4,num_units);
saccade_info = NaN(4,num_units);
saccade_info_90 = NaN(4,num_units);
saccade_info_95 = NaN(4,num_units);
saccade_shuffled_info = cell(4,num_units);

for c = 1:size(saccade_locked_firing,1);
    for unit = 1:num_units
        [fixation_info(c,unit),fixation_shuffled_info{c,unit}]...
            = estimated_mutual_information(fixation_locked_firing{c,unit},numshuffs,info_type,smval,Fs);
        [saccade_info(c,unit),saccade_shuffled_info{c,unit}]...
            = estimated_mutual_information(saccade_locked_firing{c,unit},numshuffs,info_type,smval,Fs);
        fixation_info_90(c,unit) = prctile(fixation_shuffled_info{unit},90);
        fixation_info_95(c,unit) = prctile(fixation_shuffled_info{unit},95);
        saccade_info_90(c,unit) = prctile(saccade_shuffled_info{unit},90);
        saccade_info_95(c,unit) = prctile(saccade_shuffled_info{unit},95);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Analyzing Predictive vs Reactive Spike Data---%%%
%determine which items in which sequence had a consistent number of
%predicted saccades greater than the minimum threshold
predicted_items = zeros(2,4);
percent_predicted = NaN(2,4);
for seq = 1:2
    percent_predicted(seq,:) = sum(time_to_fixation(which_sequence == seq,:) < predicted_rt)...
        ./sum(~isnan(time_to_fixation(which_sequence == seq,:)))*100;
    predicted = find(percent_predicted(seq,:) > predicted_thresh);
    predicted_items(seq,predicted) = 1;
end
percent_predicted = round(percent_predicted);

% determine  which trials and items were predicted
predicted = zeros(num_trials,4);
for seq = 1:2;
    for item = 1:4
        if predicted_items(seq,item)
            predicted(which_sequence == seq ,item) = 1;
            predicted(time_to_fixation(:,item) > predicted_rt,item) = 0;
        end
    end
end
predicted(1,:) = 0; predicted(11,:) = 0; %dont ever want 1st trial to accidentally skew results

t = -twin:twin-1;
unit_names = cfg.channel(1:num_units);
for unit = 1:num_units
    figure
    ylims = NaN(1,8);
    for c = 1:4
        subplot(2,2,c)
        hold on
        if ~all(all(isnan(saccade_locked_firing{c,unit}(which_sequence == 1 & ~predicted(:,c)',:))))
            dofill(t,saccade_locked_firing{c,unit}(which_sequence == 1 & ~predicted(:,c)',:),'red',1,smval); %reactive sequence 1
        end
        if ~all(all(isnan(saccade_locked_firing{c,unit}(which_sequence == 2 & ~predicted(:,c)',:))))
            dofill(t,saccade_locked_firing{c,unit}(which_sequence == 2 & ~predicted(:,c)',:),'blue',1,smval);%reactive sequence 2
        end
        if predicted_items(1,c);
            if ~all(all(isnan(saccade_locked_firing{c,unit}(which_sequence == 1 & predicted(:,c)'))))
                dofill(t,saccade_locked_firing{c,unit}(which_sequence == 1 & predicted(:,c)',:),'magenta',1,smval); %predictive sequence 1
            end
        end
        if predicted_items(2,c)
            if ~all(all(isnan(saccade_locked_firing{c,unit}(which_sequence == 2 & predicted(:,c)',:))))
                dofill(t,saccade_locked_firing{c,unit}(which_sequence == 2 & predicted(:,c)',:),'cyan',1,smval);%predictive sequence 2
            end
        end
        hold off
        xlabel('Time from Eye Movement (ms)')
        ylabel('Firing Rate (Hz)')
        yl = ylim;
        ylims(c) = yl(2);
        
        bitstr = [];
        if fixation_info(c,unit) > fixation_info_95(c,unit)
            bitstr = ['fix_{95} = ' num2str(fixation_info(c,unit))];
        elseif fixation_info(c,unit) > fixation_info_90(c,unit)
            bitstr = ['fix_{90} = ' num2str(fixation_info(c,unit))];
        end
        if saccade_info(c,unit) > saccade_info_95(c,unit)
            bitstr = [bitstr 'sac_{95} = ' num2str(saccade_info(c,unit))];
        elseif saccade_info(c,unit) > saccade_info_90(c,unit)
            bitstr = [bitstr 'sac_{90} = ' num2str(saccade_info(c,unit))];
        end
        
        if overlap == c
            bitstr = [bitstr ' (Overlapping)'];
        end
        
        if sum(predicted_items(:,c)) == 0 %no predictived saccades
            title(['Item # ' num2str(c) ' ' bitstr])
        elseif sum(predicted_items(:,c)) == 2 % both sequences predicted
            title(['Item # ' num2str(c) ' Seq #1 ' num2str(percent_predicted(1,c)) '% Seq #2 ' num2str(num2str(percent_predicted(2,c))) '%' ' ' bitstr])
        elseif predicted_items(1,c)
            title(['Item # ' num2str(c) ' Seq #1 ' num2str(percent_predicted(1,c)) '%' ' ' bitstr])
        elseif predicted_items(2,c)
            title(['Item # ' num2str(c) ' Seq #2 ' num2str(num2str(percent_predicted(2,c))) '%' ' ' bitstr])
        end
    end
    for c = 1:4;
        subplot(2,2,c)
        ylim([0 max(ylims)]);
    end
    if multiunit(unit)
        subtitle(['Saccade-Locked Multiunit ' unit_names{unit} ' n =' num2str(num_trials)]);
    else
        subtitle(['Saccade-Locked ' unit_names{unit} ' n =' num2str(num_trials)]);
    end
    save_and_close_fig(figure_dir,[preprocessed_data_file(1:end-11) '_' unit_names{unit} '_Sequence_Saccade_Locked_analysis']);
end

t = -twin:twin-1;
unit_names = cfg.channel(1:num_units);
for unit = 1:num_units
    figure
    ylims = NaN(1,8);
    for c = 1:4
        subplot(2,2,c)
        hold on
        if ~all(all(isnan(fixation_locked_firing{c,unit}(which_sequence == 1 & ~predicted(:,c)',:))))
            dofill(t,fixation_locked_firing{c,unit}(which_sequence == 1 & ~predicted(:,c)',:),'red',1,smval); %reactive sequence 1
        end
        if ~all(all(isnan(fixation_locked_firing{c,unit}(which_sequence == 2 & ~predicted(:,c)',:))))
            dofill(t,fixation_locked_firing{c,unit}(which_sequence == 2 & ~predicted(:,c)',:),'blue',1,smval);%reactive sequence 2
        end
        if predicted_items(1,c)
            if ~all(all(isnan(fixation_locked_firing{c,unit}(which_sequence == 1 & predicted(:,c)',:))))
                dofill(t,fixation_locked_firing{c,unit}(which_sequence == 1 & predicted(:,c)',:),'magenta',1,smval); %predictive sequence 1
            end
        end
        if predicted_items(2,c)
            if ~all(all(isnan(fixation_locked_firing{c,unit}(which_sequence == 2 & predicted(:,c)',:))))
                dofill(t,fixation_locked_firing{c,unit}(which_sequence == 2 & predicted(:,c)',:),'cyan',1,smval);%predictive sequence 2
            end
        end
        hold off
        xlabel('Time from Eye Movement (ms)')
        ylabel('Firing Rate (Hz)')
        yl = ylim;
        ylims(c) = yl(2);
        if sum(predicted_items(:,c)) == 0 %no predictived saccades
            title(['Item # ' num2str(c)])
        elseif sum(predicted_items(:,c)) == 2 % both sequences predicted
            title(['Item # ' num2str(c) ' Seq #1 ' num2str(percent_predicted(1,c)) '% Seq #2 ' num2str(num2str(percent_predicted(2,c))) '%'])
        elseif predicted_items(1,c)
            title(['Item # ' num2str(c) ' Seq #1 ' num2str(percent_predicted(1,c)) '%'])
        elseif predicted_items(2,c)
            title(['Item # ' num2str(c) ' Seq #2 ' num2str(num2str(percent_predicted(2,c))) '%'])
        end
    end
    for c = 1:4;
        subplot(2,2,c)
        ylim([0 max(ylims)]);
    end
    if multiunit(unit)
        subtitle(['Fixation-Locked Multiunit ' unit_names{unit} ' n =' num2str(num_trials)]);
    else
        subtitle(['Fixation-Locked ' unit_names{unit} ' n =' num2str(num_trials)]);
    end
    save_and_close_fig(figure_dir,[preprocessed_data_file(1:end-11) '_' unit_names{unit} '_Sequence_fixation_locked_analysis']);
end


save([data_dir preprocessed_data_file(1:8) '-Eyemovement_Locked_Sequence_results.mat'],...
    'predicted_thresh','twin','smval','predicted_items','percent_predicted',...
    'saccade_locked_firing','fixation_locked_firing','fixation_start_time',...
    'saccade_start_time','reaction_time','time_to_fixation','fixation_accuracy',...
    'fixation_duration','extrafixations','fixation_info','fixation_info_90',...
    'fixation_info_95','fixation_shuffled_info','saccade_info','saccade_info_90',...
    'saccade_info_95','saccade_shuffled_info');
end