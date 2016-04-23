function Sequence_Saccade_AnalysisV2(data_dir,figure_dir,session_data,task,predicted_rt)
% written by Seth Konig September, 2014
% updated by Seth Konig April 7, 2015 to include estimating temporal
% updated SDK 1/11/17 to handlde new format and partial session data for
% vaild trials only.
% information locked to saccades and fixations on each item.
% Significantly modified original version (Find in OLd and Misc Code) to
% include newer main function for anlalyzing all sequence trial data. SDK
% 3/3/15. SDK Removed LFP from analysis 1/6/2016
%
% Function analyizes spike times correlated with eye movements in the
% sequence task.
%
% Inputs:
%   1) data_dir: direction where preprocessed data is located
%   2) figure_dir: location of where to put save figures
%   3) session_data: data containing relevant session data
%   4) task: what task was this data come from with i.e. 'ListSQ' or 'Sequence'
%   5) predicted_rt: maximum "reaction time" for what constitutes as predictive, everything else is reactive

%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-sequence_saccade_locked_results'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set Default Parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

twin = 750;% how much time to take before and after saccade.
long_window = 4000; %how long to look from "trial start"
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
fixwin = 5;%size of fixation window on each crosshair
event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
trial_start_code = 15;
smval = 60;%gaussian 1/2 width for smoothing
predicted_thresh = 10;% percent of saccades that must be < predict_rt to item to constitue as a learned-predicted item
% for Vivian 5% was the false positive rate for 150 ms so wanted at least double the False positive rate
numshuffs = 500; %recommend this is between 100 & 1000, for bootstrapping to
% dtermine significance
info_type = 'temporal';
Fs = 1000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end

%load trial data
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials',...
    'hdr','fixationstats');
num_trials = length(cfg.trl);

%NaNs are for start and end trials otherwise cut
valid_trials(1,isnan(valid_trials(1,:))) = 1;
valid_trials(2,isnan(valid_trials(2,:))) = num_trials;

%get unit data
unit_channels = find_desired_channels(cfg,'sig');
num_units = length(unit_channels);
unit_stats = cell(1,num_units);
mu = [];
for unit = 1:num_units
    unit_stats{1,unit} = hdr.label{data(unit).whichchannel};
    for n = 1:length(unit_names)
        if strcmpi(unit_stats{1,unit},unit_names{n})
            if  multiunit(n) > 3 %likley single unit
                mu = [mu 0];
            else %likely single unit
                mu = [mu 1];
            end
            break
        end
    end
end
if any(any(cellfun(@isempty,unit_stats)))
    warning('Unit data may not have been not imported correctly')
end
multiunit = mu;
unit_names = unit_stats;
clear mu unit_stats


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Get successful trials Information by Task---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch task
    case 'ListSQ'
        
        %get important task specific information
        [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        overlap = find((sequence_locations{1}(1,:) ==  sequence_locations{2}(1,:)) & ...
            (sequence_locations{1}(2,:) ==  sequence_locations{2}(2,:)));
        
    case 'Sequence'
        
        %get important task specific information
        [itmlist,sequence_items,sequence_locations] = read_Sequence_itm_and_cnd_files(item_file,cnd_file);
        overlap = find((sequence_locations{1}(1,:) ==  sequence_locations{2}(1,:)) & ...
            (sequence_locations{1}(2,:) ==  sequence_locations{2}(2,:)));
        
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---process eye data locked to trial events---%%%
fixstats = fixationstats;
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
    for item = 1:4
        if ~isnan(fixation_numbers(item))
            
            fixation_start_time(trial,item) = fixationtimes(1,fixation_numbers(item));
            saccadeind = find(saccadetimes(2,:)+1 ==  fixation_start_time(trial,item));
            
            if ~isempty(saccadeind)
                saccade_start_time(trial,item) = saccadetimes(1,saccadeind);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Firing Rate Locked to Eye Data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%also locked to "trial start" item 1 on or saccade to item 1
saccade_locked_firing = cell(4,num_units);
fixation_locked_firing = cell(4,num_units);
trial_nums = cell(1,num_units);
for unit = 1:num_units
    for c = 1:4
        saccade_locked_firing{c,unit} = NaN(num_trials,twin*2);
        fixation_locked_firing{c,unit} = NaN(num_trials,twin*2);
        trial_nums{c,unit} = NaN(1,num_trials);
    end
end

for trial = 1:num_trials
    for unit = 1:num_units;
        if successful_sequence_trials(trial) >= valid_trials(1,unit) && successful_sequence_trials(trial) <= valid_trials(2,unit)
            spikes = find(data(unit).values{successful_sequence_trials(trial)});
            for c = 1:4;
                fixt = fixation_start_time(trial,c);
                sact = saccade_start_time(trial,c);
                
                %if no saccade detected on item want to check if previous
                %trial has saccade/begining of fixation
                if isnan(sact) && isnan(fixt)
                    continue
                elseif isnan(sact) && ~isnan(fixt)
                    if fixt == 1 %only fixation is already started when the trial starts
                        disp('Fixation already on item when trial started. Looking back to previous trial')
                        if successful_sequence_trials(trial) > 1 %can't look back if it's on the first trial there's "no" eye data
                            
                            if ~isnan(fixstats{successful_sequence_trials(trial)-1}.XY(1,end)) %is the monkey looking at the screen at the end of the trial
                                
                                fixations1 = fixstats{successful_sequence_trials(trial)-1}.fixations(:,end);%last fixation previous trial
                                fixations2 =  fixstats{successful_sequence_trials(trial)}.fixations(:,1);%first fixation this trial
                                fixationtimes1 = fixstats{successful_sequence_trials(trial)-1}.fixationtimes(:,end); %fixation times of last fixation previous trial
                                saccadetimes = fixstats{successful_sequence_trials(trial)-1}.saccadetimes(:,end); %saccade times of last fixation previous trial
                                last_location = fixstats{successful_sequence_trials(trial)-1}.XY(:,end); %last eye location from previous trial
                                
                                if fixationtimes(2,end) > saccadetimes(2,end) %Yes, then did the trial end on fixation and already there
                                    %if it did are the eyes in the same location
                                    if sqrt(sum((fixations1-fixations2).^2)) < 36 %within 1.5 dva so approx yes
                                        spikes2 = find(data(unit).values{successful_sequence_trials(trial)-1});
                                        fixt = fixationtimes1(1);
                                        fix_spikes = spikes2(spikes2 > fixt-twin & spikes2 <= fixt+twin)-fixt+twin;
                                        fixt = fixationtimes1(1)-length(data(unit).values{successful_sequence_trials(trial)-1});%should be negative
                                        fix_spikes = [fix_spikes find(spikes < fixt+twin)-fixt+twin];
                                        temp = zeros(1,twin*2);
                                        temp(fix_spikes) = 1;
                                        fixation_locked_firing{c,unit}(trial,:) = temp;
                                        
                                        sact = saccadetimes(1);
                                        sac_spikes = spikes2(spikes2 > sact-twin & spikes2 <= sact+twin)-sact+twin;
                                        sact = saccadetimes(1)-length(data(unit).values{successful_sequence_trials(trial)-1});%should be negative
                                        sac_spikes = [sac_spikes find(spikes < sact+twin)-sact+twin];
                                        temp = zeros(1,twin*2);
                                        temp(sac_spikes) = 1;
                                        saccade_locked_firing{c,unit}(trial,:) = temp;
                                        
                                        continue %so we don't overwrite this trial in the next section
                                    else %if not possibly blinked at end of trial no good estimate of eye movements possible
                                        continue
                                    end
                                else %ended on a saccade %then fixation start is real then will continue on to the regular code
                                    
                                    if sqrt(sum((last_location-fixations2).^2)) < 36 %within 1.5 dva so approx yes
                                        spikes2 = find(data(unit).values{successful_sequence_trials(trial)-1});
                                        
                                        sact = saccadetimes(1);
                                        sac_spikes = spikes2(spikes2 > sact-twin & spikes2 <= sact+twin)-sact+twin;
                                        sact = saccadetimes(1)-length(data(unit).values{successful_sequence_trials(trial)-1});%should be negative
                                        sac_spikes = [sac_spikes find(spikes < sact+twin)-sact+twin];
                                        temp = zeros(1,twin*2);
                                        temp(sac_spikes) = 1;
                                        saccade_locked_firing{c,unit}(trial,:) = temp;
                                        
                                    else %if not possibly blinked at end of trial no good estimate of eye movements possible
                                        continue
                                    end
                                end
                            else %possibly blinked at end of trial no good estimate of eye movements possible
                                continue
                            end
                        else
                            continue %just ignore this trial then for now
                        end
                    else
                        warning(['Could not identify saccade but could locate fixation! Trial#:' num2str(trial)])
                    end
                elseif ~isnan(sact) && isnan(fixt)
                    error(['Could not identify fixation but could locate saccade! Trial#:' num2str(trial)])
                end
                
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
                trial_nums{c,unit}(trial) = trial;
            end
        end
    end
end

saccade_locked_firing = laundry(saccade_locked_firing);
fixation_locked_firing = laundry(fixation_locked_firing);
trial_nums = laundry(trial_nums);


fixation_info = NaN(4,num_units,3);
fixation_info_90 = NaN(4,num_units,3);
fixation_info_95 = NaN(4,num_units,3);
fixation_shuffled_info = cell(4,num_units,3);
saccade_info = NaN(4,num_units,3);
saccade_info_90 = NaN(4,num_units,3);
saccade_info_95 = NaN(4,num_units,3);
saccade_shuffled_info = cell(4,num_units,3);

for seq = 1:3
    for c = 1:size(saccade_locked_firing,1);
        for unit = 1:num_units
            if seq == 3 %both sequences
                [fixation_info(c,unit,seq),fixation_shuffled_info{c,unit}]...
                    = estimated_mutual_information(fixation_locked_firing{c,unit},numshuffs,info_type,smval,Fs);
                [saccade_info(c,unit,seq),saccade_shuffled_info{c,unit}]...
                    = estimated_mutual_information(saccade_locked_firing{c,unit},numshuffs,info_type,smval,Fs);
            else %do each sequence individually
                [fixation_info(c,unit,seq),fixation_shuffled_info{c,unit,seq}] = estimated_mutual_information(...
                    fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == seq,:),numshuffs,info_type,smval,Fs);
                [saccade_info(c,unit,seq),saccade_shuffled_info{c,unit}] = estimated_mutual_information(...
                    saccade_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == seq,:),numshuffs,info_type,smval,Fs);
            end
            fixation_info_90(c,unit,seq) = prctile(fixation_shuffled_info{c,unit,seq},90);
            fixation_info_95(c,unit,seq) = prctile(fixation_shuffled_info{c,unit,seq},95);
            saccade_info_90(c,unit,seq) = prctile(saccade_shuffled_info{c,unit,seq},90);
            saccade_info_95(c,unit,seq) = prctile(saccade_shuffled_info{c,unit,seq},95);
        end
    end
end

%average data across all items, row 1 sequence 1, row 2 sequence 2
all_saccade_locked_firing = cell(2,num_units);
all_fixation_locked_firing = cell(2,num_units);

for unit = 1:num_units
    for c = 1:4
        all_saccade_locked_firing{1,unit} = [all_saccade_locked_firing{1,unit};...
            saccade_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 1,:)];
        all_saccade_locked_firing{2,unit} = [all_saccade_locked_firing{2,unit};...
            saccade_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 2,:)];
        
        all_fixation_locked_firing{1,unit} = [all_fixation_locked_firing{1,unit};...
            fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 1,:)];
        all_fixation_locked_firing{2,unit} = [all_fixation_locked_firing{2,unit};...
            fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 2,:)];
    end
end

all_fixation_info = NaN(3,num_units); %row 1 seq 1, row 2 seq 2, row 3 seq1 & seq2
all_fixation_info_90 = NaN(3,num_units);
all_fixation_info_95 = NaN(3,num_units);
all_fixation_shuffled_info = cell(4,num_units);
all_saccade_info = NaN(3,num_units);
all_saccade_info_90 = NaN(3,num_units);
all_saccade_info_95 = NaN(3,num_units);
all_saccade_shuffled_info = cell(4,num_units);

for unit = 1:num_units
    for c = 1:3
        if c == 3 %for seq 1 and sequence 2
            [all_fixation_info(c,unit),all_fixation_shuffled_info{c,unit}]...
                = estimated_mutual_information([all_fixation_locked_firing{1,unit}; all_fixation_locked_firing{2,unit}]...
                ,numshuffs,info_type,smval,Fs);
            [all_saccade_info(c,unit),all_saccade_shuffled_info{c,unit}]...
                = estimated_mutual_information([all_saccade_locked_firing{1,unit}; all_saccade_locked_firing{2,unit}]...
                ,numshuffs,info_type,smval,Fs);
        else
            [all_fixation_info(c,unit),all_fixation_shuffled_info{c,unit}]...
                = estimated_mutual_information(all_fixation_locked_firing{c,unit},numshuffs,info_type,smval,Fs);
            [all_saccade_info(c,unit),all_saccade_shuffled_info{c,unit}]...
                = estimated_mutual_information(all_saccade_locked_firing{c,unit},numshuffs,info_type,smval,Fs);
        end
        
        all_fixation_info_90(c,unit) = prctile(all_fixation_shuffled_info{c,unit},90);
        all_fixation_info_95(c,unit) = prctile(all_fixation_shuffled_info{c,unit},95);
        all_saccade_info_90(c,unit) = prctile(all_saccade_shuffled_info{c,unit},90);
        all_saccade_info_95(c,unit) = prctile(all_saccade_shuffled_info{c,unit},95);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Analyzing Predictive vs Reactive Spike Data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
predicted(successful_sequence_trials == 1,:) = 0; predicted(successful_sequence_trials == 11,:) = 0; %dont ever want 1st trial to accidentally skew results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plots locked to saccades---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = -twin:twin-1;
for unit = 1:num_units
    figure
    ylims = NaN(1,8);
    for c = 1:4
        subplot(2,2,c)
        hold on
        if ~all(all(isnan(saccade_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 1 & ~predicted(trial_nums{c,unit},c)',:))))
            dofill(t,saccade_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 1 & ~predicted(trial_nums{c,unit},c)',:),'blue',1,smval); %reactive sequence 1
        end
        if ~all(all(isnan(saccade_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 2 & ~predicted(trial_nums{c,unit},c)',:))))
            dofill(t,saccade_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 2 & ~predicted(trial_nums{c,unit},c)',:),'red',1,smval);%reactive sequence 2
        end
        if predicted_items(1,c);
            if ~all(all(isnan(saccade_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 1 & predicted(trial_nums{c,unit},c)'))))
                dofill(t,saccade_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 1 & predicted(trial_nums{c,unit},c)',:),'cyan',1,smval); %predictive sequence 1
            end
        end
        if predicted_items(2,c)
            if ~all(all(isnan(saccade_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 2 & predicted(trial_nums{c,unit},c)',:))))
                dofill(t,saccade_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 2 & predicted(trial_nums{c,unit},c)',:),'magenta',1,smval);%predictive sequence 2
            end
        end
        hold off
        xlabel('Time from Saccade (ms)')
        ylabel('Firing Rate (Hz)')
        yl = ylim;
        ylims(c) = yl(2);
        
        title_str = [];
        if overlap == c
            title_str = [title_str '(Overlapping) '];
        end
        
        if any(fixation_info(c,unit,:) > fixation_info_95(c,unit,:))
            title_str = ['fix_{95} = ' num2str(fixation_info(c,unit,end))];
        elseif any(fixation_info(c,unit,:) > fixation_info_90(c,unit,:))
            title_str = ['fix_{90} = ' num2str(fixation_info(c,unit,end))];
        end
        
        if sum(predicted_items(:,c)) == 0 %no predictived saccades
            title(['Item # ' num2str(c) ' ' title_str])
        elseif sum(predicted_items(:,c)) == 2 % both sequences predicted
            n1 = sum(time_to_fixation(which_sequence(trial_nums{c,unit}) == 1,c) < predicted_rt);
            n2 = sum(time_to_fixation(which_sequence(trial_nums{c,unit}) == 2,c) < predicted_rt);
            title(['Item # ' num2str(c) ' Seq #1 n_p=' num2str(n1) ' Seq #2 n_p' num2str(n2) ' ' title_str])
        elseif predicted_items(1,c)
            n1 = sum(time_to_fixation(which_sequence(trial_nums{c,unit}) == 1,c) < predicted_rt);
            title(['Item # ' num2str(c) ' Seq #1 n_p=' num2str(n1) ' ' title_str])
        elseif predicted_items(2,c)
            n2 = sum(time_to_fixation(which_sequence(trial_nums{c,unit}) == 2,c) < predicted_rt);
            title(['Item # ' num2str(c) ' Seq #1 n_p=' num2str(n2) ' ' title_str])
        end
    end
    for c = 1:4;
        subplot(2,2,c)
        ylim([0 max(ylims)]);
    end
    
    num_str = [' n_1=' num2str(sum(which_sequence(trial_nums{c,unit}) == 1)) ' ' ...
        'n_21=' num2str(sum(which_sequence(trial_nums{c,unit}) == 2))];
    if multiunit(unit)
        subtitle(['Saccade-Locked Multiunit ' unit_names{unit} num_str]);
    else
        subtitle(['Saccade-Locked ' unit_names{unit} num_str]);
    end
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_Sequence_Saccade_Locked_analysis']);
    
    figure
    hold on
    dofill(t,[all_saccade_locked_firing{1,unit}; all_saccade_locked_firing{2,unit}],'black',1,smval);%all items both sequences
    dofill(t,all_saccade_locked_firing{1,unit},'blue',1,smval);%all items sequence 1
    dofill(t,all_saccade_locked_firing{2,unit},'red',1,smval);%all items sequence 2
    hold off
    xlabel('Time from Saccade (ms)')
    ylabel('Firing Rate (Hz)')
    
    title_str = ['All items Saccade Locked n = ' num2str(size(all_saccade_locked_firing{1,unit},1)+...
        size(all_saccade_locked_firing{2,unit},1))];
    if all_saccade_info(3,unit) > all_saccade_info_95(3,unit)
        title_str = [title_str ' ' num2str(all_saccade_info(3,unit)) ' bits_{95}'];
    elseif all_saccade_info(3,unit) > all_saccade_info_90(3,unit)
        title_str = [title_str ' ' num2str(all_saccade_info(3,unit)) ' bits_{90}'];
    end
    if multiunit(unit)
        title(['Saccade-Locked Multiunit ' unit_names{unit} title_str]);
    else
        title(['Saccade-Locked ' unit_names{unit} title_str]);
    end
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_Sequence_Saccade_Locked_analysis_All_items']);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plots locked to fixations---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = -twin:twin-1;
for unit = 1:num_units
    figure
    ylims = NaN(1,8);
    for c = 1:4
        subplot(2,2,c)
        hold on
        if ~all(all(isnan(fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 1 & ~predicted(trial_nums{c,unit},c)',:))))
            dofill(t,fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 1 & ~predicted(trial_nums{c,unit},c)',:),'blue',1,smval); %reactive sequence 1
        end
        if ~all(all(isnan(fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 2 & ~predicted(trial_nums{c,unit},c)',:))))
            dofill(t,fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 2 & ~predicted(trial_nums{c,unit},c)',:),'red',1,smval);%reactive sequence 2
        end
        if predicted_items(1,c)
            if ~all(all(isnan(fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 1 & predicted(trial_nums{c,unit},c)',:))))
                dofill(t,fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 1 & predicted(trial_nums{c,unit},c)',:),'cyan',1,smval); %predictive sequence 1
            end
        end
        if predicted_items(2,c)
            if ~all(all(isnan(fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 2 & predicted(trial_nums{c,unit},c)',:))))
                dofill(t,fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == 2 & predicted(trial_nums{c,unit},c)',:),'magenta',1,smval);%predictive sequence 2
            end
        end
        hold off
        xlabel('Time from Fixation (ms)')
        ylabel('Firing Rate (Hz)')
        yl = ylim;
        ylims(c) = yl(2);
        
        
        title_str = [];
        if overlap == c
            title_str = [title_str '(Overlapping) '];
        end
        
        if any(saccade_info(c,unit,:) > saccade_info_95(c,unit,:))
            title_str = [title_str 'sac_{95} = ' num2str(saccade_info(c,unit,end))];
        elseif any(saccade_info(c,unit,:) > saccade_info_90(c,unit,:))
            title_str = [title_str 'sac_{90} = ' num2str(saccade_info(c,unit,end))];
        end
        
        
        if sum(predicted_items(:,c)) == 0 %no predictived saccades
            title(['Item # ' num2str(c) ' ' title_str])
        elseif sum(predicted_items(:,c)) == 2 % both sequences predicted
            n1 = sum(time_to_fixation(which_sequence(trial_nums{c,unit}) == 1,c) < predicted_rt);
            n2 = sum(time_to_fixation(which_sequence(trial_nums{c,unit}) == 2,c) < predicted_rt);
            title(['Item # ' num2str(c) ' Seq #1 n_p=' num2str(n1) ' Seq #2 n_p' num2str(n2) ' ' title_str])
        elseif predicted_items(1,c)
            n1 = sum(time_to_fixation(which_sequence(trial_nums{c,unit}) == 1,c) < predicted_rt);
            title(['Item # ' num2str(c) ' Seq #1 n_p=' num2str(n1) ' ' title_str])
        elseif predicted_items(2,c)
            n2 = sum(time_to_fixation(which_sequence(trial_nums{c,unit}) == 2,c) < predicted_rt);
            title(['Item # ' num2str(c) ' Seq #1 n_p=' num2str(n2) ' ' title_str])
        end
        
    end
    for c = 1:4;
        subplot(2,2,c)
        ylim([0 max(ylims)]);
    end
    
    num_str = [' n_1=' num2str(sum(which_sequence(trial_nums{c,unit}) == 1)) ' ' ...
        'n_2=' num2str(sum(which_sequence(trial_nums{c,unit}) == 2))];
    if multiunit(unit)
        subtitle(['Saccade-Locked Multiunit ' unit_names{unit} num_str]);
    else
        subtitle(['Saccade-Locked ' unit_names{unit} num_str]);
    end
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_Sequence_fixation_locked_analysis']);
    
    figure
    hold on
    dofill(t,[all_fixation_locked_firing{1,unit}; all_fixation_locked_firing{2,unit}],'black',1,smval);%all items both sequences
    dofill(t,all_fixation_locked_firing{1,unit},'blue',1,smval);%all items sequence 1
    dofill(t,all_fixation_locked_firing{2,unit},'red',1,smval);%all items sequence 2
    hold off
    xlabel('Time from fixation (ms)')
    ylabel('Firing Rate (Hz)')
    
    title_str = ['All items Fixation Locked n = ' num2str(size(all_fixation_locked_firing{1,unit},1)+...
        size(all_fixation_locked_firing{2,unit},1))];
    if all_fixation_info(3,unit) > all_fixation_info_95(3,unit)
        title_str = [title_str ' ' num2str(all_fixation_info(3,unit)) ' bits_{95}'];
    elseif all_fixation_info(3,unit) > all_fixation_info_90(3,unit)
        title_str = [title_str ' ' num2str(all_fixation_info(3,unit)) ' bits_{90}'];
    end
    if multiunit(unit)
        title(['Fixation-Locked Multiunit ' unit_names{unit} title_str]);
    else
        title(['Fixation-Locked ' unit_names{unit} title_str]);
    end
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_Sequence_fixation_Locked_analysis_All_items']);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Firing Rate Locked to "Trial Start"---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_locked_firing = cell(2,num_units);%row 1 item 1 on, row 2 saccade to item
trial_nums2 = cell(1,num_units);
for unit = 1: num_units
    trial_locked_firing{1,unit} = NaN(num_trials,twin+long_window);
    trial_locked_firing{2,unit} = NaN(num_trials,twin+long_window);
    trial_nums2{unit} = NaN(1,num_trials);
end


for trial = 1:num_trials
    for unit = 1:num_units;
        if successful_sequence_trials(trial) >= valid_trials(1,unit) && successful_sequence_trials(trial) <= valid_trials(2,unit)
            spikes = find(data(unit).values{successful_sequence_trials(trial)});
            sact = saccade_start_time(trial,1);
            if  sact >= twin %only take trials where sact >= twin so no predictive stuff for trial start
                
                item_on = cfg.trl(trial).alltim(cfg.trl(trial).allval == 23)-cfg.trl(trial).alltim(1);
                sac_spikes = spikes(spikes > item_on-twin & spikes <= item_on+long_window)-item_on+twin;
                temp = zeros(1,twin+long_window);
                temp(sac_spikes) = 1;
                trial_locked_firing{1,unit}(trial,:) = temp;
                
                sac_spikes = spikes(spikes > sact-twin & spikes <= sact+long_window)-sact+twin;
                temp = zeros(1,twin+long_window);
                temp(sac_spikes) = 1;
                trial_locked_firing{2,unit}(trial,:) = temp;
                
                trial_nums2{unit}(trial) = trial;
            end
        end
    end
end

trial_nums2 = laundry(trial_nums2);
trial_locked_firing = laundry(trial_locked_firing);

trial_info = NaN(2,num_units,3);
trial_info_90 = NaN(2,num_units,3);
trial_info_95 = NaN(2,num_units,3);
trial_shuffled_info = cell(2,num_units,3);

for seq = 1:3
    for c = 1:size(trial_locked_firing,1);
        for unit = 1:num_units
            if seq == 3 %both sequences
                [trial_info(c,unit,seq),trial_shuffled_info{c,unit}]...
                    = estimated_mutual_information(trial_locked_firing{c,unit},numshuffs,info_type,smval,Fs);
            else %do each sequence individually
                [trial_info(c,unit,seq),trial_shuffled_info{c,unit,seq}] = estimated_mutual_information(...
                    trial_locked_firing{c,unit}(which_sequence(trial_nums2{unit}) == seq,:),numshuffs,info_type,smval,Fs);
            end
            trial_info_90(c,unit,seq) = prctile(trial_shuffled_info{c,unit,seq},90);
            trial_info_95(c,unit,seq) = prctile(trial_shuffled_info{c,unit,seq},95);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plots locked to "Trial Start" -%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = -twin:long_window-1;
for unit = 1:num_units
    figure
    
    subplot(2,1,1)
    hold on
    dofill(t,trial_locked_firing{1,unit}(which_sequence(trial_nums2{unit}) == 1,:),'blue',1,smval); %sequence 1
    dofill(t,trial_locked_firing{1,unit}(which_sequence(trial_nums2{unit}) == 2,:),'red',1,smval); %sequence 2
    hold off
    xlabel('Time from Saccade to Item 1 (ms)')
    ylabel('Firing Rate (Hz)')
    yl = ylim;
    ylims(1) = yl(2);
    
    if any(trial_info(1,unit,:) > trial_info_95(1,unit,:))
        title(['sac_{95} = ' num2str(trial_info(1,unit,end))]);
    elseif any(trial_info(1,unit,:) > trial_info_90(1,unit,:))
        title(['sac_{95} = ' num2str(trial_info(1,unit,end))]);
    end
    
    subplot(2,1,2)
    hold on
    dofill(t,trial_locked_firing{2,unit}(which_sequence(trial_nums2{unit}) == 1,:),'blue',1,smval); %sequence 1
    dofill(t,trial_locked_firing{2,unit}(which_sequence(trial_nums2{unit}) == 2,:),'red',1,smval); %sequence 2
    hold off
    xlabel('Time from Item 1 On (ms)')
    ylabel('Firing Rate (Hz)')
    yl = ylim;
    ylims(2) = yl(2);
    
    if any(trial_info(2,unit,:) > trial_info_95(2,unit,:))
        title(['start_{95} = ' num2str(trial_info(2,unit,end))]);
    elseif any(trial_info(2,unit,:) > trial_info_90(2,unit,:))
        title(['start_{95} = ' num2str(trial_info(2,unit,end))]);
    end
    
    
    subplot(2,1,1)
    ylim([0 max(ylims)])
    subplot(2,1,2)
    ylim([0 max(ylims)])
    
    
    num_str = [' n_1=' num2str(sum(which_sequence(trial_nums2{1,unit}) == 1)) ' ' ...
        'n_2=' num2str(sum(which_sequence(trial_nums2{1,unit}) == 2))];
    if multiunit(unit)
        subtitle(['Trial_Start-Locked Multiunit ' unit_names{unit} num_str]);
    else
        subtitle(['Trial_Start-Locked ' unit_names{unit} num_str]);
    end
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_Sequence_Trial_Start_locked_analysis']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Save Relevant Data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%consolidate for better organization
trial_locked_info.rate = trial_info;
trial_locked_info.shuffled_info_rate = trial_shuffled_info;
trial_locked_info.shuffled_95_percentile = trial_info_95;
trial_locked_info.shuffled_90_percentile = trial_info_90;

fixation_locked_info.rate = fixation_info;
fixation_locked_info.shuffled_info_rate = fixation_shuffled_info;
fixation_locked_info.shuffled_95_percentile = fixation_info_95;
fixation_locked_info.shuffled_90_percentile = fixation_info_90;

saccade_locked_info.rate = saccade_info;
saccade_locked_info.shuffled_info_rate = saccade_shuffled_info;
saccade_locked_info.shuffled_95_percentile = saccade_info_95;
saccade_locked_info.shuffled_90_percentile = saccade_info_90;

all_fixation_locked_info.rate = all_fixation_info;
all_fixation_locked_info.shuffled_info_rate = all_fixation_shuffled_info;
all_fixation_locked_info.shuffled_95_percentile = all_fixation_info_95;
all_fixation_locked_info.shuffled_90_percentile = all_fixation_info_90;

all_saccade_locked_info.rate = all_saccade_info;
all_saccade_locked_info.shuffled_info_rate = all_saccade_shuffled_info;
all_saccade_locked_info.shuffled_95_percentile = all_saccade_info_95;
all_saccade_locked_info.shuffled_90_percentile = all_saccade_info_90;


save([data_dir task_file(1:8) '-Eyemovement_Locked_Sequence_results.mat'],...
    'predicted_thresh','twin','smval','predicted_items','percent_predicted',...
    'saccade_locked_firing','fixation_locked_firing','fixation_start_time',...
    'saccade_start_time','reaction_time','time_to_fixation','fixation_accuracy',...
    'fixation_duration','extrafixations','trial_nums','trial_locked_firing',...
    'trial_nums2','trial_locked_info','fixation_locked_info','saccade_locked_info',...
    'which_sequence','all_saccade_locked_info','all_fixation_locked_info');

end