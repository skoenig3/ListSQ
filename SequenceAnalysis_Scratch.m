clc
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Figures\';

predicted_rt = 135;%maximum "reaction time" for what constitutes as predictive, everything else is reactive
twin = 500;% how much time to take before and after saccade.
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
fixwin = 5;%size of fixation window on each crosshair
event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
item_codes = [23 25 27 29;
    24 26 28 30]; %row 1 item on, row 2 item of

trial_start_code = 15;
smval = 120;%gaussian 1/2 width for smoothing
predicted_thresh = 10;% percent of saccades that must be < 150 ms to item to constitue as a learned-predicted item
% for Vivian 5% was the false positive rate for 150 ms so wanted at least double the False positive rate
numshuffs = 10; %recommend this is between 100 & 1000, for bootstrapping to
% dtermine significance
info_type = 'temporal';
Fs = 1000;

multiunits = {zeros(1,16)};

preprocessed_data_file = 'TO160217_3-preprocessed.mat';


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
    
    event_codes = cfg.trl(trial).allval;
    event_times = event_times - event_times(event_codes == 15);
    for item = 2:4
        if ~isnan(trialdata.fixationnums(item))
            
            
            fixation_start = fixationtimes(1,trialdata.fixationnums(item));
            saccadeind = find(saccadetimes(2,:)+1 ==  fixation_start);
            if isempty(saccadeind)
                continue;
            end
            item_off_times(trial,item) =  event_times(event_codes ==  item_codes(2,item-1))...
                -saccadetimes(1,saccadeind);
        end
    end
    
end

