function Sequence_Fixation_AnalysisV2(data_dir,figure_dir,session_data)
% modfied from Sequence_Saccade_AnalysisV2 by Seth Konig August 30, 2016
% code now looks at all fixations to determine if eye movements modulate
% firinig rate. Comparable to List_Saccade_AnalysisV2.
%
% Function analyizes spike times correlated with fixations in the
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
%   2) Saves processed data to data_dir tagged with '-sequence_fixation_locked_results'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set Default Parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure_dir = [figure_dir 'Sequence Fixation Analysis\'];

twin = 500;% how much time to take before and after saccade.
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
fixwin = 5;%size of fixation window on each crosshair
event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
trial_start_code = 15;
smval = 60;%gaussian 1/2 width for smoothing
numshuffs = 1000; %recommend this is between 100 & 1000, for bootstrapping to
% dtermine significance
info_type = 'temporal';
Fs = 1000;
min_blks = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task = 'ListSQ';
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
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
if num_units == 0
    return
end
valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks);
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Get successful trials Information by Task---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get important task specific information
[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
overlap = find((sequence_locations{1}(1,:) ==  sequence_locations{2}(1,:)) & ...
    (sequence_locations{1}(2,:) ==  sequence_locations{2}(2,:)));

%preallocate space and parallel structure of cfg
successful_sequence_trials = NaN(1,length(cfg.trl));
which_sequence = NaN(1,length(cfg.trl));
for t = 1:num_trials;
    if sum(cfg.trl(t).allval == 3) >= 5; %in which sequence trials were rewarded
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

fixation_start_time = NaN(length(fixationstats),4);%when did fixation on item start
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
    
    fixation_numbers = trialdata.fixationnums; %fixation number for each item
    fixationtimes = fixationstats{trial}.fixationtimes;
    for item = 1:4
        if ~isnan(fixation_numbers(item))
            fixation_start_time(trial,item) = fixationtimes(1,fixation_numbers(item));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Firing Rate Locked to Eye Data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%also locked to "trial start" item 1 on or saccade to item 1
fixation_locked_firing = cell(1,num_units);
item_nums = cell(1,num_units);
sequence_nums =  cell(1,num_units);
for unit = 1:num_units
    fixation_locked_firing{unit} = NaN(4*num_trials,twin*2);
    item_nums{unit} = NaN(1,4*num_trials);
    sequence_nums{unit} = NaN(1,4*num_trials);
end

fix_ind = ones(1,unit);
for trial = 1:num_trials
    for unit = 1:num_units;
        if successful_sequence_trials(trial) >= valid_trials(1,unit) && successful_sequence_trials(trial) <= valid_trials(2,unit)
            spikes = find(data(unit).values{successful_sequence_trials(trial)});
            for c = 1:4;
                fixt = fixation_start_time(trial,c);
                
                if isnan(fixt)
                    continue
                elseif fixt < twin %starts before trial does
                    continue
                end
                
                if ~isnan(fixt)
                    fix_spikes = spikes(spikes > fixt-twin & spikes <= fixt+twin)-fixt+twin;
                    temp = zeros(1,twin*2);
                    temp(fix_spikes) = 1;
                    fixation_locked_firing{unit}(fix_ind(unit),:) = temp;
                    
                    item_nums{unit}(fix_ind(unit)) = c;
                    sequence_nums{unit}(fix_ind(unit)) = which_sequence(trial);
                    fix_ind(unit) = fix_ind(unit)+1;
                end
            end
        end
    end
end
fixation_locked_firing = laundry(fixation_locked_firing);
item_nums = laundry(item_nums);
sequence_nums = laundry(sequence_nums);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Determine if Units Are significantly modulated by fixation---%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for all fixations
fixation_info = [];
fixation_info.rate = NaN(1,num_units);
fixation_info.temporalstability = NaN(1,num_units);
fixation_info.rate_95 = NaN(1,num_units);
fixation_info.rate_prctile = NaN(1,num_units);
fixation_info.temporalstability_95 = NaN(1,num_units);
fixation_info.temporalstability_prctile = NaN(1,num_units);
fixation_info.shuffled_rate = cell(1,num_units);
fixation_info.shuffled_temporalstability = cell(1,num_units);

%for fixations with spikes only
fixation_info2 = [];
fixation_info2.rate = NaN(1,num_units);
fixation_info2.temporalstability = NaN(1,num_units);
fixation_info2.rate_95 = NaN(1,num_units);
fixation_info2.rate_prctile = NaN(1,num_units);
fixation_info2.temporalstability_95 = NaN(1,num_units);
fixation_info2.temporalstability_prctile = NaN(1,num_units);
fixation_info2.shuffled_rate = cell(1,num_units);
fixation_info2.shuffled_temporalstability = cell(1,num_units);

for unit = 1:num_units
    if ~isempty(fixation_locked_firing{unit})
        %don't want to run on trials with the first eye movements occuring with 500 (twin) ms of image onset
        [observed_info_rate,shuffled_info_rate]...
            = estimated_mutual_information(fixation_locked_firing{unit},numshuffs,info_type,smval,Fs);
        
        fixation_info.rate(unit) = observed_info_rate.skaggs;
        fixation_info.temporalstability(unit) = observed_info_rate.temporalstability;
        fixation_info.shuffled_rate{unit} = shuffled_info_rate.skaggs;
        fixation_info.shuffled_temporalstability{unit} = shuffled_info_rate.temporalstability;
        fixation_info.rate_95(unit) = prctile(fixation_info.shuffled_rate{unit},95);
        fixation_info.temporalstability_95(unit) = prctile(fixation_info.shuffled_temporalstability{unit},95);
        fixation_info.rate_prctile(unit) = 100*sum(fixation_info.rate(unit) > fixation_info.shuffled_rate{unit})/numshuffs;
        fixation_info.temporalstability_prctile = 100*sum(fixation_info.temporalstability(unit) > fixation_info.shuffled_temporalstability{unit})/numshuffs;
        
        
        %since many units appear sparse especially spatial ones want to
        %run on eye movemetns that actually have spikes so use these only
        fixation_firing = fixation_locked_firing{unit};
        fixation_firing(nansum(fixation_firing,2) == 0,:) = [];%remove fixations without spikes
        
        [observed_info_rate,shuffled_info_rate] = estimated_mutual_information(fixation_firing,numshuffs,info_type,smval,Fs);
        fixation_info2.rate(unit) = observed_info_rate.skaggs;
        fixation_info2.temporalstability(unit) = observed_info_rate.temporalstability;
        fixation_info2.shuffled_rate{unit} = shuffled_info_rate.skaggs;
        fixation_info2.shuffled_temporalstability{unit} = shuffled_info_rate.temporalstability;
        fixation_info2.rate_95(unit) = prctile(fixation_info2.shuffled_rate{unit},95);
        fixation_info2.temporalstability_95(unit) = prctile(fixation_info2.shuffled_temporalstability{unit},95);
        fixation_info2.rate_prctile(unit) = 100*sum(fixation_info2.rate(unit) > fixation_info2.shuffled_rate{unit})/numshuffs;
        fixation_info2.temporalstability_prctile = 100*sum(fixation_info2.temporalstability(unit) > fixation_info2.shuffled_temporalstability{unit})/numshuffs;
    end
end

%%%%%%%%%%%%%%%%%%%%
%---Plot Results---%
%%%%%%%%%%%%%%%%%%%%
t = -twin+1:twin;
for unit = 1:num_units
    if ~isempty(fixation_locked_firing{unit})
        fixation_firing = fixation_locked_firing{unit};
        fixation_firing(nansum(fixation_firing,2) == 0,:) = [];%remove fixations without spikes
        
        figure
        subplot(2,2,1)
        hold on
        dofill(t,fixation_locked_firing{unit},'blue',1,smval);
        dofill(t,fixation_firing,'green',1,smval);
        yl = ylim;
        if yl(1) < 0
            yl(1) = 0;
            ylim(yl);
        end
        plot([0 0],[yl(1) yl(2)],'k--')
        hold off
        xlabel('Time from Fixation Start (ms)')
        ylabel('Firing Rate')
        legend('Fixation','Fixation+')
        xlim([-twin twin])
        title_str = [];
        if fixation_info.rate_prctile(unit) > 95
            title_str = [title_str 'fix_{bit} = ' num2str(fixation_info.rate_prctile(unit),3) '%%'];
        end
        if  fixation_info2.rate_prctile(unit) > 95
            title_str = [title_str ' fix+_{bit} = ' num2str(fixation_info2.rate_prctile(unit),3) '%%'];
        end
        if ~isempty(title_str)
            title(sprintf(title_str))
        end
        
        fixation_firing = fixation_locked_firing{unit};
        seq = sequence_nums{unit};
        items = item_nums{unit};
        
        f1 = sum( fixation_firing,2);
        [~,fi] = sort(f1);
        [trial,time] = find(fixation_firing(fi,:) == 1);
        subplot(2,2,2)
        plot(time-500,(trial),'.k')
        ylabel('Ranked by Spike Count')
        xlabel('Time from fixation Start')
        xlim([-twin twin])
        ylim([0 max(trial)+1])
        
        subplot(2,2,3)
        [trial,time] = find(fixation_firing == 1);
        plot(time-twin,(trial),'.k')
        ylim([0 max(trial)])
        ylabel('Occurence #')
        xlabel('Time from fixation Start')
        xlim([-twin twin])
        
        subplot(2,2,4)
        b4 = 0;
        hold on
        for s = 1:2
            for c = 1:4
                fr =  fixation_firing(seq == s & items == c,:);
                [trial,time] = find(fr == 1);
                plot(time-twin,(trial+b4),'.k')
                if ~isempty(trial)
                    b4 = b4+max(trial);
                end
            end
        end
        hold off
        ylabel('By Sequence and Item #')
        xlabel('Time from fixation Start')
        xlim([-twin twin])
        
               
        n_str = [' n_{fix} = ' num2str(size(fixation_locked_firing{unit},1))];
        if multiunit(unit)
            subtitle(['Fixation-Locked Multiunit ' unit_stats{1,unit} n_str]);
        else
            subtitle(['Fixation-Locked ' unit_stats{1,unit} n_str]);
        end
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Sequence_Fixation_Locked Rasters']);
    end
end

save([data_dir task_file(1:8) '-Fixation_locked_Sequence_results.mat'],...
    'twin','smval','fixation_locked_firing','fixation_info','fixation_info2',...
    'item_nums','sequence_nums')