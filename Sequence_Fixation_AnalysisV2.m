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
%   4) predicted_rt: maximum "reaction time" for what constitutes as predictive, everything else is reactive

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
item_event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
trial_start_code = 15;
smval = 50;%gaussian 1/2 width for smoothing
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

valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,'ListSQ');
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---process eye data locked to trial events---%%%
num_trials = length(successful_sequence_trials);
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
        event_codes,event_times,1);
    
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
event_locked_firing = cell(1,num_units);
event_off_locked_firing = cell(1,num_units);
item_nums = cell(1,num_units);
sequence_nums =  cell(1,num_units);
for unit = 1:num_units
    fixation_locked_firing{unit} = NaN(4*num_trials,twin*2);
    event_locked_firing{unit} =  NaN(4*num_trials,twin*2);
    event_off_locked_firing{unit} =  NaN(4*num_trials,twin*2);
    item_nums{unit} = NaN(1,4*num_trials);
    sequence_nums{unit} = NaN(1,4*num_trials);
    reaction_time{unit} = NaN(1,4*num_trials);
    trial_numbers{unit} = NaN(1,4*num_trials);
end

fix_ind = ones(1,unit);
for trial = 1:num_trials
    for unit = 1:num_units;
        if successful_sequence_trials(trial) >= valid_trials(1,unit) && successful_sequence_trials(trial) <= valid_trials(2,unit)
            spikes = find(data(unit).values{successful_sequence_trials(trial)});
            for c = 1:4;
                fixt = fixation_start_time(trial,c);
                
                item_on = cfg.trl(trial).alltim(cfg.trl(trial).allval == item_event_codes(2*c-1)); %when item turned on
                item_on = item_on-cfg.trl(trial).alltim(1);
                item_off = cfg.trl(trial).alltim(cfg.trl(trial).allval == item_event_codes(2*c)); %when item turned off
                item_off = item_off-cfg.trl(trial).alltim(1);
                
                if isnan(fixt)
                    continue
                elseif fixt < twin %starts before trial does
                    continue
                end
                
                if ~isnan(fixt)
                    
                    %---Fixation Aligned Acitivity---%
                    fix_spikes = spikes(spikes > fixt-twin & spikes <= fixt+twin)-fixt+twin;
                    temp = zeros(1,twin*2);
                    temp(fix_spikes) = 1;
                    fixation_locked_firing{unit}(fix_ind(unit),:) = temp;
                    
                    
                    %---Event Aligned Activity---%
                    event_spikes = spikes(spikes >item_on-twin & spikes <= item_on+twin)-item_on+twin;
                    temp = zeros(1,twin*2);
                    temp(event_spikes) = 1;
                    event_locked_firing{unit}(fix_ind(unit),:) = temp;
                    
                    %---Event End/Off Aligned Activity---%
                    event_spikes = spikes(spikes >item_off-twin & spikes <= item_off+twin)-item_off+twin;
                    temp = zeros(1,twin*2);
                    temp(event_spikes) = 1;
                    event_off_locked_firing{unit}(fix_ind(unit),:) = temp;
                    
                    %---Trial Paramaters---%
                    item_nums{unit}(fix_ind(unit)) = c;
                    sequence_nums{unit}(fix_ind(unit)) = which_sequence(trial);
                    reaction_time{unit}(fix_ind(unit)) = fixt-item_on;
                    trial_numbers{unit}(fix_ind(unit)) = trial;
                    fix_ind(unit) = fix_ind(unit)+1;
                end
            end
        end
    end
end
fixation_locked_firing = laundry(fixation_locked_firing);
event_locked_firing = laundry(event_locked_firing);
event_off_locked_firing = laundry(event_off_locked_firing);
item_nums = laundry(item_nums);
sequence_nums = laundry(sequence_nums);
reaction_time = laundry(reaction_time);
trial_numbers = laundry(trial_numbers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Determine if Units Are significantly modulated by fixation---%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for all fixations
fixation_info = [];
fixation_info.rate = NaN(1,num_units);
fixation_info.temporalstability = NaN(1,num_units);
fixation_info.rate_prctile = NaN(1,num_units);
fixation_info.temporalstability_prctile = NaN(1,num_units);
fixation_info.shuffled_rate = cell(1,num_units);
fixation_info.shuffled_temporalstability = cell(1,num_units);

%for appearance for all items
all_items_info = [];
all_items_info.rate = NaN(4,2,num_units);
all_items_info.temporalstability = NaN(4,2,num_units);
all_items_info.rate_prctile = NaN(4,2,num_units);
all_items_info.temporalstability_prctile = NaN(4,2,num_units);
all_items_info.shuffled_rate = cell(4,2,num_units);
all_items_info.shuffled_temporalstability = cell(4,2,num_units);


%for fixations on each item
sequence_info = [];
sequence_info.rate = NaN(4,2,num_units);
sequence_info.temporalstability = NaN(4,2,num_units);
sequence_info.rate_prctile = NaN(4,2,num_units);
sequence_info.temporalstability_prctile = NaN(4,2,num_units);
sequence_info.shuffled_rate = cell(4,2,num_units);
sequence_info.shuffled_temporalstability = cell(4,2,num_units);

%for appearance of each item
item_info = [];
item_info.rate = NaN(4,2,num_units);
item_info.temporalstability = NaN(4,2,num_units);
item_info.rate_prctile = NaN(4,2,num_units);
item_info.temporalstability_prctile = NaN(4,2,num_units);
item_info.shuffled_rate = cell(4,2,num_units);
item_info.shuffled_temporalstability = cell(4,2,num_units);

%for disappearance of items
items_off_info = [];
items_off_info.rate = NaN(4,2,num_units);
items_off_info.temporalstability = NaN(4,2,num_units);
items_off_info.rate_prctile = NaN(4,2,num_units);
items_off_info.temporalstability_prctile = NaN(4,2,num_units);
items_off_info.shuffled_rate = cell(4,2,num_units);
items_off_info.shuffled_temporalstability = cell(4,2,num_units);

for unit = 1:num_units
    if ~isempty(fixation_locked_firing{unit})
        
        %all fixation aligned activity
        if sum(sum(fixation_locked_firing{unit}))
            [observed_info_rate,shuffled_info_rate]...
                = estimated_mutual_information(fixation_locked_firing{unit},numshuffs,info_type,smval,Fs);
            fixation_info.rate(unit) = observed_info_rate.skaggs;
            fixation_info.temporalstability(unit) = observed_info_rate.temporalstability(1);
            fixation_info.shuffled_rate{unit} = shuffled_info_rate.skaggs;
            fixation_info.shuffled_temporalstability{unit} = shuffled_info_rate.temporalstability(1,:);
            fixation_info.rate_prctile(unit) = 100*sum(fixation_info.rate(unit) > fixation_info.shuffled_rate{unit})/numshuffs;
            fixation_info.temporalstability_prctile(unit) = 100*sum(fixation_info.temporalstability(unit)...
                > fixation_info.shuffled_temporalstability{unit})/numshuffs;
        else
            fixation_info.rate(unit) = 0;
            fixation_info.temporalstability(unit) = 0;
        end
        
        %all item on aligned activity
        if sum(sum(event_locked_firing{unit})) > 0
            [observed_info_rate,shuffled_info_rate]...
                = estimated_mutual_information(event_locked_firing{unit},numshuffs,info_type,smval,Fs);
            all_items_info.rate(unit) = observed_info_rate.skaggs;
            all_items_info.temporalstability(unit) = observed_info_rate.temporalstability(1);
            all_items_info.shuffled_rate{unit} = shuffled_info_rate.skaggs;
            all_items_info.shuffled_temporalstability{unit} = shuffled_info_rate.temporalstability(1,:);
            all_items_info.rate_prctile(unit) = 100*sum(all_items_info.rate(unit) > all_items_info.shuffled_rate{unit})/numshuffs;
            all_items_info.temporalstability_prctile(unit) = 100*sum(all_items_info.temporalstability(unit)...
                > all_items_info.shuffled_temporalstability{unit})/numshuffs;
        else
            all_items_info.rate(unit) = 0;
            all_items_info.temporalstability(unit) = 0;
        end
        
        %for activity aligned to fixation on each item
        for s = 1:2
            for c = 1:4
                fix_aligned = fixation_locked_firing{unit}((sequence_nums{unit} == s) & (item_nums{unit} == c),:);
                if sum(sum(fix_aligned)) > 0
                    [observed_info_rate,shuffled_info_rate]...
                        = estimated_mutual_information(fix_aligned,numshuffs,info_type,smval,Fs);
                    sequence_info.rate(c,s,unit) = observed_info_rate.skaggs;
                    sequence_info.temporalstability(c,s,unit) = observed_info_rate.temporalstability(1);
                    sequence_info.shuffled_rate{c,s,unit} = shuffled_info_rate.skaggs;
                    sequence_info.shuffled_temporalstability{c,s,unit} = shuffled_info_rate.temporalstability(1,:);
                    sequence_info.rate_prctile(c,s,unit) = 100*sum(sequence_info.rate(c,s,unit) > sequence_info.shuffled_rate{c,s,unit})/numshuffs;
                    sequence_info.temporalstability_prctile(c,s,unit) = 100*sum(sequence_info.temporalstability(c,s,unit) ...
                        > sequence_info.shuffled_temporalstability{c,s,unit})/numshuffs;
                else
                    sequence_info.rate(c,s,unit) = 0;
                    sequence_info.temporalstability(c,s,unit) = 0;
                end
                
                
                %event on aligned activity
                item_aligned = event_locked_firing{unit}((sequence_nums{unit} == s) & (item_nums{unit} == c),:);
                if sum(sum(item_aligned)) > 0
                    [observed_info_rate,shuffled_info_rate]...
                        = estimated_mutual_information(item_aligned,numshuffs,info_type,smval,Fs);
                    item_info.rate(c,s,unit) = observed_info_rate.skaggs;
                    item_info.temporalstability(c,s,unit) = observed_info_rate.temporalstability(1);
                    item_info.shuffled_rate{c,s,unit} = shuffled_info_rate.skaggs;
                    item_info.shuffled_temporalstability{c,s,unit} = shuffled_info_rate.temporalstability(1,:);
                    item_info.rate_prctile(c,s,unit) = 100*sum(item_info.rate(c,s,unit) > item_info.shuffled_rate{c,s,unit})/numshuffs;
                    item_info.temporalstability_prctile(c,s,unit) = 100*sum(item_info.temporalstability(c,s,unit)...
                        > item_info.shuffled_temporalstability{c,s,unit})/numshuffs;
                else
                    item_info.rate(c,s,unit) = 0;
                    item_info.temporalstability(c,s,unit) = 0;
                end
                
                %event off aligned activity
                item_aligned = event_off_locked_firing{unit}((sequence_nums{unit} == s) & (item_nums{unit} == c),:);
                if sum(sum(item_aligned)) > 0
                    [observed_info_rate,shuffled_info_rate]...
                        = estimated_mutual_information(item_aligned,numshuffs,info_type,smval,Fs);
                    item_off_info.rate(c,s,unit) = observed_info_rate.skaggs;
                    item_off_info.temporalstability(c,s,unit) = observed_info_rate.temporalstability(1);
                    item_off_info.shuffled_rate{c,s,unit} = shuffled_info_rate.skaggs;
                    item_off_info.shuffled_temporalstability{c,s,unit} = shuffled_info_rate.temporalstability(1,:);
                    item_off_info.rate_prctile(c,s,unit) = 100*sum(item_off_info.rate(c,s,unit) > item_off_info.shuffled_rate{c,s,unit})/numshuffs;
                    item_off_info.temporalstability_prctile(c,s,unit) = 100*sum(item_off_info.temporalstability(c,s,unit)...
                        > item_off_info.shuffled_temporalstability{c,s,unit})/numshuffs;
                else
                    item_off_info.rate(c,s,unit) = 0;
                    item_off_info.temporalstability(c,s,unit) = 0;
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%
%---Plot Results---%
%%%%%%%%%%%%%%%%%%%%
t = -twin+1:twin;
for unit = 1:num_units
    if ~isempty(fixation_locked_firing{unit})
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---All Fixation Aligned Activity---%
        fixation_firing = fixation_locked_firing{unit};
        event_firing = event_locked_firing{unit};
        
        ff = nandens(fixation_firing,50,'gauss',1000);
        
        yls = NaN(2,2);
        
        figure
        
        subplot(2,2,1)
        dofill(t,fixation_firing,'black',1,smval);
        yls(1,:) = ylim;
        xlabel('Time from Fixation Start/Item On (ms)')
        ylabel('Firing Rate (Hz)')
        xlim([-twin twin])
        title_str = 'Fixation Aligned: ';
        if fixation_info.rate_prctile(unit) > 95
            title_str = [title_str 'fix_{bit} = ' num2str(fixation_info.rate(unit),3) ...
                ' (' num2str(fixation_info.rate_prctile(unit),3) '%%)'];
        end
        if fixation_info.temporalstability_prctile(unit) > 95
            title_str = [title_str ', \\rho_{1/2} = ' ...
                num2str(fixation_info.temporalstability(unit),2) ' (' ...
                num2str(fixation_info.temporalstability_prctile(unit),3) '%%)'];
        end
        title(sprintf(title_str))
        
        subplot(2,2,3)
        [trial,time] = find(fixation_firing == 1);
        plot(time-twin,(trial),'.k')
        ylim([0 max(trial)])
        ylabel('Occurence #')
        xlabel('Time from fixation Start (ms)')
        title('Fixation Aligned')
        xlim([-twin twin])
        ylim([0 max(trial)+1])
        box off
        
        subplot(2,2,2)
        hold on
        dofill(t,event_firing,'blue',1,smval);
        yls(2,:) = ylim;
        xlabel('Time from Fixation Start/Item On (ms)')
        ylabel('Firing Rate (Hz)')
        xlim([-twin twin])
        title_str = 'Event Aligned: ';
        if all_items_info.rate_prctile(unit) > 95
            title_str = [title_str '  Event_{bit} = ' num2str(all_items_info.rate(unit),3) ...
                ' (' num2str(all_items_info.rate_prctile(unit),3) '%%)'];
        end
        if all_items_info.temporalstability_prctile(unit) > 95
            title_str = [title_str ', \\rho_{1/2} = ' ...
                num2str( all_items_info.temporalstability(unit),2) ' (' ...
                num2str( all_items_info.temporalstability_prctile(unit),3) '%%)'];
        end
        title(sprintf(title_str))
        
        subplot(2,2,4)
        [trial,time] = find(event_firing == 1);
        plot(time-twin,(trial),'.k')
        if ~isempty(trial)
            ylim([0 max(trial)+1])
        end
        ylabel('Occurence #')
        xlabel('Time from Item On (ms)')
        title('Event Aligned')
        xlim([-twin twin])
        ylim([0 max(trial)+1])
        box off
        
        ymin = min(yls(:,1));
        ymin(ymin < 0) = 0;
        ymax = max(yls(:,2));
        
        subplot(2,2,1)
        ylim([ymin ymax])
        hold on
        plot([0 0],[ymin ymax],'k--')
        hold off
        
        subplot(2,2,2)
        ylim([ymin ymax])
        hold on
        plot([0 0],[ymin ymax],'k--')
        hold off
        
        n_str = [' n = ' num2str(size(fixation_locked_firing{unit},1))];
        if multiunit(unit)
            subtitle(['Fixation-Locked Multiunit ' task_file(1:8) ' ' unit_stats{1,unit} n_str]);
        else
            subtitle(['Fixation-Locked ' task_file(1:8) ' ' unit_stats{1,unit} n_str]);
        end
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Sequence_All_Fixation_Event_Locked']);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---Fixation Aligned Activity by Item---%
        yls = NaN(8,2);
        figure
        
        ind = 1;
        for s = 1:2
            for c = 1:4
                
                fix_algined = fixation_locked_firing{unit}((sequence_nums{unit} == s) & (item_nums{unit} == c),:);
                
                
                %---Firing rate curves---%
                subplot(4,4,c+8*(s-1))
                dofill(t,fix_algined,'black',1,smval);
                xlabel('Time from Fixation Start (ms)')
                ylabel('Firing Rate (Hz)')
                xlim([-twin twin])
                title_str = [];
                if sequence_info.rate_prctile(c,s,unit) > 95
                    title_str = [title_str 'fix_{bit} = ' num2str(sequence_info.rate_prctile(c,s,unit),3) '%%'];
                end
                if sequence_info.temporalstability_prctile(c,s,unit) > 95
                    title_str = [title_str ', \\rho_{1/2} = ' ...
                        num2str(sequence_info.temporalstability(c,s,unit),2) ...
                        ' (' num2str(sequence_info.temporalstability_prctile(c,s,unit),3) '%%)'];
                end
                if ~isempty(title_str)
                    title(sprintf(title_str))
                end
                yls(ind,:) = ylim;
                ind = ind+1;
                
                %---rasters---%
                subplot(4,4,c+8*(s-1)+4)
                [trial,time] = find(fix_algined == 1);
                plot(time-twin,(trial),'.k')
                ylabel('Trial #')
                xlabel('Time from fixation Start')
                xlim([-twin twin])
                if ~isempty(trial)
                    ylim([0 max(trial)+1])
                end
                box off
            end
        end
        
        ymin = min(yls(:,1));
        ymin(ymin < 0) = 0;
        ymax = max(yls(:,2));
        
        for s = 1:2
            for c = 1:4
                subplot(4,4,c+8*(s-1))
                ylim([ymin ymax])
                hold on
                plot([0 0],[ymin ymax],'k--')
                hold off
                
            end
        end
        if multiunit(unit)
            subtitle([task_file(1:8) ' ' unit_stats{1,unit}]);
        else
            subtitle([task_file(1:8) ' ' unit_stats{1,unit}]);
        end
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Sequence_Fixation_Locked_by_item']);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---Item/Event Aligned Activity by Item---%
        yls = NaN(8,2);
        figure
        
        ind = 1;
        for s = 1:2
            for c = 1:4
                
                item_aligned = event_locked_firing{unit}((sequence_nums{unit} == s) & (item_nums{unit} == c),:);
                rts = reaction_time{unit}((sequence_nums{unit} == s) & (item_nums{unit} == c));
                
                %---Firing rate curves---%
                subplot(4,4,c+8*(s-1))
                hold on
                dofill(t,item_aligned,'blue',1,smval);
                hold off
                xlabel('Time from Item On (ms)')
                ylabel('Firing Rate (Hz)')
                xlim([-twin twin])
                title_str = [];
                if item_info.rate_prctile(c,s,unit) > 95
                    title_str = [title_str 'Event_{bit} = ' num2str(item_info.rate_prctile(c,s,unit),3) '%%'];
                end
                if item_info.temporalstability_prctile(c,s,unit) > 95
                    title_str = [title_str ', \\rho_{1/2} = ' ...
                        num2str(item_info.temporalstability(c,s,unit),2) ...
                        ' (' num2str(item_info.temporalstability_prctile(c,s,unit),3) '%%)'];
                end
                if ~isempty(title_str)
                    title(sprintf(title_str))
                end
                yls(ind,:) = ylim;
                ind = ind+1;
                
                
                %---rasters---%
                subplot(4,4,c+8*(s-1)+4)
                [~,si] = sort(rts);
                [trial,time] = find(item_aligned(si,:) == 1);
                plot(time-twin,(trial),'.k')
                ylabel('Ranked by RT')
                xlabel('Time from fixation Start')
                xlim([-twin twin])
                if ~isempty(trial)
                    ylim([0 max(trial)+1])
                end
                box off
            end
        end
        
        ymin = min(yls(:,1));
        ymin(ymin < 0) = 0;
        ymax = max(yls(:,2));
        
        for s = 1:2
            for c = 1:4
                subplot(4,4,c+8*(s-1))
                ylim([ymin ymax])
                hold on
                plot([0 0],[ymin ymax],'k--')
                hold off
            end
        end
        if multiunit(unit)
            subtitle([task_file(1:8) ' ' unit_stats{1,unit}]);
        else
            subtitle([task_file(1:8) ' ' unit_stats{1,unit}]);
        end
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Sequence_Item_Locked_by_item']);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---Item/Event Off Aligned Activity by Item---%
        yls = NaN(8,2);
        figure
        
        ind = 1;
        for s = 1:2
            for c = 1:4
                
                item_aligned = event_off_locked_firing{unit}((sequence_nums{unit} == s) & (item_nums{unit} == c),:);
                rts = reaction_time{unit}((sequence_nums{unit} == s) & (item_nums{unit} == c));
                
                %---Firing rate curves---%
                subplot(4,4,c+8*(s-1))
                hold on
                dofill(t,item_aligned,'blue',1,smval);
                hold off
                xlabel('Time from Item On (ms)')
                ylabel('Firing Rate (Hz)')
                xlim([-twin twin])
                title_str = [];
                if item_off_info.rate_prctile(c,s,unit) > 95
                    title_str = [title_str 'Event_{bit} = ' num2str(item_off_info.rate_prctile(c,s,unit),3) '%%'];
                end
                if item_off_info.temporalstability_prctile(c,s,unit) > 95
                    title_str = [title_str ', \\rho_{1/2} = ' ...
                        num2str(item_off_info.temporalstability(c,s,unit),2) ...
                        ' (' num2str(item_off_info.temporalstability_prctile(c,s,unit),3) '%%)'];
                end
                if ~isempty(title_str)
                    title(sprintf(title_str))
                end
                yls(ind,:) = ylim;
                ind = ind+1;
                
                
                %---rasters---%
                subplot(4,4,c+8*(s-1)+4)
                [~,si] = sort(rts);
                [trial,time] = find(item_aligned(si,:) == 1);
                plot(time-twin,(trial),'.k')
                ylabel('Ranked by RT')
                xlabel('Time from fixation Start')
                xlim([-twin twin])
                if ~isempty(trial)
                    ylim([0 max(trial)+1])
                end
                box off
            end
        end
        
        ymin = min(yls(:,1));
        ymin(ymin < 0) = 0;
        ymax = max(yls(:,2));
        
        for s = 1:2
            for c = 1:4
                subplot(4,4,c+8*(s-1))
                ylim([ymin ymax])
                hold on
                plot([0 0],[ymin ymax],'k--')
                hold off
            end
        end
        if multiunit(unit)
            subtitle([task_file(1:8) ' ' unit_stats{1,unit}]);
        else
            subtitle([task_file(1:8) ' ' unit_stats{1,unit}]);
        end
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Sequence_Item_Locked_by_item_Off']);
    end
end

save([data_dir task_file(1:8) '-Fixation_locked_Sequence_results.mat'],...
    'twin','smval','fixation_locked_firing','fixation_info','sequence_info',...
    'item_nums','sequence_nums','event_locked_firing','reaction_time',....
    'all_items_info','sequence_info','unit_names','overlap','trial_numbers',...
    'event_off_locked_firing','item_off_info','item_info')