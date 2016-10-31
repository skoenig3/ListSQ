function Reward_analysis(data_dir,figure_dir,session_data)
% written by Seth Konig 9/29/16. Revised from ITI_analysis
% Modified from time_locked_analysisV2 to
% focus more specifically on the Reward period since many neurons seem to like
% this period. Analysis just for ListSQ since time_locked_analysisV2 is
% sufficient for other tasks.
%
% Inputs:
%   1) data_dir: direction where preprocessed data is located
%   2) figure_dir: location of where to put save figures
%   3) session_data: data containing relevant session data
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-time_locked_results'


figure_dir = [figure_dir 'Reward analysis\'];
task = 'ListSQ';
min_blks = 2;
twin = 150;% how much time to take prior to ITT start and end
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
reward_pulse_dur = 115;%ms for each pulse sound, audio file confirms sound last ~115 ms
IRI = 166;%ms inter-reward-interval
reward_dur = 945; %116 ms per pulse *6;
%takes into account the fact last pulse occured at 830 ms after first and
%115 pulse duration
trial_start_code = 15; %ITI code
imgon_code = 23;
reward_code = 3;
smval = 60 ;%gaussian 1/2 width for smoothing, 30 ms std
numshuffs = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials','hdr');
%grab unit data
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

if num_units == 0
    return
end
num_trials = length(cfg.trl);

valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks);
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Align Spikes to Reward Period Start---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Aligning spike times to trial events')

%get important task specific information
[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[~,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,imgon_code);

%preallocate space and parallel structure of cfg
time_locked_firing = cell(1,num_units);
trial_type = cell(1,num_units);%which sequence preceeded reward
for unit = 1:num_units
    time_locked_firing{unit} = NaN(num_trials,reward_dur+2*twin);
    trial_type{unit} = NaN(1,num_trials);
end
%trial type 1: sequence 1
%trial type 2: sequence 2

for t = 1:num_trials
    if any(cfg.trl(t).allval == imgon_code); %in which image was displayed or 1st item in sequence was displayed
        if itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) || isempty(find(cfg.trl(t).allval == reward_code,1))% sequence trials and only want rewarded ones
            continue % go to the next trial
        end
        
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
        trial_end = cfg.trl(t).alltim(end)-trial_start;
        rewards = cfg.trl(t).alltim(cfg.trl(t).allval == reward_code)-trial_start;
        reward_start = rewards(1);
        reward_end = rewards(end);
        
%         if trial_end-reward_end < twin
%             disp('what?')
%         end
        
        for unit = 1:num_units
            if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only put data in for valid trials
                
                %spikes from current trial
                spikes = find(data(unit).values{t} == 1);
                event_spikes = spikes(spikes > reward_start-twin & ...
                    spikes <= reward_start+reward_dur+twin)-reward_start+twin;
                
                temp = zeros(1,reward_dur+2*twin);
                if trial_end-(reward_start+reward_dur+twin) < 0 %trial "ends" before window
                    leftover = (reward_start+reward_dur+twin)-trial_end;
                    if t ~= num_trials
                        spks = find(data(unit).values{t+1}(1:leftover) == 1)+(trial_end-reward_start);
                        event_spikes = [event_spikes spks];
                    else
                        temp(end-leftover+1:end) = NaN;
                    end
                end
                temp(event_spikes) = 1;
                
                time_locked_firing{unit}(t,:) = temp;
                trial_type{unit}(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
            end
        end
    end
end

%remove excess NaNs associated with error trials
time_locked_firing = laundry(time_locked_firing);
trial_type = laundry(trial_type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Perform Temporal analysis and signifance testing---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Determining if spikes are locked to trial events')
Fs = data(1).fsample; %should be 1000

info_type = 'temporal';
% perform mutual information theory analysis to determine
% if neurons encode info about a particular time within the ITI
temporal_info.rate = NaN(1,num_units);
temporal_info.shuffled_rate = cell(1,num_units);
temporal_info.rate_prctile =  NaN(1,num_units);
temporal_info.temporalstability = NaN(2,num_units);
temporal_info.shuffled_temporalstability = cell(1,num_units);
temporal_info.temporalstability_prctile =  NaN(2,num_units);
for unit = 1:num_units
    if ~isempty(time_locked_firing{unit})
        
        %for previous trial type
        [observed_info_rate,shuffled_info_rate]...
            = estimated_mutual_information(time_locked_firing{unit},numshuffs,info_type,smval,Fs);
        
        temporal_info.rate(unit) = observed_info_rate.skaggs;
        temporal_info.shuffled_rate{unit} = shuffled_info_rate.skaggs;
        temporal_info.rate_prctile(unit) = 100*sum(temporal_info.rate(unit)...
            > temporal_info.shuffled_rate{unit})/numshuffs;
        temporal_info.temporalstability(:,unit) = observed_info_rate.temporalstability;
        temporal_info.shuffled_temporalstability{unit} = shuffled_info_rate.temporalstability;
        temporal_info.temporalstability_prctile(1,unit) = 100*sum(...
            observed_info_rate.temporalstability(1) > ...
            shuffled_info_rate.temporalstability(1,:))/numshuffs;
        temporal_info.temporalstability_prctile(2,unit) = 100*sum(...
            observed_info_rate.temporalstability(2) > ...
            shuffled_info_rate.temporalstability(2,:))/numshuffs;
       
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot and Save Figures of Results---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = -twin+1:reward_dur+twin;
unit_names = unit_stats(1,:);
for unit = 1:num_units
    if ~isempty(time_locked_firing{unit})
        
        
        figure
        
        %firing rate locked to Reward start
        subplot(2,2,1)
        dofill(t,time_locked_firing{unit},'k',1,smval);
        hold on
        dofill(t,time_locked_firing{unit}(trial_type{unit} == 1,:),'r',1,smval);
        dofill(t,time_locked_firing{unit}(trial_type{unit} == 2,:),'b',1,smval);
        ylims = ylim;
        if ylims(1) < 0
            ylims(1) = 0;
            ylim([0 ylims(2)])
        end
        for i = 1:6
            plot([(i-1)*IRI (i-1)*IRI],[ylims(1) ylims(2)],'k--')
        end
        xlabel('Time from  Reward start (ms)')
        ylabel('Firing Rate (Hz)')
        title_str = ['n = ' num2str(size(time_locked_firing{unit},1))];
        if temporal_info.rate_prctile(unit) > 95 || temporal_info.temporalstability_prctile(unit) > 95
            title_str = [title_str ', bit = ' num2str(temporal_info.rate_prctile(unit),3) ...
                '%, \rho_{1/2} = ' num2str(temporal_info.temporalstability(1,unit),2) ...
                ' ' num2str( temporal_info.temporalstability_prctile(1,unit),3) '%'];
        end
        xlim([-twin twin+reward_dur])
        title(title_str)
        legend('All','Seq #1','Seq #2')
        
        %raster plot locked to Reward start
        subplot(2,2,3)
        [trial,time] = find(time_locked_firing{unit} == 1);
        hold on
        for i = 1:6
            plot([(i-1)*IRI (i-1)*IRI],[0 max(trial)],'k--')
        end
        if ~isempty(trial)
            plot(time-twin,trial,'.k')
            xlim([-twin reward_dur+twin])
            ylim([0 max(trial)+1]);
        end
        ylabel('Trial #')
        xlabel('Time from  Reward start (ms)')

        
        %---Plot Rasters for Sequence 1 and Sequence 2---%
        seq1 = time_locked_firing{unit}(trial_type{unit} == 1,:);
        seq2 = time_locked_firing{unit}(trial_type{unit} == 2,:);
        
        [trial1,time1] = find(seq1 == 1);
        [trial2,time2] = find(seq2 == 1);
        if ~isempty(trial1)
            trial2 = max(trial1)+trial2+1;
        else
        end
        
        subplot(2,2,4)
        hold on
        plot(time1-twin,trial1,'.r')
        plot(time2-twin,trial2,'.b')
        hold off
        ylabel('Trial #')
        xlabel('Time from  Reward start (ms)')
        ylim([0 max(trial2)+1])
        
        if multiunit(unit)
            subtitle(['MultiUnit ' unit_names{unit}]);
        else
            subtitle(unit_names{unit});
        end
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_Reward_analysis']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Finally save all the data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([data_dir task_file(1:8) '-ListSQ-Reward_results.mat'],'time_locked_firing',...
    'trial_type','smval','temporal_info','unit_names','twin','reward_dur');
disp(['Reward Locked Data Analyis for ' task_file(1:8) ' saved']);
end