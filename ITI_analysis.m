function ITI_analysis(data_dir,figure_dir,session_data)
% written by Seth Konig 4/6/16. Updated September 5, 2016
% Modified from time_locked_analysisV2 to
% focus more specifically on the ITI period since many neurons seem to like
% this period. Analysis jsut for ListSQ since time_locked_analysisV2 is
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


figure_dir = [figure_dir 'ITI analysis\'];
task = 'ListSQ';
min_blks = 2;
twin = 200;% how much time to take prior to ITT start and end
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
trial_start_code = 15; %ITI code
ITI_dur = 1000;
trial_end_code = 20;
eye_off_code = 101;
imgon_code = 23;
imgoff_code = 24;
reward_code = 3;
smval = 60 ;%gaussian 1/2 width for smoothing
numshuffs = 1000;
p_thresh = 0.01;

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

%determine which trials units are stable for
valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Align Spikes to ITI---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Aligning spike times to trial events')

%get important task specific information
[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[~,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,imgon_code);

%preallocate space and parallel structure of cfg
time_locked_firing = cell(1,num_units);
previous_trial_type = cell(1,num_units);%what kind of trial happened previous the trial
upcomming_trial_type = cell(1,num_units);%what kind of trial happened current this trial
for unit = 1:num_units
    time_locked_firing{unit} = NaN(num_trials,ITI_dur+twin);
    previous_trial_type{unit} = zeros(1,num_trials);
    upcomming_trial_type{unit} = zeros(1,num_trials);
end
%trial type 1: sequence 1
%trial type 2: sequence 2
%trial type 3: novel image
%trial type 4: repeat image

trial_type_previous = 0;
for t = 2:num_trials
    if any(cfg.trl(t).allval == imgon_code); %in which image was displayed or 1st item in sequence was displayed

        %---this trial---%
        this_trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
        
        %---This Trial---%
        if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(end)%last trial was sequence
            if ~isempty(find(cfg.trl(t).allval == reward_code,1))%%last trial was successful
                trial_type_current = find(sequence_items == itmlist(cfg.trl(t).cnd-1000)); %store for next trial
            else
                trial_type_current = -find(sequence_items == itmlist(cfg.trl(t).cnd-1000)); %store for next trial
            end
        elseif any(cfg.trl(t).allval == imgon_code)%succesful image trial
            if  novel_vs_repeat(cfg.trl(t).cnd == img_cnd) == 1
                trial_type_current = 3;%store for next trial
            else
                trial_type_current = 4;%store for next trial
            end
        else %not a succsefful image trial
            if  novel_vs_repeat(cfg.trl(t).cnd == img_cnd) == 1
                trial_type_current = -3;%store for next trial
            else
                trial_type_current = -4;%store for next trial
            end
        end
        
        %---Last Trial---%
        last_trial_start = cfg.trl(t-1).alltim(cfg.trl(t-1).allval == trial_start_code);
        if itmlist(cfg.trl(t-1).cnd-1000) <= sequence_items(end)%last trial was sequence
            if ~isempty(find(cfg.trl(t-1).allval == reward_code,1))%%last trial was successful
                last_trial_end = cfg.trl(t-1).alltim(cfg.trl(t-1).allval == reward_code);
                last_trial_end = last_trial_end(end);
                trial_type_previous = find(sequence_items == itmlist(cfg.trl(t-1).cnd-1000)); %store for next trial
            else
                last_trial_end = cfg.trl(t-1).alltim(cfg.trl(t-1).allval == eye_off_code);
                trial_type_previous = -find(sequence_items == itmlist(cfg.trl(t-1).cnd-1000)); %store for next trial
            end
        elseif any(cfg.trl(t-1).allval == imgon_code)%succesful image trial
            last_trial_end = cfg.trl(t-1).alltim(cfg.trl(t-1).allval == imgoff_code);
            if  novel_vs_repeat(cfg.trl(t-1).cnd == img_cnd) == 1
                trial_type_previous = 3;%store for next trial
            else
                trial_type_previous = 4;%store for next trial
            end
        else %not a succsefful image trial
            last_trial_end = cfg.trl(t-1).alltim(cfg.trl(t-1).allval == eye_off_code);
            if  novel_vs_repeat(cfg.trl(t-1).cnd == img_cnd) == 1
                trial_type_previous = -3;%store for next trial
            else
                trial_type_previous = -4;%store for next trial
            end
        end

        event_start = last_trial_end-last_trial_start; %so same time axis as spikes
        for unit = 1:num_units
            if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only put data in for valid trials
                
                %spikes from last trial
                spikes = find(data(unit).values{t-1} == 1);
                event_spikes = spikes(spikes > event_start-twin & ...
                    spikes <= event_start+ITI_dur)-event_start+twin;
                
                %get spikes from this trial
                temp = zeros(1,ITI_dur+twin);
                if this_trial_start-(last_trial_end+ITI_dur) < 0 %trial "ends" before window
                    leftover = (last_trial_end+ITI_dur)-this_trial_start;
                    if t ~= num_trials
                        spks = find(data(unit).values{t}(1:leftover) == 1)+...
                            (this_trial_start-last_trial_end)+twin;
                        event_spikes = [event_spikes spks];
                    else
                        temp(end-leftover+1:end) = NaN;
                    end
                end
                temp(event_spikes) = 1;
                
                time_locked_firing{unit}(t,:) = temp;
                previous_trial_type{unit}(t) = trial_type_previous;
                upcomming_trial_type{unit}(t) = trial_type_current;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Perform Temporal analysis and signifance testing---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Determining if spikes are locked to trial events')
Fs = data(1).fsample; %should be 1000

info_type = 'temporal';
% perform mutual information theory analysis to determine
% if neurons encode info about a particular time within the ITI
ITI_temporal_info.rate = NaN(1,num_units);
ITI_temporal_info.shuffled_rate = cell(1,num_units);
ITI_temporal_info.rate_prctile =  NaN(1,num_units);
ITI_temporal_info.temporalstability = NaN(2,num_units);
ITI_temporal_info.shuffled_temporalstability = cell(1,num_units);
ITI_temporal_info.temporalstability_prctile =  NaN(2,num_units);
for unit = 1:num_units
    if  nansum(nansum(time_locked_firing{unit})) > 0
        
        %for ITI following successful previous trial
           [observed_info_rate,shuffled_info_rate]...
            = estimated_mutual_information(time_locked_firing{unit}(previous_trial_type{unit} > 0,:),numshuffs,info_type,smval,Fs);
        
        ITI_temporal_info.rate(unit) = observed_info_rate.skaggs;
        ITI_temporal_info.shuffled_rate{unit} = shuffled_info_rate.skaggs;
        ITI_temporal_info.rate_prctile(unit) = 100*sum(ITI_temporal_info.rate(unit)...
            > ITI_temporal_info.shuffled_rate{unit})/numshuffs;
        ITI_temporal_info.temporalstability(:,unit) = observed_info_rate.temporalstability;
        ITI_temporal_info.shuffled_temporalstability{unit} = shuffled_info_rate.temporalstability;
        ITI_temporal_info.temporalstability_prctile(1,unit) = 100*sum(...
            observed_info_rate.temporalstability(1) > ...
            shuffled_info_rate.temporalstability(1,:))/numshuffs;
        ITI_temporal_info.temporalstability_prctile(2,unit) = 100*sum(...
            observed_info_rate.temporalstability(2) > ...
            shuffled_info_rate.temporalstability(2,:))/numshuffs;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot and Save Figures of Results---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = -twin+1:ITI_dur;
unit_names = unit_stats(1,:);
for unit = 1:num_units
    if nansum(nansum(time_locked_firing{unit})) > 0
        
        ylims = NaN(3,2);
        
        figure
        
        %firing rate locked to ITI start
        subplot(2,2,1)
        dofill(t,time_locked_firing{unit}(previous_trial_type{unit} > 0,:),'k',1,smval);
        xlabel('Time from  ITI start (ms)')
        ylabel('Firing Rate (Hz)')
        title_str = ['n = ' num2str(size(time_locked_firing{unit},1))];
        if ITI_temporal_info.rate_prctile(unit) > 95 || ITI_temporal_info.temporalstability_prctile(unit) > 95
            title_str = [title_str ', bit = ' num2str(ITI_temporal_info.rate_prctile(unit),3) ...
                '%, \rho_{1/2} = ' num2str(ITI_temporal_info.temporalstability(1,unit),2) ...
                ' ' num2str( ITI_temporal_info.temporalstability_prctile(1,unit),3) '%'];
        end
        xlim([-twin ITI_dur])
        title(title_str)
        ylims(1,:) = ylim;
        
        %raster plot locked to ITI start
        subplot(2,2,3)
        [trial,time] = find(time_locked_firing{unit}(previous_trial_type{unit} > 0,:) == 1);
        if ~isempty(trial)
            plot(time-twin,trial,'.k')
            xlim([-twin ITI_dur])
            ylim([0 max(trial)+1]);
        end
        ylabel('Trial #')
        xlabel('Time from ITI start (ms)')
        box off
        
        
        %previous trial incluence of ITI firing rate
        subplot(2,2,2)
        hold on
        dofill(t,time_locked_firing{unit}(previous_trial_type{unit} == 1,:),'r',1,smval);
        dofill(t,time_locked_firing{unit}(previous_trial_type{unit} == 2,:),'g',1,smval);
        dofill(t,time_locked_firing{unit}(previous_trial_type{unit} == 3,:),'b',1,smval);
        dofill(t,time_locked_firing{unit}(previous_trial_type{unit} == 4,:),'m',1,smval);
        dofill(t,time_locked_firing{unit},'k',1,smval);
        hold off
        xlim([-twin twin+ITI_dur])
        xlabel('Time from  ITI start (ms)')
        ylabel('Firing Rate (Hz)')
        xlim([-twin ITI_dur])
        title('Divided By Previous Trial Type')
        
        ylims(2,:) = ylim;
        
        %current trial incluence of ITI firing rate
        subplot(2,2,4)
        hold on
        dofill(t,time_locked_firing{unit}(upcomming_trial_type{unit} == 1,:),'r',1,smval);
        dofill(t,time_locked_firing{unit}(upcomming_trial_type{unit} == 2,:),'g',1,smval);
        dofill(t,time_locked_firing{unit}(upcomming_trial_type{unit} == 3,:),'b',1,smval);
        dofill(t,time_locked_firing{unit}(upcomming_trial_type{unit} == 4,:),'m',1,smval);
        dofill(t,time_locked_firing{unit},'k',1,smval);
        hold off
        xlim([-twin ITI_dur])
        xlabel('Time from  ITI start (ms)')
        ylabel('Firing Rate (Hz)')
        title('Divided By Upcomming Trial Type')
        legend({'Sequence 1','Sequence 2','Novel Images','Repeat Images','All'},...
            'Location','NorthEast')
        
        ylims(3,:) = ylim;
        
        %set y-limits to be the same
        ymin = min(ylims(:,1));
        ymin(ymin < 0) = 0;
        ymax = max(ylims(:,2));
        
        for sb = [1 2 4]
            subplot(2,2,sb)
            hold on
            plot([0 0],[ymin ymax],'k--')
            hold off
            ylim([ymin ymax])
        end
        
        if multiunit(unit)
            subtitle(['MultiUnit' task_file(1:8) ' ' unit_names{unit}]);
        else
            subtitle([ task_file(1:8) ' ' unit_names{unit}]);
        end
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_ITI_analysis']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Finally save all the data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([data_dir task_file(1:8) '-ListSQ-ITI_results.mat'],'time_locked_firing',...
    'previous_trial_type','upcomming_trial_type','smval','ITI_temporal_info',...
    'unit_names','twin','ITI_dur');
disp(['ITI Locked Data Analyis for ' task_file(1:8) ' saved']);
end