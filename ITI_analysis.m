function ITI_analysis(data_dir,figure_dir,session_data,task)
% written by Seth Konig 4/6/16. Modified from time_locked_analysisV2 to
% focus more specifically on the ITI period since many neurons seem to like
% this period. Analysis jsut for ListSQ since time_locked_analysisV2 is
% sufficient for other tasks.
%
% Inputs:
%   1) data_dir: direction where preprocessed data is located
%   2) figure_dir: location of where to put save figures
%   3) session_data: data containing relevant session data
%   4) task: what task was this data come from with i.e. 'cvtnew','ListSQ'
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-time_locked_results'

twin = 750;% how much time to take previous and current saccade.
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
trial_start_code = 15; %ITI code
ITI_dur = 1000;%ms
trial_end_code = 20;
imgon_code = 23;
imgoff_code = 24;
reward_code = 3;
smval = 60 ;%gaussian 1/2 width for smoothing
numshuffs = 500;

%ANVOA analysis sliding window parameters
ANOVA_window = 200; %ms
ANOVA_step = 50;%ms 

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import & reformat data so that spikes are locked to ITI---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_trials = length(cfg.trl);

%NaNs are for start and end trials otherwise cut
valid_trials(1,isnan(valid_trials(1,:))) = 1;
valid_trials(2,isnan(valid_trials(2,:))) = num_trials;

disp('Aligning spike times to trial events')

%get important task specific information
[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[~,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,imgon_code);

%preallocate space and parallel structure of cfg
time_locked_firing = cell(1,num_units);
previous_trial_type = cell(1,num_units);%what kind of trial happened previous the trial
current_trial_type = cell(1,num_units);%what kind of trial happened current this trial
for unit = 1:num_units
    time_locked_firing{unit} = NaN(num_trials,ITI_dur+2*twin);
    previous_trial_type{unit} = NaN(1,num_trials);
    current_trial_type{unit} = NaN(1,num_trials);
end
%trial type 1: sequence 1
%trial type 2: sequence 2
%trial type 3: novel image
%trial type 4: repeat image

trial_type_previous = NaN;
for t = 1:num_trials
    if any(cfg.trl(t).allval == imgon_code); %in which image was displayed or 1st item in sequence was displayed
        if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(end) && isempty(find(cfg.trl(t).allval == reward_code,1))% sequence trials and only want rewarded ones
            continue % go to the next trial
        end
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
  
        for unit = 1:num_units
            if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only put data in for valid trials
                
                %spikes from current trial
                spikes = find(data(unit).values{t}(1:ITI_dur+twin))+twin;
                temp = zeros(1,ITI_dur+2*twin);
                temp(spikes) = 1;
                
                %spikes from previous trial's end
                if t > 1
                    spikes = find(data(unit).values{t-1}(end-twin:end));
                    temp(spikes) = 1;
                else
                    temp(1:twin) = NaN; 
                end
                time_locked_firing{unit}(t,:) = temp; 
                
                previous_trial_type{unit}(t) = trial_type_previous; 
                if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(end) %sequence trial
                   current_trial_type{unit}(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
                   trial_type_previous = find(sequence_items == itmlist(cfg.trl(t).cnd-1000)); %store for next trial
                else %image trial
                    if  novel_vs_repeat(cfg.trl(t).cnd == img_cnd) == 1
                          current_trial_type{unit}(t) = 3; %novel image
                          trial_type_previous = 3;%store for next trial
                    else
                         current_trial_type{unit}(t) = 4; %repeat image
                         trial_type_previous = 4;%store for next trial
                    end
                end
            end
        end
    end
end

%remove excess NaNs associated with error trials
time_locked_firing = laundry(time_locked_firing);
current_trial_type = laundry(current_trial_type);
previous_trial_type = laundry(previous_trial_type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Perform Temporal analysis and signifance testing---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Determining if spikes are locked to trial events')
Fs = data(1).fsample; %should be 1000

info_type = 'temporal';
% perform mutual information theory analysis to determine
% if neurons encode info about a particular time within the ITI
temporal_info.previous_rate = NaN(4,num_units);
temporal_info.previous_shuffled_rate = cell(4,num_units);
temporal_info.previous_rate_prctile =  NaN(4,num_units);
temporal_info.current_rate = NaN(4,num_units);
temporal_info.current_shuffled_rate = cell(4,num_units);
temporal_info.current_rate_prctile =  NaN(4,num_units);
temporal_info.rate = NaN(1,num_units);
temporal_info.shuffled_rate = cell(1,num_units);
temporal_info.rate_prctile = NaN(1,num_units);;
for unit = 1:num_units
    for type = 1:4
        if ~isempty(time_locked_firing{unit})
            
            %for previous trial type
            temp_time_lock_previous = time_locked_firing{unit}(previous_trial_type{unit} == type,:);
            [temporal_info.previous_rate(type,unit),temporal_info.previous_shuffled_info_rate{type,unit}]...
                = estimated_mutual_information(temp_time_lock_previous,numshuffs,info_type,smval,Fs);
            temporal_info.previous_rate_prctile(type,unit) = ...
                100*sum(temporal_info.previous_shuffled_info_rate{type,unit}...
                < temporal_info.previous_rate(type,unit))/numshuffs;
            
            %for current trial type
            temp_time_lock_current =  time_locked_firing{unit}(current_trial_type{unit} == type,:);
            [temporal_info.current_rate(type,unit),temporal_info.current_shuffled_info_rate{type,unit}]...
                = estimated_mutual_information(temp_time_lock_current,numshuffs,info_type,smval,Fs);
            temporal_info.current_rate_prctile(type,unit) = ...
                100*sum(temporal_info.current_shuffled_info_rate{type,unit}...
                < temporal_info.current_rate(type,unit))/numshuffs;
        end
    end
    [temporal_info.rate(unit),temporal_info.shuffled_info_rate{unit}]...
        = estimated_mutual_information(time_locked_firing{unit},numshuffs,info_type,smval,Fs);
    
    temporal_info.rate_prctile(unit) = 100*sum(temporal_info.shuffled_info_rate{unit} < ...
        temporal_info.rate(unit))/numshuffs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---ANOVA anlysis to Determine if Firing Rate is Modulated by TrialType---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_steps = (ITI_dur+2*twin)/ANOVA_step-ANOVA_window/ANOVA_step+1;%

pvals_previous = NaN(num_units,num_steps);
pvals_current = NaN(num_units,num_steps);
for unit = 1:num_units
    for step =  1:num_steps
        time_window = (step-1)*ANOVA_step+1:((step-1)*ANOVA_step+ANOVA_window);%time window to grab data
        
        %get spikes/trial for each of the conditions varied over previous
        %trial type
        temp_firing1 = sum(time_locked_firing{unit}(previous_trial_type{unit} == 1,time_window),2);
        temp_firing2 = sum(time_locked_firing{unit}(previous_trial_type{unit} == 2,time_window),2);
        temp_firing3 = sum(time_locked_firing{unit}(previous_trial_type{unit} == 3,time_window),2);
        temp_firing4 = sum(time_locked_firing{unit}(previous_trial_type{unit} == 4,time_window),2);
        
        pvals_previous(unit,step) =anova1([temp_firing1;temp_firing2;temp_firing3;temp_firing4],...
            [ones(length(temp_firing1),1);2*ones(length(temp_firing2),1);...
            3*ones(length(temp_firing3),1);4*ones(length(temp_firing4),1)],'off');
        
        %get spikes/trial for each of the conditon varied over current
        %trial type
        temp_firing1 = sum(time_locked_firing{unit}(current_trial_type{unit} == 1,time_window),2);
        temp_firing2 = sum(time_locked_firing{unit}(current_trial_type{unit} == 2,time_window),2);
        temp_firing3 = sum(time_locked_firing{unit}(current_trial_type{unit} == 3,time_window),2);
        temp_firing4 = sum(time_locked_firing{unit}(current_trial_type{unit} == 4,time_window),2);
        
         pvals_current(unit,step) =anova1([temp_firing1;temp_firing2;temp_firing3;temp_firing4],...
            [ones(length(temp_firing1),1);2*ones(length(temp_firing2),1);...
            3*ones(length(temp_firing3),1);4*ones(length(temp_firing4),1)],'off');
   end    
end

ANOVA_stats.pvals_previous = pvals_previous;
ANOVA_stats.pvals_current = pvals_current;

%use a simple bonferroni correction
significant_times_previous = zeros(num_steps,num_units);
significant_times_previous(pvals_previous < 0.05/num_steps) = 1;
significant_times_current =zeros(num_steps,num_units);
significant_times_current(pvals_current < 0.05/num_steps) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot and Save Figures of Results---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = -twin+1:ITI_dur+twin;
unit_names = unit_stats(1,:);
for unit = 1:num_units
    %%
    ylims = NaN(2,2);
    
    figure
    
    %previous trial incluence of ITI firing rate
    subplot(1,2,1)
    hold on
    dofill(t,time_locked_firing{unit}(previous_trial_type{unit} == 1,:),'r',1,smval);
    dofill(t,time_locked_firing{unit}(previous_trial_type{unit} == 2,:),'g',1,smval);
    dofill(t,time_locked_firing{unit}(previous_trial_type{unit} == 3,:),'b',1,smval);
    dofill(t,time_locked_firing{unit}(previous_trial_type{unit} == 4,:),'m',1,smval);
    dofill(t,time_locked_firing{unit},'k',1,smval);
    hold off
    xlabel('Time from  ITI start (ms)')
    ylabel('Firing Rate (Hz)')
    
    title_str =['Divided By Previous Trial Type'];
    if temporal_info.rate_prctile(unit) > 90;
        title_str = [title_str '\n ' num2str(temporal_info.rate(unit),3) ' Bits ' num2str(temporal_info.rate_prctile(unit),3) '%%'];
    end
    title(sprintf(title_str))
    
    ylims(1,:) = ylim;
    
    %current trial incluence of ITI firing rate
    subplot(1,2,2)
    hold on
    dofill(t,time_locked_firing{unit}(current_trial_type{unit} == 1,:),'r',1,smval);
    dofill(t,time_locked_firing{unit}(current_trial_type{unit} == 2,:),'g',1,smval);
    dofill(t,time_locked_firing{unit}(current_trial_type{unit} == 3,:),'b',1,smval);
    dofill(t,time_locked_firing{unit}(current_trial_type{unit} == 4,:),'m',1,smval);
    dofill(t,time_locked_firing{unit},'k',1,smval);
    hold off
    xlabel('Time from  ITI start (ms)')
    ylabel('Firing Rate (Hz)')
    title('Divided By Current Trial Type')
    legend({'Sequence 1','Sequence 2','Novel Images','Repeat Images','All'},...
        'Location','NorthEast')
        
    ylims(2,:) = ylim;

    %set y-limits to be the same
    ymin = min(ylims(:,1));
    ymin(ymin < 0) = 0;
    ymax = max(ylims(:,2));
      
    %---shade in times in which firing rate significantly varies by trial type---%
    %by prvious trial type
    subplot(1,2,1)
    ylim([ymin ymax])
    hold on
    %not efficient will want to rewrite and may depend on multiple
    %comparisons corrections method
    for step = 1:num_steps
        if significant_times_previous(step,unit) == 1
            time_window = [(step-1)*ANOVA_step+1 ((step-1)*ANOVA_step+ANOVA_window)]-twin;%time window to grab data
            h = fill([time_window(1) time_window(2) time_window(2) time_window(1) time_window(1)],...
                [ymin ymin ymax ymax 0],'k');
            uistack(h,'bottom')
            set(h,'facealpha',.25,'EdgeColor','None')
        end
    end
    hold off
    
    subplot(1,2,2)
    ylim([ymin ymax])
    hold on
    %not efficient will want to rewrite and may depend on multiple
    %comparisons corrections method
    for step = 1:num_steps
        if significant_times_current(step,unit) == 1
            time_window = [(step-1)*ANOVA_step+1 ((step-1)*ANOVA_step+ANOVA_window)]-twin;%time window to grab data
            h = fill([time_window(1) time_window(2) time_window(2) time_window(1) time_window(1)],...
                [ymin ymin ymax ymax 0],'k');
            uistack(h,'bottom')
            set(h,'facealpha',.25,'EdgeColor','None')
        end
    end
    hold off
    %%
    
    if multiunit(unit)
        subtitle(['Saccade-Locked Multiunit ' unit_names{unit} ' n =' num2str(size(time_locked_firing{unit},1))]);
    else
        subtitle(['Saccade-Locked ' unit_names{unit} ' n ='  num2str(size(time_locked_firing{unit},1))]);
    end
    
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_ITI_analysis']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Finally save all the data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([data_dir task_file(1:8) '-ListSQ-ITI_results.mat'],'time_locked_firing',...
    'previous_trial_type','current_trial_type','smval','temporal_info',...
    'unit_names','ANOVA_stats');
disp(['ITI Locked Data Analyis for ' task_file(1:8) ' saved']);
end