function Visual_Response_Analysis(data_dir,figure_dir,session_data,task)
% written by Seth Konig July 11, 2016
% Code looks at whether neurons show a visual response after image onset.
% Code parallels time_locked_analysisV2.m and is somewhat redundant.


twin = 750;% how much time to take before and after an event.
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
numshuffs = 1; %number of shuffles to do for bootstrapping
fr_threshold = 1; %peak rate must be greater than 1 Hz to process. Don't want to waste
%processing/shuffling time on "silent neurons"
smval = 60; %temporal 1/2 width of gaussian smoothing filter


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import & reformat data so that spikes are locked to events---%%%

num_trials = length(cfg.trl);

%NaNs are for start and end trials otherwise cut
valid_trials(1,isnan(valid_trials(1,:))) = 1;
valid_trials(2,isnan(valid_trials(2,:))) = num_trials;


disp('Aligning spike times to trial events')

%get important task specific information
[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,23);

%same epoch numbers as time_locked_analysisV2.m
% epoch 14 fixation cross is displayed
% epoch 15 fixation on crosshair to initialize image trials
% epoch 16 image viewing period start


% time_locked_firing since ITI period is defined by 2 events (15 & 16)
event_names = {'Fix Spot On','Fixation on Fix Spot','Image On'};
event_codes = [35 8 23];

%image duration changed part way through from 7 secs to 5 secs of
%cumulative looking times
if task_file(1:8) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end
event_durs = [zeros(1,2) imgdur]; %event durations
event_t0 = twin*ones(1,3); %in the matrix what index is equivalent to the time the event occured


%preallocate space and parallel structure of cfg
time_lock_firing = make_time_lock_cell(event_durs,num_units,num_trials,twin);

which_images = cell(1,num_units);
nvr = cell(1,num_units);
for unit = 1:num_units
    nvr{unit} = NaN(1,length(cfg.trl));
    which_images{unit} = NaN(1,length(cfg.trl));
end

for t = 1:num_trials
    if any(cfg.trl(t).allval == event_codes(2)); %in which image was displayed or 1st item in sequence was displayed
        if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(end) %then sequence trial
            continue % go to the next trial
        end
        if  any(cfg.trl(t).allval == 23)%image turns on
            trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == 15);
            event_times = get_List_event_times(cfg.trl(t),event_codes,trial_start);
            for unit = 1:num_units
                if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only put data in for valid trials
                    time_lock_firing = put_in_time_matrix(data,unit,t,event_times,time_lock_firing,twin,event_names);
                    img_index = find(cfg.trl(t).cnd == img_cnd);
                    nvr{unit}(img_index) = novel_vs_repeat(img_index);
                    which_images{unit}(img_index) = which_img(img_index);
                end
            end
        end
    end
end

%remove excess NaNs associated with error trials
time_lock_firing = laundry(time_lock_firing);
nvr = laundry(nvr);
which_images = laundry(which_images);

%for comparing novel to repeat only want the trials in which both
%images were shown
for unit = 1:num_units
    rmv = [];
    for img = 1:96
        ind = find(which_images{unit} == img);
        if length(ind) == 1 %so either novel or repeat but not both
            rmv = [rmv ind];
        end
    end
    for event = 3%4:5;
        time_lock_firing{event,unit}(rmv,:) = [];
    end
    which_images{unit}(rmv) = [];
    nvr{unit}(rmv) = [];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Perform Temporal analysis and signifance testing---%%%
disp('Determining if spikes are locked to trial events')
Fs = data(1).fsample; %should be 1000

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

%---Calculate Peak firing rate for all epochs
peak_firing_rate = NaN(size(time_lock_firing));
for unit = 1:num_units
    for event = 1:size(time_lock_firing,1)
        if ~isempty(time_lock_firing{event,unit})
            [firing_rate,~]= nandens(time_lock_firing{event,unit},smval,'gauss',Fs,'nanflt');%nandens made by nathan killian
            peak_firing_rate(event,unit) = max(firing_rate);
        end
    end
end

% perform mutual information theory analysis to determine
% if neurons encode info about a particular time within an epoch
for event = 1:size(time_lock_firing,1);
    for unit = 1:num_units
        if ~isempty(time_lock_firing{event,unit}) && any(peak_firing_rate(:,unit) > fr_threshold) %only processs active neurons
            temp_time_lock = time_lock_firing{event,unit};
            [temporal_info.rate(event,unit),temporal_info.shuffled_info_rate{event,unit}]...
                = estimated_mutual_information(temp_time_lock,numshuffs,info_type,smval,Fs);
            
            epoch_data.firing_rates{event,unit} = nansum(temp_time_lock(:,twin:end),2);
            epoch_data.dur(event,unit) = (size(temp_time_lock,2)-twin)/1000; %duration in seconds
            epoch_data.num_trials(event,unit) = sum(~isnan(epoch_data.firing_rates{event,unit}));
        end
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%---Plot and Save Figures of Results---%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for unit = 1:num_units
    if ~isempty(time_lock_firing{1,unit})
        figure
        
        ylims = NaN(1,7);
        s = [];
        
        %plot from Fixation Spot on
        [s(1),ylims(1)]=make_subplot([2,3],1,time_lock_firing{1,unit},...
            epoch_data.event_t0(1),250,epoch_data.event_names{1},....
            temporal_info.rate(1,unit),temporal_info.shuffled_95_percentile(1,unit),...
            temporal_info.shuffled_90_percentile(1,unit),smval);
        
        %plot from (cortex) Fixation on Fixation Spot on
        [s(2),ylims(2)]=make_subplot([2,3],4,time_lock_firing{2,unit},...
            epoch_data.event_t0(2),250,epoch_data.event_names{2},....
            temporal_info.rate(2,unit),temporal_info.shuffled_95_percentile(2,unit),...
            temporal_info.shuffled_90_percentile(2,unit),smval);
        
        %plot from image onset
        [s(3),ylims(3)]=make_subplot([2,3],[2 3 5 6],time_lock_firing{3,unit},...
            epoch_data.event_t0(3),250,epoch_data.event_names{3},....
            temporal_info.rate(3,unit),temporal_info.shuffled_95_percentile(3,unit),...
            temporal_info.shuffled_90_percentile(3,unit),smval);
        
        %plot from image onset for novel vs repeat
        novel_time_lock_firing = time_lock_firing{3,unit}(nvr{unit} == 1,:);
        repeat_time_lock_firing = time_lock_firing{3,unit}(nvr{unit} == 2,:);
        [s(3),ylims(3)]=make_subplot_2CND([2,3],[2 3 5 6],novel_time_lock_firing,repeat_time_lock_firing,...
            epoch_data.event_t0(3),250,epoch_data.event_names{3},...
            temporal_info.rate(3,unit),temporal_info.shuffled_95_percentile(3,unit),...
            temporal_info.shuffled_90_percentile(3,unit),smval);
        
        
        
        max_y = 5*round(max(ylims)/5);
        max_y(max_y < 1) = 1;
        for i = 1:3;
            set(s(i),'ylim',[0 max_y])
        end
        
        if multiunit(unit)
            subtitle(['Multiunit ' unit_stats{1,unit} ]);
        else
            subtitle(unit_stats{1,unit});
        end
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Visual_Response_analysis']);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Finally save all the data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat'],'time_lock_firing',...
     'smval','temporal_info','epoch_data','unit_stats','peak_firing_rate');
disp(['Time Locked Data Analyis for ' task_file(1:8) ' saved']);


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

function event_times = get_List_event_times(cfgtrl,event_codes,trial_start)
% function determines event times in a given trial
event_times = cell(1,size(event_codes,2));
for events = 1:size(event_codes,2)
    temp = cfgtrl.alltim(cfgtrl.allval == event_codes(events))-trial_start;
    temp = temp(1);
    event_times{events} = temp;
end
end

function [time_lock_firing] = put_in_time_matrix(data,unit,t,event_times,time_lock_firing,twin,event_names)
% function finds spikes that are time locked to events
spikes = find(data(unit).values{t});
for event = 1:length(event_times)
    if ~isempty(event_times{event})
        if strcmpi(event_names{event},'Image On')
            event_spikes = spikes(spikes > event_times{event}-twin & spikes <= event_times{event}+2*twin)-event_times{event}+twin;
            timevec = [event_times{event}-twin+1 event_times{event}+2*twin]-event_times{event}+twin;
            tempvec = zeros(1,timevec(2)-timevec(1)+1);
            tempvec(event_spikes) = 1;
            time_lock_firing{event,unit}(t,timevec(1):timevec(2)) = tempvec;
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


function [s,ylimit]=make_subplot(src,plot_nums,time_matrix,t0,timestep,epoch_name,inforate,info95,info90,smval)
% src: subplot row and column
% plot_nums: plot number or numbers
% inforate: observed inforamtion rate
% info95: bootstrapped 95 percentile

s = subplot(src(1),src(2),plot_nums);
t = 1:size(time_matrix,2);
if ~isempty(t)
    [~,~,~,y] =dofill(t,time_matrix,'black',1,smval);
    set(gca,'Xtick',0:timestep:size(time_matrix,2))
    set(gca,'XtickLabel',num2cell((0:timestep:size(time_matrix,2))-t0));
    ylabel('Firing Rate (Hz)')
    xlabel(['Time from ' epoch_name ' (ms)'])
    num_trial_str = ['n = ' num2str(size(time_matrix,1))];
    if any(inforate > info95)
        title([num_trial_str ', Bits_{95} = ' num2str(inforate(end))]);
    elseif any(inforate > info90)
        title([num_trial_str ', Bits_{90} = ' num2str(inforate(end))]);
    else
        title(num_trial_str);
    end
    ylimit = max(y)*1.2;
else
    ylimit = NaN;
end
end

function [s,ylimit]= make_subplot_2CND(src,plot_nums,time_matrix1,time_matrix2,t0,timestep,epoch_name,inforate,info95,info90,smval)
% essentailly same as above just 2 time_matrices
s = subplot(src(1),src(2),plot_nums);
t = 1:size(time_matrix1,2);
if ~isempty(t)
    dofill(t,time_matrix1,'blue',1,smval);
    dofill(t,time_matrix2,'red',1,smval);
end
set(gca,'Xtick',0:timestep:size(time_matrix1,2))
set(gca,'XtickLabel',num2cell((0:timestep:size(time_matrix1,2))-t0));
ylabel('Firing Rate (Hz)')
xlabel(['Time from ' epoch_name ' (ms)'])

num_trial_str = ['n_1 = ' num2str(size(time_matrix1,1)) ' n_2 = ' num2str(size(time_matrix1,1))];
if any(inforate > info95)
    title([num_trial_str ', Bits_{95} = ' num2str(inforate(end))]);
elseif any(inforate > info90)
    title([num_trial_str ', Bits_{90} = ' num2str(inforate(end))]);
else
    title(num_trial_str)
end
ylimit = ylim;
ylimit = ylimit(2);
end