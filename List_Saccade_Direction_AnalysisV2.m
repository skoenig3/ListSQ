function List_Saccade_Direction_AnalysisV2(data_dir,figure_dir,session_data)
% written by Seth Konig January 18, 2016 based on List_Saccade_Analysis
% information locked to saccades and fixations on each item.
% updated SDK 4/4/16 to handlde new format and partial session data for
% vaild trials only.CVTNEW section on 1/19/16
%
% Function analyizes spike times correlated with saccade direction in the
% List image portion of the task. Saccade direction is computed as the
% angle between fixations not the actual saccade directio leaving a
% fixation. Saccades with ampltiudes < 2 dva are ignored
%
% Inputs:
%   1) data_dir: direction where preprocessed data is located
%   2) preprocessed_data_file: preprocessed ListSQ data file with LFPs
%   divided by channel and trial.
%   3) figure_dir: location of where to put save figures
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-Eyemovement_Locked_List_results'

task = 'ListSQ';
twin = 150;% how much time to take before and after saccade.
tstart = 500; %how much after the image turns on to ignore, trying to ignore visual response
bin_deg =4; %number of degrees per bin
trial_start_code = 15;
img_on_code = 23;
img_off_code = 24;
smval = 4;% moving average filter width for saccade angle
smval2 = 12;% gaussian window std for dofill smoothign of saccade locked activity
saccade_ampltidue_threshold = 2*24; %miniumum saccade amplitude to use

degrees = [0:bin_deg:360]-180; %binned degrees
amplitude = saccade_ampltidue_threshold/24:2:25; %binned amplitude

Fs = 1000;
numshuffs = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import data & get successful trials---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
%grab unit data.
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials',...
    'hdr','fixationstats');
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

[itmlist,sequence_items] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);

num_trials = length(cfg.trl);
%NaNs are for start and end trials otherwise cut
valid_trials(1,isnan(valid_trials(1,:))) = 1;
valid_trials(2,isnan(valid_trials(2,:))) = num_trials;

%set the image duration
if str2double(task_file(3:8)) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Preallocate spaces for some variables---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just want to do this in one spot
spike_times = cell(1,num_units); %spike times organized by trial (row) and time (column)

sac_direction = cell(1,num_units);%direction from current fixation to future fixation
sac_distance = cell(1,num_units);%amplitude from current fixation to future fixation
saccade_start_time = cell(1,num_units);%when did saccade to item start
fixation_start_time = cell(1,num_units);%when did fixation on item start

for unit = 1:num_units
    spike_times{unit} = NaN(num_trials,imgdur*1.5);
    
    sac_direction{unit} = NaN(num_trials,50);
    sac_distance{unit} = NaN(num_trials,50);
    saccade_start_time{unit} = NaN(num_trials,50);
    fixation_start_time{unit} = NaN(num_trials,50);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---process eye data locked to Fixations---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for t = 1:num_trials
    if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
        fixationtimes = fixationstats{t}.fixationtimes;
        saccadetimes = fixationstats{t}.saccadetimes;
        fixations = fixationstats{t}.fixations;
        
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
        imgon = cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start+tstart;
        imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start;
        
        % if monkey isn't paying attention data isn't probably
        % worth much plus have to cut off somewhere
        if imgoff-imgon > 1.5*imgdur-1
            imgoff = imgon+1.5*imgdur-1;
        end
        
        %find fixations and saccades that did not occur during the image period;
        %should also take care of the 1st fixation on the crosshair
        
        %fixation started before image turned on
        invalid= find(fixationtimes(1,:) < imgon);
        fixationtimes(:,invalid) = [];
        
        %fixation started after the image turned off and/or firing rate could corrupted by image turning off
        invalid= find(fixationtimes(1,:) > imgoff+twin);
        fixationtimes(:,invalid) = [];
        
        %saccade started before image turned on
        invalid= find(saccadetimes(1,:) < imgon);
        saccadetimes(:,invalid) = [];
        
        %saccade started after the image turned off and/or firing rate could corrupted by image turning off
        invalid= find(saccadetimes(1,:) > imgoff+twin);
        saccadetimes(:,invalid) = [];
        %no fixation to follow so cut
        
        if ~isempty(fixationtimes)
            if saccadetimes(1,end) > fixationtimes(2,end)
                saccadetimes(:,end) = [];
            end
        end
        
        for unit = 1:num_units
            if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
                for s = 1:size(saccadetimes,2);
                    
                    %doing pre direction/amplitude locked to saccades upcomming movement direction
                    
                    sacstart = saccadetimes(1,s);
                    sacend = saccadetimes(2,s);
                    
                    %not straight forward since could have been looking
                    %off screen before saccade was detected or blinked
                    prior_fix = find(sacstart-1 ==  fixationtimes(2,:));
                    next_fix = find(sacend+1 ==  fixationtimes(1,:));
                    
                    if ~isempty(prior_fix) && ~isempty(next_fix);
                        if (next_fix-prior_fix) ~= 1
                            error('fixaiton indeces should be consecutive')
                        end
                        
                        sac_direction{unit}(t,s) = atan2d(fixations(2,next_fix)...
                            -fixations(2,prior_fix),fixations(1,next_fix)-fixations(1,prior_fix));
                        sac_distance{unit}(t,s) = sqrt((fixations(2,next_fix)...
                            -fixations(2,prior_fix)).^2+(fixations(1,next_fix)...
                            -fixations(1,prior_fix)).^2);
                        
                        saccade_start_time{unit}(t,s) = saccadetimes(1,s)-imgon; %t0 zeroed with spike_times for pre-times
                        fixation_start_time{unit}(t,s) = fixationtimes(1,next_fix)-imgon; %t0 zeroed with spike_times for post-times
                    end
                end
                
                if  all(isnan(saccade_start_time{unit}(t,:))) %add fillers so can maintain structure during NaN cleanup these
                    %data will not be processed since ampltiude is set to -1
                    sac_distance{unit}(t,1) = -1;
                    sac_direction{unit}(t,1) = -1;
                    saccade_start_time{unit}(t,1) = -1;
                    fixation_start_time{unit}(t,1) = -1;
                end
            end
        end
    end
end

%remove extra NaNs/trials
sac_direction = laundry(sac_direction);
sac_distance = laundry(sac_distance);
saccade_start_time = laundry(saccade_start_time);
fixation_start_time = laundry(fixation_start_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Get spike times matrix---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:num_trials
    if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
        
        imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start+tstart; %want to avoid visual response
        imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start;
        
        % if monkey isn't paying attention data isn't probably
        % worth much plus have to cut off somewhere
        if imgoff-imgon > 1.5*imgdur
            imgoff = imgon+1.5*imgdur;
        end
        
        for unit = 1:num_units
            if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
                spikes = find(data(unit).values{t});
                if ~isempty(spikes)
                    spikeind = spikes(spikes >= imgon & spikes <= imgoff)-imgon;
                    spikeind(spikeind == 0) = []; %happened exactly at the same time
                    temp = zeros(1,(imgoff-imgon));
                    temp(spikeind) = 1;
                    spike_times{unit}(t,1:length(temp)) = temp;
                else
                    temp = zeros(1,(imgoff-imgon));
                    spike_times{unit}(t,1:length(temp)) = temp;
                end
            end
        end
    end
end
spike_times = laundry(spike_times,1); %remove excess NaNs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Firing Rate Locked to Eye Data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%preallocate space for variables
saccade_direction_firing_pre = cell(1,num_units);
saccade_direction_firing_post = cell(1,num_units);
saccade_amplitude_firing_pre = cell(1,num_units);
saccade_amplitude_firing_post = cell(1,num_units);
saccade_locked_firing = cell(1,num_units);
fixation_locked_firing = cell(1,num_units);
for unit = 1:num_units
    saccade_direction_firing_pre{unit} = NaN(size(saccade_start_time{unit}));
    saccade_direction_firing_post{unit} = NaN(size(saccade_start_time{unit}));
    saccade_amplitude_firing_pre{unit} = NaN(size(saccade_start_time{unit}));
    saccade_amplitude_firing_post{unit} = NaN(size(saccade_start_time{unit}));
    
    saccade_locked_firing{unit} = NaN(50*num_trials,twin*2);
    fixation_locked_firing{unit} = NaN(50*num_trials,twin*2);
end

sacindex = ones(1,num_units);
for unit = 1:num_units;
    for trial = 1:size(saccade_start_time{unit},1)
        spikes = find(spike_times{unit}(trial,:) == 1);
        for s = 1:size(saccade_start_time{unit},2);
            
            sact = saccade_start_time{unit}(trial,s);
            fixt = fixation_start_time{unit}(trial,s);
            
            %hard to conistently detect < 1.5, and within fovea and want
            %saccades that may facilitate foraging/exploring not
            %exploitation and information/detail gathering at same location
            if sac_distance{unit}(trial,s) < saccade_ampltidue_threshold;
                continue
            end
            
            if ~isnan(sact) %should also be fixt as well
                sac_spikes = find(spikes >= sact-twin & spikes < sact);
                fix_spikes = find(spikes > fixt & spikes <= fixt+twin);
                
                %pre firing
                saccade_amplitude_firing_pre{unit}(trial,s) = length(sac_spikes)*(Fs/twin); %firing rate
                saccade_direction_firing_pre{unit}(trial,s) = length(sac_spikes)*(Fs/twin);%firing rate
                
                %post firing
                saccade_amplitude_firing_post{unit}(trial,s) = length(fix_spikes)*(Fs/twin); %firing rate
                saccade_direction_firing_post{unit}(trial,s) = length(fix_spikes)*(Fs/twin);%firing rate
                
                
                %to compute firing rate around time of saccades
                sac_spikes = spikes(spikes > sact-twin & spikes <= sact+twin)-sact+twin;
                fix_spikes = spikes(spikes > fixt-twin & spikes <= fixt+twin)-fixt+twin;
                
                temp = zeros(1,twin*2);
                temp(fix_spikes) = 1;
                fixation_locked_firing{unit}(sacindex(unit),:) = temp;
                
                temp = zeros(1,twin*2);
                temp(sac_spikes) = 1;
                saccade_locked_firing{unit}(sacindex(unit),:) = temp;
                
                sacindex(unit) = sacindex(unit)+1;
                
            end
        end
    end
end
saccade_locked_firing = laundry(saccade_locked_firing);
fixation_locked_firing = laundry(fixation_locked_firing);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Shuffled Firing Rate Locked to Eye Data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shuffled_mean_vectors = NaN(num_units,numshuffs,2); %peak firing rate for shuffled data
%unit by row, shuffled number by col, and fix vs sac zcol
for unit = 1:num_units;
    for shuff = 1:numshuffs
        shuffled_spikes = circshift_acrosstrials(spike_times{unit}); %shuffle spike times
        
        saccade_firing = NaN(size(saccade_start_time{unit}));
        fixation_firing = NaN(size(saccade_start_time{unit}));
        for trial = 1:size(saccade_start_time{unit},1)
            spikes = find(shuffled_spikes(trial,:) == 1);
            for s = 1:size(saccade_start_time{unit},2);
                
                %use same eye movement times
                sact = saccade_start_time{unit}(trial,s);
                fixt = fixation_start_time{unit}(trial,s);
                
                %hard to conistently detect < 1.5, and within fovea and want
                %saccades that may facilitate foraging/exploring not
                %exploitation and information/detail gathering at same location
                if sac_distance{unit}(trial,s) < saccade_ampltidue_threshold;
                    continue
                end
                
                if ~isnan(sact) %should also be fixt as well
                    sac_spikes = find(spikes >= sact-twin & spikes < sact);
                    fix_spikes = find(spikes > fixt & spikes <= fixt+twin);
                    
                    saccade_firing(trial,s) = length(sac_spikes)*(Fs/twin);%firing rate
                    fixation_firing(trial,s) = length(fix_spikes)*(Fs/twin);%firing rate
                end
            end
        end
        
        directional_firing_pre = cell(1,length(degrees));
        directional_firing_post = cell(1,length(degrees));
        for bin = 2:length(degrees)
            pre_ind = find(sac_direction{unit} < degrees(bin) & sac_direction{unit} >= degrees(bin-1));
            directional_firing_pre{bin} = [directional_firing_pre{bin}; saccade_firing(pre_ind)];
            
            post_ind = find(sac_direction{unit} < degrees(bin) & sac_direction{unit} >= degrees(bin-1));
            directional_firing_post{bin} = [directional_firing_post{bin}; fixation_firing(post_ind)];
        end
        
        means = cellfun(@nanmean,directional_firing_pre(2:end));
        means =  [means(end-6:end) means means(1:7)];
        means = filtfilt(1/smval*ones(1,smval),1,means);
        means = means(8:end-7);
        
        shuffled_mean_vectors(unit,shuff,1) = max(means);
        
        
        means = cellfun(@nanmean,directional_firing_post(2:end));
        means =  [means(end-6:end) means means(1:7)];
        means = filtfilt(1/smval*ones(1,smval),1,means);
        means = means(8:end-7);
        
        shuffled_mean_vectors(unit,shuff,2) = max(means);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Plots for Saccade Direction/Ampltiude Activity---%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directional_firing_pre = cell(length(degrees),num_units);
directional_firing_post = cell(length(degrees),num_units);
amplitude_firing_pre = cell(length(amplitude)-1,num_units);
amplitude_firing_post = cell(length(amplitude)-1,num_units);

for unit = 1:num_units
    
    for bin = 1:length(amplitude)-1
        pre_ind = find(sac_distance{unit}/24 >= amplitude(bin) & sac_distance{unit}/24 < amplitude(bin+1));
        amplitude_firing_pre{bin,unit} = [amplitude_firing_pre{bin,unit}; saccade_amplitude_firing_pre{unit}(pre_ind)];
        
        post_ind = find(sac_distance{unit}/24 >= amplitude(bin) & sac_distance{unit}/24 < amplitude(bin+1));
        amplitude_firing_post{bin,unit} = [amplitude_firing_post{bin,unit}; saccade_amplitude_firing_post{unit}(post_ind)];
        
    end
    
    for bin = 2:length(degrees)
        pre_ind = find(sac_direction{unit} < degrees(bin) & sac_direction{unit} >= degrees(bin-1));
        directional_firing_pre{bin,unit} = [directional_firing_pre{bin,unit}; saccade_direction_firing_pre{unit}(pre_ind)];
        
        post_ind = find(sac_direction{unit} < degrees(bin) & sac_direction{unit} >= degrees(bin-1));
        directional_firing_post{bin,unit} = [directional_firing_post{bin,unit}; saccade_direction_firing_post{unit}(post_ind)];
    end
end

degrees = degrees(2:end);
degrees = degrees*pi/180;
degrees = [degrees degrees(1)];
polar_firing = cell(2,num_units);
direction_info.uniformity_pval = NaN(2,num_units);
direction_info.peak_firing_rate = NaN(2,num_units);
direction_info.prefered_direction = NaN(2,num_units);
for unit = 1:num_units;
    means = cellfun(@nanmean,directional_firing_pre(2:end,unit))';
    means =  [means(end-6:end) means means(1:7)];
    means = filtfilt(1/smval*ones(1,smval),1,means);
    means = means(8:end-7);
    
    polar_firing{1,unit} = [means means(1)];
    
    direction_info.uniformity_pval(1,unit) = circ_rtest(degrees(1:end-1),means,bin_deg);
    direction_info.prefered_direction(1,unit) = circ_mean(degrees(1:end-1),means,2);
    direction_info.peak_firing_rate(1,unit) = max(means);
    direction_info.shuffled_peak(1,unit) = {shuffled_mean_vectors(unit,:,1)};
    direction_info.shuffled_peak_prctile(1,unit) = 100*sum(shuffled_mean_vectors(unit,:,1) < ...
        direction_info.peak_firing_rate(1,unit))/numshuffs;
    
    
    means = cellfun(@nanmean,directional_firing_post(2:end,unit))';
    means =  [means(end-6:end) means means(1:7)];
    means = filtfilt(1/smval*ones(1,smval),1,means);
    means = means(8:end-7);
    
    
    direction_info.uniformity_pval(2,unit) = circ_rtest(degrees(1:end-1),means,bin_deg);
    direction_info.prefered_direction(2,unit) = circ_mean(degrees(1:end-1),means,2);
    direction_info.peak_firing_rate(2,unit) = max(means);
    direction_info.shuffled_peak(2,unit) = {shuffled_mean_vectors(unit,:,2)};
    direction_info.shuffled_peak_prctile(2,unit) = 100*sum(shuffled_mean_vectors(unit,:,2) < ...
        direction_info.peak_firing_rate(2,unit))/numshuffs;
    
    polar_firing{2,unit} = [means means(1)];
end


unit_names = cfg.channel(1:num_units);
t = -twin:twin-1;
for unit = 1:num_units
    
    plot_lim = 1.25*max(direction_info.peak_firing_rate(:,unit));
    plot_lim(plot_lim < .25) = 0.25; %scale to minumum of 0.25 Hz
    
    figure
    
    %for pre-saccade firing period
    subplot(2,3,1)
    bar(amplitude(1:end-1),cellfun(@nanmean,amplitude_firing_pre(:,unit)));
    xlabel('Saccade Amplitude (dva)')
    ylabel('Firing Rate (Hz)')
    title('Pre-Saccade: Saccade Amplitude')
    
    subplot(2,3,2)
    polar(pi,plot_lim,'w') %invisible point so scales are the same
    hold on
    polar(degrees,polar_firing{1,unit},'b')
    direction_info.prefered_direction(1,unit) = circ_mean(degrees(1:end-1),means,2);
    direction_info.peak_firing_rate
    plot([0 cos(direction_info.prefered_direction(1,unit))*direction_info.peak_firing_rate(1,unit)],...
        [0 sin(direction_info.prefered_direction(1,unit))*direction_info.peak_firing_rate(1,unit)],'r') %indicate prefered direction
    hold off
    
    title_str = 'Pre-Saccade: Saccade Direction';
    direction_info.shuffled_peak_prctile
    if direction_info.uniformity_pval(1,unit) < 0.06
        title_str = [title_str '\n p = ' num2str(direction_info.uniformity_pval(1,unit),3)];
    end
    if direction_info.shuffled_peak_prctile(1,unit) > 90
        title_str = [title_str '\n Shuffled_{peak} ' num2str(direction_info.shuffled_peak_prctile(1,unit),3)];
    end
    title(sprintf(title_str))
    
    
    subplot(2,3,3)
    [~,~,~,y1] = dofill(t,saccade_locked_firing{unit},'k',1,smval2);
    xlabel('Time from Saccade Start (ms)')
    ylabel('Firing Rate (Hz)')
    
    %for post saccade diring period
    subplot(2,3,4)
    bar(amplitude(1:end-1),cellfun(@nanmean,amplitude_firing_pre(:,unit)));
    xlabel('Saccade Amplitude (dva)')
    ylabel('Firing Rate (Hz)')
    title('Post-Saccade: Saccade Amplitude')
    
    subplot(2,3,5)
    polar(pi,plot_lim,'w') %invisible point so scales are the same
    hold on
    polar(degrees,polar_firing{2,unit},'b')
    plot([0 cos(direction_info.prefered_direction(2,unit))*direction_info.peak_firing_rate(2,unit)],...
        [0 sin(direction_info.prefered_direction(2,unit))*direction_info.peak_firing_rate(2,unit)],'r') %indicate prefered direction
    hold off
    
    title_str = 'Post-Saccade: Saccade Direction';
    direction_info.shuffled_peak_prctile
    if direction_info.uniformity_pval(2,unit) < 0.06
        title_str = [title_str '\n p = ' num2str(direction_info.uniformity_pval(2,unit),3)];
    end
    if direction_info.shuffled_peak_prctile(2,unit) > 90
        title_str = [title_str '\n Shuffled_{peak} ' num2str(direction_info.shuffled_peak_prctile(2,unit),3)];
    end
    title(sprintf(title_str))
    
    subplot(2,3,6)
    [~,~,~,y2] = dofill(t,fixation_locked_firing{unit},'k',1,smval2);
    xlabel('Time from Fixation Start (ms)')
    ylabel('Firing Rate (Hz)')
    
    %scale eye locked firing rate to be the same
    ymin = 0.8*min([min(y1) min(y2)]);
    ymax = 1.25*min([max(y1) max(y2)]);
    subplot(2,3,3)
    ylim([ymin ymax])
    subplot(2,3,6)
    ylim([ymin ymax])
    
    bitstr = [];
    if multiunit(unit)
        subtitle(['Multiunit ' unit_names{unit} ' n =' num2str(sum(cellfun(@numel,amplitude_firing_post(:,unit)))) ' ' bitstr]);
    else
        subtitle([unit_names{unit} ' n =' num2str(sum(cellfun(@numel,amplitude_firing_post(:,unit)))) ' ' bitstr]);
    end
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_List_Direction_analysis']);
    
end

save([data_dir task_file(1:end-11) '-List_Direction_results.mat'],...
    'twin','smval','bin_deg','amplitude_firing_pre','amplitude_firing_post',...
    'directional_firing_pre','directional_firing_post','direction_info',...
    'fixation_locked_firing','saccade_locked_firing');