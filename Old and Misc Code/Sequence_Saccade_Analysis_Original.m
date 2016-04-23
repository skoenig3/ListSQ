function Sequence_Saccade_Analysis(data_dir,preprocessed_data_file,figure_dir,task)
% written by Seth Konig September, 2014
% Function analyses spike times correlated with eye movements in the
% sequence task.
%
% Inputs:
%   1) data_dir: direction where preprocessed data is located
%   2) preprocessed_data_file: preprocessed ListSQ data file with LFPs
%   divided by channel and trial.
%   3) figure_dir: location of where to put save figures
%   4) task: what task was this data come from with i.e. 'ListSQ' or 'Sequence'
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-sequence_saccade_locked_results'


twin = 500;% how much time to take before and after saccade.
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
fixwin = 2.5*24;%size of fixation window on each crosshair
event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
trial_start_code = 15;
smval = 120;%gaussian 1/2 width for smoothing
predicted_thresh = 10;% percent of saccades that must be < 150 ms to item to constitue as a learned-predicted item
% for Vivian 5% was the false positive rate for 150 ms so wanted at least double the False positive rate
predicted_rt = 150;%maximum "reaction time" for what constitutes as predictive, everything else is reactive

switch task
    case 'ListSQ'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import data & get successful trials---%%%
        
        load([data_dir preprocessed_data_file]);
        
        %get important task specific information
        [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
        
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

saccade_start_time = NaN(length(fixationstats),4);%when did saccade to item start
fixation_start_time = NaN(length(fixationstats),4);%when did fixation on item start
reaction_time = NaN(length(fixationstats),4); %how long after a stimulus disappears do they saccade
time_to_fixation = NaN(length(fixationstats),4); %how fast did they get to it the first time around
fixation_accuracy = NaN(length(fixationstats),4); %how far off
fixation_duration = NaN(length(fixationstats),4); %fixation duration
fixation_on_test = NaN(length(fixationstats),4); %how long the stimulus was on while fixating
extrafixations = NaN(length(fixationstats),4); %how many noncross hair fixations do they mak
for trial = 1:num_trials
    fixations = fixationstats{trial}.fixations;
    fixationtimes = fixationstats{trial}.fixationtimes;
    saccadetimes = fixationstats{trial}.saccadetimes;
    xy = fixationstats{trial}.XY;
    
    % Find saccades that were small in amplitude and determine if these
    % smaller saccades were corrective and I should essentially count
    % the pre and post saccade fixations as both on the cross hair.
    saccade_amplitudes = NaN(1,size(saccadetimes,2));
    for sacc = 1:size(saccadetimes,2);
        sacx = xy(1,saccadetimes(2,sacc))-xy(1,saccadetimes(1,sacc));
        sacy = xy(2,saccadetimes(2,sacc))-xy(2,saccadetimes(1,sacc));
        saccade_amplitudes(sacc) = sqrt(sacx.^2+sacy.^2);
    end
    if any(saccade_amplitudes(2:end) <= fixwin/2) %going to generalize  these as corrective saccades. Ignore 1st saccade
        [fixationtimes,fixations,saccadetimes] = remove_corrective_saccades(...
            xy,saccadetimes,saccade_amplitudes,fixwin);
    end
    
    locs = sequence_locations{which_sequence(trial)};
    
%         figure
%         hold on
%         plot(xy(1,:),xy(2,:),'g');
%         %         for f = 1:length(fixationtimes);
%         %             plot(fixations(1,f),fixations(2,f),'k*','markersize',5);
%         %         end
%         for c = 1:size(locs,2);
%             plot(locs(1,c),locs(2,c),'kx','markersize',12)
%         end
%         xlim([0 800])
%         ylim([0 600])
    
    valid_fixations = ones(1,size(fixations,2));
    for c = 1:size(locs,2)
        if c == 1
            valid_window = [event_times(trial,2*c-1)-twin/2 event_times(trial,2*c)];
        else
            valid_window = [event_times(trial,2*c-1) event_times(trial,2*c)];
        end
        valid_ind = valid_window(1):valid_window(2);
        
        potential_fix = [];
        for f = 1:size(fixationtimes,2);
            if valid_fixations(f)
                fixtimes = fixationtimes(1,f):fixationtimes(2,f);
                C = intersect(fixtimes,valid_ind);
                if length(C) >= 50;
                    potential_fix = [potential_fix f];
                end
            end
        end
        if ~isempty(potential_fix)
            if length(potential_fix)  == 1;
                fixstart = fixationtimes(1,potential_fix);
                valid_fixations(1:potential_fix) = 0;
                dist = sqrt((fixations(1,potential_fix)-locs(1,c)).^2+...
                    (fixations(2,potential_fix)-locs(2,c)).^2);
            elseif length(potential_fix) > 1;
                dist = sqrt((fixations(1,potential_fix)-locs(1,c)).^2+...
                    (fixations(2,potential_fix)-locs(2,c)).^2);
                [~,thefix] = min(dist);
                fixstart = fixationtimes(1,potential_fix(thefix));
                valid_fixations(1:potential_fix(thefix)) = 0;
                potential_fix = potential_fix(thefix);
                dist = dist(thefix);
                if thefix ~= 1
                    extrafixations(trial,c) = thefix-1;
                end
            end
            

            time_to_fixation(trial,c) = fixstart-event_times(trial,2*c-1);
            
            if (time_to_fixation(trial,c) < -250) && (c > 1)
                disp('error')
            end
            
            fixation_duration(trial,c) = fixationtimes(2,potential_fix)-fixationtimes(1,potential_fix)+1;
            if c < size(locs,2)
                reaction_time(trial,c) = fixationtimes(2,potential_fix)-event_times(trial,2*c);
            end
            timeON = event_times(trial,2*c-1):event_times(trial,2*c);
            fixON = fixationtimes(1,potential_fix):fixationtimes(2,potential_fix);
            fixation_on_test(trial,c) = length(intersect(timeON,fixON));
            fixation_accuracy(trial,c) = dist;
            
            fixation_start_time(trial,c) = fixationtimes(1,potential_fix);
            if saccadetimes(1,1) < fixationtimes(1,1)
                saccade_start_time(trial,c) = saccadetimes(1,potential_fix);
            elseif potential_fix ~= 1
                saccade_start_time(trial,c) = saccadetimes(1,potential_fix-1); %else remains a NaN
            end
            
            
%                         plot(xy(1,fixationtimes(1,potential_fix):fixationtimes(2,potential_fix)),...
%                             xy(2,fixationtimes(1,potential_fix):fixationtimes(2,potential_fix)),'c');
%                         plot(xy(1,fixstart),xy(2,fixstart),'m*','markersize',10)
            
        end
    end
%         pause(0.5)
%         close all
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
        dofill(t,saccade_locked_firing{c,unit}(which_sequence == 1 & ~predicted(:,c)',:),'red',1,smval); %reactive sequence 1
        dofill(t,saccade_locked_firing{c,unit}(which_sequence == 2 & ~predicted(:,c)',:),'blue',1,smval);%reactive sequence 2
        if predicted_items(1,c)
            dofill(t,saccade_locked_firing{c,unit}(which_sequence == 1 & predicted(:,c)',:),'magenta',1,smval); %predictive sequence 1
        end
        if predicted_items(2,c)
            dofill(t,saccade_locked_firing{c,unit}(which_sequence == 2 & predicted(:,c)',:),'cyan',1,smval);%predictive sequence 2
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
            title(['Item # ' num2str(c) ' Seq #1 ' num2str(percent_predicted(1,c)) '%' ])
        elseif predicted_items(2,c)
            title(['Item # ' num2str(c) ' Seq #2 ' num2str(num2str(percent_predicted(2,c))) '%'])
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
        dofill(t,fixation_locked_firing{c,unit}(which_sequence == 1 & ~predicted(:,c)',:),'red',1,smval); %reactive sequence 1
        dofill(t,fixation_locked_firing{c,unit}(which_sequence == 2 & ~predicted(:,c)',:),'blue',1,smval);%reactive sequence 2
        if predicted_items(1,c)
            dofill(t,fixation_locked_firing{c,unit}(which_sequence == 1 & predicted(:,c)',:),'magenta',1,smval); %predictive sequence 1
        end
        if predicted_items(2,c)
            dofill(t,fixation_locked_firing{c,unit}(which_sequence == 2 & predicted(:,c)',:),'cyan',1,smval);%predictive sequence 2
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Analyzing Predictive vs Reactive LFP Data---%%%

%%%---First Remove line noise & it's harmonics---%%%
Fline       = 60; %removes harmonics so 120, 180, etc.
LFPchannels = find_desired_channels(hdr,data,'LFP');
Fs = hdr.Fs;

sequence_LFPdata = cell(4,num_trials);
for trial = 1:num_trials
    
    LFPdata = [];
    for channel = 1:length(LFPchannels)
        LFPdata = [LFPdata; data(LFPchannels(channel)).values{successful_sequence_trials(trial)}];
    end
    
    % taken directly from fieldtrip ft_preproc_dftfilter.m
    % determine the size of the data
    [Nchans, Nsamples] = size(LFPdata);
    
    % determine the largest integer number of line-noise cycles that fits in the data
    sel = 1:round(floor(Nsamples * Fline/Fs) * Fs/Fline);
    
    % temporarily remove mean to avoid leakage
    mdat = mean(LFPdata(:,sel),2);
    zeroedLFPdat  = LFPdata - mdat(:,ones(1,Nsamples));
    
    % fit a sin and cos to the signal and subtract them
    time  = (0:Nsamples-1)/Fs;
    tmp  = exp(j*2*pi*Fline*time);                    % complex sin and cos
    % ampl = 2*dat*tmp'/Nsamples;                  % estimated amplitude of complex sin and cos
    ampl = 2*zeroedLFPdat(:,sel)*tmp(sel)'/length(sel);     % estimated amplitude of complex sin and cos on integer number of cycles
    est  = ampl*tmp;                               % estimated signal at this frequency
    %filt = dat - est;                              % subtract estimated signal
    LFPdata = LFPdata - est; %do want to remove the mean
    LFPdata = real(LFPdata);
    LFPdata = LFPdata./(std(LFPdata')'*ones(1,length(LFPdata))); %normalize since impedance changes absolute power
    
    for channel = 1:4;
        sequence_LFPdata{channel,trial} = LFPdata(channel,:);
    end
end

%%%---Align LFPs to Saccades and Fixations---%%%
%channel data alternate every 4th row
fixation_aligned_LFP = cell(1,4);
saccade_aligned_LFP = cell(1,4);
predicted_vs_reactive = NaN(4,4*num_trials);
sequence = NaN(4,4*num_trials);
for c = 1:4
    fixation_aligned_LFP{c} = NaN(4*num_trials,1024); %take LFP 512 ms before and after fixation
    saccade_aligned_LFP{c} = NaN(4*num_trials,1024); %take LFP 512 ms before and after saccade
    for trial = 1:num_trials
        fixt = fixation_start_time(trial,c);
        sact = saccade_start_time(trial,c);
        for channel = 1:4
            if fixt > 512
                fixation_aligned_LFP{c}(4*(trial-1)+channel,:) = sequence_LFPdata{channel,trial}(fixt-512:fixt+511);
            end
            if sact > 512 && sact+511 <= length(sequence_LFPdata{channel,trial})
                saccade_aligned_LFP{c}(4*(trial-1)+channel,:) = sequence_LFPdata{channel,trial}(sact-512:sact+511);
            end
            predicted_vs_reactive(c,4*(trial-1)+channel) = predicted(trial,c);
            sequence(c,4*(trial-1)+channel) = which_sequence(trial);
        end
    end
end

saccade_predicted_vs_reactive_contrasts = cell(2,4);
for item = 1:4
    for seq = 1:2
        if predicted_items(seq,item)
            temp = saccade_aligned_LFP{item};
            temp = temp(sequence(item,:) == seq,:);
            temp_factor = predicted_vs_reactive(item,:);
            temp_factor = temp_factor(:,sequence(item,:) == seq)';
            saccade_predicted_vs_reactive_contrasts{seq,item} = LFP_wavelet_ANOVA(temp,temp_factor);
            figure(101)
            subplot(2,2,item)
            hold on
            if seq == 1
                plot(nanmean(temp(temp_factor == 1,:))-nanmean(temp(temp_factor == 0,:)),'blue')
                plot(saccade_predicted_vs_reactive_contrasts{seq,item},'red')
            else
                plot(nanmean(temp(temp_factor == 1,:))-nanmean(temp(temp_factor == 0,:)),'cyan')
                plot(saccade_predicted_vs_reactive_contrasts{seq,item},'magenta')
            end
        end
    end
    if sum(predicted_items(:,item)) > 0 %if at least 1 item was predicted
        if sum(predicted_items(:,item)) == 2 % both sequences predicted
            title(['Item # ' num2str(item) ' Seq #1 ' num2str(percent_predicted(1,item)) '% Seq #2 ' num2str(num2str(percent_predicted(2,item))) '%'])
        elseif predicted_items(1,item)
            title(['Item # ' num2str(item) ' Seq #1 ' num2str(percent_predicted(1,item)) '%'])
        elseif predicted_items(2,item)
            title(['Item # ' num2str(item) ' Seq #2 ' num2str(num2str(percent_predicted(2,item))) '%'])
        end
        xlabel('Time from Fixation (ms)')
        xlim([0 1024])
        set(gca,'Xtick',0:256:1024)
        set(gca,'XtickLabel',num2cell(-512:256:512))
        ylabel('Normalized LFP (uV)')
    else
        axis off
    end
end
subtitle(['Saccade-Locked LFP n =' num2str(num_trials)]);
save_and_close_fig(figure_dir,[preprocessed_data_file(1:end-11) '_Sequence_Saccade_locked_wfANOVA_LFP']);


fixation_predicted_vs_reactive_contrasts = cell(2,4);
for item = 1:4
    for seq = 1:2
        if predicted_items(seq,item)
            temp = fixation_aligned_LFP{item};
            temp = temp(sequence(item,:) == seq,:);
            temp_factor = predicted_vs_reactive(item,:);
            temp_factor = temp_factor(:,sequence(item,:) == seq)';
            fixation_predicted_vs_reactive_contrasts{seq,item} = LFP_wavelet_ANOVA(temp,temp_factor);
            figure(101)
            subplot(2,2,item)
            hold on
            if seq == 1
                plot(nanmean(temp(temp_factor == 1,:))-nanmean(temp(temp_factor == 0,:)),'blue')
                plot(fixation_predicted_vs_reactive_contrasts{seq,item},'red')
            else
                plot(nanmean(temp(temp_factor == 1,:))-nanmean(temp(temp_factor == 0,:)),'cyan')
                plot(fixation_predicted_vs_reactive_contrasts{seq,item},'magenta')
            end
        end
    end
    if sum(predicted_items(:,item)) > 0 %if at least 1 item was predicted
        if sum(predicted_items(:,item)) == 2 % both sequences predicted
            title(['Item # ' num2str(item) ' Seq #1 ' num2str(percent_predicted(1,item)) '% Seq #2 ' num2str(num2str(percent_predicted(2,c))) '%'])
        elseif predicted_items(1,item)
            title(['Item # ' num2str(item) ' Seq #1 ' num2str(percent_predicted(1,item)) '%'])
        elseif predicted_items(2,item)
            title(['Item # ' num2str(item) ' Seq #2 ' num2str(num2str(percent_predicted(2,item))) '%'])
        end
        xlabel('Time from Fixation (ms)')
        set(gca,'Xtick',0:256:1024)
        xlim([0 1024])
        set(gca,'XtickLabel',num2cell(-512:256:512))
        ylabel('Normalized LFP (uV)')
    else
        axis off
    end
end
subtitle(['Fixation-Locked LFP n =' num2str(num_trials)]);
save_and_close_fig(figure_dir,[preprocessed_data_file(1:end-11) '_Sequence_fixation_locked_wfANOVA_LFP']);

save([data_dir preprocessed_data_file(1:8) '-Eyemovement_Locked_Sequence_results.mat'],...
    'predicted_thresh','twin','smval','predicted_items','percent_predicted',...
    'saccade_locked_firing','fixation_locked_firing','fixation_start_time',...
    'saccade_start_time','reaction_time','time_to_fixation','fixation_accuracy',...
    'fixation_duration','fixation_on_test','extrafixations',...
    'saccade_predicted_vs_reactive_contrasts','fixation_predicted_vs_reactive_contrasts');
end

function [fixationtimes,fixations,saccadetimes] = remove_corrective_saccades(...
    xy,saccadetimes,saccade_amplitudes,fixwin)
% written by Seth Konig on June 26, 2014
% function removes "corrective saccades" determine as saccades (other than
% the 1st recorded saccade) that are shorter than a defined threshold.
% This threshold is set by the "saccade_amplitude" variable.

corrective_saccades = find(saccade_amplitudes <= fixwin/2);
corrective_saccades(corrective_saccades == 1) = []; %ignore these
saccadetimes(:,corrective_saccades) =[];%remove corrective saccades

saccade_indexes = [];
for sac = 1:size(saccadetimes,2)
    saccade_indexes = [saccade_indexes saccadetimes(1,sac):saccadetimes(2,sac)];
end

fixation_indexes = 1:size(xy,2);
[~ , ia, ~] = intersect(fixation_indexes,saccade_indexes);
fixation_indexes(ia) = [];

[fixationtimes,fixations] = BehavioralIndexXY(fixation_indexes,xy(1,:),xy(2,:));
end

function [behaviortime, behaviormean] = BehavioralIndexXY(behavind,x,y)
%function is the same as above but also calculates mean fixation position
dind = diff(behavind);
gaps =find(dind > 1);
behaveind = zeros(length(gaps),50);
if ~isempty(gaps)
    for gapind = 1:length(gaps)+1;
        if gapind == 1;
            temp = behavind(1:gaps(gapind));
        elseif gapind == length(gaps)+1
            temp = behavind(gaps(gapind-1)+1:end);
        else
            temp = behavind(gaps(gapind-1)+1:gaps(gapind));
        end
        behaveind(gapind,1:length(temp)) = temp;
    end
else
    behaveind =  behavind;
end
behaviortime = zeros(2,size(behaveind,1));
behaviormean = zeros(2,size(behaveind,1));
for index=1:size(behaveind,1)
    rowfixind = behaveind(index,:);
    rowfixind(rowfixind == 0) = [];
    behaviortime(:,index) = [rowfixind(1);rowfixind(end)];
    behaviormean(:,index) = [mean(x(rowfixind));...
        mean(y(rowfixind))];
end
end