function cvtnew_time_locked_analysis(data_dir,figure_dir,session_data)
% originally written as part of time_locked_analysisV2 by Seth Konig August, 2014
% and updated SDK 1/11/17 to handlde new format and partial session data for
% vaild trials only for CVTNEW section on 1/19/16.
% Still too clunky so seperated again to have own cvtnew section.
% Function analyses spike times aligned to events that occur on the monitor
% at the time dicated by cortex event codes. Analysis does not analyze eye
% movements directly but when cortex says the eye had entered the fixation
% window. Updates to come on this.
%
% Inputs:
%   1) data_dir: direction where preprocessed data is located
%   2) figure_dir: location of where to put save figures
%   3) session_data: data containing relevant session data
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-time_locked_results'

task = 'cvtnew';
figure_dir = [figure_dir 'Cvtnew Time Analysis\'];
numshuffs = 1000; %number of shuffles to do for bootstrapping
min_blks = 1;%we will at least look at 1+ block data may later require 2
smval = 30; %15 ms standard deviation for shorter event durations
smval2 = 60; %30 ms standard deviation for longer event durations
twin1 = 200;% how much time to take before an event
twin2 = 1000;%how much time to look at after reward and ITI start
twin3 = 500;%minimum fixation on cross hair
twin4 = 2200; %for long window on image on,should max at 2100 but cortex could have some lag
twin5 = 500;%aligned to appearance of crosshair
twin6 = 700;%maximum response duration

%CVTNEW variable trial length
%time is + 300 ms
trial_lengths = [700 1133; ... %short
    1134 1567; ... %medium
    1568 2200]; %long

clrs = ['rgb'];

%Cortex codes
trial_start_code = 15;
cross_on_code = 35;
fixation_code = 8;
dot_on_code = 25;
dot_clrchng_code = 27;
response_code = 4;
reward_code = 3;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[task_file,~,~,multiunit,unit_names,unit_confidence,...
    sorting_quality,~,lfp_quality,~] = get_task_data(session_data,task);
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
    return %if no units exit function
end

%get trials in which spikes are valid i.e. trials for which neuron is
%stable for
valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import & reformat data so that spikes are locked to events---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([task_file(1:8) ' Aligning spike times to trial events'])

num_trials = length(cfg.trl);

ITI_start = [];
%preallocate space and parallel structure of cfg
time_lock_firing = cell(num_units,7);%event aligned spike trains
trial_duration = cell(1,num_units);
for unit = 1:num_units
    time_lock_firing{unit,1} = NaN(num_trials,twin1+twin5);%cross hair on
    time_lock_firing{unit,2} = NaN(num_trials,twin1+twin3);%fix crosshair
    time_lock_firing{unit,3} = NaN(num_trials,twin1+twin4);%dot on
    time_lock_firing{unit,4} = NaN(num_trials,twin1+twin6);%dot color change
    time_lock_firing{unit,5} = NaN(num_trials,twin5+twin1);%bar release
    time_lock_firing{unit,6} = NaN(num_trials,twin1+twin2);%reward period
    time_lock_firing{unit,7} = NaN(num_trials,twin1+twin2);%ITI period
    trial_duration{unit} = NaN(1,num_trials);
end

for t = 1:num_trials
    if any(cfg.trl(t).allval == reward_code); %in which trial was successful
        
        
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code); %start of ITI
        trial_end = cfg.trl(t).alltim(end)-trial_start; %end of trial after reward
        crosson = cfg.trl(t).alltim(cfg.trl(t).allval == cross_on_code)-trial_start; %cross hair on
        fix_cross = cfg.trl(t).alltim(cfg.trl(t).allval == fixation_code)-trial_start; %fixation on cross hair according to cortex
        doton = cfg.trl(t).alltim(cfg.trl(t).allval == dot_on_code)-trial_start; %when image turns on
        dotchange = cfg.trl(t).alltim(cfg.trl(t).allval == dot_clrchng_code)-trial_start; %when dot changes color
        monkey_response =  cfg.trl(t).alltim(cfg.trl(t).allval == response_code)-trial_start; %when release bar
        reward = cfg.trl(t).alltim(cfg.trl(t).allval == reward_code)-trial_start; %reward pulses
        reward = [reward reward(end)+200]; %200 ms inter-pulse interval, pulse duration actually ~115 ms
        
        for unit = 1:num_units
            if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only put data in for valid trials
                
                spikes = find(data(unit).values{t}); %spike trains for this trial
                trial_duration{unit}(t) = dotchange-doton;
                
                %---crosson locked firing---%
                event_spikes = spikes(spikes > crosson-twin1 & spikes <= crosson+twin5)-crosson+twin1;
                tempvec = zeros(1,twin1+twin5);
                tempvec(event_spikes) = 1;
                time_lock_firing{unit,1}(t,:) = tempvec;
                
                %---fixation on cross locked firing---%
                event_spikes = spikes(spikes > fix_cross-twin1 & spikes <= fix_cross+twin3)-fix_cross+twin1;
                tempvec = zeros(1,twin1+twin3);
                tempvec(event_spikes) = 1;
                time_lock_firing{unit,2}(t,:) = tempvec;
                
                
                %---dot on---%
                %variable duration event
                event_spikes = spikes(spikes > doton-twin1 & spikes <= dotchange)-doton+twin1;
                timevec = [doton-twin1+1 dotchange]-doton+twin1;
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                time_lock_firing{unit,3}(t,timevec(1):timevec(2)) = tempvec;
                
                %---dot color change---%
                %variable duration event
                if monkey_response-dotchange > twin6; %slight "lag" in cortex up to 10-20 ms depending on event
                    if monkey_response-dotchange  > 730
                        disp('what!?! "lag" is too much')
                    end
                    monkey_response = dotchange+twin6;
                end
                
                event_spikes = spikes(spikes > dotchange-twin1 & spikes <= monkey_response)-dotchange+twin1;
                timevec = [dotchange-twin1+1 monkey_response]-dotchange+twin1;
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                time_lock_firing{unit,4}(t,timevec(1):timevec(2)) = tempvec;
                
                %----Bar release---%
                %so can accurately look bach in time for anticipatory
                %signal since dot color changeevent period is variable
                event_spikes = spikes(spikes > monkey_response-twin5 & spikes <= monkey_response+twin1)-monkey_response+twin5;
                tempvec = zeros(1,twin5+twin1);
                tempvec(event_spikes) = 1;
                time_lock_firing{unit,5}(t,:) = tempvec;
                
                
                %----Reward Period---%
                event_spikes = spikes(spikes > reward(1)-twin1 & spikes <= reward(1)+twin2)-reward(1)+twin1;
                tempvec = zeros(1,twin1+twin2);
                tempvec(event_spikes) = 1;
                time_lock_firing{unit,6}(t,:) = tempvec;
                
                %---ITI period---%
                event_spikes = spikes(spikes > reward(end)-twin1 & spikes <= reward(end)+twin2)-reward(end)+twin1;%this trials spikes
                temp = zeros(1,twin1+twin2);
                leftover = twin2-(trial_end-reward(end));
                if t ~= num_trials
                    spks = find(data(unit).values{t+1}(1:leftover) == 1)+(trial_end-reward(end))+twin1;
                    event_spikes = [event_spikes spks];
                else
                    temp(end-leftover+1:end) = NaN;
                end
                temp(event_spikes) = 1;
                time_lock_firing{unit,7}(t,:) = temp;
                
                
                %
                %                 if t ~= length(cfg.trl)
                %                     %will need to look in next trial to get rest of data
                %                     %since usually only get ~200 ms of data
                %                     end_trial_gap = trial_end-reward(end);
                %                     rest_time = twin1+twin2-end_trial_gap;
                %                     spikes2 = find(data(unit).values{t+1}); %next trial's spike trains
                %                     spikes2(spikes2 > rest_time) = [];%remove next trial's spike trains spikes up to rest of time
                %                     spikes2 = spikes2+end_trial_gap;%reindex
                %                     event_spikes = [event_spikes spikes2]; %concatenate this plus next trial's spikes
                %                     tempvec =
                %                     tempvec(event_spikes) = 1;
                %                 end
            end
        end
    end
end

%remove excess NaNs associated with error trials
time_lock_firing = laundry(time_lock_firing);
trial_duration = laundry(trial_duration);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Perform Temporal analysis and signifance testing---%%%
disp([task_file(1:8) ': Determining if spikes are locked to trial events'])
Fs = data(1).fsample; %should be 1000

info_type = 'temporal_variable_trial_len';

epoch_data = [];
epoch_data.rate = NaN(num_units,size(time_lock_firing,2));
epoch_data.temporalstability = NaN(num_units,size(time_lock_firing,2));
epoch_data.shuffled_rate = cell(num_units,size(time_lock_firing,2));
epoch_data.shuffled_temporalstability = cell(num_units,size(time_lock_firing,2));
epoch_data.rate_prctile = NaN(num_units,size(time_lock_firing,2));
epoch_data.temporalstability_prctile = NaN(num_units,size(time_lock_firing,2));
for unit = 1:num_units
    for event = 1:size(time_lock_firing,2)
        if isempty(time_lock_firing{unit,event})
            continue
        elseif nansum(nansum(time_lock_firing{unit,event})) == 0; %no spikes
            epoch_data.rate(unit,event) = 0;
            epoch_data.rate_prctile(unit,event) = 0;
            epoch_data.temporalstability_prctile(unit,event) = 0;
        else
            if event == 3
                num_not_nans = sum(~isnan(time_lock_firing{unit,event}));
                indeces = find(num_not_nans > 25);
                temp_firing = time_lock_firing{unit,3}(:,indeces);
                
                [temporal_info,shuffled_info_rate] = ...
                    estimated_mutual_information(temp_firing,numshuffs,info_type,smval2,Fs);
                
            elseif event == 4 %aligned to do on or dot color change so variable duration
                num_not_nans = sum(~isnan(time_lock_firing{unit,event}));
                indeces = find(num_not_nans > 25);
                temp_firing = time_lock_firing{unit,3}(:,indeces);
                
                [temporal_info,shuffled_info_rate] = ...
                    estimated_mutual_information(temp_firing,numshuffs,info_type,smval,Fs);
            elseif event == 6 || event == 7 %ITI and reward
                [temporal_info,shuffled_info_rate] = ...
                    estimated_mutual_information(time_lock_firing{unit,event},numshuffs,info_type,smval2,Fs);
            else
                [temporal_info,shuffled_info_rate] = ...
                    estimated_mutual_information(time_lock_firing{unit,event},numshuffs,info_type,smval,Fs);
            end
            
            epoch_data.rate(unit,event) = temporal_info.skaggs;
            %only want 1st vs 2nd half no e/o since not using
            epoch_data.temporalstability(unit,event) = temporal_info.temporalstability(1);
            epoch_data.shuffled_rate{unit,event} = shuffled_info_rate.skaggs;
            epoch_data.shuffled_temporalstability{unit,event} = shuffled_info_rate.temporalstability(1,:);
            epoch_data.rate_prctile(unit,event) = 100*sum(temporal_info.skaggs > shuffled_info_rate.skaggs)/numshuffs;
            epoch_data.temporalstability_prctile(unit,event) = 100*sum(temporal_info.temporalstability(1) >...
                shuffled_info_rate.temporalstability(1,:))/numshuffs;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot and Save Figures of Results---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t15 = -twin1:twin5-1;
t13 = -twin1:twin3-1;
t14 = -twin1:twin4-1;
t16 = -twin1:twin6-1;
t51 = -twin5:twin1-1;
t12 = -twin1:twin2-1;
for unit = 1:num_units
    if isempty(time_lock_firing{unit,6})
        continue
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Reward and ITI figure---%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    yls = NaN(2,2);
    
    figure
    
    %---Plot Firing rate curve Aligned to Reward---%
    subplot(2,2,1)
    dofill(t12,time_lock_firing{unit,6},'black',1,smval2);
    xlim([-twin1 twin2])
    ylabel('Firing Rate (Hz)')
    xlabel('Time from ITI start (ms)')
    if epoch_data.rate_prctile(unit,6) > 90 || epoch_data.temporalstability_prctile(unit,6) > 90
        title(['bit ' num2str(epoch_data.rate_prctile(unit,6),3) '% ' ...
            '\rho_{1/2} = ' num2str(epoch_data.temporalstability(unit,6),2) ' ' ...
            num2str(epoch_data.temporalstability_prctile(unit,6),3)  '%' ])
    end
    box off
    yls(1,:) = ylim;
    
    %reward raster
    subplot(2,2,3)
    [trial,time] = find(time_lock_firing{unit,6} == 1);
    if ~isempty(trial)
        plot(time-twin1,trial,'.k')
        xlim([-twin1 twin2])
        ylim([0 max(trial)+1]);
    end
    ylabel('Trial #')
    xlabel('Time from Reward start (ms)')
    box off
    
    %---Plot Firing rate curve Aligned to ITI---%
    subplot(2,2,2)
    dofill(t12,time_lock_firing{unit,7},'black',1,smval2);
    xlim([-twin1 twin2])
    ylabel('Firing Rate (Hz)')
    xlabel('Time from ITI start (ms)')
    if epoch_data.rate_prctile(unit,7) > 90 || epoch_data.temporalstability_prctile(unit,7) > 90
        title(['bit ' num2str(epoch_data.rate_prctile(unit,7),3) '% ' ...
            '\rho_{1/2} = ' num2str(epoch_data.temporalstability(unit,7),2) ' ' ...
            num2str(epoch_data.temporalstability_prctile(unit,7),3)  '%' ])
    end
    box off
    yls(2,:) = ylim;
    
    %ITI raster
    subplot(2,2,4)
    [trial,time] = find(time_lock_firing{unit,7} == 1);
    if ~isempty(trial)
        plot(time-twin1,trial,'.k')
        xlim([-twin1 twin2])
        ylim([0 max(trial)+1]);
    end
    ylabel('Trial #')
    xlabel('Time from ITI start (ms)')
    box off
    
    %scale to same axis
    ymin = min(yls(:,1));
    ymin(ymin < 0) = 0;
    ymax = max(yls(:,2));
    
    subplot(2,2,1)
    ylim([ymin ymax])
    subplot(2,2,2)
    ylim([ymin ymax])
    
    n_str = [' n = ' num2str(size(time_lock_firing{unit,7},1))];
    if multiunit(unit)
        subtitle([' ' task_file(1:8) ' Multiunit ' unit_stats{1,unit} n_str]);
    else
        subtitle([' ' task_file(1:8) ' ' unit_stats{1,unit} n_str]);
    end
    save_and_close_fig(figure_dir,[task_file(1:8) '_' unit_stats{1,unit} '_cvtnew-time_locked_analysis_Reward_ITI']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Peri-Dot Aligned Figures---%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    yls = NaN(3,2);
    
    figure
    
    %---Plot Firing rate Aligned to Cross hair and Fixation on Cross---%
    subplot(2,3,1)
    hold on
    dofill(t15,time_lock_firing{unit,1},'black',1,smval);%cross on
    dofill(t13,time_lock_firing{unit,2},'blue',1,smval);%fix on cross
    hold off
    xlim([-twin1 twin5])
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Event (ms)')
    if epoch_data.rate_prctile(unit,2) > 90 || epoch_data.temporalstability_prctile(unit,2) > 90
        title(['bit_{fix} ' num2str(epoch_data.rate_prctile(unit,2),3) '% ' ...
            '\rho_{fix} = ' num2str(epoch_data.temporalstability(unit,2),2) ' ' ...
            num2str(epoch_data.temporalstability_prctile(unit,2),3)  '%' ])
    end
    legend('Cross On','Fix on Cross')
    box off
    yls(1,:) = ylim;
    
    %cross appearannce raster
    subplot(2,3,4)
    [trial,time] = find(time_lock_firing{unit,1} == 1);
    if ~isempty(trial)
        plot(time-twin1,trial,'.k')
        xlim([-twin1 twin5])
        ylim([0 max(trial)+1]);
    end
    ylabel('Trial #')
    xlabel('Time from Cross on (ms)')
    box off
    
    
    %---Plot Firing rate Aligned to Dot Color Change---%
    num_not_nans = sum(~isnan(time_lock_firing{unit,4}));
    indeces = find(num_not_nans > 25);
    
    subplot(2,3,2)
    hold on
    dofill2(indeces-twin1,time_lock_firing{unit,4}(:,indeces),'black',1,smval);
    hold off
    xlim([-twin1 max(indeces)-twin1])
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Dot Color change (ms)')
    if epoch_data.rate_prctile(unit,4) > 90 || epoch_data.temporalstability_prctile(unit,4) > 90
        title(['bit ' num2str(epoch_data.rate_prctile(unit,4),3) '% ' ...
            '\rho_{1/2} = ' num2str(epoch_data.temporalstability(unit,4),2) ' ' ...
            num2str(epoch_data.temporalstability_prctile(unit,4),3)  '%' ])
    end
    box off
    yls(3,:) = ylim;
    
    %color change raster
    subplot(2,3,5)
    [trial,time] = find(time_lock_firing{unit,4} == 1);
    if ~isempty(trial)
        plot(time-twin1,trial,'.k')
        xlim([-twin1 max(indeces)-twin1])
        ylim([0 max(trial)+1]);
    end
    ylabel('Trial #')
    xlabel('Time from Cross on (ms)')
    box off
    
    %---Plot Firing rate Aligned to Bar Release---%
    subplot(2,3,3)
    hold on
    dofill(t51,time_lock_firing{unit,5},'black',1,smval);
    hold off
    xlim([-twin5 twin1])
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Bar Release (ms)')
    if epoch_data.rate_prctile(unit,5) > 90 || epoch_data.temporalstability_prctile(unit,5) > 90
        title(['bit ' num2str(epoch_data.rate_prctile(unit,5),3) '% ' ...
            '\rho_{1/2} = ' num2str(epoch_data.temporalstability(unit,5),2) ' ' ...
            num2str(epoch_data.temporalstability_prctile(unit,5),3)  '%' ])
    end
    box off
    yls(4,:) = ylim;
    
    %color change raster
    subplot(2,3,6)
    [trial,time] = find(time_lock_firing{unit,5} == 1);
    if ~isempty(trial)
        plot(time-twin5,trial,'.k')
        xlim([-twin5 twin1])
        ylim([0 max(trial)+1]);
    end
    ylabel('Trial #')
    xlabel('Time from Bar Release (ms)')
    box off
    
    %scale to same axis
    ymin = min(yls(:,1));
    ymin(ymin < 0) = 0;
    ymax = max(yls(:,2));
    
    subplot(2,3,1)
    ylim([ymin ymax])
    subplot(2,3,2)
    ylim([ymin ymax])
    subplot(2,3,3)
    ylim([ymin ymax])
    
    n_str = [' n = ' num2str(size(time_lock_firing{unit,7},1))];
    if multiunit(unit)
        subtitle([' ' task_file(1:8) ' Multiunit ' unit_stats{1,unit} n_str]);
    else
        subtitle([' ' task_file(1:8) ' ' unit_stats{1,unit} n_str]);
    end
    
    save_and_close_fig(figure_dir,[task_file(1:8) '_' unit_stats{1,unit} '_cvtnew-peridot_analysis']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Dot Onset and by Trial Duration---%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure
    
    %---Plot Firing rate Aligned to Dot On  Across All trials---%
    num_not_nans = sum(~isnan(time_lock_firing{unit,3}));
    indeces = find(num_not_nans > 25);
    
    subplot(2,2,1)
    hold on
    dofill2(indeces-twin1,time_lock_firing{unit,3}(:,indeces),'black',1,smval2);
    hold off
    xlim([-twin1 twin4-twin1])
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Dot on (ms)')
    if epoch_data.rate_prctile(unit,3) > 90 || epoch_data.temporalstability_prctile(unit,3) > 90
        title(['bit ' num2str(epoch_data.rate_prctile(unit,3),3) '% ' ...
            '\rho_{1/2} = ' num2str(epoch_data.temporalstability(unit,3),2) ' ' ...
            num2str(epoch_data.temporalstability_prctile(unit,3),3)  '%' ])
    end
    box off
    yls(2,:) = ylim;
    
    %cross appearannce raster
    subplot(2,2,3)
    [trial,time] = find(time_lock_firing{unit,3} == 1);
    if ~isempty(trial)
        plot(time-twin1,trial,'.k')
        xlim([-twin1 twin4-twin1])
        ylim([0 max(trial)+1]);
    end
    ylabel('Trial #')
    xlabel('Time from Dot on (ms)')
    box off
    
    %---Plot Firing rate Aligned to Dot On  by Trial Duration---%
    [sorted_durs,sorted_index] = sort(trial_duration{unit});
    sorted_spikes = time_lock_firing{unit,3}(sorted_index,:);
    
    subplot(2,2,2)
    hold on
    for dur = 1:3
        these_trials = find(sorted_durs < trial_lengths(dur,2) & sorted_durs >= trial_lengths(dur,1));
        num_not_nans = sum(~isnan(sorted_spikes(these_trials,:)));
        indeces = find(num_not_nans > 25);
        dofill2(indeces-twin1,sorted_spikes(these_trials,indeces),clrs(dur),1,smval2);
    end
    hold off
    xlim([-twin1 twin4])
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Dot on (ms)')
    if epoch_data.rate_prctile(unit,3) > 90 || epoch_data.temporalstability_prctile(unit,3) > 90
        title(['bit ' num2str(epoch_data.rate_prctile(unit,3),3) '% ' ...
            '\rho_{1/2} = ' num2str(epoch_data.temporalstability(unit,3),2) ' ' ...
            num2str(epoch_data.temporalstability_prctile(unit,3),3)  '%' ])
    end
    box off
    
    %cross appearannce raster
    subplot(2,2,4)
    b4 = 0;
    hold on
    for dur = 1:3
        these_trials = find(sorted_durs < trial_lengths(dur,2) & sorted_durs >= trial_lengths(dur,1));
        [trial,time] = find(sorted_spikes(these_trials,:) == 1);
        if ~isempty(trial)
            trial = trial+b4;
            plot(time-twin1,trial,[clrs(dur) '.'])
            b4 = max(trial);
            ylim([0 max(trial)+1]);
        end
    end
    hold off
    xlim([-twin1 twin4])
    ylabel('Trial #')
    xlabel('Time from Dot on (ms)')
    box off
    
    n_str = [' n = ' num2str(size(time_lock_firing{unit,7},1))];
    if multiunit(unit)
        subtitle([' ' task_file(1:8) ' Multiunit ' unit_stats{1,unit} n_str]);
    else
        subtitle([' ' task_file(1:8) ' ' unit_stats{1,unit} n_str]);
    end
    
    save_and_close_fig(figure_dir,[task_file(1:8) '_' unit_stats{1,unit} '_cvtnew-dot_on_analysis']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Finally save all the data---%%%
save([data_dir task_file(1:8) '-cvtnew-time_locked_results.mat'],...
    'time_lock_firing','smval','smval2','temporal_info','epoch_data','unit_stats',...
    'twin1','twin2','twin3','twin4','twin5','twin6');
disp(['Time Locked Data Analyis for ' task_file(1:8) ' saved']);

end
