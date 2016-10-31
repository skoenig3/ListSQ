% Population Saccadic Modulation
% written Seth Konig 8/12/16
%
% code imports data from List_Saccade_AnalysisV2 results and looks at how
% modulated the whole population is by saccades. Analysis includes...
% 1) Average saccade/fixation-locked firing rate
% 2) AP Differences
%   a) % of neurons significantly modulated by eye movements
%   b) degree of modulation e.g. % change in firing rate?
% 3) Relative timing of modulation i.e. lag
%   a) on average
%   b) distribution
clar
task = 'ListSQ';
Fs = 1000;%Hz

avg_saccade_firing_rates = []; %significant firing rates
all_firing_rates = []; %significant + non-significant
avg_saccade_firing_rates_limited =[];%for significant ones limted to duration of saccade+fixation
all_not_normalized = []; %significant firing rate not-normalized
all_mean_not_normalized = [];%significant firing rate but only mean subtracted
all_abs_normalized = [];%significant firing rate without flipping

twin1 = 200;% how much time to take before eye movement starts, in case neurons are encoding upcomming eye movement
twin2 = 400;%how much time to take after eye movement has started

sig_count = zeros(1,3);


max_firing = NaN(350,2);
monkey_count = zeros(2,2);
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Population Figures\Place Cell Eye Movements\';
for monkey = 2:-1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for session = 1:length(session_data)
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{session},task);
        if isempty(task_file)
            continue
        end
        disp(num2str(session))
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials','cfg','hdr','data','fixationstats');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        
        if exist([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat']);
            load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat']);
        else
            continue
        end
        
        for unit = 1:num_units
            if ~isempty(saccade_locked_firing{unit})
                if all(temporal_info.saccade.shuffled_temporalstability_prctile(:,unit) > 95) %significant stability and reliability
                    col = 1;%significant modulation
                    sig_count(1) = sig_count(1)+1;
                elseif temporal_info.saccade.shuffled_temporalstability_prctile(1,unit) > 95 %significant stability
                    col = 2;%significant modulation
                    sig_count(2) = sig_count(2)+1;
                else %not saccade modulated
                    sig_count(3) = sig_count(3)+1;
                    col = 3;%not saccade modulated
                end
                
                %spatial signficant
                %                 if ((spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.shuffled_spatialstability_prctile(unit) > 95)) ...
                %
                %                 if (fixation_info.rate(unit) >  fixation_info_95.rate(unit)) &&  (saccade_info.rate(unit) >  saccade_info_95.rate(unit))
                %                     place_sig_count(1,1) = place_sig_count(1,1)+1;
                %                 elseif (fixation_info2.rate(unit) >  fixation_info2_95.rate(unit)) &&  (saccade_info2.rate(unit) >  saccade_info2_95.rate(unit))
                %                     place_sig_count(2,1) = place_sig_count(2,1)+1;
                %                 elseif (fixation_info.rate(unit) >  fixation_info_95.rate(unit)) || (saccade_info.rate(unit) >  saccade_info_95.rate(unit))
                %                     place_sig_count(1,3) = place_sig_count(1,3)+1;
                %                 elseif (fixation_info2.rate(unit) >  fixation_info2_95.rate(unit)) || (saccade_info2.rate(unit) >  saccade_info2_95.rate(unit))
                %                     place_sig_count(2,3) = place_sig_count(2,3)+1;
                %                 else
                %                     place_sig_count(1,2) = place_sig_count(1,2)+1;
                %                 end
                %                 end
                
                saccade_firing = saccade_locked_firing{unit}(saccade_information{unit}(:,4) > image_on_twin,:);
                saccade_firing(nansum(saccade_firing,2) == 0,:) = [];%remove saccades without spikes
                [firing_rate,smoothed_trials]= nandens2(saccade_firing,30,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                
                
                if col == 1 || col == 2
                    %want to look only looking at window around saccade/fixation
                    info = saccade_information{unit}((saccade_information{unit}(:,4) > image_on_twin),:);
                    limited_firing = NaN(size(saccade_firing));
                    info(info(:,7) > twin2,7) = twin2; %if next fixation duration is > twin set to twin
                    for f = 1:size(limited_firing,1)
                        ind = twin1:twin1+info(f,7);
                        limited_firing(f,ind) = smoothed_trials(f,ind);
                    end
                    limited_firing(nansum(limited_firing,2) == 0,:) = [];
                    num_not_nans = sum(~isnan(limited_firing));%will be fewer sikes
                    not_enough_samples = find(num_not_nans < .5*size(limited_firing,1)); %median duration
                     limited_firing(:,not_enough_samples) = NaN;  %remove any time points with less than 1/4 of the samples on avg ~< 500 trials

                    limited_firing_rate = nanmean(limited_firing);
%                    num_not_nans = sum(~isnan(limited_firing));%will be fewer sikes

% %                     %plot firing rate curve over time for spikes limited time period of 1 fixation around saccade

% %                     limited_firing(:,not_enough_samples) = NaN;  %remove any time points with less than 1/4 of the samples on avg ~< 500 trials
% %                     [limited_firing_rate,~]= nandens2(limited_firing,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
%                     limited_firing_rate = nandens3(limited_firing,smval,Fs);
% %                     limited_firing_rate(not_enough_samples) = NaN;
% 
                    limited_firing_rate = limited_firing_rate-nanmean(limited_firing_rate);
                    if abs(nanmin(limited_firing_rate)) > nanmax(limited_firing_rate); %trough is larger than peak
                        limited_firing_rate = limited_firing_rate./nanmin(limited_firing_rate);
                    else
                        limited_firing_rate = limited_firing_rate/nanmax(limited_firing_rate);
                    end
                    avg_saccade_firing_rates_limited = [avg_saccade_firing_rates_limited; limited_firing_rate];
                    
                    
                    all_not_normalized = [all_not_normalized; firing_rate]; %significant firing rate not-normalizd
                    firing_rate = firing_rate-mean(firing_rate); %remove mean firing rate
                    all_mean_not_normalized = [all_mean_not_normalized; firing_rate];%significant firing rate but only mean subtracted
                    
                    all_abs_normalized = [all_abs_normalized; firing_rate/max(abs(firing_rate))];%significant firing rate without flipping
                    
                    firing_rate = firing_rate-mean(firing_rate);
                    if abs(min(firing_rate)) > max(firing_rate) %trough is larger than peak
                        firing_rate = firing_rate/min(firing_rate);
                    else
                        firing_rate = firing_rate/max(firing_rate);
                    end
                    avg_saccade_firing_rates = [avg_saccade_firing_rates; firing_rate];
                    monkey_count(1,monkey) = monkey_count(1,monkey)+1;
                else
                    firing_rate = firing_rate-mean(firing_rate);
                    if abs(min(firing_rate)) > max(firing_rate) %trough is larger than peak
                        firing_rate = firing_rate/min(firing_rate);
                    else
                        firing_rate = firing_rate/max(firing_rate);
                    end
                    all_firing_rates = [all_firing_rates; firing_rate];
                end
                monkey_count(2,monkey) = monkey_count(2,monkey)+1;
            end
        end
    end
end
%%
t = -twin1:twin2-1;
figure
subplot(2,3,1)
hold on
plot([t(1) t(end)],[0 0],'k')
plot([0 0],[-1 1],'k--')
plot(t,mean([all_firing_rates; avg_saccade_firing_rates]))
hold off
xlabel('Time from Saccade Start (ms)')
ylabel('Normalized Firing Rate')
title('All Units')

subplot(2,3,4)
hold on
plot([t(1) t(end)],[0 0],'k')
plot([0 0],[-1 1],'k--')
plot(t,mean(avg_saccade_firing_rates))
hold off
xlabel('Time from Saccade Start (ms)')
ylabel('Normalized Firing Rate')
title('All Significant Units-Normalized to Max Deviation')

subplot(2,3,2)
plot(t,mean(all_not_normalized))
yl = ylim;
hold on
plot([t(1) t(end)],[mean(mean(all_not_normalized)) mean(mean(all_not_normalized))],'k')
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlabel('Time from Saccade Start (ms)')
ylabel('Firing Rate (Hz)')
title('All Significant Units-Not Normalized')

subplot(2,3,3)
plot(t,mean(all_mean_not_normalized))
yl = ylim;
hold on
plot([t(1) t(end)],[0 0],'k')
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlabel('Time from Saccade Start (ms)')
ylabel('Firing Rate (Hz) from Mean')
title('All Significant Units-Only Mean Subtracted')

subplot(2,3,5)
hold on
plot([t(1) t(end)],[0 0],'k')
plot([0 0],[-1 1],'k--')
plot(t,mean(all_abs_normalized))
hold off
xlabel('Time from Saccade Start (ms)')
ylabel('Normalized Firing Rate')
title('All Significant Units-Normalized to Abs Max')


%clustered firing rate
[U,S,V] = pca(avg_saccade_firing_rates);
T = kmeans(U,4);
subplot(2,3,6)
hold all
plot([0 0],[-1 1],'k--')
plot([-twin1 twin2],[0 0],'k')
for i = 1:max(T)
    if sum(T == i) > 1
        plot(t,mean(avg_saccade_firing_rates(T == i,:)))
    else
        plot(t,avg_saccade_firing_rates(T == i,:));
    end
end
hold off

subtitle('Population Averages')
%% Significant cells full time period
figure
[m,i] = max(avg_saccade_firing_rates,[],2);
[mm,ii] = sort(i);
imagesc(t,[1:size(avg_saccade_firing_rates,1)],avg_saccade_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(avg_saccade_firing_rates,1)],'w--');
hold off
xlim([-twin1 twin2])
xlabel('Time from Saccade Start')
ylabel('Neuron #')
%% Significant Absolute Normalization
figure
[m,i] = max(all_abs_normalized,[],2);
[mm,ii] = sort(i);
imagesc(t,[1:size(all_abs_normalized,1)],all_abs_normalized(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_abs_normalized,1)],'w--');
hold off
xlim([-twin1 twin2])
xlabel('Time from Saccade Start')
ylabel('Neuron #')
%%
%% Significant Absolute Normalization
figure
[m,i] = max(all_abs_normalized,[],2);
[mm,ii] = sort(i);
imagesc(t,[1:size(all_abs_normalized,1)],all_abs_normalized(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_abs_normalized,1)],'w--');
hold off
xlim([-twin1 twin2])
xlabel('Time from Saccade Start')
ylabel('Neuron #')
%% Significant Limited
figure
[m,i] = max(avg_saccade_firing_rates_limited,[],2);
[mm,ii] = sort(i);
imagesc(t,[1:size(avg_saccade_firing_rates_limited,1)],avg_saccade_firing_rates_limited(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(avg_saccade_firing_rates_limited,1)],'w--');
hold off
xlim([-twin1 twin2])
xlabel('Time from Saccade Start')
ylabel('Neuron #')
caxis([-1 1])