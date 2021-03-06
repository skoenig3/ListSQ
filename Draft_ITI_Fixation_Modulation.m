% Population Fixation Modulation
% written Seth Konig 8/12/16
%
% code imports data from List_fixation_AnalysisV results and looks at how
% modulated the whole population is by fixations.
% Code does the following
% 1) Summarizes fixation modulation for for place cell and non-place cells
% 2) Tracks AP location, unit counts, and which monkey (not currently used)
% 3) Asks how modulate whole population of all units is
% 6) Copies relevant figures for eye movment modulated cells to summary directory

clar %clear, clc
task = 'ListSQ';
Fs = 1000;%Hz
min_num_fix = 250; %at least 250 fixatoins with a certain duration to analyze for time limited to fixation duration
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
ITI_Codes = [15 16];
twin1 = 200;
twin2 = 400;

summary_directory = 'C:\Users\seth.koenig\Desktop\Significant Units\Eye Movement\'; %where to store significant units
if ~isdir(summary_directory)
    mkdir(summary_directory)
end

%---Misc. Parameters (Mostly Eye Movement Modulated Cells)---%
monkey_all_unit_count = zeros(2,2);%row 1 modulated row 2 non-modulated, column by monkey
all_multi_unit_count = zeros(1,2); %all_multi_unit_countunits
all_eye_cell_unit_names = {}; %eye modulated cell unit names
all_eye_cell_monkeys = []; %1s and 2s
eye_cell_AP_location = []; %AP location of recorded place cell
eye_cell_place_cell_status = [];%whether unit is place cell or not, 1 for place, 0 for non-place
eye_cell_direction_cell_status = []; %whether unit is directionally modulated or not, 1 for direction, 0 for non-direction
eye_cell_amplitude_cell_status = [];%whether unit is amplitude modulated or not, 1 for amplitude, 0 for non-amplitude


%---Fixation Algined Firing Rate Curves---%
avg_fixation_firing_rates = []; %significant firing rates
flipped_avg_fixation_firing_rates = []; %significant firing rates +/- max
all_firing_rates = []; %significant + non-significant
all_not_normalized = []; %significant firing rate not-normalized
all_mean_not_normalized = [];%significant firing rate but only mean subtracted
all_peaks = [];%significant
all_peaks2 = [];%not significant


all_ITI_firing = [];

monkeys = {'Vivian','Tobii'};
figure_dir = {};
for monk = 2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for sess data---%%%
    %only need to run when somethings changed or sesss have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir{1} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir{2} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    end
    
    for sess = 1:length(session_data)
        %read in task file data
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        %load task file data
        load([data_dir task_file(1:end-11) '-preprocessed.mat']);
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        if num_units == 0
            continue
        end
        
        num_trials = length(cfg.trl); %number of trials
        valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
        if all(isnan(valid_trials(1,:)))
            continue
        end
        
        %load spatial analysis data
        disp(task_file(1:8))
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'spatial_info')
        
        %load ege movement modulation analysis data
        load([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat']);
        
        %should also load direction modulation
        load([data_dir task_file(1:8) '-Saccade_Direction_and_Amplitude_Analysis.mat'])
        
        if smval ~= 30
            error('Smoothing Value (2xStd) does not match expectations!')
        end
        
        for unit = 1:num_units
            if isempty(temporal_info.fixation.shuffled_rate{unit})
                continue
            elseif ~isempty(temporal_info.fixation.shuffled_rate{unit}) && (length(temporal_info.fixation.shuffled_rate{unit}) < 1000)
                error('Should have 1000 shuffles') %make sure file has right number of shuffles
            end
            if ~isempty(spatial_info.shuffled_info_rate{unit})&& length(spatial_info.shuffled_info_rate{unit}) < 1000
                error('Should have 1000 shuffles') %make sure file has right number of shuffles
            elseif isempty(spatial_info.shuffled_info_rate{unit})
                continue
            end
            
            if (temporal_info.fixation.shuffled_temporalstability_prctile(1,unit) > 95) ... %significant stability
                    && (temporal_info.fixation.shuffled_rate_prctile(unit) > 95) % %skagg 95%+
                
                %---Misc. Parameters (Mostly Eye Movement Modulated Cells)---%
                monkey_all_unit_count(1,monk) = monkey_all_unit_count(1,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                all_eye_cell_unit_names = [all_eye_cell_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}]; %unit name
                all_eye_cell_monkeys = [all_eye_cell_monkeys monk]; %1s and 2s
                eye_cell_AP_location = [eye_cell_AP_location chamber_zero(1)+ session_data{sess}.location(1)]; %AP location of recorded cell
                
                %is unit spatially modulated
                if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                        && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                    eye_cell_place_cell_status = [eye_cell_place_cell_status 1]; %place cell
                else
                    eye_cell_place_cell_status = [eye_cell_place_cell_status 0]; %non-place cell
                end
                
                %is unit directionally modulated
                if mrls.all_saccades_shuffled_prctile(unit) > 95
                    eye_cell_direction_cell_status = [eye_cell_direction_cell_status 1];
                else
                    eye_cell_direction_cell_status = [eye_cell_direction_cell_status 0];
                end
                
                %is unit amplitude modulated
                if amplitude_correlations_percentile(unit) > 97.5
                    eye_cell_amplitude_cell_status = [eye_cell_amplitude_cell_status 1];
                else
                    eye_cell_amplitude_cell_status = [eye_cell_amplitude_cell_status 0];
                end
                
                %---Fixation Algined Firing Rate Curves---%
                info = fixation_information{unit}((fixation_information{unit}(:,4) > image_on_twin),:); %fixaion info after image onset
                fixation_firing = fixation_locked_firing{unit}(fixation_information{unit}(:,4) > image_on_twin,:); %fixation aligned spikes after image onset
                firing_rate = nandens(fixation_firing,30,'gauss',Fs,'nanflt');%fixation aligned firing rate
                
                all_not_normalized = [all_not_normalized; firing_rate]; %significant firing rate not-normalized
                all_mean_not_normalized = [all_mean_not_normalized; firing_rate-nanmean(firing_rate)];%significant firing rate but only mean subtracted
                
                firing_rate = firing_rate-nanmean(firing_rate);
                fr2 = firing_rate;
                if nanmax(fr2) < abs(nanmin(fr2))
                   fr2 = fr2/nanmin(fr2);
                else
                    fr2 = fr2/nanmax(fr2);
                end
                flipped_avg_fixation_firing_rates = [flipped_avg_fixation_firing_rates; fr2];
                
                firing_rate = firing_rate/nanmax(abs(firing_rate));
                avg_fixation_firing_rates = [avg_fixation_firing_rates; firing_rate]; %significant firing rates
                all_firing_rates = [all_firing_rates; firing_rate]; %significant + non-significant
                
            else
                %---Misc Parameters---%
                monkey_all_unit_count(2,monk) = monkey_all_unit_count(2,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                
                %---Fixation Algined Firing Rate Curves---%
                info = fixation_information{unit}((fixation_information{unit}(:,4) > image_on_twin),:); %fixaion info after image onset
                fixation_firing = fixation_locked_firing{unit}(fixation_information{unit}(:,4) > image_on_twin,:); %fixation aligned spikes after image onset
                firing_rate = nandens(fixation_firing,30,'gauss',Fs,'nanflt');%fixation aligned firing rate
                firing_rate = firing_rate-nanmean(firing_rate);
                firing_rate = firing_rate/nanmax(abs(firing_rate));
                all_firing_rates = [all_firing_rates; firing_rate]; %significant + non-significant
            end
            
            fixind = 1;
            ITI_firing = NaN(1000,twin1+twin2);
            for t=1:length(cfg.trl)
                ITI_end = cfg.trl(t).alltim(cfg.trl(t).allval == ITI_Codes(2))- cfg.trl(t).alltim(cfg.trl(t).allval == ITI_Codes(1));
                if isempty(fixationstats{t}.fixationtimes)
                    continue
                end
                if  length(data(unit).values{t}) < 2000
                   continue 
                end
                    
                ITI_fixations = find(fixationstats{t}.fixationtimes(2,:) <= ITI_end);
                spikes = data(unit).values{t};
               
                for itf = 1:length(ITI_fixations)
                    fixtimes = fixationstats{t}.fixationtimes(:,ITI_fixations(itf));
                    if fixtimes(1) <= twin1
                        if t == 1
                            continue
                        else
                             spikespre = data(unit).values{t-1};
                             leftover = twin1-fixtimes(1);
                             ITI_firing(fixind,:) = [spikespre(end-leftover+1:end) spikes(1:fixtimes(1)+twin2)];
                        end
                    else
                        ITI_firing(fixind,:) = spikes(fixtimes(1)-twin1:fixtimes(1)+twin2-1); 
                    end 
                    fixind = fixind+1;
                end
            end
            ITI_firing = laundry(ITI_firing);
            
            figure
            subplot(1,2,1)
            [trial,times] = find(ITI_firing == 1);
            plot(times-twin1,trial,'k.');
            xlabel('Time From Fixation Start (ms)')
            ylabel('Fixation #')
            if ~isempty(trial)
            ylim([0 max(trial)])
            end
            
            subplot(1,2,2)
            iti_firing_rate = nandens(ITI_firing,30,'gauss',Fs,'nanflt');%fixation aligned firing rate
            hold on
            plot(-twin1:twin2-1,iti_firing_rate);
            firing_rate = nandens(fixation_firing,30,'gauss',Fs,'nanflt');%fixation aligned firing rate
            plot(-twin1:twin2-1,firing_rate);
            hold off
            xlabel('Time From Fixation Start (ms)')
            ylabel('Firing Rate (Hz)')
            legend('ITI Eye Locked','Image Eye Locked')
            title(['Significant Correlation %: ' num2str(temporal_info.fixation.shuffled_temporalstability_prctile(1,unit),3) '%'])
            
            save_and_close_fig('C:\Users\seth.koenig\Desktop\ITIEye\',[task_file unit_stats{1,unit} '_ITI_Eye'])
            
            iti_firing_rate = iti_firing_rate-nanmean(iti_firing_rate(1:twin1));
            iti_firing_rate = iti_firing_rate/max(iti_firing_rate);
            all_ITI_firing = [all_ITI_firing; iti_firing_rate];
        end
    end
end
%%
disp(['Found ' num2str(length(all_eye_cell_unit_names)) ' eye movement modulated neurons!'])
disp([num2str(sum(eye_cell_place_cell_status == 1)) '/' num2str(length(all_eye_cell_unit_names)) ' are also place cells'])
disp([num2str(sum(eye_cell_direction_cell_status == 1)) '/' num2str(length(all_eye_cell_unit_names)) ' are also saccade direction cells'])
disp([num2str(sum(eye_cell_amplitude_cell_status == 1)) '/' num2str(length(all_eye_cell_unit_names)) ' are also saccade amplitude cells'])
disp([num2str(sum(eye_cell_place_cell_status == 0 & eye_cell_direction_cell_status == 0 & eye_cell_amplitude_cell_status == 0)) '/' num2str(length(all_eye_cell_unit_names)) ' are non-spatial'])

%%
figure
subplot(2,2,1)
plot(-twin1:twin2-1,nanmean(all_ITI_firing))
xlabel('Time from Fixation Start(ms)')
ylabel('Normalized Firing Rate')
title('Average of All Neurons')

subplot(2,2,2)
[m,i] = max(all_ITI_firing,[],2);
[mm,ii] = sort(i);
imagesc(-twin1:twin2-1,[1:size(all_ITI_firing,1)],all_ITI_firing(ii,:))
xlabel('Time from Fixation Start(ms)')
ylabel('Normalized Firing Rate')
title('Average of All Neurons')
vals = all_ITI_firing(:);
caxis([-nanstd(vals) 1])
colorbar


data = all_ITI_firing;
data(isnan(data(:,1)),:) = [];
[U,S,V] = pca(data,3);
T = kmeans(U,3);
t = -twin1:twin2-1;
subplot(2,2,3)
hold all
for i = 1:max(T)
    if sum(T == i) > 1
        plot(t,mean(data(T == i,:)))
    else
        plot(t,data(T == i,:));
    end
end
plot([0 0],[-1 1],'k--')
plot([-twin1 twin2],[0 0],'k')
xlabel('Time from Fixation Start (ms)')
ylabel('Normalized Firing')
title('PCA all neurons')
%%
t = -twin1:twin2-1;

figure

%---Plot Average Firing Rate For all Neurons---%
subplot(2,3,1)
hold on
plot([t(1) t(end)],[0 0],'k')
plot([0 0],[-1 1],'k--')
plot(t,mean([all_firing_rates; avg_fixation_firing_rates]))
hold off
xlabel('Time from fixation Start (ms)')
ylabel('Normalized Firing Rate')
title('All Units')

%---Plot Normalized Average Firing Rate for Only Eye Movement Modulated Neurons---%
subplot(2,3,4)
hold on
plot([t(1) t(end)],[0 0],'k')
plot([0 0],[-1 1],'k--')
plot(t,mean(avg_fixation_firing_rates),'b')
plot(t,mean(flipped_avg_fixation_firing_rates),'r')
hold off
xlabel('Time from fixation Start (ms)')
ylabel('Normalized Firing Rate')
title('Normalized Significant Units')

%---Plot Not-Normalized Average Firing Rate for Only Eye Movement Modulated Neurons---%
subplot(2,3,2)
plot(t,mean(all_not_normalized))
yl = ylim;
hold on
plot([t(1) t(end)],[mean(mean(all_not_normalized)) mean(mean(all_not_normalized))],'k')
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlabel('Time from fixation Start (ms)')
ylabel('Firing Rate (Hz)')
title('All Significant Units-Not Normalized')

%---Plot Not-Normalized Only Mean subtracted Average Firing Rate for Only Eye Movement Modulated Neurons---%
subplot(2,3,3)
plot(t,mean(all_mean_not_normalized))
yl = ylim;
hold on
plot([t(1) t(end)],[0 0],'k')
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlabel('Time from fixation Start (ms)')
ylabel('Firing Rate (Hz) from Mean')
title('All Significant Units-Only Mean Subtracted')


%clustered firing rate
[U,S,V] = pca(avg_fixation_firing_rates);
T = kmeans(U,4);
subplot(2,3,5)
hold all
plot([0 0],[-1 1],'k--')
plot([-twin1 twin2],[0 0],'k')
for i = 1:max(T)
    if sum(T == i) > 1
        plot(t,mean(avg_fixation_firing_rates(T == i,:)))
    else
        plot(t,avg_fixation_firing_rates(T == i,:));
    end
end
hold off
title('Clustered Curves')

subtitle('Population Averages')
%% All units
% all_all_peaks = [all_peak all_peaks2];
all = [all_firing_rates; avg_fixation_firing_rates];
figure
[m,i] = max(all,[],2);
%[mm,ii] = sort(all_all_peaks);
[mm,ii] = sort(i);
imagesc(t,[1:size(all,1)],all(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(all,1)],'w--');
hold off
xlim([-twin1 twin2])
xlabel('Time from fixation Start (ms)')
ylabel('Neuron #')

%% Significant cells full time period
figure
[mm,ii] = sort(all_peaks);
imagesc(t,[1:size(avg_fixation_firing_rates,1)],avg_fixation_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(avg_fixation_firing_rates,1)],'w--');
hold off
xlim([-twin1 twin2])
xlabel('Time from fixation Start')
ylabel('Neuron #')
caxis([-1 1])
%%
figure
[m,i] = max(avg_fixation_firing_rates(:,200:end),[],2);
[mm,ii] = sort(i);
imagesc(t,[1:size(avg_fixation_firing_rates,1)],avg_fixation_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(avg_fixation_firing_rates,1)],'w--');
hold off
xlim([-twin1 twin2])
xlabel('Time from fixation Start')
ylabel('Neuron #')
%% Significant Limited
%%
figure
[m,i] = max(avg_fixation_firing_rates_limited,[],2);
[mm,ii] = sort(i);
imagesc(t,[1:size(avg_fixation_firing_rates_limited,1)],avg_fixation_firing_rates_limited(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(avg_fixation_firing_rates_limited,1)],'w--');
hold off
xlim([-twin1 twin2])
xlabel('Time from fixation Start')
ylabel('Neuron #')
caxis([-1 1])
%%
%Non-spatial neurons only

fire_rates_non_spatial = avg_fixation_firing_rates(eye_cell_place_cell_status == 0 & ...
    eye_cell_direction_cell_status == 0 & eye_cell_amplitude_cell_status == 0,:);


flipped_rates_non_spatial = flipped_avg_fixation_firing_rates(eye_cell_place_cell_status == 0 & ...
    eye_cell_direction_cell_status == 0 & eye_cell_amplitude_cell_status == 0,:);

figure
subplot(1,2,1)
ylim([-1 1]);
yl = ylim;
hold on
plot([t(1) t(end)],[0 0],'k')
plot([0 0],[yl(1) yl(2)],'k--')
plot(t,mean(fire_rates_non_spatial),'b')
plot(t,mean(flipped_rates_non_spatial),'r')
hold off
xlabel('Time from fixation Start (ms)')
ylabel('Normalized Firing Rate (')
title('Average Firing Rate')
legend(' ',' ','Max Norm','Max/Min Norm')
axis square

subplot(1,2,2)
[m,i] = max(fire_rates_non_spatial(:,200:end),[],2);
[mm,ii] = sort(i);
imagesc(t,[1:size(fire_rates_non_spatial,1)],fire_rates_non_spatial(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(fire_rates_non_spatial,1)],'w--');
hold off
xlim([-twin1 twin2])
xlabel('Time from fixation Start')
ylabel('Neuron #')
subtitle([num2str(size(fire_rates_non_spatial,1)) ' Non-spatial Neurons'])
axis square