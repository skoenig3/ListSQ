function List_Saccade_AnalysisV2(data_dir,figure_dir,session_data)
% written by Seth Konig December, 2014
% updated by Seth Konig April 7, 2015 to include estimating temporal
% information locked to saccades and fixations on each item.
% updated SDK 1/11/17 to handlde new format and partial session data for
% vaild trials only.
%
% Function analyizes spike times correlated with eye movements in the
% List image portion of the ListSQ task.
%
% Inputs:
%   1) data_dir: directory where preprocessed data is located
%   2) figure_dir: location of where to put save figures
%   3) session_data: contains relavent session information
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-Eyemovement_Locked_List_results'

figure_dir = [figure_dir 'List Saccade Analysis\'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set important Analysis parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%looking at all potential visual place cells has mean delay to peak of ~150 ms
%with very few happening befre fixation (so start of saccade) and ending by 400 ms for sure
twin1 = 100;% how much time to take before eye movement starts, in case neurons are encoding upcomming eye movement
twin2 = 400;%how much time to take after eye movement has started
image_on_twin = 500;%how much time to ignore eye movements for to remove strong visual response though some may last longer
trial_start_code = 15;
trial_end_code = 20;
imgon_code = 23;
imgoff_code = 24;
smval = 30 ;%gaussian 1/2 width for smoothing, for the moment want to smooth at high frequency modulations
numshuffs = 1000; %recommend this is between 100 & 1000, for bootstrapping to
% dtermine significance
info_type = 'temporal';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)

min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Load important Session Data and Information---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task = 'ListSQ';
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
if isempty(task_file)
    disp('No ListSQ file could be found. Exiting function...')
    return
end
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials','hdr','fixationstats');

%grab unit data
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

if num_units == 0
    disp([task_file(1:8) ': no units could be found. Exiting function...'])
    return %since no units skip analysis
end

%---determine number of blocks unit was stable for
%remove units with too few trials
%these are the absolute minimum data required to do data analysis
valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end

%set the image duration
if str2double(task_file(3:8)) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end
Fs = data(1).fsample; %should be 1000

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Get successful trials---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([task_file ': running saccade modulation anlaysis..'])

%get important task specific information
[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[~,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,imgon_code);


%preallocate space and parallel structure of cfg
num_trials = length(cfg.trl);
image_trials = zeros(1,num_trials);
for t = 1:num_trials %only take trials in which image was actually shown
    if any(cfg.trl(t).allval == 23) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end);
        image_trials(t) = 1;
    end
end

%remove excess data associated with non image trials
image_trials = find(image_trials);
num_trials = length(image_trials);

%only take the eye data from successful trials
fixationstats = fixationstats(image_trials);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Process eye data locked to Eye Movements---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%trial #, when did eye movement start within trial, ordinal #, eye movement start relative to image onset
saccade_info = NaN(10,25*length(fixationstats));

sac_ind = 1;
for t = 1:num_trials %will use all trials since may be useful to look at eye movements later for all trials
    fixationtimes = fixationstats{t}.fixationtimes;
    saccadetimes = fixationstats{t}.saccadetimes;
    xy = fixationstats{t}.XY;
    
    trial_start = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == trial_start_code);
    trial_end = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == trial_end_code);
    imgon = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == imgon_code)-trial_start;
    imgoff = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == imgoff_code)-trial_start;
    nvr = novel_vs_repeat(img_cnd == cfg.trl(image_trials(t)).cnd); %novel or repeat image
    
    if any(isnan(nvr))%presentation error so skip trial, see get_image_numbers.m
        continue
    end

    % if monkey isn't paying attention data isn't probably
    % worth much plus have to cut off somewhere
    if imgoff-imgon > 1.5*imgdur-1
        imgoff = imgon+1.5*imgdur-1;
    end
    
    %find fiations and saccades that did not occur during the image period;
    %should also take care of the 1st fixation on the crosshair
    
    %fixation started before image turned on
    invalid= find(fixationtimes(1,:) < imgon);
    fixationtimes(:,invalid) = [];
    
    %fixation started after the image turned off and/or firing rate could corrupted by image turning off
    invalid= find(fixationtimes(1,:) > imgoff-twin2);
    fixationtimes(:,invalid) = [];
    
    %saccade started before image turned on
    invalid= find(saccadetimes(1,:) < imgon);
    saccadetimes(:,invalid) = [];
    
    %saccade started after the image turned off and/or firing rate could corrupted by image turning off
    invalid= find(saccadetimes(1,:) > imgoff-twin2);
    saccadetimes(:,invalid) = [];
    
    for s = 1:size(saccadetimes,2);
        next_fix = find(fixationtimes(1,:) == saccadetimes(2,s)+1);%next fixation should start immediately after
        if isempty(next_fix) %trial ended or eye outside of image
            continue %try next one
        end
        sacamp = sqrt(sum((xy(:,saccadetimes(2,s))-xy(:,saccadetimes(1,s))).^2)); %saccade amplitude
        next_fix_dur = fixationtimes(2,next_fix)-fixationtimes(1,next_fix)+1;%next fixation duration
        if sacamp >= min_sac_amp && next_fix_dur >= min_fix_dur %next fixation has to be long enough & saccade large enough
            saccade_info(1,sac_ind) = t; %trial #
            saccade_info(2,sac_ind) = saccadetimes(1,s); %saccade start time from trial start
            saccade_info(3,sac_ind) = s; %ordinal saccade #
            saccade_info(4,sac_ind) = saccadetimes(1,s)-imgon; %saccade start relative to image onset
            saccade_info(5,sac_ind) = sacamp; %saccade amplitude
            saccade_info(6,sac_ind) = next_fix_dur;%duration of next fixation
            saccade_info(7,sac_ind) = fixationtimes(2,next_fix)-saccadetimes(1,s);%end of next fixation from saccade start
            saccade_info(8,sac_ind) = atan2d(xy(2,saccadetimes(2,s))-xy(2,saccadetimes(1,s)),xy(1,saccadetimes(2,s))-xy(1,saccadetimes(1,s)));%saccade_direction
            prior_fix = find(fixationtimes(2,:) == saccadetimes(1,s)-1);%find fixation before in case want to know
            if ~isempty(prior_fix)
                saccade_info(9,sac_ind) = saccadetimes(1,s)-fixationtimes(1,prior_fix);%how long can look back from saccade start may not use
            end
            saccade_info(10,sac_ind) = nvr; %novel vs repeat images
            sac_ind = sac_ind+1;
        end
    end
end

%Remove excess NaNs;
saccade_info = laundry(saccade_info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Firing Rate Locked to Eye Movements---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pre-allocate space for variables
saccade_locked_firing = cell(1,num_units); %firing rate locked to saccade
saccade_information = cell(1,num_units); %start time relative to image onset and ordinal #
for unit = 1:num_units
    saccade_locked_firing{unit}= NaN(size(saccade_info,2),twin1+twin2);
    saccade_information{unit} = NaN(size(saccade_info,2),10);
end

sac_ind = ones(1,num_units);
for trial = 1:num_trials
    for unit = 1:num_units;
        if image_trials(trial) >= valid_trials(1,unit) && image_trials(trial) <= valid_trials(2,unit) %only valid trials
            saccades = find(saccade_info(1,:) == trial);
            spikes = find(data(unit).values{image_trials(trial)});
            
            %collect spike times relative to saccade start
            for s = 1:length(saccades);
                sact = saccade_info(2,saccades(s));
                sac_spikes = spikes(spikes > sact-twin1 & spikes <= sact+twin2)-sact+twin1;
                temp = zeros(1,twin1+twin2);
                temp(sac_spikes) = 1;
                saccade_locked_firing{unit}(sac_ind(unit),:) = temp;
                saccade_information{unit}(sac_ind(unit),:) = saccade_info(:,saccades(s))';
                sac_ind(unit) = sac_ind(unit)+1;
            end
        end
    end
end
%remove excess NaNs
saccade_locked_firing = laundry(saccade_locked_firing);
saccade_information = laundry(saccade_information);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Mutual Info for Eye Movements and Firing Rate---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temporal_info.saccade = [];
temporal_info.saccade.rate = NaN(1,num_units);
temporal_info.saccade.shuffled_rate = cell(1,num_units);
temporal_info.saccade.shuffled_rate_prctile = NaN(1,num_units);
temporal_info.saccade.temporalstability = NaN(2,num_units);  %row 1 by half, row 2 even/odd trials
temporal_info.saccade.shuffled_temporalstability =cell(1,num_units);
temporal_info.saccade.shuffled_temporalstability_prctile = NaN(2,num_units);%row 1 by half, row 2 even/odd trials


for unit = 1:num_units
    if ~isempty(saccade_locked_firing{unit})
              
        %don't want to run on trials with the first eye movements occuring with 500 (twin) ms of image onset
        saccade_firing = saccade_locked_firing{unit}(saccade_information{unit}(:,4) > image_on_twin,:);
        
        %since many units appear sparse especially spatial ones want to
        %run on eye movemetns that actually have spikes so use these only
        %doesn't change statistics just decreases run time
        saccade_firing(nansum(saccade_firing,2) == 0,:) = [];%remove saccades without spikes
        
        [observed_info_rate,shuffled_info_rate]...
            = estimated_mutual_information(saccade_firing,numshuffs,info_type,smval,Fs);
        temporal_info.saccade.rate(unit) = observed_info_rate.skaggs;
        temporal_info.saccade_shuffled_info.rate{unit} = shuffled_info_rate.skaggs;
        temporal_info.saccade.shuffled_rate_prctile(unit) = 100*sum(...
            observed_info_rate.skaggs > shuffled_info_rate.skaggs)/numshuffs;
        temporal_info.saccade.temporalstability(:,unit) = observed_info_rate.temporalstability;
        temporal_info.saccade_shuffled_info.temporalstability{unit} = shuffled_info_rate.temporalstability;
        temporal_info.saccade.shuffled_temporalstability_prctile(1,unit) = ...
            100*sum(observed_info_rate.temporalstability(1) > ...
            shuffled_info_rate.temporalstability(1,:))/numshuffs;
        temporal_info.saccade.shuffled_temporalstability_prctile(2,unit) = ...
            100*sum(observed_info_rate.temporalstability(2) > ...
            shuffled_info_rate.temporalstability(2,:))/numshuffs;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot Rasters by Various Variables of Interst---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = -twin1:twin2-1;
unit_names = unit_stats(1,:);
for unit = 1:num_units
    if ~isempty(saccade_locked_firing{unit})
        
        %---For Saccades---%
        %ignore anything within 500 ms of image onset
        info = saccade_information{unit}((saccade_information{unit}(:,4) > image_on_twin),:);
        saccade_firing = saccade_locked_firing{unit}(saccade_information{unit}(:,4) > image_on_twin,:);
        
        %remove trials without spikes, used for shuffling anlaysis above
        saccade_firing_plus = saccade_firing;
        saccade_firing_plus(nansum(saccade_firing,2) == 0,:) = [];
        
        %want to look only looking at window around saccade/fixation
        limited_firing = NaN(size(saccade_firing));
        info(info(:,7) > twin2,7) = twin2; %if next fixation duration is > twin set to twin
        for f = 1:size(limited_firing,1)
            ind = twin1:twin1+info(f,7);
            limited_firing(f,ind) = saccade_firing(f,ind);
        end
        limited_firing(nansum(limited_firing,2) == 0,:) = [];
        
        figure
        
        %plot firing rate curve over time
        hold on
        subplot(3,3,1)
        dofill(t,saccade_firing_plus(2:2:end,:),'blue',1,smval);%even trials
        dofill(t,saccade_firing_plus(1:2:end,:),'red',1,smval);%odd trials
        dofill(t,saccade_firing_plus,'black',1,smval); %all trials with spikes
        yl = ylim;
        if yl(1) < 0;
            yl(1) =0;
        end
        plot([0 0],[yl(1) yl(2)],'k--')
        xlim([-twin1 twin2])
        hold off
        legend('Even','Odd','All','Location','Best')
        ylim(yl);
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Saccade Start')
        set(gca,'Xtick',[-100 0 100 200 300 400])
        set(gca,'XMinorTick','on','YMinorTick','on')
        grid on
        ylim(yl);
        title(['Bit ' num2str(temporal_info.saccade.shuffled_rate_prctile(unit),3) '% '...
            '\rho_{1/2} = ' num2str(temporal_info.saccade.temporalstability(1,unit),2) ...
            ' (' num2str(temporal_info.saccade.shuffled_temporalstability_prctile(1,unit),3) ...
            '%) \rho_{e/o} = ' num2str(temporal_info.saccade.temporalstability(2,unit),2) ...
            ' (' num2str(temporal_info.saccade.shuffled_temporalstability_prctile(2,unit),3) '%)'])
        
        %plot raster over time by saccade occurence with session
        subplot(3,3,2)
        [trial,time] = find(saccade_firing_plus == 1);
        plot(time-twin1,(trial),'.k')
        ylim([0 max(trial)])
        ylabel('Occurence #')
        set(gca,'Xtick',[-100 0 100 200 300 400])
        xlabel('Time from Saccade Start')
        title('Raster whole period')
        
        %plot firing rate curve over time for spikes limited time period of 1 fixation around saccade
        num_not_nans = sum(~isnan(limited_firing));%will be fewer sikes
        not_enough_samples = find(num_not_nans < .5*size(limited_firing,1)); %median duration
%         limited_firing(:,not_enough_samples) = NaN;  %remove any time points with less than 1/4 of the samples on avg ~< 500 trials
%         [firing_rate,~]= nandens2(limited_firing,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
%         firing_rate(not_enough_samples) = NaN;  %remove any time points with less than 1/4 of the samples on avg ~< 500 trials
        limited_firing_rate = nandens3(limited_firing,smval,Fs);

        
        subplot(3,3,4)
        plot(t,limited_firing_rate,'k');
        hold on
        yl = ylim;
        if yl(1) < 0;
            yl(1) =0;
        end
        plot([0 0],[yl(1) yl(2)],'k--')
        hold off
        xlim([-twin1 twin2])
        ylim(yl);
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Saccade Start')
        set(gca,'Xtick',[-100 0 100 200 300 400])
        set(gca,'XMinorTick','on','YMinorTick','on')
        grid on
        ylim(yl);
        title('Firing Rate: 1 fixation & 1 saccade before')
        
        %make raster plot for spikes limited time period of 1 fixation around saccade
        subplot(3,3,5)
        [trial,time] = find(limited_firing == 1);
        plot(time-twin1,trial,'.k')
        xlim([-twin1 twin2])
        ylim([0 max(trial)])
        set(gca,'Xtick',[-100 0 100 200 300 400])
        ylabel('Occurence #')
        xlabel('Time from Saccade Start')
        title('Raster: 1 fixation & 1 saccade before')
        
        
        %plot raster over time by spikes/saccade epoch aka by average firing rate
        f1 = sum(saccade_firing,2);
        [~,fi] = sort(f1);
        [trial,time] = find(saccade_firing(fi,:) == 1);
        subplot(3,3,3)
        plot(time-twin1,(trial),'.k')
        ylim([0 max(trial)])
        set(gca,'Xtick',[-100 0 100 200 300 400])
        ylabel('Ranked Spike Count')
        xlabel('Time from Saccade Start')
        title('For all saccades including zero spike counts')
        
        %plot raster over time by saccade order within image presentation
        [~,fi] = sort(info(:,3));
        subplot(3,3,6)
        [trial,time] = find(saccade_firing(fi,:) == 1);
        plot(time-twin1,trial,'.k')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-100 0 100 200 300 400])
        ylabel('Ranked Ordinal Saccade #')
        xlabel('Time from Saccade Start')
        title('Discrete Time within Image Period')
        
        %plot raster by saccade amplitude
        [~,fi] = sort(info(:,5));
        subplot(3,3,7)
        [trial,time] = find(saccade_firing(fi,:) == 1);
        plot(time-twin1,trial,'.k')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-100 0 100 200 300 400])
        ylabel('Ranked Saccade Amplitude')
        xlabel('Time from Saccade Start')
        title('Saccade Amplitude')
        
        %plot raster by saccade direction
        [~,fi] = sort(info(:,8));
        subplot(3,3,8)
        [trial,time] = find(saccade_firing(fi,:) == 1);
        plot(time-twin1,trial,'.k')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-100 0 100 200 300 400])
        ylabel('Ranked Saccade Direction')
        xlabel('Time from Saccade Start')
        title('Saccade Direction')
        
        %plot raster by fixation duration
        [~,fi] = sort(info(:,7));
        subplot(3,3,9)
        [trial,time] = find(saccade_firing(fi,:) == 1);
        dur = info(fi,7);
        dur(dur > twin2) = twin2;
        plot(time-twin1,trial,'.k')
        hold on
        plot(dur,1:length(fi),'r.')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-100 0 100 200 300 400])
        ylabel('Ranked Fixation Duration')
        xlabel('Time from Saccade Start')
        title('Fixation Duration')
        
        n_str = ['   n_{sac} =' num2str(size(saccade_locked_firing{1,unit},1))];
        if multiunit(unit)
            subtitle(['Saccade-Locked Multiunit ' unit_names{unit} n_str]);
        else
            subtitle(['Saccade-Locked ' unit_names{unit} n_str]);
        end
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} 'Eye_Locked_analysis_Saccade_Rasters']);
    end
end

save([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat'],...
    'twin1','twin2','image_on_twin','smval','saccade_locked_firing',...
    'saccade_information','temporal_info','unit_names');
end