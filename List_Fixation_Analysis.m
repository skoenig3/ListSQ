function List_Fixation_Analysis(data_dir,figure_dir,session_data)
% modified from Finalish Version of List_Fixation_AnalysisV2 to only look at
% fixation aligned data instead of Fixation aligned.
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

figure_dir = [figure_dir 'List Fixation Analysis\'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set important Analysis parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%looking at all potential visual place cells has mean delay to peak of ~150 ms
%with very few happening befre fixation (so start of Fixation) and ending by 400 ms for sure
twin1 = 200;% how much time to take before eye movement starts, in case neurons are encoding upcomming eye movement
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
% info_type2 = 'temporal_eye';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)

min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro Fixations too small and hard to detect
min_num_fix = 250; %at least 250 fixatoins with a certain dureation to analyze

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
disp([task_file ': running Fixation modulation anlaysis..'])

%get important task specific information
[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[~,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,imgon_code);


%preallocate space and parallel structure of cfg
num_trials = length(cfg.trl);
image_trials = zeros(1,num_trials);
original_trial_num = zeros(1,num_trials);
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
fixation_info = NaN(10,25*length(fixationstats));

spike_times = cell(1,num_units);%image onset aligned data
for unit = 1:num_units
    spike_times{unit} = NaN(num_trials,imgdur*1.5);
end

fix_ind = 1;
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
    
    %find fiations and Fixations that did not occur during the image period;
    %should also take care of the 1st fixation on the crosshair
    
    %fixation started before image turned on
    invalid= find(fixationtimes(1,:) < imgon);
    fixationtimes(:,invalid) = [];
    
    %fixation started after the image turned off and/or firing rate could corrupted by image turning off
    invalid= find(fixationtimes(2,:) > imgoff);
    fixationtimes(:,invalid) = [];
    
    %Fixation started before image turned on
    invalid= find(saccadetimes(1,:) < imgon);
    saccadetimes(:,invalid) = [];
    
    %Fixation started after the image turned off and/or firing rate could corrupted by image turning off
    invalid= find(saccadetimes(2,:) > imgoff);
    saccadetimes(:,invalid) = [];
    
    for unit = 1:num_units
        if image_trials(t) >= valid_trials(1,unit) && image_trials(t) <= valid_trials(2,unit) %only valid trials
            spikes = find(data(unit).values{image_trials(t)});
            temp = zeros(1,imgoff-imgon+1);
            spikes(spikes <= imgon) = [];
            spikes(spikes > imgoff) = [];
            spikes = spikes-imgon;
            temp(spikes) = 1;
            spike_times{unit}(t,1:length(temp)) = temp;
        end
    end
    
    for f = 1:size(fixationtimes,2);
        prior_sac = find(fixationtimes(1,f) == saccadetimes(2,:)+1);%next fixation should start immediately after
        if isempty(prior_sac) %trial ended or eye outside of image
            continue %try next one
        end
        sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %Fixation amplitude
        fix_dur = fixationtimes(2,f)-fixationtimes(1,f)+1;%next fixation duration
        if sacamp >= min_sac_amp && fix_dur >= min_fix_dur %next fixation has to be long enough & Fixation large enough
            fixation_info(1,fix_ind) = t; %trial #
            fixation_info(2,fix_ind) = fixationtimes(1,f); %fixation start time from trial start
            fixation_info(3,fix_ind) = f; %ordinal fixation #
            fixation_info(4,fix_ind) = fixationtimes(1,f)-imgon; %fixation start relative to image onset
            fixation_info(5,fix_ind) = fix_dur;%duration of fixation
            fixation_info(6,fix_ind) = sacamp; %Fixation amplitude
            fixation_info(7,fix_ind) = saccadetimes(2,prior_sac)-saccadetimes(1,prior_sac)+1;%Saccade duration to know how far to go back
            fixation_info(8,fix_ind) = atan2d(xy(2,saccadetimes(2,prior_sac))-xy(2,saccadetimes(1,prior_sac)),...
                xy(1,saccadetimes(2,prior_sac))-xy(1,saccadetimes(1,prior_sac)));%Fixation_direction
            fixation_info(9,fix_ind) = nvr; %novel vs repeat images
            fix_ind = fix_ind+1;
        end
    end
end

%Remove excess NaNs;
fixation_info = laundry(fixation_info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Firing Rate Locked to Eye Movements---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pre-allocate space for variables
fixation_locked_firing = cell(1,num_units); %firing rate locked to Fixation
fixation_information = cell(1,num_units); %start time relative to image onset and ordinal #
for unit = 1:num_units
    fixation_locked_firing{unit}= NaN(size(fixation_info,2),twin1+twin2);
    fixation_information{unit} = NaN(size(fixation_info,2),9);
end

fix_ind = ones(1,num_units);
for trial = 1:num_trials
    for unit = 1:num_units;
        if image_trials(trial) >= valid_trials(1,unit) && image_trials(trial) <= valid_trials(2,unit) %only valid trials
            fixations = find(fixation_info(1,:) == trial);
            spikes = find(data(unit).values{image_trials(trial)});
            
            %collect spike times relative to Fixation start
            for f = 1:length(fixations);
                fixt = fixation_info(2,fixations(f));
                sac_spikes = spikes(spikes > fixt-twin1 & spikes <= fixt+twin2)-fixt+twin1;
                temp = zeros(1,twin1+twin2);
                temp(sac_spikes) = 1;
                fixation_locked_firing{unit}(fix_ind(unit),:) = temp;
                fixation_information{unit}(fix_ind(unit),:) = fixation_info(:,fixations(f))';
                fix_ind(unit) = fix_ind(unit)+1;
            end
        end
    end
end
%remove excess NaNs
fixation_locked_firing = laundry(fixation_locked_firing);
fixation_information = laundry(fixation_information);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Mutual Info for Eye Movements and Firing Rate---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temporal_info.fixation = [];
temporal_info.fixation.rate = NaN(1,num_units);
temporal_info.fixation.shuffled_rate = cell(1,num_units);
temporal_info.fixation.shuffled_rate_prctile = NaN(1,num_units);
temporal_info.fixation.temporalstability = NaN(2,num_units);  %row 1 by half, row 2 even/odd trials
temporal_info.fixation.shuffled_temporalstability =cell(1,num_units);
temporal_info.fixation.shuffled_temporalstability_prctile = NaN(2,num_units);%row 1 by half, row 2 even/odd trials

%for limited time period of ~1 saccade+ ~1 fixation
temporal_info_limited.fixation = [];
temporal_info_limited.fixation.rate = NaN(1,num_units);
temporal_info_limited.fixation.shuffled_rate = cell(1,num_units);
temporal_info_limited.fixation.shuffled_rate_prctile = NaN(1,num_units);
temporal_info_limited.fixation.temporalstability = NaN(2,num_units);  %row 1 by half, row 2 even/odd trials
temporal_info_limited.fixation.shuffled_temporalstability =cell(1,num_units);
temporal_info_limited.fixation.shuffled_temporalstability_prctile = NaN(2,num_units);%row 1 by half, row 2 even/odd trials
for unit = 1:num_units
    if ~isempty(fixation_locked_firing{unit})
        
        %don't want to run on trials with the first eye movements occuring with 500 (twin) ms of image onset
        fixation_firing = fixation_locked_firing{unit}(fixation_information{unit}(:,4) > image_on_twin,:);
        
        %since many units appear sparse especially spatial ones want to
        %run on eye movemetns that actually have spikes so use these only
        %doesn't change statistics just decreases run time
        fixation_firing(nansum(fixation_firing,2) == 0,:) = [];%remove Fixations without spikes
        
        [observed_info_rate,shuffled_info_rate]...
            = estimated_mutual_information(fixation_firing,numshuffs,info_type,smval,Fs);
        temporal_info.fixation.rate(unit) = observed_info_rate.skaggs;
        temporal_info.fixation_shuffled_info.rate{unit} = shuffled_info_rate.skaggs;
        temporal_info.fixation.shuffled_rate_prctile(unit) = 100*sum(...
            observed_info_rate.skaggs > shuffled_info_rate.skaggs)/numshuffs;
        temporal_info.fixation.temporalstability(:,unit) = observed_info_rate.temporalstability;
        temporal_info.fixation_shuffled_info.temporalstability{unit} = shuffled_info_rate.temporalstability;
        temporal_info.fixation.shuffled_temporalstability_prctile(1,unit) = ...
            100*sum(observed_info_rate.temporalstability(1) > ...
            shuffled_info_rate.temporalstability(1,:))/numshuffs;
        temporal_info.fixation.shuffled_temporalstability_prctile(2,unit) = ...
            100*sum(observed_info_rate.temporalstability(2) > ...
            shuffled_info_rate.temporalstability(2,:))/numshuffs;
        
        %---limited to 1 eye movement---%
        fix_times = NaN(size(fixation_info,2),twin1+twin2);
        for f = 1:size(fixation_info,2)
            if fixation_info(4,f) > image_on_twin
                ind = fixation_info(4,f)- fixation_info(7,f): fixation_info(4,f)+ fixation_info(5,f);
                ind(ind < 1) = [];
                spikes = spike_times{unit}(fixation_info(1,f),ind);
                ind = ind-fixation_info(4,f);
                spikes(ind < -45) = [];
                ind(ind < -45) = [];
                spikes(ind > twin2) = [];
                ind(ind > twin2) = [];
                if length(ind) < 100
                    disp('what')
                end
                if any(isnan(spikes))
                    error('WTF')
                end
                fix_times(f,ind+twin1) = spikes;
            end
        end
        fix_times = laundry(fix_times,1);
        n_trials = floor(size(fix_times,1)/2);
        fix1 = fix_times(1:n_trials,:);
        fix2 = fix_times(n_trials+1:end,:);
        num_not_nans1 = sum(~isnan(fix1),1);%will be fewer sikes
        num_not_nans1 = num_not_nans1 < min_num_fix;%.5*size(limited_firing,1)); %median duration
        num_not_nans2 = sum(~isnan(fix2),1);%will be fewer sikes
        num_not_nans2 = num_not_nans2 < min_num_fix;%.5*size(limited_firing,1)); %median duration
        num_not_nans = find(num_not_nans1 | num_not_nans2);
        num_not_nans = unique([1:twin1-45 num_not_nans]);
        lmfr1 = nandens3(fix1,smval,Fs);
        lmfr2 = nandens3(fix2,smval,Fs);
        lmfr1(num_not_nans) = [];
        lmfr2(num_not_nans) = [];
        temporal_info_limited.fixation.temporalstability(1,unit) = ...
            corr(lmfr1(:),lmfr2(:),'row','pairwise','type','Spearman');
        
        lmfr = nandens3(fix_times,smval,Fs);
        lmfr(num_not_nans) = [];
        total_time = 1:length(lmfr);
        p_x = total_time/sum(total_time);
        lambda = nansum(nansum(lmfr.*p_x));
        temporal_info_limited.fixation.rate(unit) = p_log_p(lambda,lmfr,p_x);
        
        
        shuffled_corrs = NaN(1,numshuffs);
        shuffled_rate = NaN(1,numshuffs);
        parfor shuff = 1:numshuffs
            shuffled_spike_times  = circshift_row_long(spike_times{unit});
            
            fix_times_shuff = NaN(size(fixation_info,2),twin1+twin2);
            for f = 1:size(fixation_info,2)
                if fixation_info(4,f) > image_on_twin
                    ind = fixation_info(4,f)- fixation_info(7,f): fixation_info(4,f)+ fixation_info(5,f);
                    ind(ind < 1) = [];
                    spikes =  shuffled_spike_times(fixation_info(1,f),ind);
                    ind = ind-fixation_info(4,f);
                    spikes(ind <= -twin1) = [];
                    ind(ind <= -twin1) = [];
                    spikes(ind > twin2) = [];
                    ind(ind > twin2) = [];
                    if any(isnan(spikes))
                        error('WTF')
                    end
                    fix_times_shuff(f,ind+twin1) = spikes;
                end
            end
            fix_times_shuff = laundry(fix_times_shuff,1);
            n_trials = floor(size(fix_times_shuff,1)/2);
            fix1 = fix_times_shuff(1:n_trials,:);
            fix2 = fix_times_shuff(n_trials+1:end,:);
            lmfr1 = nandens3(fix1,smval,Fs);
            lmfr2 = nandens3(fix2,smval,Fs);
            lmfr1(num_not_nans) = [];
            lmfr2(num_not_nans) = [];
            shuffled_corrs(shuff) = corr(lmfr1(:),lmfr2(:),'row','pairwise','type','Spearman');
            
            lmfr = nandens3(fix_times_shuff,smval,Fs);
            lmfr(num_not_nans) = [];
            lambda = nansum(nansum(lmfr.*p_x));
            shuffled_rate(shuff) = p_log_p(lambda,lmfr,p_x);
        end
        temporal_info_limited.fixation_shuffled_info.temporalstability{unit} = shuffled_corrs;
        temporal_info_limited.fixation.shuffled_temporalstability_prctile(1,unit) = ...
            100*sum(temporal_info_limited.fixation.temporalstability(1,unit) > ...
            shuffled_corrs)/numshuffs;
        
        temporal_info_limited.fixation.shuffled_rate{unit} = shuffled_rate;
        temporal_info_limited.fixation.shuffled_rate_prctile(unit) = 100*sum(...
            temporal_info_limited.fixation.rate(unit) > shuffled_rate)/numshuffs;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot Rasters by Various Variables of Interst---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = -twin1:twin2-1;
unit_names = unit_stats(1,:);
for unit = 1:num_units
    if ~isempty(fixation_locked_firing{unit})
        
        %---For Fixations---%
        %ignore anything within 500 ms of image onset
        info = fixation_information{unit}((fixation_information{unit}(:,4) > image_on_twin),:);
        fixation_firing = fixation_locked_firing{unit}(fixation_information{unit}(:,4) > image_on_twin,:);
        n_trials = round(size(fixation_firing,1)/2);
        
        
        figure
        
        %plot firing rate curve over time
        hold on
        subplot(3,3,1)
        dofill(t,fixation_firing(1:n_trials,:),'blue',1,smval);%even trials
        dofill(t,fixation_firing(n_trials+1:end,:),'red',1,smval);%odd trials
        dofill(t,fixation_firing,'black',1,smval); %all trials with spikes
        yl = ylim;
        if yl(1) < 0;
            yl(1) =0;
        end
        plot([0 0],[yl(1) yl(2)],'k--')
        xlim([-twin1 twin2])
        hold off
        legend('1st 1/2','2nd 1/2','All','Location','Best')
        ylim(yl);
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Fixation Start')
        set(gca,'Xtick',[-100 0 100 200 300 400])
        set(gca,'XMinorTick','on','YMinorTick','on')
        grid on
        ylim(yl);
        title(['Bit ' num2str(temporal_info.fixation.shuffled_rate_prctile(unit),3) '% '...
            '\rho_{1/2} = ' num2str(temporal_info.fixation.temporalstability(1,unit),2) ...
            ' (' num2str(temporal_info.fixation.shuffled_temporalstability_prctile(1,unit),3) ...
            '%) \rho_{e/o} = ' num2str(temporal_info.fixation.temporalstability(2,unit),2) ...
            ' (' num2str(temporal_info.fixation.shuffled_temporalstability_prctile(2,unit),3) '%)'])
        
        %plot raster over time by Fixation occurence with session
        subplot(3,3,2)
        [trial,time] = find(fixation_firing == 1);
        plot(time-twin1,(trial),'.k')
        ylim([0 max(trial)])
        ylabel('Occurence #')
        set(gca,'Xtick',[-100 0 100 200 300 400])
        xlabel('Time from Fixation Start')
        title('Raster whole period')
        
        info(info(:,5) > twin2,5) = twin2; %set max fixation duration to <= twin2
        info(info(:,7) >= twin1,7) = twin1-1; %if prior Saccade duration to <= twin1
        
        %want to look only looking at window around saccade/fixation
        limited_firing = NaN(size(fixation_firing));
        info(info(:,5) > twin2,5) = twin2; %if next fixation duration is > twin set to twin !!!!!!!!!!!shouldn't this be 5?
        info(info(:,5) > twin2,5) = twin2; %if next fixation duration is > twin set to twin
        for f = 1:size(limited_firing,1)
            ind = twin1-info(f,7):twin1+info(f,5);
            limited_firing(f,ind) = fixation_firing(f,ind);
        end
        
        %plot firing rate curve over time for spikes limited time period of 1 fixation around saccade
        num_not_nans = sum(~isnan(limited_firing));%will be fewer sikes
        not_enough_samples = find(num_not_nans < min_num_fix);%.5*size(limited_firing,1)); %median duration
        %         limited_firing(:,not_enough_samples) = NaN;  %remove any time points with less than 1/4 of the samples on avg ~< 500 trials
        %         [firing_rate,~]= nandens2(limited_firing,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
        %         firing_rate(not_enough_samples) = NaN;  %remove any time points with less than 1/4 of the samples on avg ~< 500 trials
        [~,limited_firing_rate] = nandens3(limited_firing,smval,Fs);
        limited_firing_rate(:,not_enough_samples) = NaN;
        limited_firing_rate(:,1:twin1-45) = NaN; %don't want to go back more than medianish saccade duration
        limited_firing(:,1:twin1-45) = NaN; %don't want to go back more than medianish saccade duration
        
        subplot(3,3,4)
        plot(t,nanmean(limited_firing_rate),'k');
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
        xlabel('Time from Fixation Start')
        set(gca,'Xtick',[-100 0 100 200 300 400])
        set(gca,'XMinorTick','on','YMinorTick','on')
        grid on
        ylim(yl);
        title(['1 Eye Movement, \rho_{1/2} = ' num2str(temporal_info_limited.fixation.temporalstability(1,unit),2) ...
            ' (' num2str(temporal_info_limited.fixation.shuffled_temporalstability_prctile(1,unit),3) '%), ' ...
            'Bit ' num2str(temporal_info_limited.fixation.shuffled_rate_prctile(unit),3) ' %'])
        
        fix_durs = info(:,5);
        [~,sorted_fix_durs] = sort(fix_durs);
        %make raster plot for spikes limited time period of 1 fixation around Fixation
        lmf = limited_firing(sorted_fix_durs,:);
        subplot(3,3,5)
        [trial,time] = find(lmf == 1);
        plot(time-twin1,trial,'.k')
        xlim([-twin1 twin2])
        ylim([0 max(trial)])
        set(gca,'Xtick',[-100 0 100 200 300 400])
        ylabel('Sorted by Fixation Duration')
        xlabel('Time from Fixation Start')
        title('Raster Limited to 1 Eye Movement')
        
        
        %plot raster over time by spikes/Fixation epoch aka by average firing rate
        f1 = sum(fixation_firing,2);
        [~,fi] = sort(f1);
        [trial,time] = find(fixation_firing(fi,:) == 1);
        subplot(3,3,3)
        plot(time-twin1,(trial),'.k')
        ylim([0 max(trial)])
        set(gca,'Xtick',[-100 0 100 200 300 400])
        ylabel('Ranked Spike Count')
        xlabel('Time from Fixation Start')
        title('For all Fixations including zero spike counts')
        
        %plot raster over time by Fixation order within image presentation
        [~,fi] = sort(info(:,3));
        subplot(3,3,6)
        [trial,time] = find(fixation_firing(fi,:) == 1);
        plot(time-twin1,trial,'.k')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-100 0 100 200 300 400])
        ylabel('Ranked Ordinal Fixation #')
        xlabel('Time from Fixation Start')
        title('Discrete Time within Image Period')
        
        %         %plot raster by Fixation amplitude
        %         [~,fi] = sort(info(:,5));
        %         subplot(3,3,7)
        %         [trial,time] = find(fixation_firing(fi,:) == 1);
        %         plot(time-twin1,trial,'.k')
        %         hold off
        %         ylim([0 max(trial)])
        %         set(gca,'Xtick',[-100 0 100 200 300 400])
        %         ylabel('Ranked Fixation Amplitude')
        %         xlabel('Time from Fixation Start')
        %         title('Fixation Amplitude')
        
        %plot raster by novel vs repeat
        subplot(3,3,7)
        [trial,time] = find(fixation_firing(info(:,9) == 1,:) == 1);
        plot(time-twin1,trial,'.b')
        if ~isempty(trial)
            b4 = max(trial);
        else
            b4 = 0;
        end
        hold on
        [trial,time] = find(fixation_firing(info(:,9) == 2,:) == 1);
        trial= trial+b4;
        plot(time-twin1,trial,'.r')
        hold off
        if ~isempty(trial)
            ylim([0 max(trial)])
        end
        set(gca,'Xtick',[-100 0 100 200 300 400])
        ylabel('Ranked Fixation Duration')
        xlabel('Time from Fixation Start')
        title('Novel (blue) vs Repeat (Red)')
        
        %plot raster by Saccade direction
        [~,fi] = sort(info(:,8));
        subplot(3,3,8)
        [trial,time] = find(fixation_firing(fi,:) == 1);
        plot(time-twin1,trial,'.k')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-100 0 100 200 300 400])
        ylabel('Ranked Saccade Direction')
        xlabel('Time from Fixation Start')
        title('Saccade Direction')
        
        %plot raster by fixation duration
        [~,fi] = sort(info(:,5));
        subplot(3,3,9)
        [trial,time] = find(fixation_firing(fi,:) == 1);
        dur = info(fi,5);
        dur(dur > twin2) = twin2;
        plot(time-twin1,trial,'.k')
        hold on
        plot(dur,1:length(fi),'r.')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-100 0 100 200 300 400])
        ylabel('Ranked Fixation Duration')
        xlabel('Time from Fixation Start')
        title('Fixation Duration')
        
        n_str = ['   n_{sac} =' num2str(size(fixation_locked_firing{1,unit},1))];
        if multiunit(unit)
            subtitle(['Fixation-Locked Multiunit ' unit_names{unit} n_str]);
        else
            subtitle(['Fixation-Locked ' unit_names{unit} n_str]);
        end
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} 'Eye_Locked_analysis_Fixation_Rasters']);
    end
end

save([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat'],...
    'twin1','twin2','image_on_twin','smval','fixation_locked_firing',...
    'fixation_information','temporal_info','unit_names','spike_times',...
    'temporal_info_limited');
end

function [info] = p_log_p(lambda,lambda_x,p_x)
%mutual info equation for Skagg Score
%from Skaggs, McNaughton, and Gothard 1993
plogp = lambda_x.*log2(lambda_x/lambda);
info = nansum(nansum(plogp.*p_x)); %bits/second
end