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
%   3) session_data
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-Eyemovement_Locked_List_results'

figure_dir = [figure_dir 'List Saccade Analysis\'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set important Analysis parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
twin = 500;% how much time to take before and after saccade as well as how much time to ignore after image onset
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
trial_start_code = 15;
trial_end_code = 20;
imgon_code = 23;
imgoff_code = 24;
smval = 60 ;%gaussian 1/2 width for smoothing, for the moment want to smooth at high frequency modulations
numshuffs = 1000; %recommend this is between 100 & 1000, for bootstrapping to
% dtermine significance
info_type = 'temporal';
Fs = 1000;
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
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);

if num_units == 0
    disp([task_file(1:8) ': no units could be found. Exiting function...'])
    return %since no units skip analysis
end

num_trials = length(cfg.trl);
%NaNs are for start and end trials otherwise cut
valid_trials(1,isnan(valid_trials(1,:))) = 1;
valid_trials(2,isnan(valid_trials(2,:))) = num_trials;

if str2double(task_file(3:8)) < 140805 %7 second viewing with fewere sequence trials
    minimum_trials_1 = 101; %for session start minimum number of trials: includes fam block and 16 novel/repeat images + sequence trials
    minimum_trials_2 =  80;%for rest of session: includes 16 novel/repeat images + sequence trials
else
    minimum_trials_1 = 117; %for session start minimum number of trials: includes fam block and 16 novel/repeat images + sequence trials
    minimum_trials_2 =  96;%for rest of session: includes 16 novel/repeat images + sequence trials
end

for unit = 1:num_units
    start_end = valid_trials(:,unit);
    if isnan(start_end(1))
        start_end(1) = 1;
    end
    if isnan(start_end(2))
        start_end(2) = length(cfg.trl);
    end
    start_end(start_end == 0) = 1;
    min_trial = cfg.trl(start_end(1)).cnd-1000; %get the first condition number
    max_trial = cfg.trl(start_end(2)).cnd-1000; %get the last condition number
    
    if min_trial < 22 %includes fam block
        min_trial = 22; %remove then count from there
    end
    num_blks = floor((max_trial-min_trial+1)/minimum_trials_2);
    
    if num_blks < min_blks
        valid_trials(:,unit) = NaN;
    end
end

%set the image duration
if str2double(task_file(3:8)) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end

if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Get successful trials---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([task_file ': running saccade modulation anlaysis..'])

%get important task specific information
[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);

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
saccade_info = NaN(9,25*length(fixationstats));
fixation_info = NaN(8,25*length(fixationstats));

sac_ind = 1;
fix_ind = 1;
for t = 1:num_trials %will use all trials since may be useful to look at eye movements later for all trials
    fixationtimes = fixationstats{t}.fixationtimes;
    saccadetimes = fixationstats{t}.saccadetimes;
    xy = fixationstats{t}.XY;
    
    trial_start = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == trial_start_code);
    trial_end = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == trial_end_code);
    imgon = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == imgon_code)-trial_start;
    imgoff = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == imgoff_code)-trial_start;
    
    
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
    invalid= find(fixationtimes(1,:) > imgoff+twin);
    fixationtimes(:,invalid) = [];
    
    %saccade started before image turned on
    invalid= find(saccadetimes(1,:) < imgon);
    saccadetimes(:,invalid) = [];
    
    %saccade started after the image turned off and/or firing rate could corrupted by image turning off
    invalid= find(saccadetimes(1,:) > imgoff+twin);
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
            sac_ind = sac_ind+1;
        end
    end
    for f = 1:size(fixationtimes,2);
        if  (fixationtimes(2,f)-fixationtimes(1,f)+1) >= min_fix_dur %if fixation is long enough
            fix_dur = fixationtimes(2,next_fix)-fixationtimes(1,next_fix)+1;%next fixation duration
            fixation_info(1,fix_ind) = t; %trial #
            fixation_info(2,fix_ind) = fixationtimes(1,f); %fixation start time from trial start
            fixation_info(3,fix_ind) = f; %ordinal fixation #
            fixation_info(4,fix_ind) = fixationtimes(1,f)-imgon; %fixation start time relative to image onset
            fixation_info(5,fix_ind) = fixationtimes(2,f)-fixationtimes(1,f)+1; %fixation duration so how far forward can we look
            prior_sac = find(fixationtimes(1,f) == saccadetimes(2,:)+1);%find fixaiton before in case want to know
            if ~isempty(prior_sac);
                fixation_info(6,fix_ind) = saccadetimes(2,prior_sac)-saccadetimes(1,prior_sac)+1; %how long before aka saccade duration
                fixation_info(7,fix_ind) = sqrt(sum((xy(:,saccadetimes(2,prior_sac))...
                    -xy(:,saccadetimes(1,prior_sac))).^2)); %prior saccade amplitude%saccade amplitude
                fixation_info(8,fix_ind) = atan2d(xy(2,saccadetimes(2,prior_sac))...
                    -xy(2,saccadetimes(1,prior_sac)),xy(1,saccadetimes(2,prior_sac))-xy(1,saccadetimes(1,prior_sac)));%prior saccade direction
            end
            fix_ind = fix_ind+1;
        end
    end
end

%Remove excess NaNs;
saccade_info = laundry(saccade_info);
fixation_info = laundry(fixation_info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Firing Rate Locked to Eye Movements---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pre-allocate space for variables
fixation_locked_firing = cell(1,unit); %firing rate locked to fixation
saccade_locked_firing = cell(1,unit); %firing rate locked to saccade
saccade_information = cell(1,unit); %start time relative to image onset and ordinal #
fixation_information = cell(1,unit);%start time relative to image onset and ordinal #
for unit = 1:num_units
    fixation_locked_firing{unit}= NaN(size(fixation_info,2),twin*2);
    saccade_locked_firing{unit}= NaN(size(saccade_info,2),twin*2);
    saccade_information{unit} = NaN(size(saccade_info,2),9);
    fixation_information{unit} = NaN(size(fixation_info,2),8);
end

fix_ind = ones(1,num_units);
sac_ind = ones(1,num_units);
for trial = 1:num_trials
    for unit = 1:num_units;
        if image_trials(trial) >= valid_trials(1,unit) && image_trials(trial) <= valid_trials(2,unit) %only valid trials
            fixations = find(fixation_info(1,:) == trial);
            saccades = find(saccade_info(1,:) == trial);
            spikes = find(data(unit).values{image_trials(trial)});
            
            %collect spike times relative to fixation onset
            for f = 1:length(fixations);
                fixt = fixation_info(2,fixations(f));
                fix_spikes = spikes(spikes > fixt-twin & spikes <= fixt+twin)-fixt+twin;
                temp = zeros(1,twin*2);
                temp(fix_spikes) = 1;
                fixation_locked_firing{unit}(fix_ind(unit),:) = temp;
                fixation_information{unit}(fix_ind(unit),:) = fixation_info(:,fixations(f))';
                fix_ind(unit) = fix_ind(unit)+1;
            end
            
            %collect spike times relative to saccade start
            for s = 1:length(saccades);
                sact = saccade_info(2,saccades(s));
                sac_spikes = spikes(spikes > sact-twin & spikes <= sact+twin)-sact+twin;
                temp = zeros(1,twin*2);
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
fixation_locked_firing = laundry(fixation_locked_firing);
fixation_information = laundry(fixation_information);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot Rasters by Various Variables of Interst---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clrs = ['rgbkrgbkrgbkrgbkrgbkrgbk'];
t = -500:499;
for unit = 1:num_units
    if ~isempty(saccade_locked_firing{unit})
        
        %---For Saccades---%
        %ignore anything within 500 ms of image onset
        info = saccade_information{unit}((saccade_information{unit}(:,4) > twin),:);
        saccade_firing = saccade_locked_firing{unit}(saccade_information{unit}(:,4) > twin,:);
        [firing_rate,~]= nandens(saccade_firing,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
        
        
        %plot firing rate curve over time
        figure
        subplot(3,3,1)
        plot(t,firing_rate);
        yl = ylim;
        if yl(1) < 0;
            yl(1) =0;
        end
        hold on
        plot([0 0],[yl(1) yl(2)],'k--')
        hold off
        ylim(yl);
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Saccade Start')
        set(gca,'Xtick',[-500 -250 0 250 500])
        set(gca,'XMinorTick','on','YMinorTick','on')
        grid on
        ylim(yl);
        
        %plot raster over time by saccade occurence with session
        subplot(3,3,2)
        [trial,time] = find(saccade_firing == 1);
        plot(time-500,(trial),'.k')
        ylim([0 max(trial)])
        ylabel('Occurence #')
        set(gca,'Xtick',[-500 -250 0 250 500])
        xlabel('Time from Saccade Start')
        
        %plot raster over time by spikes/saccade epoch aka by average firing rate
        f1 = sum( saccade_firing,2);
        [~,fi] = sort(f1);
        [trial,time] = find(saccade_firing(fi,:) == 1);
        subplot(3,3,3)
        plot(time-500,(trial),'.k')
        ylim([0 max(trial)])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Ranked Spike Count')
        xlabel('Time from Saccade Start')
        
        %plot raster over time by saccade order within image presentation
        [~,fi] = sort(info(:,3));
        subplot(3,3,4)
        [trial,time] = find(saccade_firing(fi,:) == 1);
        plot(time-500,trial,'.k')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Ranked Ordinal Saccade #')
        xlabel('Time from Saccade Start')
        
        %plot raster by saccade amplitude
        [~,fi] = sort(info(:,5));
        subplot(3,3,5)
        [trial,time] = find(saccade_firing(fi,:) == 1);
        plot(time-500,trial,'.k')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Ranked Saccade Amplitude')
        xlabel('Time from Saccade Start')
        
        %plot raster by saccade direction
        [~,fi] = sort(info(:,8));
        subplot(3,3,6)
        [trial,time] = find(saccade_firing(fi,:) == 1);
        plot(time-500,trial,'.k')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Ranked Saccade Direction')
        xlabel('Time from Saccade Start')
        
        %want to look only looking at window around saccade/fixation
        limited_firing = NaN(size(saccade_firing));
        info(isnan(info(:,9)),9) = round(nanmean(info(:,9)));%go back on average if no fixation before
        info(info(:,9) >= twin,9) = twin-1; %if prior fixation duration is > twin set to twin
        info(info(:,7) > twin,7) = twin; %if next fixation duration is > twin set to twin
        for f = 1:size(limited_firing,1)
            ind = twin-info(f,9):twin+info(f,7);
            limited_firing(f,ind) = saccade_firing(f,ind);
        end
        
        %make raster plot for spikes limited time period of 1 fixation around saccade
        subplot(3,3,8)
        [trial,time] = find(limited_firing == 1);
        plot(time-500,trial,'.k')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Occurence #')
        xlabel('Time from Saccade Start')
        title('Time Period 1 fixation arround saccade')
        
        
        
        %plot firing rate curve over time for spikes limited time period of 1 fixation around saccade
        num_not_nans = sum(~isnan(limited_firing));%will be fewer sikes
        not_enough_samples = find(num_not_nans < .5*size(limited_firing,1));
        limited_firing(:,not_enough_samples) = NaN;  %remove any time points with less than 1/4 of the samples on avg ~< 500 trials
        
        [firing_rate,~]= nandens(limited_firing,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
        
        subplot(3,3,7)
        plot(t,firing_rate);
        yl = ylim;
        if yl(1) < 0;
            yl(1) =0;
        end
        hold on
        plot([0 0],[yl(1) yl(2)],'k--')
        hold off
        ylim(yl);
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Saccade Start')
        set(gca,'Xtick',[-500 -250 0 250 500])
        set(gca,'XMinorTick','on','YMinorTick','on')
        grid on
        ylim(yl);
        title('Time Period 1 fixation arround saccade')
        
        n_str = [' n_{sac} =' num2str(size(saccade_locked_firing{1,unit},1))];
        if multiunit(unit)
            subtitle(['Saccade-Locked Multiunit ' unit_names{unit} n_str]);
        else
            subtitle(['Saccade-Locked ' unit_names{unit} n_str]);
        end
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} 'Eye_Locked_analysis_Saccade_Rasters']);
        
        
        
        %---For fixations---%
        %ignore anything within 500 ms of image onset
        info = fixation_information{unit}((fixation_information{unit}(:,4) > twin),:);
        fixation_firing = fixation_locked_firing{unit}(fixation_information{unit}(:,4) > twin,:);
        [firing_rate,~]= nandens(fixation_firing,30,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
        
        
        %plot firing rate curve over time
        figure
        subplot(3,3,1)
        plot(t,firing_rate);
        yl = ylim;
        if yl(1) < 0;
            yl(1) =0;
        end
        hold on
        plot([0 0],[yl(1) yl(2)],'k--')
        hold off
        ylim(yl);
        ylabel('Firing Rate (Hz)')
        xlabel('Time from fixation Start')
        set(gca,'Xtick',[-500 -250 0 250 500])
        set(gca,'XMinorTick','on','YMinorTick','on')
        grid on
        ylim(yl);
        
        %plot raster over time by fixation occurence with session
        subplot(3,3,2)
        [trial,time] = find(fixation_firing == 1);
        plot(time-500,(trial),'.k')
        ylim([0 max(trial)])
        ylabel('Occurence #')
        set(gca,'Xtick',[-500 -250 0 250 500])
        xlabel('Time from fixation Start')
        
        %plot raster over time by spikes/fixation epoch aka by average firing rate
        f1 = sum( fixation_firing,2);
        [~,fi] = sort(f1);
        [trial,time] = find(fixation_firing(fi,:) == 1);
        subplot(3,3,3)
        plot(time-500,(trial),'.k')
        ylim([0 max(trial)])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Ranked Spike Count')
        xlabel('Time from fixation Start')
        
        %plot raster over time by fixation order within image presentation
        [~,fi] = sort(info(:,3));
        subplot(3,3,4)
        [trial,time] = find(fixation_firing(fi,:) == 1);
        plot(time-500,trial,'.k')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Ranked Ordinal fixation #')
        xlabel('Time from fixation Start')
        
        %plot raster by Fixation Duration
        [~,fi] = sort(info(:,5));
        subplot(3,3,5)
        [trial,time] = find(fixation_firing(fi,:) == 1);
        plot(time-500,trial,'.k')
        hold on
        dur = info(fi,5);
        dur(dur > 500) = 500;
        plot(dur,1:length(fi),'r.')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Ranked Fixation Duration')
        xlabel('Time from fixation Start')
        
        %plot raster by prior saccade's direction
        [~,fi] = sort(info(:,8));
        subplot(3,3,6)
        [trial,time] = find(fixation_firing(fi,:) == 1);
        plot(time-500,trial,'.k')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Ranked Saccade Direction')
        xlabel('Time from fixation Start')
        
        %plot raster by prior saccade's amplitude
        [~,fi] = sort(info(:,7));
        subplot(3,3,7)
        [trial,time] = find(fixation_firing(fi,:) == 1);
        plot(time-500,trial,'.k')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Ranked Saccade Amplitude')
        xlabel('Time from fixation Start')
        
        %want to look only looking at window around fixation/fixation
        limited_firing = NaN(size(fixation_firing));
        if sum(isnan(info(:,6))) ~= 0
            info(isnan(info(:,6)),6) = round(nanmean(info(:,6)));%go back on average if no saccade before
        end
        info(info(:,5) > twin,5) = twin; %fixations longer than twin should be set to twin
        for f = 1:size(limited_firing,1)
            ind = twin-info(f,6):twin+info(f,5);
            limited_firing(f,ind) = fixation_firing(f,ind);
        end
        
        %make raster plot for spikes limited time period of 1 fixation around fixation
        subplot(3,3,9)
        [trial,time] = find(limited_firing == 1);
        plot(time-500,trial,'.k')
        hold off
        ylim([0 max(trial)])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Occurence #')
        xlabel('Time from fixation Start')
        title('Time Period 1 saccade before fixation')
        
        %plot firing rate curve over time for spikes limited time period of 1 fixation around fixation
        num_not_nans = sum(~isnan(limited_firing));%will be fewer sikes
        not_enough_samples = find(num_not_nans < .5*size(limited_firing,1));
        limited_firing(:,not_enough_samples) = NaN;  %remove any time points with less than 1/4 of the samples on avg ~< 500 trials
        
        [firing_rate,~]= nandens(limited_firing,30,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
        
        subplot(3,3,8)
        plot(t,firing_rate);
        yl = ylim;
        if yl(1) < 0;
            yl(1) =0;
        end
        hold on
        plot([0 0],[yl(1) yl(2)],'k--')
        hold off
        ylim(yl);
        ylabel('Firing Rate (Hz)')
        xlabel('Time from fixation Start')
        set(gca,'Xtick',[-500 -250 0 250 500])
        set(gca,'XMinorTick','on','YMinorTick','on')
        grid on
        ylim(yl);
        title('Time Period 1 saccade before fixation')
        
        n_str = [' n_{fix} = ' num2str(size(fixation_locked_firing{1,unit},1))];
        if multiunit(unit)
            subtitle(['Fixation-Locked Multiunit ' unit_names{unit} n_str]);
        else
            subtitle(['Fixation-Locked ' unit_names{unit} n_str]);
        end
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} 'Eye_Locked_analysis_Fixation_Rasters']);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Mutual Info for Eye Movements and Firing Rate---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fixation_info = [];
fixation_info.rate = NaN(1,num_units);
fixation_info.temporalstability = NaN(1,num_units);
fixation_info_95.rate = NaN(1,num_units);
fixation_info_95.temporalstability = NaN(1,num_units);
fixation_shuffled_info.rate = cell(1,num_units);
fixation_shuffled_info.temporalstability = cell(1,num_units);
saccade_info = [];
saccade_info.rate = NaN(1,num_units);
saccade_info.temporalstability = NaN(1,num_units);
saccade_info_95.rate = NaN(1,num_units);
saccade_info_95.temporalstability = NaN(1,num_units);
saccade_shuffled_info.rate = cell(1,num_units);
saccade_shuffled_info.temporalstability = cell(1,num_units);

fixation_info2 = [];
fixation_info2.rate = NaN(1,num_units);
fixation_info2.temporalstability = NaN(1,num_units);
fixation_info2_95.rate = NaN(1,num_units);
fixation_info2_95.temporalstability = NaN(1,num_units);
fixation_shuffled_info2.rate = cell(1,num_units);
fixation_shuffled_info2.temporalstability = cell(1,num_units);
saccade_info2 = [];
saccade_info2.rate = NaN(1,num_units);
saccade_info2.temporalstability = NaN(1,num_units);
saccade_info2_95.rate = NaN(1,num_units);
saccade_info2_95.temporalstability = NaN(1,num_units);
saccade_shuffled_info2.rate = cell(1,num_units);
saccade_shuffled_info2.temporalstability = cell(1,num_units);

for unit = 1:num_units
    if ~isempty(fixation_locked_firing{unit})
        %don't want to run on trials with the first eye movements occuring with 500 (twin) ms of image onset
        [observed_info_rate,shuffled_info_rate]...
            = estimated_mutual_information(fixation_locked_firing{unit}(fixation_information{unit}(:,4) > twin,:),numshuffs,info_type,smval,Fs);
        
        fixation_info.rate(unit) = observed_info_rate.skaggs;
        fixation_info.temporalstability(unit) = observed_info_rate.temporalstability;
        fixation_shuffled_info.rate{unit} = shuffled_info_rate.skaggs;
        fixation_shuffled_info.temporalstability{unit} = shuffled_info_rate.temporalstability;
        fixation_info_95.rate(unit) = prctile(fixation_shuffled_info.rate{unit},95);
        fixation_info_95.temporalstability(unit) = prctile(fixation_shuffled_info.temporalstability{unit},95);
        
        [observed_info_rate,shuffled_info_rate]...
            = estimated_mutual_information(saccade_locked_firing{unit}(saccade_information{unit}(:,4) > twin,:),numshuffs,info_type,smval,Fs);
        saccade_info.rate(unit) = observed_info_rate.skaggs;
        saccade_info.temporalstability(unit) = observed_info_rate.temporalstability;
        saccade_shuffled_info.rate{unit} = shuffled_info_rate.skaggs;
        saccade_shuffled_info.temporalstability{unit} = shuffled_info_rate.temporalstability;
        saccade_info_95.rate(unit) = prctile(saccade_shuffled_info.rate{unit},95);
        saccade_info_95.temporalstability(unit) = prctile(saccade_shuffled_info.temporalstability{unit},95);
        
        %since many units appear sparse especially spatial ones want to
        %run on eye movemetns that actually have spikes so use these only
        saccade_firing = saccade_locked_firing{unit}(saccade_information{unit}(:,4) > twin,:);
        fixation_firing = fixation_locked_firing{unit}(fixation_information{unit}(:,4) > twin,:);
        saccade_firing(nansum(saccade_firing,2) == 0,:) = [];%remove saccades without spikes
        fixation_firing(nansum(fixation_firing,2) == 0,:) = [];%remove fixations without spikes
        
        [observed_info_rate,shuffled_info_rate] = estimated_mutual_information(fixation_firing,numshuffs,info_type,smval,Fs);
        fixation_info2.rate(unit) = observed_info_rate.skaggs;
        fixation_info2.temporalstability(unit) = observed_info_rate.temporalstability;
        fixation_shuffled_info2.rate{unit} = shuffled_info_rate.skaggs;
        fixation_shuffled_info2.temporalstability{unit} = shuffled_info_rate.temporalstability;
        fixation_info2_95.rate(unit) = prctile(fixation_shuffled_info.rate{unit},95);
        fixation_info2_95.temporalstability(unit) = prctile(fixation_shuffled_info.temporalstability{unit},95);
        
        [observed_info_rate,shuffled_info_rate] = estimated_mutual_information(saccade_firing,numshuffs,info_type,smval,Fs);
        saccade_info2.rate(unit) = observed_info_rate.skaggs;
        saccade_info2.temporalstability(unit) = observed_info_rate.temporalstability;
        saccade_shuffled_info2.rate{unit} = shuffled_info_rate.skaggs;
        saccade_shuffled_info2.temporalstability{unit} = shuffled_info_rate.temporalstability;
        saccade_info2_95.rate(unit) = prctile(saccade_shuffled_info.rate{unit},95);
        saccade_info2_95.temporalstability(unit) = prctile(saccade_shuffled_info.temporalstability{unit},95);
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot Firing Rates Locked to Eye Movements---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%only plot for fixations 1-16 though
t = -twin:twin-1;
unit_names = unit_stats(1,:);
peak_firing_rate = NaN(2,unit);
sb = [1:10 13 14];
for unit = 1:num_units
    if isempty(fixation_locked_firing{unit})
        continue
    end
    
    figure
    ylims = NaN(2,13);
    
    %---Plots by Ordinal Eye Movement #---%
    for c = 1:12
        subplot(4,4,sb(c))
        hold on
        [~,~,~,y1] = dofill(t,fixation_locked_firing{unit}(fixation_information{unit}(:,3) == c,:),'blue',1,smval);
        [~,~,~,y2] = dofill(t,saccade_locked_firing{unit}(saccade_information{unit}(:,3) == c,:),'red',1,smval);
        hold off
        ylims(1,c) = 0.9*min(min(y1),min(y2));
        ylims(2,c) = 1.1*max(max(y1),max(y2));
        
        title(['Time from Eye Movement ' num2str(c) ' (ms)'])
        ylabel('Firing Rate (Hz)')
    end
    
    
    %---Plots for all fixations and Saccades---%
    %will ignore first fixations/saccades within twin (500 ms) of image onset
    subplot(4,4,[11 12 15 16])
    hold on
    [~,~,~,y1] = dofill(t,fixation_locked_firing{unit}(fixation_information{unit}(:,4) > twin,:),'blue',1,smval);
    [~,~,~,y2] =dofill(t,saccade_locked_firing{unit}(saccade_information{unit}(:,4) > twin,:),'red',1,smval);
    hold off
    ylims(1,13) = 0.9*min(min(y1),min(y2));
    ylims(2,13) = 1.1*max(max(y1),max(y2));
    peak_firing_rate(1,unit) = prctile(y1,99);
    peak_firing_rate(2,unit) = prctile(y2,99);
    
    
    xlabel('Time from Eye Movement (ms)')
    ylabel('Firing Rate (Hz)')
    legend('Fixations','Saccades')
    
    title_str =  'All Fixations/Saccades';
    if  (fixation_info.rate(unit) > fixation_info_95.rate(unit)) && (fixation_info.temporalstability(unit) > fixation_info_95.temporalstability(unit))
        title_str = [title_str ' fix_{95}'];
    end
    if(saccade_info.rate(unit) > saccade_info_95.rate(unit)) && (saccade_info.temporalstability(unit) > saccade_info_95.temporalstability(unit))
        title_str = [title_str ' sac_{95}'];
    end
    title(title_str);
    
    
    ymin = min(ylims(1,:));
    ymin(ymin < 0) = 0;
    ymax = max(ylims(2,:));
    for c = 1:13;
        if c == 13
            subplot(4,4,[11 12 15 16])
            ylim([ymin ymax]);
            hold on
            plot([0 0],[ymin ymax],'k')
            hold off
        else
            subplot(4,4,sb(c))
            ylim([ymin ymax]);
            hold on
            plot([0 0],[ymin ymax],'k')
            hold off
        end
    end
    
    n_str = [' n_{sac} =' num2str(size(saccade_locked_firing{1,unit},1)) ...
        ' n_{fix} = ' num2str(size(fixation_locked_firing{1,unit},1))];
    if multiunit(unit)
        subtitle(['Saccade-Locked Multiunit ' unit_names{unit} n_str]);
    else
        subtitle(['Saccade-Locked ' unit_names{unit} n_str]);
    end
    
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_List_All_Eye_Locked_analysis']);
end

save([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat'],...
    'twin','smval','saccade_locked_firing','fixation_locked_firing',...
    'fixation_information','saccade_information','peak_firing_rate',...
    'fixation_info','fixation_info_95','fixation_shuffled_info',...
    'saccade_info','saccade_info_95','saccade_shuffled_info',...
    'fixation_info2','fixation_info2_95','fixation_shuffled_info2',...
    'saccade_info2','saccade_info2_95','saccade_shuffled_info2');
end