function Spike_Eye_Frequency_Modulation_Analysis(data_dir,figure_dir,session_data)
%written by Seth Konig 3/21/17
%Code analysis Spike data autocorrelations, eye data correlations, and
%spike-eye crosscorrelations for significant frequency modulation. Code
%does this for List and Sequences parts of ListSQ seperately as the eye
%movement rate during List is ~4 Hz but is ~2-3 Hz during Sequences. A
%particular interest of this anlaysis is whether we see 4 Hz modulation (or
%other theta frequencies) during Sequeences and List images or if the 4 Hz
%modulation is just due to eye movements.

%---Misc Parameters---%
figure_dir = [figure_dir 'Freq Modulation\']; %mother dir
clrs = ['rgbmc'];%colors
task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
twin = 500;%maximum lag on correlations
tm = -twin:twin;%time vector
min_spike_count = 250;%don't analyze epoch if not enough spikes, can be both or either List/Sequence
image_on_twin = 500;%ignore first 500 ms
img_on_code = 23;%also first item on code in sequence task
img_off_code = 24;
reward_code = 3;
trial_start_code = 15;
pthresh = 99;%p < 0.01, bonferroni corrected 0.05/5 since testing 5 frequencies

%---PSD Autocorrelation Parmaeters---%
PSD_numshuffs = 1000;%10000;
NFFT = 2^12;%number of frequency points from 0-Fs/2
f2 = Fs/2*linspace(0,1,NFFT/2+1);%frequency values
lower_bound_frequencies = [4 10 20 30 100];

%find indicese in f2 that correspond to start of these bands
lower_bound_frequencies_ind = NaN(1,length(lower_bound_frequencies));
for f = 1:length(lower_bound_frequencies);
    ind = find(f2 >= lower_bound_frequencies(f));
    ind = ind(1);
    lower_bound_frequencies_ind(f) =  ind;
end

%---Low Frequency Filter---%
smFreqVal = 2;%how much smoothing is done in Hz e.g. +/- 1.0 Hz
freq_filt = ones(1,round(smFreqVal/mode(diff(f2)))); %square window/moving average
freq_filt = freq_filt/sum(freq_filt);%zero sum filter
Hz1 = round(smFreqVal/2/mode(diff(f2)));
Hz4 = Hz1*4;

%---High Frequency Filter---%
smFreqVal2 = 10;%how much smoothing is done in Hz e.g. +/- 5.0 Hz
freq_filt2 = ones(1,round(smFreqVal2/mode(diff(f2)))); %square window/moving average
freq_filt2 = freq_filt2/sum(freq_filt2);%zero sum filter
Hz400 = find(f2 >= 400);
Hz400 = Hz400(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Import Task Data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task_file = get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
disp(task_file(1:8))
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'cfg','hdr','data',...
    'fixationstats','valid_trials','excitatory_inhibitory');

[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
[~,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

if num_units == 0; %if no units exit function
    disp([task_file(1:8) ': no units could be found. Exiting function...'])
    return;
end

valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end

%remove units with too few trials
%these are the absolute minimum data required to do data analysis
%set the image duration
if str2double(task_file(3:8)) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end
imgdur = imgdur-image_on_twin;%how much data to take

%get important task specific information
[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[which_img,~,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);

%load spatial modulation info too
load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'spatial_info')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Variables to be Saved---%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---Whole session spike autocorrelation Analysis---%
whole_spike_times = cell(1,num_units);
whole_spike_count = zeros(1,num_units);
whole_spike_ac_observed_ac = cell(1,num_units);
whole_spike_ac_observed_FFT = cell(1,num_units);
whole_spike_ac_shuffled_ac_mean = cell(num_units,length(lower_bound_frequencies));
whole_spike_ac_shuffled_FFT_99 = cell(num_units,length(lower_bound_frequencies));
whole_spike_ac_significant_freq = cell(1,num_units);

%---List Spike Autocorrelation Analysis---%
list_spike_times = cell(1,num_units);
list_spike_ac_observed_ac = cell(1,num_units);
list_spike_ac_observed_FFT = cell(1,num_units);

%---Sequence Spike Autocorrelation Analysis---%
seq_spike_times = cell(1,num_units);
seq_spike_ac_observed_ac = cell(1,num_units);
seq_spike_ac_observed_FFT = cell(1,num_units);

%---List Eye Autocorrelation Analysis---%
list_eye_times = cell(1,num_units);
list_eye_ac_observed_ac = cell(1,num_units);
list_eye_ac_observed_FFT = cell(1,num_units);

%---Sequence Eye Autocorrelation Analysis---%
seq_eye_times = cell(1,num_units);
seq_eye_ac_observed_ac = cell(1,num_units);
seq_eye_ac_observed_FFT = cell(1,num_units);

%---List Spike-Eye Crosscorrelation Analysis---%
list_spike_eye_xc_observed_xc = cell(1,num_units);
list_spike_eye_xc_observed_max_xc = NaN(1,num_units);
list_spike_eye_xc_observed_FFT = cell(1,num_units);

%---Sequence Spike-Eye Crosscorrelation Analysis---%
seq_spike_eye_xc_observed_xc = cell(1,num_units);
seq_spike_eye_xc_observed_max_xc = NaN(1,num_units);
seq_spike_eye_xc_observed_FFT = cell(1,num_units);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Get Sequence Task Successful Trial Info---%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%faster to run outside loop only 1x
seq_fixationstats = fixationstats; %set as new variable because written over below
seq_cfg = cfg; %set as new variable because written over below

%preallocate space and parallel structure of cfg
successful_sequence_trials = NaN(1,length(cfg.trl)); %sequence trial number
which_sequence = NaN(1,length(cfg.trl)); %which sequence 1 or 2
seq_trial_dur = NaN(1,length(cfg.trl));%trial duration
for t = 1:length(cfg.trl);
    if sum(cfg.trl(t).allval == reward_code) == 6; %in which sequence trials were rewarded
        which_sequence(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
        item1on =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code); %when image turned on
        reward = cfg.trl(t).alltim(cfg.trl(t).allval == reward_code); %when image turned off
        reward = reward(1); %start of reward period
        successful_sequence_trials(t) = t;
        seq_trial_dur(t) = reward-item1on+1;
    end
end

%clean up excess nans
successful_sequence_trials = laundry(successful_sequence_trials);
which_sequence = laundry(which_sequence);
seq_trial_dur = laundry(seq_trial_dur);

%find out trial durations
max_seq_trial_dur = max(seq_trial_dur);
min_seq_trial_dur = min(seq_trial_dur);
if min_seq_trial_dur < 3750 %most trials are much longer than this so not removing many trials
    number_of_short_trials = sum(seq_trial_dur < 3750);
    disp([num2str(number_of_short_trials) ' trial(s) too short'])
    seq_trial_dur(seq_trial_dur < 3750) = NaN;
    min_seq_trial_dur = min(seq_trial_dur);
end
num_seq_trials = length(successful_sequence_trials); %number of sequence trials

%now only take successful sequence trials
seq_fixationstats = seq_fixationstats(successful_sequence_trials);
seq_cfg.trl = seq_cfg.trl(successful_sequence_trials);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Grab Eye and Spike Times for Both Tasks---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for unit =1:num_units
    if all(isnan(valid_trials(:,unit)))
        continue %no data for this neuron
    end
    
    whole_spike_times{unit} = cell2mat(data(unit).values(valid_trials(1,unit):valid_trials(2,unit)));
    whole_spike_count(unit)= sum(whole_spike_times{unit});
    if whole_spike_count(unit) < min_spike_count %oo few spikes
        continue
    end
    
    %---List Image Task---%
    img_ind = 1;
    list_eye_times{unit} = NaN(192,imgdur);
    list_spike_times{unit} = NaN(192,imgdur);
    num_trials = length(cfg.trl);
    for trial = 1:num_trials
        if trial >= valid_trials(1,unit) && trial <= valid_trials(2,unit) %only valid trials for this unit
            if any(cfg.trl(trial).allval == img_on_code) && itmlist(cfg.trl(trial).cnd-1000) > sequence_items(end) %only want image trials
                %---Get Trial Info---%
                trial_start = cfg.trl(trial).alltim(cfg.trl(trial).allval == trial_start_code);%start of ITI
                imgon = cfg.trl(trial).alltim(cfg.trl(trial).allval == img_on_code)-trial_start; %when image turned on
                imgon = imgon+image_on_twin;%don't use first 500 ms due to large visual response
                %imgoff = cfg.trl(trial).alltim(cfg.trl(trial).allval == img_off_code)-trial_start;%when image turned off
                imgoff = imgon+imgdur;%minimum image viewing period
                
                %---Get Spike Times---%
                spikes = data(unit).values{trial}; %trial spike trains
                list_spike_times{unit}(img_ind,:) = spikes(imgon:imgon+imgdur-1);
                
                %---Grab Fixation Start Times---%
                fixationtimes = fixationstats{trial}.fixationtimes; %fixtaion start and end times
                
                %remove fixations that started before image turned on
                invalid= find(fixationtimes(1,:) <= imgon);
                fixationtimes(:,invalid) = [];
                
                %remove fixations started after image turned off
                invalid= find(fixationtimes(1,:) > imgoff);
                fixationtimes(:,invalid) = [];
                
                fixationtimes = fixationtimes-imgon;%zero to item1 on
                
                temp_eye = zeros(1,imgdur);
                temp_eye(fixationtimes(1,:)) = 1;
                list_eye_times{unit}(img_ind,:) = temp_eye;
                
                img_ind = img_ind+1;
            end
        end
    end
    
    %---Sequence Task---%
    seq_eye_times{unit} = NaN(length(seq_cfg.trl),min_seq_trial_dur);
    seq_spike_times{unit} = NaN(length(seq_cfg.trl),min_seq_trial_dur);
    num_trials = length(seq_cfg.trl);
    for trial = 1:num_trials
        if isnan(seq_trial_dur(trial))
            continue
        end
        if successful_sequence_trials(trial) >= valid_trials(1,unit) && successful_sequence_trials(trial) <= valid_trials(2,unit)
            %---Get Trial Info---%
            trial_start = seq_cfg.trl(trial).alltim(seq_cfg.trl(trial).allval == trial_start_code);%start of ITI
            item1on = seq_cfg.trl(trial).alltim(seq_cfg.trl(trial).allval == img_on_code)-trial_start; %when image turned on
            %reward = seq_cfg.trl(trial).alltim(seq_cfg.trl(trial).allval == reward_code)-trial_start;
            %reward = reward(1);%start of reward period
            reward = item1on+min_seq_trial_dur;
            
            %---Get Spike Times---%
            spikes = data(unit).values{successful_sequence_trials(trial)}; %trial spike trains
            seq_spike_times{unit}(trial,:) = spikes(item1on:item1on+min_seq_trial_dur-1);
            
            %---Grab Fixation Start Times---%
            fixationtimes = seq_fixationstats{trial}.fixationtimes; %fixtaion start and end times
            
            %remove fixations that started before image turned on
            invalid= find(fixationtimes(1,:) <= item1on);
            fixationtimes(:,invalid) = [];
            
            %remove fixations started after image turned off
            invalid= find(fixationtimes(1,:) > reward);
            fixationtimes(:,invalid) = [];
            
            fixationtimes = fixationtimes-item1on;%zero to item1 on
            
            temp_eye = zeros(1,min_seq_trial_dur);
            temp_eye(fixationtimes(1,:)) = 1;
            seq_eye_times{unit}(trial,:) = temp_eye;
        end
    end
end
%---Remove Excess NaNs---%
list_spike_times = laundry(list_spike_times);
list_eye_times = laundry(list_eye_times);
seq_spike_times = laundry(seq_spike_times);
seq_eye_times = laundry(seq_eye_times);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Observed Values---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for unit = 1:num_units
    if ~isempty(list_spike_times{unit})
        list_spike_count = sum(list_spike_times{unit}(:));
        seq_spike_count = sum(seq_spike_times{unit}(:));
        
        %---List Image Period Only---%
        if list_spike_count > min_spike_count %at least ~100 spikes
            [list_spike_ac_observed_ac{unit},list_spike_ac_observed_FFT{unit}] = ...
                calculate_trial_averaged_ac(list_spike_times{unit},twin,NFFT,freq_filt);
            [list_eye_ac_observed_ac{unit},list_eye_ac_observed_FFT{unit}] = ...
                calculate_trial_averaged_ac(list_eye_times{unit},twin,NFFT,freq_filt);
            [list_spike_eye_xc_observed_xc{unit},list_spike_eye_xc_observed_FFT{unit},...
                list_spike_eye_xc_observed_max_xc(unit)] =calculate_trial_averaged_xc...
                (list_spike_times{unit},list_eye_times{unit},twin,NFFT,freq_filt);
        end
        
        %---Sequence Trial Period Only---%
        if seq_spike_count > min_spike_count %at least ~100 spikes
            [seq_spike_ac_observed_ac{unit},seq_spike_ac_observed_FFT{unit}] = ...
                calculate_trial_averaged_ac(seq_spike_times{unit},twin,NFFT,freq_filt);
            [seq_eye_ac_observed_ac{unit},seq_eye_ac_observed_FFT{unit}] = ...
                calculate_trial_averaged_ac(seq_eye_times{unit},twin,NFFT,freq_filt);
            [seq_spike_eye_xc_observed_xc{unit},seq_spike_eye_xc_observed_FFT{unit},...
                seq_spike_eye_xc_observed_max_xc(unit)] =calculate_trial_averaged_xc...
                (seq_spike_times{unit},seq_eye_times{unit},twin,NFFT,freq_filt);
        end
        
        %---Whole Session---%
        if whole_spike_count(unit) > min_spike_count
            [whole_spike_ac_observed_ac{unit}, whole_spike_ac_observed_FFT{unit}] =...
                calculate_whole_ac(whole_spike_times{unit},twin,NFFT,freq_filt,freq_filt2);
        end
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Shuffled Values---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for unit = 1:num_units
    if whole_spike_count(unit) > min_spike_count
        whole_spike_ac_significant_freq{unit} = NaN(2,5);
        last_loc = 0;
        for freq = 1:length(lower_bound_frequencies)
            if freq == 1%lowest frequency band may need whitening
                %frequencies = f2(lower_bound_frequencies_ind(freq)-Hz1:lower_bound_frequencies_ind(freq+1)+Hz1);
                PSD =  whole_spike_ac_observed_FFT{unit}(1,lower_bound_frequencies_ind(freq)-Hz1:lower_bound_frequencies_ind(freq+1)+Hz1);
                %PSD = PSD.*frequencies;%whiten ~1/f, really only for lower frequencies
                [pks,locs] = findpeaks(PSD,'MinPeakDistance',Hz4);
                pks(locs < Hz1) = [];
                locs(locs < Hz1) = [];%before band of interest
                pks(locs < last_loc) = [];%too close to last loc
                locs(locs < last_loc) = [];%too close to last loc
                pks(locs > length(PSD)-Hz1) = []; %after band of interest
                locs(locs > length(PSD)-Hz1) = []; %after band of interest
                if isempty(locs)%no peaks don't even run
                    continue
                else
%                     real_locs = NaN(1,length(locs));
%                     PSD =  whole_spike_ac_observed_FFT{unit}(1,lower_bound_frequencies_ind(freq)-Hz1:lower_bound_frequencies_ind(freq+1)+Hz1);
%                     for lc = 1:length(locs)%since whitened need to correct for shifted peak
%                         [~,maxind] = max(PSD(locs(lc)-2*Hz1:locs(lc)+2*Hz1));
%                         real_locs(lc) = lower_bound_frequencies_ind(freq)+maxind-1;
%                     end
                    [shuffled_ac,shuffled_FFT]...
                        = whole_jitter_window_shuffling(PSD_numshuffs,whole_spike_times{unit},...
                        NFFT,twin,freq_filt,lower_bound_frequencies(freq),smFreqVal);
                end
            elseif freq == 5 %highest frequency band with larger smoothing window
                PSD =  whole_spike_ac_observed_FFT{unit}(2,lower_bound_frequencies_ind(freq):Hz400);
                [pks,locs] = findpeaks(PSD,'MinPeakDistance',Hz4);
                pks(locs < last_loc) = [];%too close to last loc
                locs(locs < last_loc) = [];%too close to last loc
                if isempty(locs)%no peaks don't even run
                    continue
                else
                    [shuffled_ac,shuffled_FFT]...
                        = whole_jitter_window_shuffling(PSD_numshuffs,whole_spike_times{unit},...
                        NFFT,twin,freq_filt,lower_bound_frequencies(freq),smFreqVal);
                end
                
            else
                PSD =  whole_spike_ac_observed_FFT{unit}(1,lower_bound_frequencies_ind(freq)-Hz1:lower_bound_frequencies_ind(freq+1)+Hz1);
                [pks,locs] = findpeaks(PSD,'MinPeakDistance',Hz4);
                pks(locs < Hz1) = [];
                locs(locs < Hz1) = [];%before band of interest
                pks(locs < last_loc) = [];%too close to last loc
                locs(locs < last_loc) = [];%too close to last loc
                pks(locs > length(PSD)-Hz1) = []; %after band of interest
                locs(locs > length(PSD)-Hz1) = []; %after band of interest
                if isempty(locs)%no peaks don't even run
                    continue
                else
                    [shuffled_ac,shuffled_FFT]...
                        = whole_jitter_window_shuffling(PSD_numshuffs,whole_spike_times{unit},...
                        NFFT,twin,freq_filt,lower_bound_frequencies(freq),smFreqVal);
                end
            end
            
            whole_spike_ac_shuffled_ac_mean{unit,freq} = mean(shuffled_ac);
            whole_spike_ac_shuffled_FFT_99{unit,freq} = prctile(shuffled_FFT,99,1);
            
            if length(locs) > 1
                %take largest peak
                [~,sp] = sort(pks);
                locs = locs(sp);
                locs = locs(end);
            end
            last_loc = locs-(length(PSD)-2*Hz1)+Hz4;
            if freq == 5
                freq_ind = locs+lower_bound_frequencies_ind(freq)-1;
                cluster_observed_FFT = sum(whole_spike_ac_observed_FFT{unit}(2,freq_ind-Hz1:freq_ind+Hz1),2);
            else
                freq_ind = locs+lower_bound_frequencies_ind(freq)-Hz1-1;
                cluster_observed_FFT = sum(whole_spike_ac_observed_FFT{unit}(1,freq_ind-Hz1:freq_ind+Hz1),2);
                
            end
            
            cluster_shuffled_FFT = sum(shuffled_FFT(:,freq_ind-Hz1:freq_ind+Hz1),2);
            if 100*sum(cluster_observed_FFT > cluster_shuffled_FFT)/PSD_numshuffs > pthresh
                whole_spike_ac_significant_freq{unit}(:,freq) = [f2(freq_ind),1];
            else
                whole_spike_ac_significant_freq{unit}(:,freq) = [f2(freq_ind),0];
            end
            
        end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot all of the Results---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq_titles = {'Theta Range','Alpha Range','Beta Range','Gamma Range','HFO/Burst Range'};
xlimits = [3 11;9 21;19 31; 30 100;100 400];
clrs = lines(5);
for unit = 1:num_units
    if whole_spike_count(unit) < min_spike_count
        continue
    end
    
    
    %---Plot All Observed AutoCorrelations---%
    figure
    %spike autocorrelation
    subplot(2,3,1)
    leg = [];
    hold on
    if ~isempty(list_spike_ac_observed_ac{unit})
        dofill(tm,list_spike_ac_observed_ac{unit}-mean(list_spike_ac_observed_ac{unit}(:)),'r',1,20);
        leg = {'Images'};
    end
    if ~isempty(seq_spike_ac_observed_ac{unit})
        dofill(tm,seq_spike_ac_observed_ac{unit}-mean(seq_spike_ac_observed_ac{unit}(:)),'b',1,20);
        leg = [leg {'Sequence'}];
    end
    dofill(tm,whole_spike_ac_observed_ac{unit},'k',1,20);
    hold off
    xlim([-twin twin])
    set(gca,'Xtick',[-500 -250 0 250 500])
    ylabel('Correlation (a.u.)')
    xlabel('lag (ms)')
    title('Spike AutoCorrelation, 10 ms smooth')
    legend([leg,'Whole'])
    %
    %     %spike x eye cross correlation
    %     subplot(2,3,4)
    %     hold on
    %     leg = [];
    %     if ~isempty(list_spike_ac_observed_ac{unit})
    %         dofill(tm,list_spike_eye_xc_observed_xc{unit},'g',1,20);
    %         leg = {'Images'};
    %     end
    %     if ~isempty(seq_spike_ac_observed_ac{unit})
    %         dofill(tm,seq_spike_eye_xc_observed_xc{unit},'m',1,20);
    %         leg = [leg {'Sequence'}];
    %     end
    %     hold off
    %     xlim([-twin twin])
    %     set(gca,'Xtick',[-500 -250 0 250 500])
    %     ylabel('Correlation (a.u.)')
    %     xlabel('lag (ms)')
    %     title('Spike/Eye CrossCorrelation, 10 ms smooth')
    %     legend(leg)
    
    for sb = 2:6
        subplot(2,3,sb)
        hold on
        if ~isempty(list_spike_ac_observed_ac{unit})
            %plot(f2,list_spike_eye_xc_observed_FFT{unit},'g');%has too much power
            plot(f2,list_spike_ac_observed_FFT{unit},'r');
        end
        if ~isempty(seq_spike_ac_observed_ac{unit})
            %plot(f2,seq_spike_eye_xc_observed_FFT{unit},'m');%has too much power
            plot(f2,seq_spike_ac_observed_FFT{unit},'b');
        end
        if sb == 6
            plot(f2,whole_spike_ac_observed_FFT{unit}(2,:),'k')
        else
            plot(f2,whole_spike_ac_observed_FFT{unit}(1,:),'k')
        end
        xlim([xlimits(sb-1,1) xlimits(sb-1,2)])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        hold off
        title(freq_titles(freq))
        
    end
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Observed_AutoCorrelations'])
    %%
%     figure
%     subplot(2,3,1)
%     hold on
%     for freq = 1:5
%         if  ~isnan(whole_spike_ac_significant_freq{unit}(1,freq))
%             dofill(tm,whole_spike_ac_shuffled_ac_mean{unit,freq},clrs(freq,:),1,20);
%         end
%     end
%     dofill(tm,whole_spike_ac_observed_ac{unit},'k',1,20);
%     hold off
%     xlim([-twin twin])
%     set(gca,'Xtick',[-500 -250 0 250 500])
%     ylabel('Correlation (a.u.)')
%     xlabel('lag (ms)')
%     title('Spike AutoCorrelation, 10 ms smooth')
%  
%     for sb = 2:6
%         subplot(2,3,sb)
%         %for freq = 1:5
%         freq = sb-1;
%         if  ~isnan(whole_spike_ac_significant_freq{unit}(1,freq))
%             p = plot(f2,whole_spike_ac_shuffled_FFT_99{unit,freq});
%             set(p,'color',clrs(freq,:))
%             hold on
%             freq_ind = find(f2 == whole_spike_ac_significant_freq{unit}(1,freq));
%             if freq == 5
%                 pow = whole_spike_ac_observed_FFT{unit}(2,freq_ind);
%             else
%                 pow = whole_spike_ac_observed_FFT{unit}(1,freq_ind);
%             end
%             if whole_spike_ac_significant_freq{unit}(2,freq) == 1%significant
%                 plot(f2(freq_ind),pow,['rx'],'Markersize',15)
%             else
%                 plot(f2(freq_ind),pow,['ro'],'Markersize',15)
%             end
%         end
%         %end
%         
%         if sb == 6
%             plot(f2,whole_spike_ac_observed_FFT{unit}(2,:),'k','linewidth',3)
%         else
%             plot(f2,whole_spike_ac_observed_FFT{unit}(1,:),'k','linewidth',3)
%         end
%         xlim([xlimits(freq,1) xlimits(freq,2)])
%         yl = ylim;
%         ylim([0 yl(2)])
%         xlabel('Frequency (Hz)')
%         ylabel('Relative Power')
%         box off
%         hold off
%         
%         if whole_spike_ac_significant_freq{unit}(2,freq) == 1%significant
%              freq_ind = find(f2 == whole_spike_ac_significant_freq{unit}(1,freq));
%             title(sprintf([freq_titles{freq} ', peak @ ' num2str(f2(freq_ind),3) ' Hz']))
%         else
%             title(freq_titles{freq})
%         end
%     end
% 
%     %%
%    
%     save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Whole_AutoCorrelations_Analysis'])
    
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Save all of the Results---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save([data_dir task_file(1:8) '-Freq_Mod_Analysis.mat'],...
%     'twin','smFreqVal','smFreqVal','NFFT','f2','min_spike_count','PSD_numshuffs','unit_stats',...
%     'spatial_info','pthresh','lower_bound_frequencies','lower_bound_frequencies_ind',...
%     'freq_filt','freq_filt2','whole_spike_times','whole_spike_count','whole_spike_ac_observed_ac',...
%     'whole_spike_ac_observed_FFT','whole_spike_ac_shuffled_ac_mean','whole_spike_ac_shuffled_FFT_99',...
%     'whole_spike_ac_significant_freq','list_spike_times','list_spike_ac_observed_ac',...
%     'list_spike_ac_observed_FFT','seq_spike_times','seq_spike_ac_observed_ac',...
%     'seq_spike_ac_observed_FFT','list_eye_times','list_eye_ac_observed_ac',...
%     'list_eye_ac_observed_FFT','seq_eye_times','seq_eye_ac_observed_ac',...
%     'seq_eye_ac_observed_FFT','list_spike_eye_xc_observed_xc','list_spike_eye_xc_observed_max_xc',...
%     'list_spike_eye_xc_observed_FFT','seq_spike_eye_xc_observed_xc',...
%     'seq_spike_eye_xc_observed_max_xc','seq_spike_eye_xc_observed_FFT',...
%     'Hz1','Hz4','Hz400');
%%
end