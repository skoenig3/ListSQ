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
min_spike_count = 100;%don't analyze epoch if not enough spikes, can be both or either List/Sequence
image_on_twin = 500;%ignore first 500 ms
img_on_code = 23;%also first item on code in sequence task
img_off_code = 24;
reward_code = 3;
trial_start_code = 15;
%load whitening filter used to detect peaks
load(['C:\Users\seth.koenig\Documents\MATLAB\ListSQ\AutoCorrWhitenFilter.mat']);
whitening_filter(1) = 10*max(whitening_filter);
whitening_filter = whitening_filter(1:16:end);%FFT was originally sampled at 2^16 

%---PSD Autocorrelation Parmaeters---%
PSD_numshuffs = 1000;%10000;
xc_PSD_numshuffs =1000;%cross correlation between eye movements and spikes
NFFT = 2^12;%number of frequency points from 0-Fs/2
f2 = Fs/2*linspace(0,1,NFFT/2+1);%frequency values
smFreqVal = 2;%how much smoothing is done in Hz e.g. +/- 1.0 Hz
freq_filt = ones(1,round(smFreqVal/mode(diff(f2)))); %square window/moving average
freq_filt = freq_filt/sum(freq_filt);%zero sum filter


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

%---List Spike Autocorrelation Analysis---%
list_spike_times = cell(1,num_units);
list_spike_ac_observed_ac = cell(1,num_units);
list_spike_ac_observed_FFT = cell(1,num_units);
list_spike_shuffled_ac = cell(5,num_units);
list_spike_ac_shuffled_FFT = cell(5,num_units);
list_spike_ac_shuffled_FFT_99 = cell(5,num_units);
list_spike_ac_freq_of_interest = cell(5,num_units);

%---Sequence Spike Autocorrelation Analysis---%
seq_spike_times = cell(1,num_units);
seq_spike_ac_observed_ac = cell(1,num_units);
seq_spike_ac_observed_FFT = cell(1,num_units);
seq_spike_shuffled_ac = cell(5,num_units);
seq_spike_ac_shuffled_FFT = cell(5,num_units);
seq_spike_ac_shuffled_FFT_99 = cell(5,num_units);
seq_spike_ac_freq_of_interest = cell(5,num_units);

%---List Eye Autocorrelation Analysis---%
list_eye_times = cell(1,num_units);
list_eye_ac_observed_ac = cell(1,num_units);
list_eye_ac_observed_FFT = cell(1,num_units);
list_eye_shuffled_ac = cell(3,num_units);
list_eye_ac_shuffled_FFT = cell(3,num_units);
list_eye_ac_shuffled_FFT_99 = cell(3,num_units);
list_eye_ac_freq_of_interest = cell(3,num_units);

%---Sequence Eye Autocorrelation Analysis---%
seq_eye_times = cell(1,num_units);
seq_eye_ac_observed_ac = cell(1,num_units);
seq_eye_ac_observed_FFT = cell(1,num_units);
seq_eye_shuffled_ac = cell(3,num_units);
seq_eye_ac_shuffled_FFT = cell(3,num_units);
seq_eye_ac_shuffled_FFT_99 = cell(3,num_units);
seq_eye_ac_freq_of_interest = cell(3,num_units);

%---List Spike-Eye Crosscorrelation Analysis---%
list_spike_eye_xc_observed_xc = cell(1,num_units);
list_spike_eye_xc_observed_max_xc = NaN(1,num_units);
list_spike_eye_xc_observed_FFT = cell(1,num_units);
list_spike_eye_xc_shuffled_xc = cell(1,num_units);
list_spike_eye_xc_shuffled_max_xc = cell(1,num_units);
list_spike_eye_xc_shuffled_FFT = cell(1,num_units);
list_spike_eye_xc_shuffled_FFT_99 = cell(1,num_units);

%---Sequence Spike-Eye Crosscorrelation Analysis---%
seq_spike_eye_xc_observed_xc = cell(1,num_units);
seq_spike_eye_xc_observed_max_xc = NaN(1,num_units);
seq_spike_eye_xc_observed_FFT = cell(1,num_units);
seq_spike_eye_xc_shuffled_xc = cell(1,num_units);
seq_spike_eye_xc_shuffled_max_xc = cell(1,num_units);
seq_spike_eye_xc_shuffled_FFT = cell(1,num_units);
seq_spike_eye_xc_shuffled_FFT_99 = cell(1,num_units);

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
        
        if list_spike_count > min_spike_count %at least ~100 spikes
            [list_spike_ac_observed_ac{unit},list_spike_ac_observed_FFT{unit}] = ...
                calculate_trial_averaged_ac(list_spike_times{unit},twin,NFFT,freq_filt);
            [list_eye_ac_observed_ac{unit},list_eye_ac_observed_FFT{unit}] = ...
                calculate_trial_averaged_ac(list_eye_times{unit},twin,NFFT,freq_filt);
            [list_spike_eye_xc_observed_xc{unit},list_spike_eye_xc_observed_FFT{unit},...
                list_spike_eye_xc_observed_max_xc(unit)] =calculate_trial_averaged_xc...
                (list_spike_times{unit},list_eye_times{unit},twin,NFFT,freq_filt);
        end
        if seq_spike_count > min_spike_count %at least ~100 spikes
            [seq_spike_ac_observed_ac{unit},seq_spike_ac_observed_FFT{unit}] = ...
                calculate_trial_averaged_ac(seq_spike_times{unit},twin,NFFT,freq_filt);
            [seq_eye_ac_observed_ac{unit},seq_eye_ac_observed_FFT{unit}] = ...
                calculate_trial_averaged_ac(seq_eye_times{unit},twin,NFFT,freq_filt);
            [seq_spike_eye_xc_observed_xc{unit},seq_spike_eye_xc_observed_FFT{unit},...
                seq_spike_eye_xc_observed_max_xc(unit)] =calculate_trial_averaged_xc...
                (seq_spike_times{unit},seq_eye_times{unit},twin,NFFT,freq_filt);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Shuffled Values---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f100 = find(f2 >= 100);
f100 = f100(1);
f450 = find(f2 >= 450);
f450 = f450(1);
Hz4 = round(4/mode(diff(f2)));%peaks with minimum distance of 4 Hz away

for unit = 1:num_units
    if ~isempty(list_spike_ac_observed_ac{unit})
        
        %Spike AutoCorrelation
        list_spike_ac_freq_of_interest{unit} = find_frequencies_of_interest(...
            list_spike_ac_observed_FFT{unit},whitening_filter,f2,f100,f450,Hz4);
        for freq = 1:length(list_spike_ac_freq_of_interest{unit})
            [list_spike_shuffled_ac{freq,unit},list_spike_ac_shuffled_FFT{freq,unit}]...
                = jitter_window_shuffling(PSD_numshuffs,list_spike_times{unit},NFFT,twin,...
                1,freq_filt,list_spike_ac_freq_of_interest{unit}(freq),smFreqVal);
            list_spike_ac_shuffled_FFT_99{freq,unit} = prctile(list_spike_ac_shuffled_FFT{freq,unit},99,1);
        end
        
        %Eye AutoCorrelation
        list_eye_ac_freq_of_interest{unit} = find_frequencies_of_interest_eye(...
            list_eye_ac_observed_FFT{unit},f2,f100,Hz4);
        for freq = 1:length(list_eye_ac_freq_of_interest{unit})
            [list_eye_shuffled_ac{freq,unit},list_eye_ac_shuffled_FFT{freq,unit}]...
                = jitter_window_shuffling(PSD_numshuffs,list_eye_times{unit},NFFT,twin,...
                1,freq_filt,list_eye_ac_freq_of_interest{unit}(freq),smFreqVal);
            list_eye_ac_shuffled_FFT_99{freq,unit} = prctile(list_eye_ac_shuffled_FFT{freq,unit},99,1);
        end
        
        
        %Spike/Eye Cross Correlation
        [list_spike_eye_xc_shuffled_xc{unit},list_spike_eye_xc_shuffled_FFT{unit},...
            list_spike_eye_xc_shuffled_max_xc{unit}] = shuffle_cross_correlation(...
            list_spike_times{unit},list_eye_times{unit},twin,NFFT,freq_filt,xc_PSD_numshuffs);
        list_spike_eye_xc_shuffled_FFT_99{unit} = ...
            prctile(list_spike_eye_xc_shuffled_FFT{unit},95,1);
    end
    if ~isempty(seq_spike_ac_observed_ac{unit})
        
        %Spike AutoCorrelation
        seq_spike_ac_freq_of_interest{unit} = find_frequencies_of_interest(...
            seq_spike_ac_observed_FFT{unit},whitening_filter,f2,f100,f450,Hz4);
        for freq = 1:length(seq_spike_ac_freq_of_interest{unit})
            [seq_spike_shuffled_ac{freq,unit},seq_spike_ac_shuffled_FFT{freq,unit}]...
                = jitter_window_shuffling(PSD_numshuffs,seq_spike_times{unit},NFFT,twin,...
                1,freq_filt,seq_spike_ac_freq_of_interest{unit}(freq),smFreqVal);
            seq_spike_ac_shuffled_FFT_99{freq,unit} = prctile(seq_spike_ac_shuffled_FFT{freq,unit},99,1);
        end
        
        %Eye AutoCorrelation
        seq_eye_ac_freq_of_interest{unit} = find_frequencies_of_interest_eye(...
            seq_eye_ac_observed_FFT{unit},f2,f100,Hz4);
        for freq = 1:length(seq_eye_ac_freq_of_interest{unit})
            [seq_eye_shuffled_ac{freq,unit},seq_eye_ac_shuffled_FFT{freq,unit}]...
                = jitter_window_shuffling(PSD_numshuffs,seq_eye_times{unit},NFFT,twin,...
                1,freq_filt,seq_eye_ac_freq_of_interest{unit}(freq),smFreqVal);
            seq_eye_ac_shuffled_FFT_99{freq,unit} = prctile(seq_eye_ac_shuffled_FFT{freq,unit},99,1);
        end
        
        
        %Spike/Eye Cross Correlation
        [seq_spike_eye_xc_shuffled_xc{unit},seq_spike_eye_xc_shuffled_FFT{unit},...
            seq_spike_eye_xc_shuffled_max_xc{unit}] = shuffle_cross_correlation(...
            seq_spike_times{unit},seq_eye_times{unit},twin,NFFT,freq_filt,xc_PSD_numshuffs);
        seq_spike_eye_xc_shuffled_FFT_99{unit} = ...
            prctile(seq_spike_eye_xc_shuffled_FFT{unit},95,1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot all of the Results---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for unit = 1:num_units
    if isempty(list_spike_times{unit})
        continue
    end
    %%
    %---Plot Spike AutoCorrelation Results---%
    %for list
    if ~isempty(list_spike_shuffled_ac{1,unit})
        figure
        subplot(2,3,1)
        hold on
        for freq = 5:-1:1
            dofill(tm,list_spike_shuffled_ac{freq,unit},clrs(freq),1,20);
        end
        dofill(tm,list_spike_ac_observed_ac{unit},'k',1,20);
        hold off
        xlim([-twin twin])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Correlation (a.u.)')
        xlabel('lag (ms)')
        title('Spike AutoCorr, 10 ms smooth')
        
        subplot(2,3,4)
        hold on
        for freq = 5:-1:1
            dofill(tm,list_spike_shuffled_ac{freq,unit},clrs(freq),1,4);
        end
        dofill(tm,list_spike_ac_observed_ac{unit},'k',1,4);
        hold off
        xlim([-25 25])
        ylabel('Correlation (a.u.)')
        xlabel('lag (ms)')
        title('Spike AutoCorr, 2 ms smooth')
        
        subplot(2,3,2)
        hold on
        for freq = 5:-1:1
            plot(f2,list_spike_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == list_spike_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),list_spike_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,list_spike_ac_observed_FFT{unit},'k')
        hold off
        xlim([0 12])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Delta/Theta Range')
        
        subplot(2,3,3)
        hold on
        for freq = 5:-1:1
            plot(f2,list_spike_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == list_spike_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),list_spike_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,list_spike_ac_observed_FFT{unit},'k')
        hold off
        xlim([8 30])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Alpha/Beta Range')
        
        subplot(2,3,5)
        hold on
        for freq = 5:-1:1
            plot(f2,list_spike_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == list_spike_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),list_spike_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,list_spike_ac_observed_FFT{unit},'k')
        hold off
        xlim([30 100])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Gamma Range')
        
        subplot(2,3,6)
        hold on
        for freq = 5:-1:1
            plot(f2,list_spike_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == list_spike_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),list_spike_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,list_spike_ac_observed_FFT{unit},'k')
        hold off
        xlim([100 450])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('HFO/Burst Range')
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_List_Spike_AutoCorr_analysis'])
    end
    
    %for Sequence
    if ~isempty(seq_spike_shuffled_ac{1,unit})
        
        figure
        subplot(2,3,1)
        hold on
        for freq = 5:-1:1
            dofill(tm,seq_spike_shuffled_ac{freq,unit},clrs(freq),1,20);
        end
        dofill(tm,seq_spike_ac_observed_ac{unit},'k',1,20);
        hold off
        xlim([-twin twin])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Correlation (a.u.)')
        xlabel('lag (ms)')
        title('Spike AutoCorr, 10 ms smooth')
        
        subplot(2,3,4)
        hold on
        for freq = 5:-1:1
            dofill(tm,seq_spike_shuffled_ac{freq,unit},clrs(freq),1,4);
        end
        dofill(tm,seq_spike_ac_observed_ac{unit},'k',1,4);
        hold off
        xlim([-25 25])
        ylabel('Correlation (a.u.)')
        xlabel('lag (ms)')
        title('Spike AutoCorr, 2 ms smooth')
        
        subplot(2,3,2)
        hold on
        for freq = 5:-1:1
            plot(f2,seq_spike_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == seq_spike_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),seq_spike_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,seq_spike_ac_observed_FFT{unit},'k')
        hold off
        xlim([0 12])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Delta/Theta Range')
        
        subplot(2,3,3)
        hold on
        for freq = 5:-1:1
            plot(f2,seq_spike_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == seq_spike_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),seq_spike_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,seq_spike_ac_observed_FFT{unit},'k')
        hold off
        xlim([8 30])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Alpha/Beta Range')
        
        subplot(2,3,5)
        hold on
        for freq = 5:-1:1
            plot(f2,seq_spike_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == seq_spike_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),seq_spike_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,seq_spike_ac_observed_FFT{unit},'k')
        hold off
        xlim([30 100])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Gamma Range')
        
        subplot(2,3,6)
        hold on
        for freq = 5:-1:1
            plot(f2,seq_spike_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == seq_spike_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),seq_spike_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,seq_spike_ac_observed_FFT{unit},'k')
        hold off
        xlim([100 450])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('HFO/Burst Range')
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Sequence_Spike_AutoCorr_analysis'])
    end
    %%
    %---Plot eye AutoCorrelation Results---%
    if ~isempty(list_eye_shuffled_ac{1,unit})
        %for list
        figure
        subplot(2,3,1)
        hold on
        for freq = length(list_eye_ac_freq_of_interest{unit}):-1:1
            dofill(tm,list_eye_shuffled_ac{freq,unit},clrs(freq),1,20);
        end
        dofill(tm,list_eye_ac_observed_ac{unit},'k',1,20);
        hold off
        xlim([-twin twin])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Correlation (a.u.)')
        xlabel('lag (ms)')
        title('eye AutoCorr, 10 ms smooth')
        
        subplot(2,3,4)
        hold on
        for freq = length(list_eye_ac_freq_of_interest{unit}):-1:1
            dofill(tm,list_eye_shuffled_ac{freq,unit},clrs(freq),1,4);
        end
        dofill(tm,list_eye_ac_observed_ac{unit},'k',1,4);
        hold off
        xlim([-25 25])
        ylabel('Correlation (a.u.)')
        xlabel('lag (ms)')
        title('eye AutoCorr, 2 ms smooth')
        
        subplot(2,3,2)
        hold on
        for freq = length(list_eye_ac_freq_of_interest{unit}):-1:1
            plot(f2,list_eye_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == list_eye_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),list_eye_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,list_eye_ac_observed_FFT{unit},'k')
        hold off
        xlim([0 12])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Delta/Theta Range')
        
        subplot(2,3,3)
        hold on
        for freq = length(list_eye_ac_freq_of_interest{unit}):-1:1
            plot(f2,list_eye_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == list_eye_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),list_eye_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,list_eye_ac_observed_FFT{unit},'k')
        hold off
        xlim([8 30])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Alpha/Beta Range')
        
        subplot(2,3,5)
        hold on
        for freq = length(list_eye_ac_freq_of_interest{unit}):-1:1
            plot(f2,list_eye_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == list_eye_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),list_eye_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,list_eye_ac_observed_FFT{unit},'k')
        hold off
        xlim([30 100])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Gamma Range')
        
        subplot(2,3,6)
        hold on
        for freq = length(list_eye_ac_freq_of_interest{unit}):-1:1
            plot(f2,list_eye_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == list_eye_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),list_eye_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,list_eye_ac_observed_FFT{unit},'k')
        hold off
        xlim([100 450])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('HFO/Burst Range')
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_List_eye_AutoCorr_analysis'])
    end
    
    %for Sequence
    if ~isempty(seq_eye_shuffled_ac{1,unit})
        
        figure
        subplot(2,3,1)
        hold on
        for freq = length(list_eye_ac_freq_of_interest{unit}):-1:1
            dofill(tm,seq_eye_shuffled_ac{freq,unit},clrs(freq),1,20);
        end
        dofill(tm,seq_eye_ac_observed_ac{unit},'k',1,20);
        hold off
        xlim([-twin twin])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Correlation (a.u.)')
        xlabel('lag (ms)')
        title('eye AutoCorr, 10 ms smooth')
        
        subplot(2,3,4)
        hold on
        for freq = length(list_eye_ac_freq_of_interest{unit}):-1:1
            dofill(tm,seq_eye_shuffled_ac{freq,unit},clrs(freq),1,4);
        end
        dofill(tm,seq_eye_ac_observed_ac{unit},'k',1,4);
        hold off
        xlim([-25 25])
        ylabel('Correlation (a.u.)')
        xlabel('lag (ms)')
        title('eye AutoCorr, 2 ms smooth')
        
        subplot(2,3,2)
        hold on
        for freq = length(list_eye_ac_freq_of_interest{unit}):-1:1
            plot(f2,seq_eye_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == seq_eye_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),seq_eye_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,seq_eye_ac_observed_FFT{unit},'k')
        hold off
        xlim([0 12])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Delta/Theta Range')
        
        subplot(2,3,3)
        hold on
        for freq = length(list_eye_ac_freq_of_interest{unit}):-1:1
            plot(f2,seq_eye_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == seq_eye_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),seq_eye_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,seq_eye_ac_observed_FFT{unit},'k')
        hold off
        xlim([8 30])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Alpha/Beta Range')
        
        subplot(2,3,5)
        hold on
        for freq = length(list_eye_ac_freq_of_interest{unit}):-1:1
            plot(f2,seq_eye_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == seq_eye_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),seq_eye_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,seq_eye_ac_observed_FFT{unit},'k')
        hold off
        xlim([30 100])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Gamma Range')
        
        subplot(2,3,6)
        hold on
        for freq = length(list_eye_ac_freq_of_interest{unit}):-1:1
            plot(f2,seq_eye_ac_shuffled_FFT_99{freq,unit},clrs(freq));
            freq_ind = find(f2 == seq_eye_ac_freq_of_interest{unit}(freq));
            plot(f2(freq_ind),seq_eye_ac_observed_FFT{unit}(freq_ind),[clrs(freq) 'o'],'markersize',5)
        end
        plot(f2,seq_eye_ac_observed_FFT{unit},'k')
        hold off
        xlim([100 450])
        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('HFO/Burst Range')
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Sequence_eye_AutoCorr_analysis'])
    end
    %%
    %---Plot Spike/Eye Cross Correlation Results---%
    if ~isempty(list_spike_eye_xc_shuffled_xc{unit})
        %for list first
        figure
        subplot(2,3,1)
        dofill(tm,list_spike_eye_xc_shuffled_xc{unit},'red',1,20);
        hold on
        dofill(tm,list_spike_eye_xc_observed_xc{unit},'black',1,20);
        hold off
        xlim([-twin twin])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Correlation (a.u.)')
        xlabel('lag (ms)')
        title('Spike/Eye XC')
        
        subplot(2,3,4)
        hist(list_spike_eye_xc_shuffled_max_xc{unit},25)
        yl = ylim;
        hold on
        plot([list_spike_eye_xc_observed_max_xc(unit) list_spike_eye_xc_observed_max_xc(unit)],...
            [yl(1) yl(2)],'g')
        hold off
        box off
        xlabel('Maximum Correlation')
        ylabel('Count')
        corr_pct = 100*sum(list_spike_eye_xc_observed_max_xc(unit) > ...
            list_spike_eye_xc_shuffled_max_xc{unit})/xc_PSD_numshuffs;
        title(['Corr. Percentile: ' num2str(corr_pct,3)])
        
        subplot(2,3,2)
        plot(f2,list_spike_eye_xc_observed_FFT{unit},'k');
        hold on
        plot(f2,list_spike_eye_xc_shuffled_FFT_99{unit},'r');
        hold off
        xlim([0 12])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Delta/Theta Range')
        
        subplot(2,3,3)
        plot(f2,list_spike_eye_xc_observed_FFT{unit},'k');
        hold on
        plot(f2,list_spike_eye_xc_shuffled_FFT_99{unit},'r');
        hold off
        xlim([8 30])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Alpha/Beta Range')
        
        subplot(2,3,5)
        plot(f2,list_spike_eye_xc_observed_FFT{unit},'k');
        hold on
        plot(f2,list_spike_eye_xc_shuffled_FFT_99{unit},'r');
        hold off
        xlim([30 100])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Gamma Range')
        
        subplot(2,3,6)
        plot(f2,list_spike_eye_xc_observed_FFT{unit},'k');
        hold on
        plot(f2,list_spike_eye_xc_shuffled_FFT_99{unit},'r');
        hold off
        xlim([100 450])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('HFO/Burst Range')
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_List_Eye_Spike_XC_analysis'])
    end
    
    %for sequence second
    if ~isempty(seq_spike_eye_xc_shuffled_xc{unit})
        
        figure
        subplot(2,3,1)
        dofill(tm,seq_spike_eye_xc_shuffled_xc{unit},'red',1,20);
        hold on
        dofill(tm,seq_spike_eye_xc_observed_xc{unit},'black',1,20);
        hold off
        xlim([-twin twin])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Correlation (a.u.)')
        xlabel('lag (ms)')
        title('Spike/Eye XC')
        
        subplot(2,3,4)
        hist(seq_spike_eye_xc_shuffled_max_xc{unit},25)
        yl = ylim;
        hold on
        plot([seq_spike_eye_xc_observed_max_xc(unit) seq_spike_eye_xc_observed_max_xc(unit)],...
            [yl(1) yl(2)],'g')
        hold off
        box off
        xlabel('Maximum Correlation')
        ylabel('Count')
        corr_pct = 100*sum(seq_spike_eye_xc_observed_max_xc(unit) > ...
            seq_spike_eye_xc_shuffled_max_xc{unit})/xc_PSD_numshuffs;
        title(['Corr. Percentile: ' num2str(corr_pct,3)])
        
        subplot(2,3,2)
        plot(f2,seq_spike_eye_xc_observed_FFT{unit},'k');
        hold on
        plot(f2,seq_spike_eye_xc_shuffled_FFT_99{unit},'r');
        hold off
        xlim([0 12])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Delta/Theta Range')
        
        subplot(2,3,3)
        plot(f2,seq_spike_eye_xc_observed_FFT{unit},'k');
        hold on
        plot(f2,seq_spike_eye_xc_shuffled_FFT_99{unit},'r');
        hold off
        xlim([8 30])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Alpha/Beta Range')
        
        subplot(2,3,5)
        plot(f2,seq_spike_eye_xc_observed_FFT{unit},'k');
        hold on
        plot(f2,seq_spike_eye_xc_shuffled_FFT_99{unit},'r');
        hold off
        xlim([30 100])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('Gamma Range')
        
        subplot(2,3,6)
        plot(f2,seq_spike_eye_xc_observed_FFT{unit},'k');
        hold on
        plot(f2,seq_spike_eye_xc_shuffled_FFT_99{unit},'r');
        hold off
        xlim([100 450])
        xlabel('Frequency (Hz')
        ylabel('Relative Power')
        box off
        title('HFO/Burst Range')
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Sequence_Eye_Spike_XC_analysis'])
    end
    
    
    %---Plot Auto/Cross-Correlations Under Different Conditions---%
    list_spike_count = sum(list_spike_times{unit}(:));
    seq_spike_count = sum(seq_spike_times{unit}(:));

    figure
    subplot(2,3,1)
    hold on
    if list_spike_count > min_spike_count
        dofill(tm,list_eye_ac_observed_ac{unit},'red',1,20);
    end
    if seq_spike_count > min_spike_count
        dofill(tm,seq_eye_ac_observed_ac{unit},'blue',1,20);
    end
    hold off
    yl  = ylim;
    ylim([0 yl(2)])
    xlim([-twin twin])
    set(gca,'Xtick',[-500 -250 0 250 500])
    ylabel('Correlation (a.u.)')
    xlabel('lag (ms)')
    legend('List','Seq','Location','NorthWest')
    title('Eye Movement Autocorrelation')
    
    subplot(2,3,4)
    hold on
    if list_spike_count > min_spike_count
        dofill(tm,list_spike_ac_observed_ac{unit},'green',1,20);
    end
    if seq_spike_count > min_spike_count
        dofill(tm,seq_spike_ac_observed_ac{unit},'black',1,20);
    end
    hold off
    yl  = ylim;
    ylim([0 yl(2)])
    xlim([-twin twin])
    set(gca,'Xtick',[-500 -250 0 250 500])
    ylabel('Correlation (a.u.)')
    legend('List','Seq','Location','SouthWest')
    title('Spike Autocorrelation')
    
    
    subplot(2,3,2)
    hold on
    if list_spike_count > min_spike_count
        dofill(tm,list_spike_eye_xc_observed_xc{unit},'cyan',1,20);
    end
    if seq_spike_count > min_spike_count
        dofill(tm,seq_spike_eye_xc_observed_xc{unit},'magenta',1,20);
    end
    hold off
    xlim([-twin twin])
    set(gca,'Xtick',[-500 -250 0 250 500])
    ylabel('Correlation (a.u.)')
    xlabel('lag (ms)')
    legend('List','Seq')
    title('Spike Autocorrelation')
    xlabel('lag (ms)')
    legend('List','Seq','Location','SouthWest')
    title('Spike/Eye Crosscorrelation')
    
    subplot(2,3,3)
    hold on
    if list_spike_count > min_spike_count
        plot(f2,list_spike_ac_observed_FFT{unit}/mean(list_spike_ac_observed_FFT{unit}),'green');
        plot(f2,list_eye_ac_observed_FFT{unit}/mean(list_eye_ac_observed_FFT{unit}),'red');
        plot(f2,list_spike_eye_xc_observed_FFT{unit}/mean(list_spike_eye_xc_observed_FFT{unit}),'cyan');
    end
    if seq_spike_count > min_spike_count
        plot(f2,seq_spike_ac_observed_FFT{unit}/mean(seq_spike_ac_observed_FFT{unit}),'black');
        plot(f2,seq_eye_ac_observed_FFT{unit}/mean(seq_eye_ac_observed_FFT{unit}),'blue');
        plot(f2,seq_spike_eye_xc_observed_FFT{unit}/mean(seq_spike_eye_xc_observed_FFT{unit}),'magenta');
    end
    hold off
    xlim([0 12])
    ylabel('Relative Power')
    xlabel('Frequency')
    legend('List Spike','List Eye','List Spike/Eye','Seq Spike','Seq Eye','Seq Spike/Eye')
    
    subplot(2,3,5)
    hold on
    if list_spike_count > min_spike_count
        plot(f2,list_spike_ac_observed_FFT{unit}/mean(list_spike_ac_observed_FFT{unit}),'green');
        plot(f2,list_eye_ac_observed_FFT{unit}/mean(list_eye_ac_observed_FFT{unit}),'red');
        plot(f2,list_spike_eye_xc_observed_FFT{unit}/mean(list_spike_eye_xc_observed_FFT{unit}),'cyan');
    end
    if seq_spike_count > min_spike_count
        plot(f2,seq_spike_ac_observed_FFT{unit}/mean(seq_spike_ac_observed_FFT{unit}),'black');
        plot(f2,seq_eye_ac_observed_FFT{unit}/mean(seq_eye_ac_observed_FFT{unit}),'blue');
        plot(f2,seq_spike_eye_xc_observed_FFT{unit}/mean(seq_spike_eye_xc_observed_FFT{unit}),'magenta');
    end
    hold off
    xlim([8 30])
    ylabel('Relative Power')
    xlabel('Frequency')
    
    subplot(2,3,6)
    hold on
    if list_spike_count > min_spike_count
        plot(f2,list_spike_ac_observed_FFT{unit}/mean(list_spike_ac_observed_FFT{unit}),'green');
        plot(f2,list_eye_ac_observed_FFT{unit}/mean(list_eye_ac_observed_FFT{unit}),'red');
        plot(f2,list_spike_eye_xc_observed_FFT{unit}/mean(list_spike_eye_xc_observed_FFT{unit}),'cyan');
    end
    if seq_spike_count > min_spike_count
        plot(f2,seq_spike_ac_observed_FFT{unit}/mean(seq_spike_ac_observed_FFT{unit}),'black');
        plot(f2,seq_eye_ac_observed_FFT{unit}/mean(seq_eye_ac_observed_FFT{unit}),'blue');
        plot(f2,seq_spike_eye_xc_observed_FFT{unit}/mean(seq_spike_eye_xc_observed_FFT{unit}),'magenta');
    end
    hold off
    xlim([30 100])
    ylabel('Relative Power')
    xlabel('Frequency')
    
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_all_observed_freq_modulation'])
end

%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Save all of the Results---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
save([data_dir task_file(1:8) '-Freq_Mod_Analysis.mat'],...
    'twin','smFreqVal','NFFT','f2','min_spike_count','PSD_numshuffs',...
    'xc_PSD_numshuffs','unit_stats','spatial_info','list_spike_times',...
    'list_eye_times','seq_spike_times','seq_eye_times','whitening_filter',...
'list_spike_times','list_spike_ac_observed_ac','list_spike_ac_observed_FFT',...
'list_spike_shuffled_ac','list_spike_ac_shuffled_FFT_99','list_spike_ac_freq_of_interest',...
'seq_spike_times','seq_spike_ac_observed_ac','seq_spike_ac_observed_FFT',...
'seq_spike_shuffled_ac','seq_spike_ac_shuffled_FFT_99','seq_spike_ac_freq_of_interest',...
'list_eye_times','list_eye_ac_observed_ac','list_eye_ac_observed_FFT',...
'list_eye_shuffled_ac','list_eye_ac_shuffled_FFT_99','list_eye_ac_freq_of_interest',...
'seq_eye_times','seq_eye_ac_observed_ac','seq_eye_ac_observed_FFT',...
'seq_eye_shuffled_ac','seq_eye_ac_shuffled_FFT_99','seq_eye_ac_freq_of_interest',...
'list_spike_eye_xc_observed_xc','list_spike_eye_xc_observed_max_xc',...
'list_spike_eye_xc_observed_FFT','list_spike_eye_xc_shuffled_xc',...
'list_spike_eye_xc_shuffled_max_xc','list_spike_eye_xc_shuffled_FFT_99',...
'seq_spike_eye_xc_observed_xc','seq_spike_eye_xc_observed_max_xc',...
'seq_spike_eye_xc_observed_FFT','seq_spike_eye_xc_shuffled_xc',...
'seq_spike_eye_xc_shuffled_max_xc','seq_spike_eye_xc_shuffled_FFT_99');
%%
end

function FOI = find_frequencies_of_interest(PSD,whitened_PSD,f2,f100,f450,Hz4)
FOI= NaN(1,5);%frequencies of interest
[~,maxi] = max(PSD(1:f100));%orignal max below 100 Hz
FOI(1) = f2(maxi);
[~,maxi] = max(PSD(f100:f450));%max power from 100-450 Hz
FOI(5) = f2(maxi)+f2(f100);
whitened_PSD = PSD./whitened_PSD;
[pks,locs] = findpeaks(whitened_PSD(1:f100),'MinPeakDistance',Hz4);
pks(locs < Hz4/4) = [];%remove peaks less than 1 Hz
locs(locs < Hz4/4) = [];%remove peaks less than 1 Hz
[~,is] = sort(pks);%sort peaks by magnitude
locs = locs(is);%sort locs by peak
locs = locs(end-2:end);%take 3 largest peaks
FOI(2:4)= f2(locs);

%whitening shifts observed peaks
real_peak = NaN(1,5);
real_peak(1) = FOI(1);
real_peak(5) = FOI(5);
for freq = 2:4
    fmin = find(f2 <= FOI(freq)-1/(.05*FOI(freq)));%scale window to 1/frequency
    if isempty(fmin)
        fmin = 1;
    end
    fmin = fmin(end);
    fmax = find(f2 >= FOI(freq)+1/(.05*FOI(freq)));%scale window to 1/frequency
    fmax = fmax(1);
    [~,maxind] = max(PSD(fmin:fmax));
    real_peak(freq) = fmin+maxind-1;
end
real_peak(2:4) = f2(real_peak(2:4));
FOI = sort(real_peak);
end

function FOI = find_frequencies_of_interest_eye(PSD,f2,f100,Hz4)
[pks,locs] = findpeaks(PSD(1:f100),'MinPeakDistance',Hz4/2);%at least 2 Hz apart
pks(locs < Hz4/4) = [];%remove peaks less than 1 Hz
locs(locs < Hz4/4) = [];%remove peaks less than 1 Hz
[~,is] = sort(pks);%sort peaks by magnitude
locs = locs(is);%sort locs by peak
locs = locs(end-2:end);%take 3 largest peaks, probably really only 2 but just in case
FOI = f2(locs);
FOI = sort(FOI);
end