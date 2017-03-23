function [shuffled_auto_corr,shuffled_fft] = jitter_window_shuffling(PSD_numshuffs,trial_spike_times,NFFT,twin,trial_length,flt,FOI,smFreqVal)

shuffled_auto_corr = NaN(PSD_numshuffs,twin*2+1);
shuffled_fft = NaN(PSD_numshuffs,NFFT/2+1);

bin_size = round(1/(FOI-smFreqVal/2)/2*1000);%~window size of 1/2 the frequency
if bin_size < 10 %save on computational time but also no point in going lower
    bin_size = 10;
end

if all(trial_length == trial_length(1));
    num_window = ceil(size(trial_spike_times,2)/bin_size);
    
    spike_windows = cell(1,num_window);
    for win = 1:num_window
        if win == num_window
            spike_windows{win} = trial_spike_times(:,1+bin_size*(win-1):end);
        else
            spike_windows{win} = trial_spike_times(:,1+bin_size*(win-1):bin_size*win);
        end
    end
    if sum(sum(trial_spike_times)) ~= sum(cellfun(@sum,(cellfun(@sum,spike_windows,'UniformOutput',false))))
        error('Where did the spikes go')
    end
    
    parfor shuff = 1:PSD_numshuffs
        shuff_window = cell(1,num_window);
        for  win = 1:num_window
            shuff_window{win} =  shake(spike_windows{win},2);
        end
        shuffled_times = cell2mat(shuff_window);

        shuff_xc = NaN(size(shuffled_times,1),twin*2+1);
        for f = 1:size(shuffled_times,1);
            this_trial = shuffled_times(f,:);
            shuff_xc(f,:) = xcorr(this_trial,this_trial,twin);
        end
        
        %Calculated and Normalize Shuffled Average Cross Correlation
        shuff_xc(:,twin+1) = 0;
        s = mean(shuff_xc);
        s = s/max(s); %normalize to 1
        s(twin+1) = 1; %set zeroth lag to 1
        shuffled_auto_corr(shuff,:) = s;
        
        %Calculate Shuffled Power Spectrum
        Y2 = fft(s-mean(s),NFFT);
        Y2 = abs(Y2(1:NFFT/2+1)).^2;
        Y2 = filtfilt(flt,1,Y2); %smooth to compare to observed values
        shuffled_fft(shuff,:)= Y2;        
    end
    shuffled_auto_corr = mean(shuffled_auto_corr);%don't need to waste memory for visualization only
end