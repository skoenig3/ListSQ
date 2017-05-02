function [shuffled_auto_corr,smoothed_shuffled_fft] = whole_jitter_window_shuffling(PSD_numshuffs,spike_times,NFFT,twin,freq_filt,FOI,smFreqVal)

shuffled_auto_corr = NaN(PSD_numshuffs,twin*2+1);

bin_size = round(1/(FOI-smFreqVal/2)/2*1000);%~window size of 1/2 the time period of the frequency below this 1
num_wind = floor(length(spike_times)/bin_size);%number of windows to shuffl over

%remove small number of extras at end, would take way too much
%computational time to try and take the last few ms of data
spike_times = spike_times(1:num_wind*bin_size);
spike_times = gpuArray(spike_times);
spike_times= reshape(spike_times,bin_size,num_wind);

%---Calculate Shuffled AutoCorrelations---%
parfor shuff = 1:PSD_numshuffs
    xc = calculate_shuffled_xc(spike_times,twin);
    shuffled_auto_corr(shuff,:) = gather(xc);
end

%---Calculate FFT (PSD) of Autocorrelation---%
shuffled_fft = fft(shuffled_auto_corr,NFFT,2);%FFT
shuffled_fft = abs(shuffled_fft(:,1:NFFT/2+1)).^2;%power


smoothed_shuffled_fft = zeros(size(shuffled_fft));
parfor shuff = 1:PSD_numshuffs
    smoothed_shuffled_fft(shuff,:) = filtfilt(freq_filt,1,shuffled_fft(shuff,:));%smooth power
end

end

% function xc = calculate_shuffled_xc(spike_times,twin)
%     spike_prime = shake(spike_times,1);
%     spike_prime = spike_prime(1:end);
%     xc = xcorr(spike_prime,spike_prime,twin);
%     xc(:,twin+1) = 0;%set center value to zero for now
%     xc = xc/max(xc);%normalize to 1
%     xc(twin+1) = 1;%set zeroth lag to 1
%     xc = xc-mean(xc); %subtract mean
% end