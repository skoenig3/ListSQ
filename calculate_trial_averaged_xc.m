function [trial_averaged_xc,trial_averaged_xc_FFT,max_xc] = calculate_trial_averaged_xc...
(event_times1,event_times2,twin,NFFT,freq_filt)
%written by Seth Konig 3/21/17
%calculates trial averaged cross correlation by first calculating cross
%correlation for exch trial, then averages xcross trials, and finally
%normalizes

if all(size(event_times1)~= size(event_times2))
    error('Event Matrices must be same size')
end
%---Calculate Autocorrelation---%
trial_averaged_xc= NaN(size(event_times1,1),twin*2+1);
for tr = 1:size(event_times1,1);
    temp1 = event_times1(tr,:);
    temp1(isnan(temp1)) = []; %remove trailing nans if they exist

        temp2 = event_times2(tr,:);
    temp2(isnan(temp2)) = []; %remove trailing nans if they exist

    if length(temp1) ~= length(temp2)
        error('Indexing Error?... vectors must be same size')
    end
    trial_averaged_xc(tr,:) = xcorr(temp1,temp2,twin);
end
trial_averaged_xc = mean(trial_averaged_xc);
binned_average = bin1(trial_averaged_xc,5);%5ms bins, faster than smoothing
max_xc = max(abs(binned_average));%largest correlation value
trial_averaged_xc = trial_averaged_xc-mean(trial_averaged_xc);%finally remove DC power
trial_averaged_xc = trial_averaged_xc/max(abs(trial_averaged_xc));

%---Calculate FFT (PSD) of Autocorrelation---%
Y = fft(trial_averaged_xc-mean(trial_averaged_xc),NFFT);%FFT
trial_averaged_xc_FFT= abs(Y(1:NFFT/2+1)).^2;%power
trial_averaged_xc_FFT = filtfilt(freq_filt,1,trial_averaged_xc_FFT);%smooth power
end