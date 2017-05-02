function [trial_averaged_ac,trial_averaged_ac_FFT] = calculate_trial_averaged_ac(event_times,twin,NFFT,freq_filt)
%written by Seth Konig 3/21/17
%calculates trial averaged cross correlation by first calculating cross
%correlation for each trial, then averages across trials, and finally
%normalizes


%---Calculate Autocorrelation---%
trial_averaged_ac= NaN(size(event_times,1),twin*2+1);
for tr = 1:size(event_times,1);
    temp = event_times(tr,:);
    temp(isnan(temp)) = []; %remove trailing nans if they exist
    trial_averaged_ac(tr,:) = xcorr(temp,temp,twin);
end


trial_averaged_ac(:,twin+1) = 0; %set center value to zero for now
trial_averaged_ac = mean(trial_averaged_ac); %average across trials
trial_averaged_ac = trial_averaged_ac/max(trial_averaged_ac); %normalize to 1
trial_averaged_ac(twin+1) = 1; %set zeroth lag to 1
%trial_averaged_ac = trial_averaged_ac-mean(trial_averaged_ac);%finally remove DC power

%---Calculate FFT (PSD) of Autocorrelation---%
Y = fft(trial_averaged_ac-mean(trial_averaged_ac),NFFT);%FFT
trial_averaged_ac_FFT= abs(Y(1:NFFT/2+1)).^2;%power
trial_averaged_ac_FFT = filtfilt(freq_filt,1,trial_averaged_ac_FFT);%smooth power
end