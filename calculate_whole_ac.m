function [xc,ac_FFTs] = calculate_whole_ac(event_times,twin,NFFT,freq_filt,freq_filt2)
%written by Seth Konig 3/21/17
%calculates cross correlation for whole sesssion. Should be faster than 
%calculate_trial_averaged_ac since not iterating through trials

%---Calculate Autocorrelation---%
xc = xcorr(event_times,event_times,twin);
xc(twin+1)=0; %set center value to zero for now
xc = xc/max(xc);%normalize to 1
xc(twin+1) = 1;%set zeroth lag to 1
xc = xc-mean(xc);%subtract mean

%---Calculate FFT (PSD) of Autocorrelation---%
ac_FFT = fft(xc,NFFT);%FFT
ac_FFT = abs(ac_FFT(1:NFFT/2+1)).^2;%power
ac_FFT1 = filtfilt(freq_filt,1,ac_FFT);%smooth power
ac_FFT2 = filtfilt(freq_filt2,1,ac_FFT);%smooth power
ac_FFTs = [ac_FFT1; ac_FFT2];