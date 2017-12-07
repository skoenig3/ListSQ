function [xc,xc_FFTs] = calculate_whole_xc(event_times1,event_times2,...
    twin,NFFT,freq_filt,freq_filt2)
%written by Seth Konig 9/15/17 modified from calculate_whole_ac
%calculates cross correlation for whole sesssion. Should be faster than 
%calculate_trial_averaged_ac since not iterating through trials

%---Calculate Autocorrelation---%
xc = xcorr(event_times1,event_times2,twin);
xc(twin+1)=0; %set center value to zero for now
xc = xc/max(xc);%normalize to 1
xc(twin+1) = 1;%set zeroth lag to 1
xc = xc-mean(xc);%subtract mean

%---Calculate FFT (PSD) of Autocorrelation---%
xc_FFT = fft(xc,NFFT);%FFT
xc_FFT = abs(xc_FFT(1:NFFT/2+1)).^2;%power
xc_FFT1 = filtfilt(freq_filt,1,xc_FFT);%smooth power
xc_FFT2 = filtfilt(freq_filt2,1,xc_FFT);%smooth power
xc_FFTs = [xc_FFT1; xc_FFT2];