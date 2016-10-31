function [filtsignal,phi,amp] = filter_Hilbert(signal,Fs,bpass,meth)
%filter_Hilbert - Filter data and perform Hilbert transform
% 
%Jon Rueckemann 2016

nyqst=Fs/2;
bpass=bpass./nyqst; %normalizes fq by Nyquist

switch lower(meth)
    case {'butter','iir'}
        [x,y]=butter(3,bpass); %Third order Butterworth filter
    case 'fir'
        filter_order=2./bpass(1); %nsamples in one cycle of lower fq bound
        x=fir1(round(3*filter_order),bpass); %3 cycles of lower bound
        y=1;
end

filtsignal=filtfilt(x,y,signal);
h_signal=hilbert(filtsignal);
phi=angle(h_signal);
amp=abs(h_signal);
end