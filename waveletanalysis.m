function [trialpower,trialphase,wfq,meanpower,meanphase,powervar,phasevar] = waveletanalysis(signal,varargin)
%WAVELETANALYSIS - Uses Morlet wavelets to determine amplitude and phase
%information for an LFP signal.
%
%signal - TxN array of voltages, T trials with N timepoints per trial
%Fs - integer. Sampling rate
%
%optional (input as dyads):
%wfq - Mx1 array. Frequencies of wavelets
%           (min_fq, max_fq, n_fq, and logfqaxis are unused)
%min_fq - Minimum frequency of wavelets
%max_fq - Maximum frequency of wavelets
%n_fq - Number of wavelet frequencies
%logfqaxis - Binary. Logarithmically space sampled frequencies between the
%           minimum and maximum.  Otherwise, linear scale
%n_cycles - scalar.  Number of cycles defining the wavelet Gaussian
%wtime - 1xT array.  Duration of wavelet
%
%Jon Rueckemann 2016
%(adapted from Michael Cohen; Analyzing Neural Time Series Data (2014))

%Default values
Fs=1000; %1000Hz
wavdur=1;
min_fq=4;
max_fq=80;
n_fq=30;
n_cycles=[];
logfqaxis=true;

%Extract user input
extractvarargin(varargin);

wtime = -wavdur:1/Fs:wavdur; %Define wavelet duration

%Define wavelet frequencies
if ~exist('wfq','var')
    if min_fq==max_fq
        wfq=min_fq;
    elseif logfqaxis
        wfq=logspace(log10(min_fq),log10(max_fq),n_fq);
    else
        wfq=linspace(min_fq,max_fq,n_fq); %#ok<UNRCH>
    end
end

%Duration of wavelet Gaussian in seconds via the number of cycles per Hz
if isempty(n_cycles)
    %Adaptively scale number of cycles per wavelet logarithmically,
    %following the logarithmic trendline from 3 cycles for 2Hz to 10 cycles
    %for 80Hz
    s=(1.1963/pi)*wfq.^-0.6736;
else
    s=n_cycles./(2*pi*wfq);
end

%Define convolution parameters
n_wavelet=numel(wtime);
n_data=numel(signal);
n_convolution=n_wavelet+n_data-1;
n_conv_pow2=pow2(nextpow2(n_convolution));
half_of_wavelet_size=(n_wavelet-1)/2;

%Get FFT of data
[n_trials,n_pts]=size(signal);
sigfft=fft(reshape(signal',1,[]),n_conv_pow2);

%Initialize output
convmat=zeros(n_fq,n_pts,n_trials);

%Iterate through frequencies and compute convolution
for fi=1:n_fq
    %FFT of wavelet
    %(first term is the sine wave in Hz and the second term is the gaussian 
    %window defined by the number of cycles above)
    wavelet=fft(exp(2i*pi*wfq(fi).*wtime).*exp(-wtime.^2./(2*(s(fi)^2))),n_conv_pow2);
    wavelet=wavelet./max(wavelet); %scale transform
    
    %Convolve wavelet with signal
    sigconv=ifft(wavelet.*sigfft);
    sigconv=sigconv(1:n_convolution);
    sigconv=sigconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    %Restructure data into trials
    convmat(fi,:,:)=reshape(sigconv,1,n_pts,n_trials);
end

%Calculate power and phase information for each trial
trialpower=convmat.*conj(convmat);
trialphase=angle(convmat);

%Calculate the trial mean power 
if nargout>3
    meanpower=mean(trialpower,3);
end
if nargout>5
    powervar=var(trialpower,0,3);
end

%Calculate the trial phase angle and circular variance
if nargout>4
    tempsum=sum(exp(1i*trialphase),3);
    meanphase=angle(tempsum);
end
if nargout>6
    phasevar=1-(abs(tempsum)./n_trials); %circular var = 1-MRL
end
end