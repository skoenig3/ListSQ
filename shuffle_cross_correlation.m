function [shuffled_xc,shuffled_xc_FFT,shuffled_max_xc] = shuffle_cross_correlation(...
    event_times1,event_times2,twin,NFFT,freq_filt,numshuffs)
%written by Seth Konig 3/21/17
%calculates shuffled crosscorelation and FFT/PSD by rotating the event
%times of 1 of the matrices relative to the other (doesn't rotate).

shuffled_xc = NaN(numshuffs,2*twin+1);
shuffled_xc_FFT = NaN(numshuffs,NFFT/2+1);
shuffled_max_xc = NaN(1,numshuffs);
parfor shuff = 1:numshuffs
    
    %---Randomly Rotate Event Time1 Relative to Event Times 2---%
    shuffled_event_times = circshift_row(event_times1);

    %---Calculate Shuffled Cross Correlation---%
    [shuffled_xc(shuff,:),shuffled_xc_FFT(shuff,:),shuffled_max_xc(shuff)] = ...
        calculate_trial_averaged_xc...
        (shuffled_event_times,event_times2,twin,NFFT,freq_filt);
end
shuffled_xc = mean(shuffled_xc);%don't need to waste memory for visualization only