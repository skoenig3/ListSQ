function [observed_info_rate,shuffled_info_rate]= estimated_mutual_information(trial_data,numshuffs,info_type,smval,Fs)
% written by Seth Konig August, 2014
% estimates the "mutual information" rate (aka Skaggs Infomrmation) between
% spikes and variable of interest. Further the function calculates the expected
% mutual information for spikie trains rotated in time within a given trial.
% Updated 1/12/16 by adding case 'spatial_noshuff' so can compute mutual
% information without shuffling
% Updated 4/3/16 to shuffle spatial information by rotating whole spike
% train across all trials instead of within trials. Eye position can be
% heaviltiy biased within trials but less so across trials this should give
% a better representation of the data.
% Inputs:
%   1) trial_data: arranged in trials by row and time from designated event
%   by colomuns. 1's for when spikes occured within a trial and 0s for when
%   spikes didn't occur. if trial lengths are variable then NaNs should be
%   used to fill in time columns in which time points did not occur
%   2) numshuffs: number of shuffles to compute expected mutual information given
%   the statistics of the spike trains
%   3) info_type: the type of data i.e. 'temporal', 'spatial', or 'directional'
%   4) smval: smoothing parameter for temporal or spatial filter.
%       a) For temporal,s mval is gaussian 1/2 width.
%       b) For spatial smval.std is the standard deviation of the 2D
%       gaussian filter and smval.size is the size of the 2D filter.
%   5) Fs: sampling rate i.e. 1000 Hz
%
% Outputs:
%   1) observed_info_rate: mutual information rate for observed data
%   2) shuffled_info_rate: expected mutual information given the statistics
%   of the spike trains computed using bootstrapping. It is a matrix
%   containined numshuff worht of bootstrapped infomration rates

if ~isempty(trial_data)
    switch info_type
        case 'temporal'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate the observed mutual information rate---%%%
            total_time = sum(~isnan(trial_data)); %total number of time_points
            ntrials = floor(size(trial_data,1)/2); %number of trials/2
            
            %probability of time being observed at a given point across all
            %trials, not necessarily the same
            p_x = total_time/sum(total_time);
            [lambda_x,~]= nandens(trial_data,smval,'gauss',Fs,'nanflt');%nandens made by nathan killian
            lambda = nansum(nansum(lambda_x.*p_x));
            [ observed_info_rate.skaggs] = p_log_p(lambda,lambda_x,p_x);
            
            if observed_info_rate.skaggs < 0
                error('Information should be greater than 0. WTF?')
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate the observed temporal stability---%%%
            [lambda_x1,~]= nandens(trial_data(1:ntrials,:),smval,'gauss',Fs,'nanflt');%firing rate curve for first half
            [lambda_x2,~]= nandens(trial_data(ntrials+1:end,:),smval,'gauss',Fs,'nanflt');%firing rate curve for second half
            observed_info_rate.temporalstability = corr2(lambda_x1,lambda_x2);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate the expected mutual information rate using bootstrapping---%%%
            
            if observed_info_rate.skaggs ~= 0 %i.e. not zero firing rate or 100% uniform
                skaggs = NaN(1,numshuffs);
                temporalstability = NaN(1,numshuffs);
                parfor shuffled = 1:numshuffs;
                    shuffled_firing = circshift_row(trial_data);
                    
                    [lambda_x,~]= nandens(shuffled_firing,smval,'gauss',Fs,'nanflt');
                    lambda = nansum(nansum(lambda_x.*p_x));
                    skaggs(shuffled) = p_log_p(lambda,lambda_x,p_x);
                    
                    if skaggs(shuffled) <= 0
                        disp('Error: Information should be greater than 0. WTF?')
                    end
                    
                    [lambda_x1,~]= nandens(shuffled_firing(1:ntrials,:),smval,'gauss',Fs,'nanflt');%firing rate curve for first half
                    [lambda_x2,~]= nandens(shuffled_firing(ntrials+1:end,:),smval,'gauss',Fs,'nanflt');%firing rate curve for second half
                    temporalstability(shuffled) = corr2(lambda_x1,lambda_x2);
                end
                shuffled_info_rate.skaggs = skaggs;
                shuffled_info_rate.temporalstability =temporalstability;
            else
                shuffled_info_rate.skaggs = [];
                shuffled_info_rate.temporalstability = [];
            end
            
        case 'spatial'
            
            %get data and filtering parameters
            eyepos = trial_data{1}; %eye position by trial odd rows x, even rows y
            spike_times = trial_data{2};%spike times aligned to image eye data by trial
            imageX = trial_data{3}(1);%horizontal size of the image
            imageY = trial_data{3}(2);%horizontal size of the image
            min_bin_dur = trial_data{4}; %minimum amount of time in binned location to use
            
            binsize = smval(1); %number of pixels per spatial bin
            filter_width = smval(2); %std of 2D gaussian smoothing filter
            filter_size = filter_width*10;
            H = fspecial('gaussian',filter_size,filter_width);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate the observed mutual information rate---%%%
            
            [filtered_time] = filter_time(eyepos,imageX,imageY,Fs,binsize,H);
            filtered_time(filtered_time < min_bin_dur) = NaN;
            [filtered_space] = filter_space(eyepos,spike_times,imageX,imageY,binsize,H);
            
            lambda_x = filtered_space./filtered_time; %observed firing rate over space
            p_x = filtered_time/nansum(nansum(filtered_time));
            lambda = nansum(nansum(lambda_x.*p_x));
            [observed_info_rate.skaggs] = nansum(p_log_p(lambda,lambda_x,p_x));
            
            if observed_info_rate.skaggs < 0
                error('Information should be greater than 0. WTF?')
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate the observed spatial stability---%%%
            %look at correlation between firing locations for first and 2nd half of data
            
            ntrials = floor(size(spike_times,1)/2);
            
            %---get rate map for first half vs second half---%
            %first half
            filtered_time1 = filter_time(eyepos(1:ntrials*2,:),imageX,imageY,Fs,binsize,H);
            filtered_time1(filtered_time1 < min_bin_dur) = NaN; 
            filtered_space1 = filter_space(eyepos(1:ntrials*2,:),...
                spike_times(1:ntrials,:),imageX,imageY,binsize,H);
            ratemap1 = filtered_space1./filtered_time1;
            ratemap1(isnan(ratemap1)) = 0;
            
            %second half
            filtered_time2 = filter_time(eyepos(ntrials*2+1:ntrials*2*2,:),...
                imageX,imageY,Fs,binsize,H);
            filtered_time2(filtered_time2 < min_bin_dur) = NaN;
            filtered_space2 = filter_space(eyepos(ntrials*2+1:ntrials*2*2,:),...
                spike_times(ntrials+1:ntrials*2,:),imageX,imageY,binsize,H);
            ratemap2 = filtered_space2./filtered_time2;
            ratemap2(isnan(ratemap2)) = 0;
            
            observed_info_rate.spatialstability = corr2(ratemap1,ratemap2);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate the expected information rates and spatial stability using bootstrapping---%%%
            
            if observed_info_rate.skaggs ~= 0 %i.e. not zero firing rate or 100% uniform
                skaggs = NaN(1,numshuffs);
                spatialstability = NaN(1,numshuffs);
                parfor shuffled = 1:numshuffs;
                    
                    %previous method of shuffling data within trial
                    %shuffled_firing = circshift_row(spike_times);
                    shuffled_firing = circshift_acrosstrials(spike_times);
                    
                    
                    %---Calcualte Shuffled Mutual Information---%
                    %observed firing rate over space is the only thing changing
                    [filtered_space] = filter_space(eyepos,shuffled_firing,imageX,imageY,binsize,H);
                    
                    lambda_x = filtered_space./filtered_time;
                    lambda = nansum(nansum(lambda_x.*p_x));
                    skaggs(shuffled) = nansum(p_log_p(lambda,lambda_x,p_x));
                    if  skaggs(shuffled) < 0
                        error('Information should be greater than 0. WTF?')
                    end
                    
                    %---Calcualte Shuffled Spatial Stability---%
                    [filtered_space1] = filter_space(eyepos(1:ntrials*2,:),...
                        shuffled_firing(1:ntrials,:),imageX,imageY,binsize,H);
                    ratemap1 = filtered_space1./filtered_time1;
                    ratemap1(isnan(ratemap1)) = 0;
                    
                    [filtered_space2] = filter_space(eyepos(ntrials*2+1:ntrials*2*2,:),...
                        shuffled_firing(ntrials+1:ntrials*2,:),imageX,imageY,binsize,H);
                    ratemap2 = filtered_space2./filtered_time2;
                    ratemap2(isnan(ratemap2)) = 0;
                    
                    spatialstability(shuffled) = corr2(ratemap1,ratemap2);
                end
                shuffled_info_rate.skaggs = skaggs;
                shuffled_info_rate.spatialstability = spatialstability;
            else
                shuffled_info_rate.skaggs = [];
                shuffled_info_rate.spatialstability = [];
            end
            
            
        case 'spatial_cvtnew' %special case since rotating spike train is not
            %an effective measure due correlation in dot position over time
            
            %get data and filtering parameters
            eyepos = trial_data{1}; %eye position by trial odd rows x, even rows y
            spike_times = trial_data{2};%spike times aligned to image eye data by trial
            imageX = trial_data{3}(1);%horizontal size of the image
            imageY = trial_data{3}(2);%horizontal size of the image
            
            binsize = smval(1); %number of pixels per spatial bin
            filter_width = smval(2); %std of 2D gaussian smoothing filter
            filter_size = filter_width*10;
            H = fspecial('gaussian',filter_size,filter_width);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate the observed mutual information rate---%%%
            
            filtered_time = filter_time(eyepos,imageX,imageY,Fs,binsize,H);
            filtered_space = filter_space(eyepos,spike_times,imageX,imageY,binsize,H);
            
            lambda_x = filtered_space./filtered_time; %observed firing rate over space
            p_x = filtered_time/nansum(nansum(filtered_time));
            lambda = nansum(nansum(lambda_x.*p_x));
            [observed_info_rate.skaggs] = nansum(p_log_p(lambda,lambda_x,p_x));
            
            if observed_info_rate.skaggs < 0
                error('Information should be greater than 0. WTF?')
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate the observed spatial stability---%%%
            %look at correlation between firing locations for first and 2nd half of data
            
            ntrials = floor(size(spike_times,1)/2);
            
            %---get rate map for first half vs second half---%
            %first half
            [filtered_time1] = filter_time(eyepos(1:ntrials*2,:),...
                imageX,imageY,Fs,binsize,H);
            [filtered_space1] = filter_space(eyepos(1:ntrials*2,:),...
                spike_times(1:ntrials,:),imageX,imageY,binsize,H);
            ratemap1 = filtered_space1./filtered_time1;
            ratemap1(isnan(ratemap1)) = 0;
            
            %second half
            [filtered_time2] = filter_time(eyepos(ntrials*2+1:ntrials*2*2,:),...
                imageX,imageY,Fs,binsize,H);
            [filtered_space2] = filter_space(eyepos(ntrials*2+1:ntrials*2*2,:),...
                spike_times(ntrials+1:ntrials*2,:),imageX,imageY,binsize,H);
            ratemap2 = filtered_space2./filtered_time2;
            ratemap2(isnan(ratemap2)) = 0;
            
            observed_info_rate.spatialstability = corr2(ratemap1,ratemap2);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate the expected information rates and spatial stability using bootstrapping---%%%
            shuffled_info_rate.skaggs = NaN(1,numshuffs);
            shuffled_info_rate.spatialstability = NaN(1,numshuffs);
            if observed_info_rate.skaggs ~= 0 %i.e. not zero firing rate or 100% uniform
                for shuffled = 1:numshuffs;
                    
                    %previous method of shuffling data within trial
                    %shuffled_firing = circshift_row(spike_times);
                    shuffled_firing = circshift_cvtnew(spike_times);
                    
                    
                    %---Calcualte Shuffled Mutual Information---%
                    %observed firing rate over space is the only thing changing
                    [filtered_space] = filter_space(eyepos,shuffled_firing,imageX,imageY,binsize,H);
                    
                    lambda_x = filtered_space./filtered_time;
                    lambda = nansum(nansum(lambda_x.*p_x));
                    shuffled_info_rate.skaggs(shuffled) = nansum(p_log_p(lambda,lambda_x,p_x));
                    if shuffled_info_rate.skaggs(shuffled) < 0
                        error('Information should be greater than 0. WTF?')
                    end
                    
                    %---Calcualte Shuffled Spatial Stability---%
                    [filtered_space1] = filter_space(eyepos(1:ntrials*2,:),...
                        shuffled_firing(1:ntrials,:),imageX,imageY,binsize,H);
                    ratemap1 = filtered_space1./filtered_time1;
                    ratemap1(isnan(ratemap1)) = 0;
                    
                    [filtered_space2] = filter_space(eyepos(ntrials*2+1:ntrials*2*2,:),...
                        shuffled_firing(ntrials+1:ntrials*2,:),imageX,imageY,binsize,H);
                    ratemap2 = filtered_space2./filtered_time2;
                    ratemap2(isnan(ratemap2)) = 0;
                    
                    shuffled_info_rate.spatialstability(shuffled) = corr2(ratemap1,ratemap2);
                end
            end
            
        case 'spatial_noshuff'
            shuffled_info_rate = [];
            
            %get data and filtering parameters
            eyepos = trial_data{1}; %eye position by trial odd rows x, even rows y
            spike_times = trial_data{2};%spike times aligned to image eye data by trial
            imageX = trial_data{3}(1);%horizontal size of the image
            imageY = trial_data{3}(2);%horizontal size of the image
            min_bin_dur = trial_data{4}; %minimum amount of time in binned location to use
            
            binsize = smval(1); %number of pixels per spatial bin
            filter_width = smval(2); %std of 2D gaussian smoothing filter
            filter_size = filter_width*10;
            H = fspecial('gaussian',filter_size,filter_width);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate the observed mutual information rate---%%%
            filtered_time = filter_time(eyepos,imageX,imageY,Fs,binsize,H);
            filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
            filtered_space = filter_space(eyepos,spike_times,imageX,imageY,binsize,H);
            lambda_x = filtered_space./filtered_time; %observed firing rate over space
            p_x = filtered_time/nansum(nansum(filtered_time));
            lambda = nansum(nansum(lambda_x.*p_x));
            [observed_info_rate.skaggs] = sum(p_log_p(lambda,lambda_x,p_x));
            
            if observed_info_rate.skaggs < 0
                error('Information should be greater than 0. WTF?')
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate the observed spatial stability---%%%
            %look at correlation between firing locations for first and 2nd half of data
            
            ntrials = floor(size(spike_times,1)/2);
            
            %---get rate map for first half vs second half---%
            %first half
            filtered_time1 = filter_time(eyepos(1:ntrials*2,:),imageX,imageY,Fs,binsize,H);
            filtered_time1(filtered_time1 < min_bin_dur) = NaN; 
            filtered_space1 = filter_space(eyepos(1:ntrials*2,:),...
                spike_times(1:ntrials,:),imageX,imageY,binsize,H);
            ratemap1 = filtered_space1./filtered_time1;
            ratemap1(isnan(ratemap1)) = 0;
            
            %second half
            filtered_time2 = filter_time(eyepos(ntrials*2+1:ntrials*2*2,:),...
                imageX,imageY,Fs,binsize,H);
            filtered_time2(filtered_time2 < min_bin_dur) = NaN;
            filtered_space2 = filter_space(eyepos(ntrials*2+1:ntrials*2*2,:),...
                spike_times(ntrials+1:ntrials*2,:),imageX,imageY,binsize,H);
            ratemap2 = filtered_space2./filtered_time2;
            ratemap2(isnan(ratemap2)) = 0;
            
            observed_info_rate.spatialstability = corr2(ratemap1,ratemap2);
            
            
            
    end
else
    observed_info_rate = NaN;
    shuffled_info_rate = [];
end
end

function [info] = p_log_p(lambda,lambda_x,p_x)
%mutual info equation
plogp = lambda_x.*log2(lambda_x/lambda);
info = nansum(plogp.*p_x);
end