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
            
            %probability of time being observed at a given point across all
            %trials, not necessarily the same
            p_x = total_time/sum(total_time);
            [lambda_x,~]= nandens(trial_data,smval,'gauss',Fs,'nanflt');%nandens made by nathan killian
            lambda = nansum(nansum(lambda_x.*p_x));
            [observed_info_rate] = p_log_p(lambda,lambda_x,p_x);
            
            if observed_info_rate < 0
                error('Information should be greater than 0. WTF?')
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate the expected mutual information rate using bootstrapping---%%%
            shuffled_info_rate = NaN(1,numshuffs);
            if observed_info_rate ~= 0 %i.e. not zero firing rate or 100% uniform
                for shuffled = 1:numshuffs;
                    shuffled_firing = circshift_row(trial_data);
                    [lambda_x,~]= nandens(shuffled_firing,smval,'gauss',Fs,'nanflt');
                    lambda = nansum(nansum(lambda_x.*p_x));
                    shuffled_info_rate(shuffled) = p_log_p(lambda,lambda_x,p_x);
                     
                    if shuffled_info_rate(shuffled) <= 0
                        disp('Error: Information should be greater than 0. WTF?')
                    end
                end
            end
            
        case 'spatial'
                    
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
            
            [filtered_time] = get_smoothed_Time(eyepos,imageX,imageY,Fs,binsize,H);
            filtered_time(filtered_time < 0.025) = NaN;
            [filtered_space] = get_smoothed_Space(eyepos,spike_times,imageX,imageY,binsize,H);

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
            [filtered_time1] = get_smoothed_Time(eyepos(1:ntrials*2,:),...
                imageX,imageY,Fs,binsize,H);
            [filtered_space1] = get_smoothed_Space(eyepos(1:ntrials*2,:),...
                spike_times(1:ntrials,:),imageX,imageY,binsize,H);      
            ratemap1 = filtered_space1./filtered_time1;
            ratemap1(isnan(ratemap1)) = 0; 
            
            %second half
            [filtered_time2] = get_smoothed_Time(eyepos(ntrials*2+1:ntrials*2*2,:),...
                imageX,imageY,Fs,binsize,H);
            [filtered_space2] = get_smoothed_Space(eyepos(ntrials*2+1:ntrials*2*2,:),...
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
                    shuffled_firing = circshift_acrosstrials(spike_times);
                    
                    
                    %---Calcualte Shuffled Mutual Information---%
                    %observed firing rate over space is the only thing changing
                    [filtered_space] = get_smoothed_Space(eyepos,shuffled_firing,imageX,imageY,binsize,H);

                    lambda_x = filtered_space./filtered_time;
                    lambda = nansum(nansum(lambda_x.*p_x));
                    shuffled_info_rate.skaggs(shuffled) = nansum(p_log_p(lambda,lambda_x,p_x));
                    if shuffled_info_rate.skaggs(shuffled) < 0
                        error('Information should be greater than 0. WTF?')
                    end
                    
                    %---Calcualte Shuffled Spatial Stability---%
                    [filtered_space1] = get_smoothed_Space(eyepos(1:ntrials*2,:),...
                        shuffled_firing(1:ntrials,:),imageX,imageY,binsize,H);
                    ratemap1 = filtered_space1./filtered_time1;
                    ratemap1(isnan(ratemap1)) = 0;
                    
                    [filtered_space2] = get_smoothed_Space(eyepos(ntrials*2+1:ntrials*2*2,:),...
                        shuffled_firing(ntrials+1:ntrials*2,:),imageX,imageY,binsize,H);
                    ratemap2 = filtered_space2./filtered_time2;
                    ratemap2(isnan(ratemap2)) = 0;
                    
                    shuffled_info_rate.spatialstability(shuffled) = corr2(ratemap1,ratemap2);
                end
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
            
            [filtered_time] = get_smoothed_Time(eyepos,imageX,imageY,Fs,binsize,H);
            [filtered_space] = get_smoothed_Space(eyepos,spike_times,imageX,imageY,binsize,H);

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
            [filtered_time1] = get_smoothed_Time(eyepos(1:ntrials*2,:),...
                imageX,imageY,Fs,binsize,H);
            [filtered_space1] = get_smoothed_Space(eyepos(1:ntrials*2,:),...
                spike_times(1:ntrials,:),imageX,imageY,binsize,H);      
            ratemap1 = filtered_space1./filtered_time1;
            ratemap1(isnan(ratemap1)) = 0; 
            
            %second half
            [filtered_time2] = get_smoothed_Time(eyepos(ntrials*2+1:ntrials*2*2,:),...
                imageX,imageY,Fs,binsize,H);
            [filtered_space2] = get_smoothed_Space(eyepos(ntrials*2+1:ntrials*2*2,:),...
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
                    [filtered_space] = get_smoothed_Space(eyepos,shuffled_firing,imageX,imageY,binsize,H);

                    lambda_x = filtered_space./filtered_time;
                    lambda = nansum(nansum(lambda_x.*p_x));
                    shuffled_info_rate.skaggs(shuffled) = nansum(p_log_p(lambda,lambda_x,p_x));
                    if shuffled_info_rate.skaggs(shuffled) < 0
                        error('Information should be greater than 0. WTF?')
                    end
                    
                    %---Calcualte Shuffled Spatial Stability---%
                    [filtered_space1] = get_smoothed_Space(eyepos(1:ntrials*2,:),...
                        shuffled_firing(1:ntrials,:),imageX,imageY,binsize,H);
                    ratemap1 = filtered_space1./filtered_time1;
                    ratemap1(isnan(ratemap1)) = 0;
                    
                    [filtered_space2] = get_smoothed_Space(eyepos(ntrials*2+1:ntrials*2*2,:),...
                        shuffled_firing(ntrials+1:ntrials*2,:),imageX,imageY,binsize,H);
                    ratemap2 = filtered_space2./filtered_time2;
                    ratemap2(isnan(ratemap2)) = 0;
                    
                    shuffled_info_rate.spatialstability(shuffled) = corr2(ratemap1,ratemap2);
                end
            end
            
        case 'spatial_noshuff'
            %same as above without the shuffling
            eyepos = trial_data{1}; %eye position by trial odd rows x, even rows y
            spike_times = trial_data{2};%spike times aligned to image eye data by trial
            imageX = trial_data{3}(1);%horizontal size of the image
            imageY = trial_data{3}(2);%horizontal size of the image
            
            binsize = smval(1); %number of pixels per spatial bing
            filter_width = smval(2); %std of 2D gaussian smoothing filter
            H = fspecial('gaussian',[filter_width*10+5 filter_width*10+5],filter_width);
            
            %calculate the total time spent at any locaitons in binned pixels
            spatial_time = time_per_pixel(eyepos,imageX,imageY,Fs);
            filtered_time = bin2(spatial_time,binsize,binsize);
            filtered_time = imfilter(filtered_time,H);
            filtered_time(filtered_time < 0.001) = NaN;
            filtered_time = filtered_time(end:-1:1,:); %flip so rightside up
            
            %caluclate total spikes over space
            [firing_location] = pixel_spike_location(eyepos,spike_times,imageX,imageY);
            filtered_space = bin2(firing_location,binsize,binsize);
            filtered_space = imfilter(filtered_space,H);
            filtered_space(filtered_space == 0) = NaN;
            filtered_space = filtered_space(end:-1:1,:); %flip so rightside up
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate the observed mutual information rate---%%%
            lambda_x = filtered_space./filtered_time; %observed firing rate over space
            p_x = filtered_time/nansum(nansum(filtered_time));
            lambda = nansum(nansum(lambda_x.*p_x));
            [observed_info_rate] = sum(p_log_p(lambda,lambda_x,p_x));
            
            if observed_info_rate < 0
                error('Information should be greater than 0. WTF?')
            end
            
            
            shuffled_info_rate = []; 
            
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

function [filtered_time] = get_smoothed_Time(eyepos,imageX,imageY,Fs,binsize,H)

%calculate the total time spent at any locaitons in binned pixels
spatial_time = time_per_pixel(eyepos,imageX,imageY,Fs);
filtered_time = bin2(spatial_time,binsize,binsize);
filtered_time = imfilter(filtered_time,H);
filtered_time(filtered_time < 0.001) = NaN;
filtered_time = filtered_time(end:-1:1,:); %flip so rightside up
end

function [filtered_space] = get_smoothed_Space(eyepos,spike_times,imageX,imageY,binsize,H)
%caluclate total spikes over space
[firing_location] = pixel_spike_location(eyepos,spike_times,imageX,imageY);
filtered_space = bin2(firing_location,binsize,binsize);
filtered_space = imfilter(filtered_space,H);
filtered_space(filtered_space == 0) = NaN;
filtered_space = filtered_space(end:-1:1,:); %flip so rightside up
end