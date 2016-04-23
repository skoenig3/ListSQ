function spatial_analysisV2(data_dir,figure_dir,session_data,task)
% written by Seth Konig August, 2014
% Function analyizes where the monkey is attending when spikes.
% updated SDK 1/11/16 to handlde new format and partial session data for
% vaild trials only.CVTNEW section on 1/19/16
%
% Inputs:
%   1) data_dir: direction where preprocessed data is located
%   2) figure_dir: location of where to put save figures
%   3) session_data: data containing relevant session data
%   4) task: what task was this data come from with i.e. 'cvtnew','ListSQ'
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-spatial_analysis_results.mat'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set important task parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
twin = 500;% here how much after image onset to ignore data
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
%temporal_shifts = [-150 -100 -50 0 50 100 150 200];%shifts in ms of spike trains to position data
%could at -200 but easier to  plot 8x2 than 9x2 on a square grid

%filtering and estimated mutual information parameters
binsize = 25; %pixels per bin spatial bin in either dimension ~1/2 dva
filter_width = 2; %std of 2D guassian filter ~ 2 dva
numshuffs = 500; %recommend this is between 100 & 1000
fr_threshold = 1; %peak rate must be greater than 1 Hz to process. Don't want to waste
%processing/shuffling time on "silent neurons"
filter_size = filter_width*10;
H = fspecial('gaussian',filter_size,filter_width);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials','hdr');
%grab unit data
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

num_trials = length(cfg.trl);
%NaNs are for start and end trials otherwise cut
valid_trials(1,isnan(valid_trials(1,:))) = 1;
valid_trials(2,isnan(valid_trials(2,:))) = num_trials;

switch task
    case {'cvtnew','CVTNEW'}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import & reformat data so that spikes are locked to dot position---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %also need to load meta data which is the dot position
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'meta');
        
        
        disp('Collecting Spike Locations')
        
        %set/get some general important info
        num_trials = length(cfg.trl);
        
        %NaNs are for start and end trials otherwise cut
        valid_trials(1,isnan(valid_trials(1,:))) = 1;
        valid_trials(2,isnan(valid_trials(2,:))) = num_trials;
        
        Fs = data(1).fsample; %should be 1000
        ITIstart_code = 15;
        dot_on_code = 25;
        dot_clrchng_code = 27;
        bar_code_response = 4; %monkey made a move
        reward_code = 3;
        imageX = 800; %horizontal size of the screen
        imageY = 600; %horizontal size of the screen
        
        %preallocate space and parallel structure of cfg
        dotpos = cell(1,num_units);
        spike_times = cell(1,num_units);
        for unit = 1:num_units
            spike_times{unit} = NaN(length(cfg.trl),3000);
            dotpos{unit} = NaN(2*length(cfg.trl),3000);
        end
        
        for t = 1:length(cfg.trl);
            if any(cfg.trl(t).allval == reward_code); %only take correct trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                pathon =  cfg.trl(t).alltim(cfg.trl(t).allval == dot_on_code)-trial_start; %meta data starts at 1st 666 which appears before dot actually turns on in event array
                dot_clrchng = cfg.trl(t).alltim(cfg.trl(t).allval == dot_clrchng_code)-trial_start;
                responded = cfg.trl(t).alltim(cfg.trl(t).allval == bar_code_response)-trial_start;
                
                xn = round(interp1(meta(t).sample,meta(t).x,meta(t).sample(1):meta(t).sample(end),'cubic'));
                yn = round(interp1(meta(t).sample,meta(t).y,meta(t).sample(1):meta(t).sample(end),'cubic'));
                
                for unit = 1:num_units
                    if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only put data in for valid trials
                        spikes = find(data(unit).values{t});
                        spikeind = spikes(spikes > pathon & spikes <= responded)-pathon;
                        spikeind(spikeind < 1) = []; %should only happen when spikes occur at the same time as tstart
                        temp = zeros(1,length(xn));
                        temp(spikeind) = 1;
                        spike_times{unit}(t,1:length(temp)) = temp;
                        
                        dotpos{unit}(2*t-1,1:length(xn)) = xn;
                        dotpos{unit}(2*t,  1:length(xn)) = yn;
                    end
                end
            end
        end
        
        %remove excess NaNs associated with error trials
        dotpos = laundry(dotpos);
        spike_times = laundry(spike_times);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Perform Spatial analysis and signifance testing---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Determining if neurons are spatially modulated')
        
        info_type = 'spatial_cvtnew';%type of mutual information analysis to perform
        
        %will perform spatial analysis on novel images, repeat images, and all images
        spatial_info.rate = NaN(1,num_units); %the observed information rate in bits/sec
        spatial_info.spatialstability = NaN(1,num_units); %the observed spatial correlation over time
        spatial_info.shuffled_info_rate = cell(1,num_units); %bootstrapped information rate in bits/sec expected by chance from spike train
        spatial_info.shuffled_spatialstability = cell(1,num_units); %bootstrapped spatial stability expected by chance from spike train
        trial_data{3} = [imageX imageY];
        for unit = 1:num_units
            trial_data{1} = dotpos{unit};
            trial_data{2} = spike_times{unit};
            [observed_info_rate,shuffled_info_rate]...
                = estimated_mutual_information(trial_data,numshuffs,info_type,[binsize,filter_width],Fs);
            
            spatial_info.rate(unit) = observed_info_rate.skaggs;
            spatial_info.spatialstability(unit) = observed_info_rate.spatialstability;
            spatial_info.shuffled_info_rate{unit} = shuffled_info_rate.skaggs;
            spatial_info.shuffled_spatialstability{unit} = shuffled_info_rate.spatialstability;
        end
        
        %estimate the significance threshold has the 95-percentile of the
        %bootstrapped data
        for unit = 1:num_units
            if ~isempty(spatial_info.shuffled_info_rate{unit})
                spatial_info.shuffled_rate_prctile(unit) = ...
                    100*sum(spatial_info.rate(unit) > ...
                    spatial_info.shuffled_info_rate{unit})/numshuffs;
                spatial_info.shuffled_spatialstability_prctile(unit) = ...
                    100*sum(spatial_info.spatialstability(unit) > ...
                    spatial_info.shuffled_spatialstability{unit})/numshuffs;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Plot and Save Figures of Results---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        task_type = 'cvtnew_spatial';
        unit_names.multiunit = multiunit;
        unit_names.name = cfg.channel;
        spatial_analysis_plotsV2(figure_dir,task_file,dotpos,spike_times,spatial_info,task_type,...
            unit_names,[binsize,filter_width],imageX,imageY,NaN,Fs);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Finally save all the data---%%%
        save([data_dir task_file(1:10) '-spatial_analysis_results.mat'],...
            'spike_times','dotpos','spatial_info','binsize','filter_width')
        disp(['Spatial Data Analyis for ' task_file(1:10) ' saved']);
        
    case 'ListSQ'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Reformat data so that spikes are locked to eye movements---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Collecting Spike Locations')
        
        %set/get some general important info
        [eyechans] = find_desired_channels(cfg,'eye');
        num_trials = length(cfg.trl);
        Fs = data(1).fsample; %should be 1000
        ITIstart_code = 15;
        img_on_code= 23;
        img_off_code = 24;
        imageX = 800; %horizontal size of the image
        imageY = 600; %horizontal size of the image
        
        %set the image duration
        if str2double(task_file(3:8)) < 140805
            imgdur = 7000;
        else
            imgdur = 5000;
        end
        
        %get important task specific information
        [itmlist,sequence_items] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);
        
        %preallocate space and parallel structure of cfg
        eyepos = cell(1,num_units);
        spike_times = cell(1,num_units);
        which_images = cell(1,num_units);
        nvr = cell(1,num_units);
        for unit = 1:num_units
            spike_times{unit} = NaN(length(cfg.trl),imgdur*1.5);
            eyepos{unit} = NaN(2*length(cfg.trl),imgdur*1.5);
            which_images{unit} = NaN(1,length(cfg.trl));
            nvr{unit} = NaN(1,length(cfg.trl));
        end
        
        for t = 1:num_trials
            if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                
                imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start+twin; %want to avoid visual response
                imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start;
                
                % if monkey isn't paying attention data isn't probably
                % worth much plus have to cut off somewhere
                if imgoff-imgon > 1.5*imgdur-1
                    imgoff = imgon+1.5*imgdur-1;
                end
                
                xn = round(data(eyechans(1)).values{t}(imgon:imgoff));
                yn = round(data(eyechans(2)).values{t}(imgon:imgoff));
                
                img_index = find(img_cnd == cfg.trl(t).cnd);
                for unit = 1:num_units
                    if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
                        eyepos{unit}(2*t-1,1:length(xn)) = xn;
                        eyepos{unit}(2*t,1:length(yn)) = yn;
                        
                        which_images{unit}(t) = which_img(img_index);
                        nvr{unit}(t) = novel_vs_repeat(img_index);
                        
                        spikes = find(data(unit).values{t});
                        if ~isempty(spikes)
                            spikeind = spikes(spikes >= imgon & spikes <= imgoff)-imgon;
                            spikeind(spikeind == 0) = []; %happened exactly at the same time
                            temp = zeros(1,length(xn));
                            temp(spikeind) = 1;
                            spike_times{unit}(t,1:length(temp)) = temp;
                        else
                            temp = zeros(1,length(xn));
                            spike_times{unit}(t,1:length(temp)) = temp;
                        end
                    end
                end
            end
        end
        
        %remove excess NaNs associated with error trials
        eyepos = laundry(eyepos,1);
        spike_times = laundry(spike_times,1);
        which_images = laundry(which_images);
        nvr = laundry(nvr);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Perform Spatial analysis and signifance testing---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Determining if neurons are spatially modulated')
        
        info_type = 'spatial';%type of mutual information analysis to perform
        
        %--calculate peak firing rate---%
        peak_firing_rate = NaN(3,num_units);
        for unit = 1:num_units
            for condition = 1:3
                if condition == 1 %novel images
                    trial_data{1} = select_eyepos(eyepos{unit},nvr{unit} == 1);
                    trial_data{2} = spike_times{unit}(nvr{unit} == 1,:);
                elseif condition == 2 %repeat images
                    trial_data{1} = select_eyepos(eyepos{unit},nvr{unit} == 2);
                    trial_data{2} = spike_times{unit}(nvr{unit} == 2,:);
                elseif condition == 3 %all images
                    trial_data{1} = eyepos{unit};
                    trial_data{2} = spike_times{unit};
                end
                
                [filtered_time] = get_smoothed_Time(trial_data{1},imageX,imageY,Fs,binsize,H);
                filtered_time(filtered_time < 0.025) = NaN;
                [filtered_space] = get_smoothed_Space(trial_data{1},trial_data{2},imageX,imageY,binsize,H);
                
                fr = filtered_space./filtered_time; %observed firing rate over space
                peak_firing_rate(condition,unit) = prctile(fr(1:end),99);%don't want to grab outliers
            end
        end
                
        %will perform spatial analysis on novel images, repeat images, and all images
        spatial_info.rate = NaN(3,num_units); %the observed information rate in bits/sec
        spatial_info.shuffled_info_rate = cell(3,num_units); %bootstrapped information rate in bits/sec expected by chance from spike train
        trial_data{3} = [imageX imageY];
        for unit = 1:num_units
            if any(peak_firing_rate(:,unit) > fr_threshold) %don't run on firing rates that are too low
                for condition = 1:3
                    if condition == 1 %novel images
                        trial_data{1} = select_eyepos(eyepos{unit},nvr{unit} == 1);
                        trial_data{2} = spike_times{unit}(nvr{unit} == 1,:);
                    elseif condition == 2 %repeat images
                        trial_data{1} = select_eyepos(eyepos{unit},nvr{unit} == 2);
                        trial_data{2} = spike_times{unit}(nvr{unit} == 2,:);
                    elseif condition == 3 %all images
                        trial_data{1} = eyepos{unit};
                        trial_data{2} = spike_times{unit};
                    end
                    [observed_info_rate,shuffled_info_rate]...
                        = estimated_mutual_information(trial_data,numshuffs,info_type,[binsize,filter_width],Fs);
                    
                    spatial_info.rate(condition,unit) = observed_info_rate.skaggs;
                    spatial_info.spatialstability(condition,unit) = observed_info_rate.spatialstability;
                    spatial_info.shuffled_info_rate{condition,unit} = shuffled_info_rate.skaggs;
                    spatial_info.shuffled_spatialstability{condition,unit} = shuffled_info_rate.spatialstability;
                end
            end
        end
        
        %estimate the significance threshold has the 95-percentile of the
        %bootstrapped data
        for condition = 1:3
            for unit = 1:num_units
                if ~isempty(spatial_info.shuffled_info_rate{condition,unit})
                    spatial_info.shuffled_rate_prctile(condition,unit) = ...
                        100*sum(spatial_info.rate(condition,unit) > ...
                        spatial_info.shuffled_info_rate{condition,unit})/numshuffs;
                    spatial_info.shuffled_spatialstability_prctile(condition,unit) = ...
                        100*sum(spatial_info.spatialstability(condition,unit) > ...
                        spatial_info.shuffled_spatialstability{condition,unit})/numshuffs;
                else
                       spatial_info.shuffled_rate_prctile(condition,unit) = NaN;
                       spatial_info.shuffled_spatialstability_prctile(condition,unit)=NaN;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Plot and Save Figures of Results---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        task_type = 'List_spatial';
        unit_names.multiunit = multiunit;
        unit_names.name = unit_stats(1,:);
        spatial_analysis_plotsV2(figure_dir,task_file,eyepos,spike_times,spatial_info,task_type,...
            unit_names,[binsize,filter_width],imageX,imageY,nvr,Fs,peak_firing_rate);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Perform Spatial analysis on Time shifted Data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         spike_times_ts = [];
%         spatial_info_ts = [];
%         spike_times_ts = cell(length(temporal_shifts),num_units);
%         for ts = 1:length(temporal_shifts)
%             for unit = 1:num_units
%                 spike_times_ts{ts,unit} = circshift_row_non_random(spike_times{unit},temporal_shifts(ts));
%             end
%         end
%         
%         info_type = 'spatial_noshuff'; %don't need to shuffle the data again already did above
%         %will perform spatial analysis on novel images, repeat images, and all images
%         spatial_info_ts.rate = NaN(3,num_units,length(temporal_shifts)); %the observed information rate in bits/sec
%         spatial_info_ts.shuffled_95_percentile = NaN(3,num_units,length(temporal_shifts));
%         spatial_info_ts.shuffled_90_percentile = NaN(3,num_units,length(temporal_shifts));
%         trial_data{3} = [imageX imageY];
%         
%         for condition = 1:3
%             for unit = 1:num_units
%                 for ts = 1:length(temporal_shifts)
%                     if condition == 1 %novel images
%                         trial_data{1} = select_eyepos(eyepos{unit},nvr{unit} == 1);
%                         trial_data{2} = spike_times_ts{ts,unit}(nvr{unit} == 1,:);
%                     elseif condition == 2 %repeat images
%                         trial_data{1} = select_eyepos(eyepos{unit},nvr{unit} == 2);
%                         trial_data{2} = spike_times_ts{ts,unit}(nvr{unit} == 2,:);
%                     elseif condition == 3 %all images
%                         trial_data{1} = eyepos{unit};
%                         trial_data{2} = spike_times_ts{ts,unit};
%                     end
%                     [spatial_info_ts.rate(condition,unit,ts),~]...
%                         = estimated_mutual_information(trial_data,numshuffs,info_type,[binsize,filter_width],Fs);
%                 end
%                 spatial_info_ts.shuffled_95_percentile(condition,unit) = spatial_info.shuffled_95_percentile(condition,unit); %grab the 95% from above
%                 spatial_info_ts.shuffled_90_percentile(condition,unit) = spatial_info.shuffled_90_percentile(condition,unit); %grab the 90% from above
%             end
%         end
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%---Plot and Save Figures of Time Shifted Results---%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         task_type = 'List_spatial_time_shifted';
%         unit_names.multiunit = multiunit;
%         unit_names.name = unit_stats(1,:);
%         spatial_info_ts.temporal_shifts = temporal_shifts;
%         spatial_analysis_plotsV2(figure_dir,task_file,eyepos,spike_times_ts,spatial_info_ts,task_type,...
%             unit_names,[binsize,filter_width],imageX,imageY,nvr,Fs);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Finally save all the data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
            'spike_times','eyepos','spatial_info','nvr','which_images',...
            'binsize','filter_width','peak_firing_rate')
        %'spike_times_ts','spatial_info_ts')
        disp(['Spatial Data Analyis for ' task_file(1:end-11) ' saved']);
end
end


function [selected] = select_eyepos(eyepos,select_rows)
%x eye position is in odd rows and y ye position is in even rows
%select rows should be a logical index
x = eyepos(1:2:end,:);
y = eyepos(2:2:end,:);
x = x(select_rows,:);
y = y(select_rows,:);
selected = NaN(sum(select_rows)*2,size(x,2));
selected(1:2:end,:) = x;
selected(2:2:end,:) = y;
end

function B = circshift_row_non_random(A,D)
%D: for shift amount
%function is equivalent to circshift_row except D is stated instead of random
%assumes NaNs are at the end of the data
if isempty(A)
    B = A;
elseif D > 0 %so shift right
    if any(any(isnan(A))) %if their are NaNs must be fillers to speed up computation
        [m, n] = size(A);
        B = NaN(m, n);
        
        for i = (1 : m)
            n = sum(~isnan(A(i,:)));
            if n~=0
                temp = A(i,1:n);
                B(i,1:n) = [temp(n - D + 1 : n) temp(1 : n - D)];
            end
        end
    else
        [m, n] = size(A);
        B = zeros(m, n);
        
        for i = (1 : m)
            B(i,:) = [A(i,(n - D + 1 : n)) A(i,(1 : n - D ))];
        end
    end
else %D < 0 so shift left
    D = -D;
    if any(any(isnan(A))) %if their are NaNs must be fillers to speed up computation
        [m, n] = size(A);
        B = NaN(m, n);
        
        for i = (1 : m)
            n = sum(~isnan(A(i,:)));
            if n~=0
                temp = A(i,1:n);
                B(i,1:n) = [temp(D+1: end) temp(1 : D)];
            end
        end
    else
        [m, n] = size(A);
        B = zeros(m, n);
        
        for i = (1 : m)
            B(i,:) = [A(i,(D+1: end)) A(i,(1 : D))];
        end
    end
end
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