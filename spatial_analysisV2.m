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
%
% Code rechecked ListSQ section for bugs October 17-18, 2016 SDK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set important task parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_on_twin = 500;%how much time to ignore eye movements for to remove
%strong visual response (though some may last longer) combined with strong central bias

numshuffs = 1000;% # of shuffles for bootstraspping, recommend this is at least 1000
binsize = 12; %pixels per bin spatial bin in either dimension 1/2 dva
filter_width = 6; %std of 2D guassian filter ~ 3 dva, could use 2 dva (4) as well
H = define_spatial_filter(filter_width);
info_type = 'spatial';%type of mutual information analysis to perform


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials',...
    'hdr','whole_session_mean_firing_rate','excitatory_inhibitory');
%grab unit data
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

if num_units == 0; %if no units exit function
    disp([task_file(1:8) ': no units could be found. Exiting function...'])
    return;
end

%---determine number of blocks unit was stable for
%remove units with too few trials
%these are the absolute minimum data required to do data analysis
if strcmpi(task,'cvtnew')
    min_blks = 1; %~100 trials may want to change
else
    min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
    %this means at least 64 viewed images counting novel+repeat presentations
end
valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end

switch task
    case {'cvtnew','CVTNEW'}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import & reformat data so that spikes are locked to dot position---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        disp('Collecting Spike Locations')
        
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
        
        %--calculate peak firing rate---%
        peak_firing_rate = NaN(1,num_units);
        for unit = 1:num_units
            trial_data{1} = dotpos{unit};
            trial_data{2} = spike_times{unit};
            
            if nansum(nansum(trial_data{2})) == 0 %this occassiionally happens with really specifc/low firing rate neurons...
                %and we don't want a NaN
                peak_firing_rate(unit) = 0;
                continue
            end
            [filtered_time] = filter_time(trial_data{1},imageX,imageY,Fs,binsize,H);
            filtered_time(filtered_time == 0) = NaN;
            [filtered_space] = filter_space(trial_data{1},trial_data{2},imageX,imageY,binsize,H);
            
            fr = filtered_space./filtered_time; %observed firing rate over space
            peak_firing_rate(unit) = prctile(fr(1:end),99);%don't want to grab outliers
        end
        
        info_type = 'spatial_cvtnew';%type of mutual information analysis to perform
        
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
        imageX = 800; %horizontal size of the image in pixels
        imageY = 600; %horizontal size of the image in pixels
        
        %set the image duration
        if str2double(task_file(3:8)) < 140805 %first sets done by PW before this data had 7 s images
            imgdur = 7000;
        else %rest of image presentations were 5 seconds
            imgdur = 5000;
        end
        
        %get important task specific information
        [itmlist,sequence_items] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);
        
        %---preallocate space and parallel structure of cfg---%
        %Each unit gets own cell since each unit may have different amount of trials/data
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
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code); %should be first code in sequence
                imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %want to avoid visual response
                imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start;
                
                % if monkey isn't paying attention and looked away image presentation
                % is now longer than imgdur (because of cumulative looking time) 
                % so data isn't probably worth much plus have to cut off somewhere
                if imgoff-imgon > 1.5*imgdur-1 %cut off trial at 1.5x length of desired looking time
                    imgoff = imgon+1.5*imgdur-1;
                end
                imgon = imgon+image_on_twin; %to avoid visual response and strong central bias
                
                xn = round(data(eyechans(1)).values{t}(imgon:imgoff)); %horizontal eye data
                yn = round(data(eyechans(2)).values{t}(imgon:imgoff)); %vertical eye data
                
                img_index = find(img_cnd == cfg.trl(t).cnd); %get image index
                
                if any(isnan(which_img(img_index)))%presentation error so skip trial, see get_image_numbers.m
                    continue
                end

                for unit = 1:num_units
                    if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
                        eyepos{unit}(2*t-1,1:length(xn)) = xn; %horizontal eye data
                        eyepos{unit}(2*t,1:length(yn)) = yn; %vertical eye data
                       
                        which_images{unit}(t) = which_img(img_index);
                        nvr{unit}(t) = novel_vs_repeat(img_index);
                        
                        spikes = find(data(unit).values{t});
                        if ~isempty(spikes)
                            spikeind = spikes(spikes > imgon & spikes <= imgoff)-imgon;
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
        eyepos = laundry(eyepos,1); %only remove rows aka trials without data
        spike_times = laundry(spike_times,1);  %only remove rows aka trials without data
        which_images = laundry(which_images);
        nvr = laundry(nvr);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Perform Spatial analysis and signifance testing---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Determining if neurons are spatially modulated')
        
        %--calculate peak firing rate---%
        peak_firing_rate = NaN(3,num_units);%row 1 novel, row 2 repeat, row 3 all images
        for unit = 1:num_units
            if isnan(valid_trials(1,unit)); %no valid trials for unit
                continue
            end
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
                if nansum(nansum(trial_data{2})) == 0 %this occassiionally happens with really specifc/low firing rate neurons...
                    %and we don't want a NaN
                    peak_firing_rate(condition,unit) = 0;
                    continue
                end
                
                [ratemap,~] = get_firing_rate_map(trial_data,imageX,imageY,binsize,H,Fs,'all');
                peak_firing_rate(condition,unit) = prctile(ratemap(:),97.5);%don't want to grab outliers
            end
        end
        
        %---Determine if Neuron is more spatial than chance---%
        spatial_info.rate = NaN(1,num_units); %the observed information rate in bits/sec
        spatial_info.shuffled_info_rate = cell(1,num_units); %bootstrapped information rate in bits/sec expected by chance from spike train
        spatial_info.shuffled_rate_prctile = NaN(1,num_units);
        
        %spatial correlation first vs second half row 1 spearman row 2 Kendall
        spatial_info.spatialstability_halves = NaN(1,num_units);
        spatial_info.shuffled_spatialstability_halves = cell(1,num_units);
        spatial_info.spatialstability_halves_prctile = NaN(1,num_units);
        
        %spatial correlation even and odd trials row 1 spearman row 2 Kendall
        spatial_info.spatialstability_even_odd = NaN(1,num_units);
        spatial_info.shuffled_spatialstability_even_odd = cell(1,num_units);
        spatial_info.spatialstability_even_odd_prctile = NaN(1,num_units);
        
        trial_data{3} = [imageX imageY];
        for unit = 1:num_units
            if ~isempty(spike_times{unit})
                trial_data{1} = eyepos{unit};
                trial_data{2} = spike_times{unit};
                [observed_info_rate,shuffled_info_rate]...
                    = estimated_mutual_information(trial_data,numshuffs,info_type,[binsize,filter_width],Fs);
                
                %skaggs information score
                spatial_info.rate(unit) = observed_info_rate.skaggs;
                spatial_info.shuffled_info_rate{unit} = shuffled_info_rate.skaggs;
                spatial_info.shuffled_rate_prctile(unit) = 100*sum(...
                    spatial_info.rate(unit) > spatial_info.shuffled_info_rate{unit})/numshuffs;
                
                %spatial correlation first half vs second half
                spatial_info.spatialstability_halves(unit) = observed_info_rate.spatialstability_halves;
                spatial_info.shuffled_spatialstability_halves{unit} = shuffled_info_rate.spatialstability_halves;
                spatial_info.spatialstability_halves_prctile(unit) = ...
                    100*sum(spatial_info.spatialstability_halves(unit) > ...
                    spatial_info.shuffled_spatialstability_halves{unit})/numshuffs;
            
                
                %spatial correlation for even and odd trials
                spatial_info.spatialstability_even_odd(unit) = observed_info_rate.spatialstability_even_odd;
                spatial_info.shuffled_spatialstability_even_odd{unit} =...
                    shuffled_info_rate.spatialstability_even_odd;
                spatial_info.spatialstability_even_odd_prctile(unit) = ...
                    100*sum(spatial_info.spatialstability_even_odd(unit) > ...
                    spatial_info.shuffled_spatialstability_even_odd{unit})/numshuffs;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Plot and Save Figures of Results---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        task_type = 'List_spatial';
        unit_names.multiunit = multiunit;
        unit_names.name = unit_stats(1,:);
        unit_names.avg_firing_rate = whole_session_mean_firing_rate;
        unit_names.putative_EI = excitatory_inhibitory;
        spatial_analysis_plotsV2(figure_dir,task_file,eyepos,spike_times,spatial_info,task_type,...
            unit_names,[binsize,filter_width],imageX,imageY,nvr,Fs,peak_firing_rate);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Finally save all the data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
            'spike_times','eyepos','spatial_info','nvr','which_images',...
            'binsize','filter_width','peak_firing_rate','image_on_twin','Fs',...
            'whole_session_mean_firing_rate','excitatory_inhibitory','unit_names')
        disp(['ListSQ Spatial Data Analyis for ' task_file(1:end-11) ' saved']);
end
end