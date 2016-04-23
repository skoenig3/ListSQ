function spatial_analysis(data_dir,preprocessed_data_file,figure_dir,task)
% written by Seth Konig August, 2014
% Function analyses where the monkey is attending when spikes occur
%
% Inputs:
%   1) data_dir: direction where preprocessed data is located
%   2) preprocessed_data_file: preprocessed ListSQ data file with LFPs
%   divided by channel and trial.
%   3) figure_dir: location of where to put save figures
%   4) task: what task was this data come from with i.e. 'cvtnew','ListSQ'
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-spatial_analysis_results.mat'

%filtering and estimated mutual information parameters
binsize = 12; %pixels per bin spatial bin in either dimension ~1/2 dva
filter_width = 3; %std of 2D guassian filter ~ 2 dva
numshuffs = 100; %recommend this is between 100 & 1000

switch task
    case 'cvtnew'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import & reformat data so that spikes are locked to dot position---%%%
        load([data_dir preprocessed_data_file],'data','cfg','multiunit','num_units','meta','hdr');
        
        [eyechans] = find_desired_channels(cfg,'eye');
        disp('Collecting Spike Locations')
        
        %set/get some general important info
        num_trials = length(cfg.trl);
        Fs = data(1).fsample; %should be 1000
        ITIstart_code = 15;
        dot_on_code = 25;
        dot_clrchng_code = 27;
        bar_code_response = 4; %monkey made a move
        reward_code = 3;
        imageX = 800; %horizontal size of the screen
        imageY = 600; %horizontal size of the screen
        
        %preallocate space and parallel structure of cfg
        dotpos = NaN(2*length(cfg.trl),3000);
        spike_times = cell(1,num_units);
        for unit = 1:num_units
            spike_times{unit} = NaN(length(cfg.trl),3000);
        end
        
        for t = 1:length(cfg.trl);
            if any(cfg.trl(t).allval == reward_code); %only take correct trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                pathon =  cfg.trl(t).alltim(cfg.trl(t).allval == dot_on_code)-trial_start; %meta data starts at 1st 666 which appears before dot actually turns on in event array
                dot_clrchng = cfg.trl(t).alltim(cfg.trl(t).allval == dot_clrchng_code)-trial_start;
                responded = cfg.trl(t).alltim(cfg.trl(t).allval == bar_code_response)-trial_start;
                
                xn = interp1(meta(t).sample,meta(t).x,meta(t).sample(1):meta(t).sample(end),'cubic');
                yn = interp1(meta(t).sample,meta(t).y,meta(t).sample(1):meta(t).sample(end),'cubic');
                
                dotpos(2*t-1,1:length(xn)) = round(xn);
                dotpos(2*t,1:length(xn)) = round(yn);
                
                for unit = 1:num_units
                    spikes = find(data(unit).values{t});
                    if ~isempty(spikes)
                        spikeind = spikes(spikes > pathon & spikes <= responded)-pathon;
                        spikeind(spikeind < 1) = []; %should only happen when spikes occur at the same time as tstart
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
        
        %remove excess NaNs associated with error trials
        dotpos = laundry(dotpos);
        spike_times = laundry(spike_times);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Perform Spatial analysis and signifance testing---%%%
        disp('Determining if neurons are spatially modulated')
        
        info_type = 'spatial';%type of mutual information analysis to perform
        
        %will perform spatial analysis on novel images, repeat images, and all images
        spatial_info.rate = NaN(1,num_units); %the observed information rate in bits/sec
        spatial_info.shuffled_info_rate = cell(1,num_units); %bootstrapped information rate in bits/sec expected by chance from spike train
        spatial_info.shuffled_95_percentile = NaN(1,num_units);
        spatial_info.shuffled_90_percentile = NaN(1,num_units);
        trial_data{1} = dotpos;
        trial_data{3} = [imageX imageY];
        for unit = 1:num_units
            trial_data{2} = spike_times{unit};
            [spatial_info.rate(unit),spatial_info.shuffled_info_rate{unit}]...
                = estimated_mutual_information(trial_data,numshuffs,info_type,[binsize,filter_width],Fs);
        end
        
        %estimate the significance threshold has the 95-percentile of the
        %bootstrapped data
        for unit = 1:num_units
            if ~isempty(spatial_info.shuffled_info_rate{unit})
                spatial_info.shuffled_95_percentile(unit) =...
                    prctile(spatial_info.shuffled_info_rate{unit},95);
                spatial_info.shuffled_90_percentile(unit) =...
                    prctile(spatial_info.shuffled_info_rate{unit},90);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Plot and Save Figures of Results---%%%
        task_type = 'cvtnew_spatial';
        unit_names.multiunit = multiunit;
        unit_names.name = cfg.channel;
        spatial_analysis_plots(figure_dir,preprocessed_data_file,dotpos,spike_times,spatial_info,task_type,...
            unit_names,[binsize,filter_width],imageX,imageY,NaN,Fs);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Finally save all the data---%%%
        save([data_dir preprocessed_data_file(1:10) '-spatial_analysis_results.mat'],...
            'spike_times','dotpos','spatial_info','binsize','filter_width')
        disp(['Spatial Data Analyis for ' preprocessed_data_file(1:10) ' saved']);
        
    case 'ListSQ'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import & reformat data so that spikes are locked to eye movements---%%%
        load([data_dir preprocessed_data_file],'data','cfg','multiunit','num_units','item_set');
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
        if str2double(preprocessed_data_file(3:8)) < 140805
            imgdur = 7000;
        else
            imgdur = 5000;
        end
        
        %get important task specific information
        [itmlist,sequence_items] = read_ListSQ_itm_and_cnd_files(item_set);
        [which_img,novel_vs_repeat] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);
        
        %preallocate space and parallel structure of cfg
        eyepos = NaN(2*length(cfg.trl),imgdur*1.5);
        spike_times = cell(1,num_units);
        for unit = 1:num_units
            spike_times{unit} = NaN(length(cfg.trl),imgdur*1.5);
        end
        
        for t = 1:num_trials
            if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                
                imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start+500;%to ignore 1st fixation
                imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start;
                
                % if monkey isn't paying attention data isn't probably
                % worth much plus have to cut off somewhere
                if imgoff-imgon > 1.5*imgdur-1
                    imgoff = imgon+1.5*imgdur-1;
                end
                
                xn = data(eyechans(1)).values{t}(imgon:imgoff);
                yn = data(eyechans(2)).values{t}(imgon:imgoff);
                eyepos(2*t-1,1:length(xn)) = round(xn);
                eyepos(2*t,1:length(xn)) = round(yn);
                
                for unit = 1:num_units
                    spikes = find(data(unit).values{t});
                    if ~isempty(spikes)
                        spikeind = spikes(spikes > imgon & spikes <= imgoff)-imgon;
                        spikeind(spikeind < 1) = []; %should only happen when spikes occur at the same time as tstart
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
        
        %remove excess NaNs associated with error trials
        eyepos = laundry(eyepos')';
        eyepos = [eyepos NaN(size(eyepos,1),7500-size(eyepos,2))]; %sometimes remove structure from one but not all
        spike_times = laundry(spike_times);
        which_img = which_img(1:size(spike_times{1},1));
        novel_vs_repeat = novel_vs_repeat(1:size(spike_times{1},1));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Perform Spatial analysis and signifance testing---%%%
        disp('Determining if neurons are spatially modulated')
        
        info_type = 'spatial';%type of mutual information analysis to perform
        
        %will perform spatial analysis on novel images, repeat images, and all images
        spatial_info.rate = NaN(3,num_units); %the observed information rate in bits/sec
        spatial_info.shuffled_info_rate = cell(3,num_units); %bootstrapped information rate in bits/sec expected by chance from spike train
        spatial_info.shuffled_95_percentile = NaN(3,num_units);
        spatial_info.shuffled_90_percentile = NaN(3,num_units);
        trial_data{3} = [imageX imageY];
        for condition = 1:3
            for unit = 1:num_units
                if condition == 1 %novel images
                    trial_data{1} = select_eyepos(eyepos,novel_vs_repeat == 1);
                    trial_data{2} = spike_times{unit}(novel_vs_repeat == 1,:);
                elseif condition == 2 %repeat images
                    trial_data{1} = select_eyepos(eyepos,novel_vs_repeat == 1);
                    trial_data{2} = spike_times{unit}(novel_vs_repeat == 2,:);
                elseif condition == 3 %all images
                    trial_data{1} = eyepos;
                    trial_data{2} = spike_times{unit};
                end
                [spatial_info.rate(condition,unit),spatial_info.shuffled_info_rate{condition,unit}]...
                    = estimated_mutual_information(trial_data,numshuffs,info_type,[binsize,filter_width],Fs);
            end
        end
        
        %estimate the significance threshold has the 95-percentile of the
        %bootstrapped data
        for event = 1:3
            for unit = 1:num_units
                if ~isempty(spatial_info.shuffled_info_rate{event,unit})
                    spatial_info.shuffled_95_percentile(event,unit) =...
                        prctile(spatial_info.shuffled_info_rate{event,unit},95);
                    spatial_info.shuffled_90_percentile(event,unit) =...
                        prctile(spatial_info.shuffled_info_rate{event,unit},90);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Plot and Save Figures of Results---%%%
        task_type = 'List_spatial';
        unit_names.multiunit = multiunit;
        unit_names.name = cfg.channel;
        spatial_analysis_plots(figure_dir,preprocessed_data_file,eyepos,spike_times,spatial_info,task_type,...
            unit_names,[binsize,filter_width],imageX,imageY,novel_vs_repeat,Fs);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Finally save all the data---%%%
        save([data_dir preprocessed_data_file(1:10) '-spatial_analysis_results.mat'],'spike_times',...
            'eyepos','spatial_info','novel_vs_repeat','binsize','filter_width')
        disp(['Spatial Data Analyis for ' preprocessed_data_file(1:10) ' saved']);
end
end

function [which_img,novel_vs_repeat] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code)
% define which imagse are novel and which images are repeat
which_img = NaN(1,96*2);
img_count = 1;
for t = 1:length(cfg.trl);
    if itmlist(cfg.trl(t).cnd-1000) > sequence_items(end)
        if any(cfg.trl(t).allval == img_on_code);
            which_img(img_count) = itmlist(cfg.trl(t).cnd-1000);
            img_count = img_count+1;
        end
    end
end
which_img = which_img-which_img(1)+1;

novel_vs_repeat = NaN(1,96*2);
for img = 1:max(which_img)
    imgind = find(which_img == img);
    dimgind = find(diff(imgind) > 1);
    if ~isempty(imgind);
        novel_vs_repeat(imgind(1:dimgind)) = 1;
        novel_vs_repeat(imgind(dimgind+1:end)) = 2;
    end
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