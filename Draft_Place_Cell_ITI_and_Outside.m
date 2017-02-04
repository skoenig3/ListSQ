%see what place cells are doing during ITI and Outside image
%written by Seth Konig 9.23.16
%
%rechecked for bugs 1/22/2017 by SDK not really a draft anymore but
%keeping same name

clar %clear,clc
%set(0,'DefaultFigureVisible','OFF');


ITIstart_code = 15;
img_on_code= 23;
img_off_code = 24;
imageX = 800; %horizontal size of the image in pixels
imageY = 600; %horizontal size of the image in pixels

task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size

figure_dir = 'C:\Users\seth.koenig\Desktop\Outside\';
monkey_count = zeros(1,2);

monkeys = {'Vivian','Tobii'};
for monk =2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    end
    
    for sess = 46%:length(session_data)
        %read in task file data
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        %load task file data
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials',...
            'stability_attribute','cfg','hdr','data','outside_xy','fixationstats');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        if num_units == 0
            continue
        end
        
        num_trials = length(cfg.trl); %number of trials
        valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
        if all(isnan(valid_trials(1,:)))
            continue
        end
        
        %load spatial analysis data
        disp(task_file(1:8))
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
            'spatial_info','spike_times','eyepos','binsize','filter_width')
        
        
        %set the image duration
        if str2double(task_file(3:8)) < 140805 %first sets done by PW before this data had 7 s images
            max_trial_dur = 7000;
        else %rest of image presentations were 5 seconds
            max_trial_dur = 5000;
        end
        max_trial_dur = max_trial_dur+2000;%buffer for ITI, fixation on cross, etc.
        
        %set/get some general important info
        num_trials = length(cfg.trl);
        H = define_spatial_filter(filter_width);
        Fs = data(1).fsample; %should be 1000
      
        
        for unit = 1:num_units
         if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                    && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                
                monkey_count(monk) = monkey_count(monk)+1;
                
                spike_times_outside = [];
                eyepos_outside = [];
                for t = valid_trials(1,unit):valid_trials(2,unit) %for valit trial numbers
                    if length(data(unit).values{t}) > 1.5*max_trial_dur %if trial duration is longer than max trial duration
                        temp_spk = data(unit).values{t}(1:1.5*max_trial_dur); %spikes for this trial
                        temp_eye = outside_xy{t}(:,1:1.5*max_trial_dur); %eye pos outside screen border for this trial
                    else %if trial is shorter
                        temp_spk = NaN(1,1.5*max_trial_dur);
                        temp_spk(1:length(data(unit).values{t}))= data(unit).values{t}; %spike times for this trial
                        temp_eye = NaN(2,1.5*max_trial_dur); 
                        temp_eye(:,1:length(outside_xy{t})) = outside_xy{t}; %eye data outside screen for this trial
                    end
                    spike_times_outside = [spike_times_outside; temp_spk];
                    eyepos_outside = [eyepos_outside;  temp_eye];
                end
                
                [ratemap_outside,filtered_time_outside] = outside_ratemap(eyepos_outside,Fs,binsize,H,spike_times_outside); %ratemap off screen
                [ratemap,filtered_time] = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all'); %ratemap during image
                maxfr = prctile(ratemap(:),97.5);%scale to 97.5th percentile for visualization
%                 minfr = prctile(ratemap(:),2.5); 
                
                figure
                subplot(2,2,1)
                h = imagesc(ratemap);
                set(h,'alphadata',~isnan(filtered_time));
                axis off
                axis equal
                title('Image Viewing Rate Map')
                colormap('jet')
                clims = caxis;
%                 clims(1) = minfr;
%                 clims(2) = maxfr;
                caxis(clims);
                colorbar
                
                subplot(2,2,2)
                h = imagesc(ratemap_outside);
                set(h,'alphadata',~isnan(filtered_time_outside));
                axis off
                axis equal
                title('Rate Map Outside All trials')
                caxis(clims)
                
                combined_ratemap = ratemap_outside;
                combined_ratemap(13:62,18:83) = ratemap;
                combined_time = filtered_time_outside;
                combined_time(13:62,18:83) = filtered_time;
                subplot(2,2,3)
                h = imagesc(combined_ratemap);
                set(h,'alphadata',~isnan(combined_time));
                axis off
                axis equal
                title('Combined Ratemaps')
                colormap('jet')
                caxis(clims)
                
                subtitle([task_file(1:end-11) ' ' unit_names{unit}])
                
                save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_outside']);
            end
        end
    end
end
set(0,'DefaultFigureVisible','ON')