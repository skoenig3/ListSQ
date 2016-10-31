%see what place cells are doing during ITI and Outside image
%written by Seth Konig 9.23.16

%get recording locations for view cells
clar
monkeys = {'Vivian','Tobii'};

task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)

imageX = 800;
imageY = 600;

Fs = 1000;
monkey_count = zeros(1,2);
all_unit_names = {};
all_context_gain = [];
for monk = 2:-1:1
    monkey = monkeys{monk};
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
        
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    end
    
    for sess = 1:length(session_data)
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials',...
            'stability_attribute','cfg','hdr','data','fixationstats','outside_xy');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        
        if num_units == 0
            continue
        end
        
        num_trials = length(cfg.trl);
        valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks);
        if all(isnan(valid_trials(1,:)))
            continue
        end
        
        %set the image duration
        if str2double(task_file(3:8)) < 140805
            imgdur = 7000+1500; %add 1500 for ITI and period before during crosshair
        else
            imgdur = 5000+1500;%add 1500 for ITI and period before during crosshair
        end

        disp(task_file(1:8))
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'])   
        
        filter_size = filter_width*10;
        H = fspecial('gaussian',filter_size,filter_width);

        for unit = 1:num_units
            if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                && (spatial_info.spatialstability_even_odd_prctile(2,unit) > 95) ... %spatial consistency 95%+
                && (spatial_info.spatialstability_halves_prctile(2,unit) > 95) %spatial stability 95%+
                monkey_count(monk) = monkey_count(monk)+1;
                
                spike_times_outside = [];
                eyepos_outside = [];
                ITI_eye = [];
                ITI_spikes = [];
                for t = valid_trials(1,unit):valid_trials(2,unit)
                    if length(data(unit).values{t}) > 1.5*imgdur
                        temp_spk = data(unit).values{t}(1:1.5*imgdur);
                        temp_eye = outside_xy{t}(:,1:1.5*imgdur);
                    else
                        temp_spk = NaN(1,1.5*imgdur);
                        temp_spk(1:length(data(unit).values{t}))= data(unit).values{t};
                        temp_eye = NaN(2,1.5*imgdur);
                        temp_eye(:,1:length(outside_xy{t})) = outside_xy{t};
                    end
                    spike_times_outside = [spike_times_outside; temp_spk];
                    eyepos_outside = [eyepos_outside;  temp_eye];
                    
                    
                    ITI_eye = [ITI_eye; fixationstats{t}.XY(:,1:1000)];
                    ITI_spikes = [ITI_spikes data(unit).values{t}(1:1000)];
                     
                end
                ITI_eye = round(ITI_eye);
                ITI_eye(ITI_eye < 1) = 1;
                
                 [ratemap_outside,filtered_time_outside] = outside_ratemap(eyepos_outside,Fs,binsize,H,spike_times_outside);
                 
                 [ratemap,filtered_time] = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all');
                 
                 [ratemap_ITI,filtered_time_ITI] = get_firing_rate_map({ITI_eye,ITI_spikes},imageX,imageY,binsize,H,Fs,'all');
                 
                 figure
                 subplot(2,2,1)
                 h = imagesc(ratemap);
                 set(h,'alphadata',~isnan(filtered_time));
                 axis off
                 axis equal
                 title('Image Viewing Rate Map')
                 colormap('jet')
                 clims = caxis;
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
                 
                 subplot(2,2,4)
                 h = imagesc(ratemap_ITI);
                 set(h,'alphadata',~isnan(filtered_time_ITI));
                 axis off
                 axis equal
                 title('ITI Rate Map')
                 colormap('jet')
                 caxis(clims)
            end
        end
    end
end