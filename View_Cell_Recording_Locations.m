%get recording locations for view cells
clar
monkeys = {'Vivian','Tobii'};

task = 'ListSQ';

locations = [2:17];
all_unit_positions = zeros(1,length(locations));
view_cell_poisitions = zeros(1,length(locations));
for monk = 1:2
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
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials','stability_attribute','cfg','hdr','data');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
            'spatial_info','peak_firing_rate')
        
        %these are the absolute minimum data required to do data analysis may want
        %to be more strigent later but not worth doing analysis (especially
        %shuffling) on such few trials for these neurons
        if str2double(task_file(3:8)) < 140805 %7 second viewing with fewere sequence trials
            minimum_trials_1 = 101; %for session start minimum number of trials: includes fam block and 16 novel/repeat images + sequence trials
            minimum_trials_2 =  80;%for rest of session: includes 16 novel/repeat images + sequence trials
        else
            minimum_trials_1 = 117; %for session start minimum number of trials: includes fam block and 16 novel/repeat images + sequence trials
            minimum_trials_2 =  96;%for rest of session: includes 16 novel/repeat images + sequence trials
        end
        
        
        for unit = 1:size(valid_trials,2)
            %unstable/usable units
            if stability_attribute(unit) ~=  1 %not used unit
                if any(~isnan(peak_firing_rate(:,unit)))
                    error('Unit was stable and should have been processed')
                end
                continue %move to the next unit
            end
            
            %usable units
            if multiunit(unit) ==1 %then multiunit and don't want
                peak_firing_rate(:,unit) = NaN;
                spatial_info.shuffled_rate_prctile(:,unit) = NaN;
                spatial_info.shuffled_spatialstability_prctile(:,unit) = NaN;
                continue
            end
            

            if any(isnan(spatial_info.rate(:,unit)) & peak_firing_rate(:,unit) > 1)
                error('Firing rate for unit is > 1 Hz but not processed')
            end
            
            %---determine number of blocks unit was stable for
            start_end = valid_trials(:,unit);
            if isnan(start_end(1))
                start_end(1) = 1;
            end
            if isnan(start_end(2))
                start_end(2) = length(cfg.trl);
            end
            start_end(start_end == 0) = 1;
            min_trial = cfg.trl(start_end(1)).cnd-1000; %get the first condition number
            max_trial = cfg.trl(start_end(2)).cnd-1000; %get the last condition number
            
            if max_trial < minimum_trials_1
                error('Unit should have been stable for more trials to be processed')
            end
            
            if min_trial < 22 %includes fam block
                min_trial = 22; %remove then count from there
            end
            num_blks = floor((max_trial-min_trial+1)/minimum_trials_2);
            
            if num_blks < 2 %remove from count
                peak_firing_rate(:,unit) = NaN;
                spatial_info.shuffled_rate_prctile(:,unit) = NaN;
                spatial_info.shuffled_spatialstability_prctile(:,unit) = NaN;
            end
            
        end
        
        total_units = sum(sum(~isnan(peak_firing_rate)) > 0);
        if total_units == 0
            continue
        end
        total_spatial_units = sum(sum((spatial_info.shuffled_rate_prctile > 95) & (spatial_info.shuffled_spatialstability_prctile > 95)) > 0);
        AP_location = round(chamber_zero(1)+ session_data{sess}.location(1));
        location_index = find(locations == AP_location);
        
        all_unit_positions(location_index) = all_unit_positions(location_index)+total_units;
        view_cell_poisitions(location_index) = view_cell_poisitions(location_index)+total_spatial_units;
        
    end
end
%%
half_threshold = 10.5;
prop = view_cell_poisitions./all_unit_positions;
prop_posterior = nanmean(prop(locations < half_threshold));
prop_anterior =  nanmean(prop(locations > half_threshold));

figure
bar(prop)
hold on
plot([9.5 9.5],[0 0.5],'k--')
hold off
xlabel('AP Recording Location')
ylabel('Proportion of Units that are view cells')
xlim([4 15])
title(sprintf([num2str(100*prop_anterior,2) '%% anterior & ' num2str(100*prop_posterior,2) '%% of Units are view cells']))