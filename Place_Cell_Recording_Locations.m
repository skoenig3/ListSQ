%get recording locations for view cells
clar
monkeys = {'Vivian','Tobii'};

task = 'ListSQ';

locations = [2:17];
all_unit_positions = zeros(1,length(locations));
place_cell_positions = zeros(1,length(locations));
all_peak_firing_rates = [];
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
        
        if num_units == 0
            continue
        end
        
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
            
            %usable units
            if multiunit(unit) ==1 
                continue
            end
            
            if isnan(spatial_info.rate(unit)) || spatial_info.rate(unit) == 0 %didn't run on unit since not stable
                continue
            end
            
            AP_location = round(chamber_zero(1)+ session_data{sess}.location(1));
            location_index = find(locations == AP_location);
            if spatial_info.shuffled_rate_prctile(unit) > 99  %%&& spatial_info.shuffled_spatialstability_prctile(unit) > 95
                all_peak_firing_rates = [all_peak_firing_rates peak_firing_rate(3,unit)];
                place_cell_positions(location_index) = place_cell_positions(location_index)+1;
            end
            all_unit_positions(location_index) = all_unit_positions(location_index)+1;            
        end
    end
end
%%
half_threshold = 10.5;
prop = place_cell_positions./all_unit_positions;
prop_posterior = nanmean(prop(locations < half_threshold));
prop_anterior =  nanmean(prop(locations > half_threshold));

figure
bar(prop)
hold on
plot([half_threshold half_threshold],[0 0.5],'k--')
hold off
xlabel('AP Recording Location')
ylabel('Proportion of Units that are view cells')
xlim([4 15])
title(sprintf([num2str(100*prop_anterior,2) '%% anterior & ' num2str(100*prop_posterior,2) '%% of Units are view cells']))