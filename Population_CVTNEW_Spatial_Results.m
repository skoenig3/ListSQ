%Population cvtnew spatial analysis results
%draft written SDK 1/28/2017


task = 'cvtnew';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size

%---Misc. Parameters (Mostly Place Cells)---%
monkey_all_unit_count = zeros(2,2);%row 1 place row 2 non-place, column by monkey
all_multi_unit_count = zeros(1,2); %all_multi_unit_countunits
all_place_cell_unit_names = {}; %place cell unit names
all_place_cell_monkeys = []; %1s and 2s
all_place_cell_spatial_skaggs = [];
all_non_place_cell_spatial_skaggs = [];

%---Spatial Correlation Values for all cells---%
all_place_cell_spatial_corrs = []; %place cell spatial correlations
all_non_place_cell_spatial_corrs = [];%non-place cell spatial correlations

monkeys = {'Vivian','Tobii'};
figure_dir = {};
for monk =2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir{1} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir{2} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    end
    
    for sess = 1:length(session_data)
        %read in task file data
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        %load task file data
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials',...
            'stability_attribute','cfg','hdr','data');
        
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
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'spatial_info')
              
        for unit = 1:num_units
            if ~isempty(spatial_info.shuffled_info_rate{unit})&& length(spatial_info.shuffled_info_rate{unit}) < 1000
                error('Should have 1000 shuffles') %make sure file has right number of shuffles
            elseif isempty(spatial_info.shuffled_info_rate{unit})
               continue  
            end
            
            if spatial_info.shuffled_rate_prctile(unit) > 95 && spatial_info.spatialstability_halves_prctile(unit) > 95
                %---Misc. Parameters (Mostly Place Cells)---%
                monkey_all_unit_count(1,monk) = monkey_all_unit_count(1,monk)+1;%row 1 place row 2 non-place, column by monkey
                all_place_cell_unit_names = [all_place_cell_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}];
                all_place_cell_monkeys = [all_place_cell_monkeys monk]; %1s and 2s
                
                %---Spatial Correlation Values for all cells---%
                all_place_cell_spatial_corrs = [all_place_cell_spatial_corrs spatial_info.spatialstability_halves(unit)]; %place cell spatial correlations
                all_place_cell_spatial_skaggs = [all_place_cell_spatial_skaggs spatial_info.shuffled_rate_prctile(unit)];
                
            else
                
                monkey_all_unit_count(2,monk) = monkey_all_unit_count(2,monk)+1;%row 1 place row 2 non-place, column by monkey
                all_non_place_cell_spatial_corrs = [all_non_place_cell_spatial_corrs  spatial_info.spatialstability_halves(unit)];%non-place cell spatial correlations
                all_non_place_cell_spatial_skaggs = [all_non_place_cell_spatial_skaggs spatial_info.shuffled_rate_prctile(unit)];
            end
            
            
        end
    end
end