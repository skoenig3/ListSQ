% Population Saccadic Modulation
% written Seth Konig 8/12/16
%
% code imports data from List_fixation_AnalysisV results and looks at how
% modulated the whole population is by fixations.
% Code does the following
% 1) Summarizes fixation modulation for for place cell and non-place cells
% 2) Tracks AP location, unit counts, and which monkey (not currently used)
% 3) Asks how modulate whole population of all units is
% 6) Copies relevant figures for eye movment modulated cells to summary directory

clar %clear, clc
task = 'ListSQ';
Fs = 1000;%Hz
min_num_fix = 250; %at least 250 fixatoins with a certain duration to analyze for time limited to fixation duration
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)

summary_directory = 'C:\Users\seth.koenig\Desktop\Significant Units\Eye Movement\'; %where to store significant units
if ~isdir(summary_directory)
    mkdir(summary_directory)
end

%---Misc. Parameters (Mostly Eye Movement Modulated Cells)---%
monkey_all_unit_count = zeros(2,2);%row 1 modulated row 2 non-modulated, column by monkey
all_multi_unit_count = zeros(1,2); %all_multi_unit_countunits
all_eye_cell_unit_names = {}; %eye modulated cell unit names
all_eye_cell_monkeys = []; %1s and 2s
eye_cell_AP_location = []; %AP location of recorded place cell
eye_cell_place_cell_status = [];%whether unit is place cell or not, 1 for place, 0 for non-place
eye_cell_direction_cell_status = []; %whether unit is directionally modulated or not, 1 for direction, 0 for non-direction
eye_cell_amplitude_cell_status = [];%whether unit is amplitude modulated or not, 1 for amplitude, 0 for non-amplitude

%---Fixation Algined Firing Rate Curves---%
avg_fixation_firing_rates = []; %significant firing rates
all_firing_rates = []; %significant + non-significant
avg_fixation_firing_rates_limited =[];%for significant ones limted to duration of 1 preceding saccade+ 1 fixation
all_not_normalized = []; %significant firing rate not-normalized
all_mean_not_normalized = [];%significant firing rate but only mean subtracted
all_peaks = [];%significant
all_peaks2 = [];%not significant

monkeys = {'Vivian','Tobii'};
figure_dir = {};
for monk = 2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for sess data---%%%
    %only need to run when somethings changed or sesss have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir{1} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_session_Unit_Data_Vivian.mat'])
        
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
        load([data_dir 'Across_session_Unit_Data_Tobii.mat'])
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
        load([data_dir task_file(1:end-11) '-Freq_Mod_Analysis.mat'])
        
        for unit = 1:num_units
            if 1
            end
        end
    end
end

%%
iei = [];
for t = 1:281
    dat = seq_eye_times{2}(t,:);
    dat = find(dat == 1);
    iei = [iei diff(dat(1:2:end))];
end