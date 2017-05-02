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
all_unit_names = {}; %eye modulated cell unit names
all_monkeys = []; %1s and 2s
place_cell_status = [];%whether unit is place cell or not, 1 for place, 0 for non-place
eye_cell_status = []; %whether unit modulated by eye movements
seq_eye_cell_status = [];%whether unit  is modulated by eye movements in sequence task

%---Summary for Eye Auto Correlation---%
list_eye_auto_correlation_freq = cell(2,347);
seq_eye_auto_correlation_freq = cell(2,347);

%---Summary for spike Auto Correlation---%
list_spike_auto_correlation_freq = cell(2,347);
seq_spike_auto_correlation_freq = cell(2,347);
whole_spike_auto_correlation_freq = cell(2,347);
sig_peak_frequency = NaN(347,5);

%---Summary for Cross Corelation Spike x Eye---%
list_spike_eye_cross_correlation_freq = cell(2,347);
seq_spike_eye_cross_correlation_freq = cell(2,347);

monkeys = {'Vivian','Tobii'};
figure_dir = {};
count = 0;
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
        load([data_dir task_file(1:8) '-Freq_Mod_Analysis.mat'])
        load([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat'],'temporal_info');
        load([data_dir task_file(1:8) '-Fixation_locked_Sequence_results.mat'],'fixation_info');
        
        Hz1 = round(smFreqVal/mode(diff(f2))/2); %window around peak frequency to test
        for unit = 1:num_units
            if all(isnan(valid_trials(:,unit))) %not a usable unit
                continue
            end
           
            %---Determine Spike x Eye Cross Correlation Signficnant Frequencies---%
            %for list
            if ~isempty(whole_spike_ac_significant_freq{unit})
                count = count+1;
                
                            
                %---Misc Info---%
                all_monkeys = [all_monkeys monk];
                all_unit_names = [all_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}];
                
                %---Determine Modulation Status Using Original Analysis---%
                %is unit spatially modulated during free viewing?
                if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                        && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                    place_cell_status = [place_cell_status 1]; %place cell
                else
                    place_cell_status = [place_cell_status 0]; %non-place cell
                end
                
                %is unit modulated by eye movemets in list task
                if (temporal_info.fixation.shuffled_temporalstability_prctile(1,unit) > 95) ... %significant stability
                        && (temporal_info.fixation.shuffled_rate_prctile(unit) > 95) % %skagg 95%+
                    eye_cell_status = [eye_cell_status 1];
                else
                    eye_cell_status = [eye_cell_status 0];
                end
                
                
                %is unit modulated by eye movements during the sequence task
                if fixation_info.rate_prctile(unit) > 95 && fixation_info.temporalstability_prctile(unit) > 95
                    seq_eye_cell_status=[seq_eye_cell_status 1];
                else
                    seq_eye_cell_status=[seq_eye_cell_status 0];
                end
                

                
                %---Whole Session Data---%
                whole_spike_auto_correlation_freq{1,count} = whole_spike_ac_observed_ac{unit};%autocorrelation
                whole_spike_auto_correlation_freq{2,count} = whole_spike_ac_observed_FFT{unit};%frequency content of autocorrelation
                sig_peak_frequency(count,whole_spike_ac_significant_freq{unit}(2,:) == 1) = ...
                    whole_spike_ac_significant_freq{unit}(1,whole_spike_ac_significant_freq{unit}(2,:) == 1); %which peaks if they existed were signifiant
                
                %---List and Sequecne Only---%
                if ~isempty(list_spike_ac_observed_ac{unit})
                    list_eye_auto_correlation_freq{1,count} = list_eye_ac_observed_ac{unit};
                    list_eye_auto_correlation_freq{2,count} = list_eye_ac_observed_FFT{unit};
                    
                    list_spike_auto_correlation_freq{1,count} = list_spike_ac_observed_ac{unit};
                    list_spike_auto_correlation_freq{2,count} = list_spike_ac_observed_FFT{unit};
                    
                    list_spike_eye_cross_correlation_freq{1,count} = list_spike_eye_xc_observed_xc{unit};
                    list_spike_eye_cross_correlation_freq{2,count} = list_spike_eye_xc_observed_FFT{unit};
                end
                
                if ~isempty(seq_spike_ac_observed_ac{unit})
                    seq_eye_auto_correlation_freq{1,count} = seq_spike_ac_observed_ac{unit};
                    seq_eye_auto_correlation_freq{2,count} = seq_spike_eye_xc_observed_FFT{unit};
                    
                    seq_spike_auto_correlation_freq{1,count} = seq_spike_ac_observed_ac{unit};
                    seq_spike_auto_correlation_freq{2,count} = seq_spike_ac_observed_FFT{unit};
                    
                    seq_spike_eye_cross_correlation_freq{1,count} = seq_spike_eye_xc_observed_xc{unit};
                    seq_spike_eye_cross_correlation_freq{2,count} = seq_spike_eye_xc_observed_FFT{unit};
                    
                end
            end
        end
    end
end
sig_peak_frequency(count+1:end,:) = [];
%%
spf =sig_peak_frequency(:);
spf(isnan(sig_peak_frequency)) = [];

spflow = spf(spf < 30);
spfgam = spf(spf > 30 & spf < 100);

