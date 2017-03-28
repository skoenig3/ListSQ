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

%---Summary for Spike x Eye Cross Correlation---%
corr_eye_cell_status = [];
seq_corr_eye_cell_status =[];
list_cross_correlation_peak_freq = cell(1,347);
list_cross_correlation_median_freq = cell(1,347);
list_cross_correlation_all_freq = cell(1,347);
seq_cross_correlation_peak_freq = cell(1,347);
seq_cross_correlation_median_freq = cell(1,347);
seq_cross_correlation_all_freq = cell(1,347);

%---Summary for Eye Auto Correlation---%
list_eye_auto_correlation_freq = cell(1,347);
seq_eye_auto_correlation_freq = cell(1,347);

%---Summary for spike Auto Correlation---%
list_spike_auto_correlation_freq = cell(1,347);
seq_spike_auto_correlation_freq = cell(1,347);

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
            
            %---Determine Modulation Status based on Autocorrelation---%
            %is unit modulated by eye movements during list task
            if isnan(list_spike_eye_xc_observed_max_xc(unit))
                corr_eye_cell_status = [corr_eye_cell_status NaN];
            else
                if 300*sum(list_spike_eye_xc_observed_max_xc(unit) >...
                        list_spike_eye_xc_shuffled_max_xc{unit})...
                        /length(list_spike_eye_xc_shuffled_max_xc{unit}) > 95
                    corr_eye_cell_status =[corr_eye_cell_status 1];
                else
                    corr_eye_cell_status =[corr_eye_cell_status 0];
                end
            end
            
            %is unit modulate by eye movements during seq task
            if isnan(seq_spike_eye_xc_observed_max_xc(unit))
                seq_corr_eye_cell_status = [seq_corr_eye_cell_status NaN];
            else
                if 300*sum(seq_spike_eye_xc_observed_max_xc(unit) >...
                        seq_spike_eye_xc_shuffled_max_xc{unit})...
                        /length(seq_spike_eye_xc_shuffled_max_xc{unit}) > 95
                    seq_corr_eye_cell_status =[seq_corr_eye_cell_status 1];
                else
                    seq_corr_eye_cell_status =[seq_corr_eye_cell_status 0];
                end
            end
            
            
            count = count+1;
            %---Determine Spike x Eye Cross Correlation Signficnant Frequencies---%
            %for list
            if ~isempty(list_spike_eye_xc_observed_FFT{unit})
                sig_times = find(list_spike_eye_xc_observed_FFT{unit} > ...
                    list_spike_eye_xc_shuffled_FFT_99{unit});
                sig_gaps = findgaps(sig_times);
                rmv = [];
                median_freq = [];
                peak_freq = [];
                for g = 1:size(sig_gaps,1)
                    gp = sig_gaps(g,:);
                    gp(gp == 0) = [];
                    if length(gp) < 2*Hz1+1
                        rmv = [rmv g];
                    else
                        median_freq = [median_freq median(gp)];
                        [~,mxi]=max(list_spike_eye_xc_observed_FFT{unit}(gp));
                        peak_freq = [peak_freq mxi+gp(1)-1];
                    end
                end
                sig_gaps(rmv,:) = [];
                sig_gaps(sig_gaps == 0) = [];%remove zeros
                
                list_cross_correlation_peak_freq{count} = peak_freq;
                list_cross_correlation_median_freq{count} = median_freq;
                list_cross_correlation_all_freq = sig_gaps;
            end
            
            %for seq
            if ~isempty(seq_spike_eye_xc_observed_FFT{unit})
                sig_times = find(seq_spike_eye_xc_observed_FFT{unit} > ...
                    seq_spike_eye_xc_shuffled_FFT_99{unit});
                sig_gaps = findgaps(sig_times);
                rmv = [];
                median_freq = [];
                peak_freq = [];
                for g = 1:size(sig_gaps,1)
                    gp = sig_gaps(g,:);
                    gp(gp == 0) = [];
                    if length(gp) < 2*Hz1+1
                        rmv = [rmv g];
                    else
                        median_freq = [median_freq median(gp)];
                        [~,mxi]=max(seq_spike_eye_xc_observed_FFT{unit}(gp));
                        peak_freq = [peak_freq mxi+gp(1)-1];
                    end
                end
                sig_gaps(rmv,:) = [];
                sig_gaps(sig_gaps == 0) = [];%remove zeros
                
                seq_cross_correlation_peak_freq{count} = peak_freq;
                seq_cross_correlation_median_freq{count} = median_freq;
                seq_cross_correlation_all_freq = sig_gaps;
            end
            
            %--Determine Sigifnicant Frequency Modulation of Eye Movements---%
            %for list
            if ~isempty(list_eye_ac_observed_FFT{unit})
                sig_freq = NaN(2,length(list_eye_ac_freq_of_interest{unit}));
                for freq = 1:length(list_eye_ac_freq_of_interest{unit})
                    freq_ind = find(f2 == list_eye_ac_freq_of_interest{unit}(freq));
                    if sum(list_eye_ac_observed_FFT{unit}(freq_ind-Hz1:freq_ind+Hz1)) > ...
                            sum(list_eye_ac_shuffled_FFT_99{freq,unit}(freq_ind-Hz1:freq_ind+Hz1))
                        sig_freq(1,freq) = 1;
                        sig_freq(2,freq) = freq_ind;
                    else
                        sig_freq(1,freq) = 0;
                        sig_freq(2,freq) = freq_ind;
                    end
                end
                list_eye_auto_correlation_freq{count} = sig_freq;
            end

            %for seq
            if ~isempty(seq_eye_ac_observed_FFT{unit})
                sig_freq = NaN(2,length(seq_eye_ac_freq_of_interest{unit}));
                for freq = 1:length(seq_eye_ac_freq_of_interest{unit})
                    freq_ind = find(f2 == seq_eye_ac_freq_of_interest{unit}(freq));
                    if sum(seq_eye_ac_observed_FFT{unit}(freq_ind-Hz1:freq_ind+Hz1)) > ...
                            sum(seq_eye_ac_shuffled_FFT_99{freq,unit}(freq_ind-Hz1:freq_ind+Hz1))
                        sig_freq(1,freq) = 1;
                        sig_freq(2,freq) = freq_ind;
                    else
                        sig_freq(1,freq) = 0;
                        sig_freq(2,freq) = freq_ind;
                    end
                end
                seq_eye_auto_correlation_freq{count} = sig_freq;
            end
            
            %--Determine Sigifnicant Frequency Modulation of Spike Trains---%
            %for list
            if ~isempty(list_spike_ac_observed_FFT{unit})
                sig_freq = NaN(2,length(list_spike_ac_freq_of_interest{unit}));
                for freq = 1:length(list_spike_ac_freq_of_interest{unit})
                    freq_ind = find(f2 == list_spike_ac_freq_of_interest{unit}(freq));
                    if sum(list_spike_ac_observed_FFT{unit}(freq_ind-Hz1:freq_ind+Hz1)) > ...
                            sum(list_spike_ac_shuffled_FFT_99{freq,unit}(freq_ind-Hz1:freq_ind+Hz1))
                        sig_freq(1,freq) = 1;
                        sig_freq(2,freq) = freq_ind;
                    else
                        sig_freq(1,freq) = 0;
                        sig_freq(2,freq) = freq_ind;
                    end
                end
                list_spike_auto_correlation_freq{count} = sig_freq;
            end
            
            %for seq
            if ~isempty(seq_spike_ac_observed_FFT{unit})
                sig_freq = NaN(2,length(seq_spike_ac_freq_of_interest{unit}));
                for freq = 1:length(seq_spike_ac_freq_of_interest{unit})
                    freq_ind = find(f2 == seq_spike_ac_freq_of_interest{unit}(freq));
                    if sum(seq_spike_ac_observed_FFT{unit}(freq_ind-Hz1:freq_ind+Hz1)) > ...
                            sum(seq_spike_ac_shuffled_FFT_99{freq,unit}(freq_ind-Hz1:freq_ind+Hz1))
                        sig_freq(1,freq) = 1;
                        sig_freq(2,freq) = freq_ind;
                    else
                        sig_freq(1,freq) = 0;
                        sig_freq(2,freq) = freq_ind;
                    end
                end
                seq_spike_auto_correlation_freq{count} = sig_freq;
            end
        end
    end
end
%%

%---Eye Autocorrelation Frequency Content---%
list_sig_eye_ac_count = zeros(size(list_eye_auto_correlation_freq{1},2),length(f2));
list_sig_eye_ac_unit_count = 0;

seq_sig_eye_ac_count = zeros(size(seq_eye_auto_correlation_freq{1},2),length(f2));
seq_sig_eye_ac_unit_count = 0;

for unit = 1:length(list_eye_auto_correlation_freq);
    if ~isempty(list_eye_auto_correlation_freq{unit})
        dat = list_eye_auto_correlation_freq{unit};
        for freq = 1:size(dat,2)
            if dat(1,freq) == 1;
                list_sig_eye_ac_count(freq,dat(2,freq)) = list_sig_eye_ac_count(freq,dat(2,freq))+1;
            end
        end
        list_sig_eye_ac_unit_count = list_sig_eye_ac_unit_count+1;
    end
    
    if ~isempty(seq_eye_auto_correlation_freq{unit})
        dat = seq_eye_auto_correlation_freq{unit};
        for freq = 1:size(dat,2)
            if dat(1,freq) == 1;
                seq_sig_eye_ac_count(freq,dat(2,freq)) = seq_sig_eye_ac_count(freq,dat(2,freq))+1;
            end
        end
        seq_sig_eye_ac_unit_count = seq_sig_eye_ac_unit_count+1;
    end
end

%---List Task Eye Movement Autocorrelation---%
figure
subplot(2,2,1)
plot(f2,sum(list_sig_eye_ac_count./list_sig_eye_ac_unit_count))
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('All Significant Frequencies')

subplot(2,2,2)
plot(f2,list_sig_eye_ac_count(1,:)./list_sig_eye_ac_unit_count)
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('1st Significant Frequency')

subplot(2,2,3)
plot(f2,list_sig_eye_ac_count(2,:)./list_sig_eye_ac_unit_count)
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('2nd Significant Frequency')

subplot(2,2,4)
plot(f2,list_sig_eye_ac_count(3,:)./list_sig_eye_ac_unit_count)
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('3rd Significant Frequency')
subtitle('List Eye')

%---Sequence Task Eye Movement Autocorrelation---%
figure
subplot(2,2,1)
plot(f2,sum(seq_sig_eye_ac_count./seq_sig_eye_ac_unit_count))
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('All Significant Frequencies')

subplot(2,2,2)
plot(f2,seq_sig_eye_ac_count(1,:)./seq_sig_eye_ac_unit_count)
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('1st Significant Frequency')

subplot(2,2,3)
plot(f2,seq_sig_eye_ac_count(2,:)./seq_sig_eye_ac_unit_count)
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('2nd Significant Frequency')

subplot(2,2,4)
plot(f2,seq_sig_eye_ac_count(3,:)./seq_sig_eye_ac_unit_count)
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('3rd Significant Frequency')
subtitle('Sequence Eye')
%%
%---spike Autocorrelation Frequency Content---%
list_sig_spike_ac_count = zeros(size(list_spike_auto_correlation_freq{1},2),length(f2));
list_sig_spike_ac_unit_count = 0;

seq_sig_spike_ac_count =  zeros(size(seq_spike_auto_correlation_freq{1},2),length(f2));
seq_sig_spike_ac_unit_count = 0;

for unit = 1:length(list_spike_auto_correlation_freq);
    if ~isempty(list_spike_auto_correlation_freq{unit})
        dat = list_spike_auto_correlation_freq{unit};
        for freq = 1:size(dat,2)
            if dat(1,freq) == 1;
                list_sig_spike_ac_count(freq,dat(2,freq)) = list_sig_spike_ac_count(freq,dat(2,freq))+1;
            end
        end
        list_sig_spike_ac_unit_count = list_sig_spike_ac_unit_count+1;
    end
    
    if ~isempty(seq_spike_auto_correlation_freq{unit})
        dat = seq_spike_auto_correlation_freq{unit};
        for freq = 1:size(dat,2)
            if dat(1,freq) == 1;
                seq_sig_spike_ac_count(freq,dat(2,freq)) = seq_sig_spike_ac_count(freq,dat(2,freq))+1;
            end
        end
        seq_sig_spike_ac_unit_count = seq_sig_spike_ac_unit_count+1;
    end
end
%%
%---List Task spike Movement Autocorrelation---%
figure
subplot(2,3,1)
dofill(f2,sum(list_sig_spike_ac_count./list_sig_spike_ac_unit_count),'black',1,5); %all trials
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('All Significant Frequencies')

subplot(2,3,2)
dofill(f2,list_sig_spike_ac_count(1,:)./list_sig_spike_ac_unit_count,'black',1,5); %all trials
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('1st Significant Frequency')

subplot(2,3,3)
dofill(f2,list_sig_spike_ac_count(2,:)./list_sig_spike_ac_unit_count,'black',1,5); %all trials
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('2nd Significant Frequency')

subplot(2,3,4)
dofill(f2,list_sig_spike_ac_count(3,:)./list_sig_spike_ac_unit_count,'black',1,5); %all trials
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('3rd Significant Frequency')
subtitle('List Spike Trains')

subplot(2,3,5)
dofill(f2,list_sig_spike_ac_count(4,:)./list_sig_spike_ac_unit_count,'black',1,5); %all trials
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('4rd Significant Frequency')

subplot(2,3,6)
dofill(f2,list_sig_spike_ac_count(5,:)./list_sig_spike_ac_unit_count,'black',1,20); %all trials
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('5rd Significant Frequency')
xlim([100 450])

subtitle('List Spike Trains')
%%
%---Sequence Task spike Movement Autocorrelation---%
figure
subplot(2,3,1)
dofill(f2,sum(seq_sig_spike_ac_count./seq_sig_spike_ac_unit_count),'black',1,5); %all trials
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('All Significant Frequencies')

subplot(2,3,2)
dofill(f2,seq_sig_spike_ac_count(1,:)./seq_sig_spike_ac_unit_count,'black',1,5); %all trials
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('1st Significant Frequency')

subplot(2,3,3)
dofill(f2,seq_sig_spike_ac_count(2,:)./seq_sig_spike_ac_unit_count,'black',1,5); %all trials
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('2nd Significant Frequency')

subplot(2,3,4)
dofill(f2,seq_sig_spike_ac_count(3,:)./seq_sig_spike_ac_unit_count,'black',1,5); %all trials
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('3rd Significant Frequency')
subtitle('seq Spike Trains')

subplot(2,3,5)
dofill(f2,seq_sig_spike_ac_count(4,:)./seq_sig_spike_ac_unit_count,'black',1,5); %all trials
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('4rd Significant Frequency')

subplot(2,3,6)
dofill(f2,seq_sig_spike_ac_count(5,:)./seq_sig_spike_ac_unit_count,'black',1,20); %all trials
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
title('5rd Significant Frequency')
xlim([100 450])

subtitle('Sequence Spike Trains')

%%
%---Sequence Task and List Spike Autocorrelation---%
figure
hold on
plot(f2,300*sum(seq_sig_spike_ac_count./seq_sig_spike_ac_unit_count),'k')
plot(f2,300*sum(list_sig_spike_ac_count./list_sig_spike_ac_unit_count),'r')
dofill(f2,sum(seq_sig_spike_ac_count./seq_sig_spike_ac_unit_count),'black',1,5); %all trials
dofill(f2,sum(list_sig_spike_ac_count./list_sig_spike_ac_unit_count),'red',1,5); %all trials
hold off
xlim([0 30])
xlabel('Frequency (Hz)')
ylabel('% of Units')
legend('Sequence','List')
title('All Significant Frequencies')

figure
hold on
plot(f2,300*sum(seq_sig_spike_ac_count./seq_sig_spike_ac_unit_count),'k')
plot(f2,300*sum(list_sig_spike_ac_count./list_sig_spike_ac_unit_count),'r')
dofill(f2,sum(seq_sig_spike_ac_count./seq_sig_spike_ac_unit_count),'black',1,5); %all trials
dofill(f2,sum(list_sig_spike_ac_count./list_sig_spike_ac_unit_count),'red',1,5); %all trials
hold off
xlim([30 100])
xlabel('Frequency (Hz)')
ylabel('% of Units')
legend('Sequence','List')
title('All Significant Frequencies')


figure
hold on
plot(f2,300*sum(seq_sig_spike_ac_count./seq_sig_spike_ac_unit_count),'k')
plot(f2,300*sum(list_sig_spike_ac_count./list_sig_spike_ac_unit_count),'r')
dofill(f2,sum(seq_sig_spike_ac_count./seq_sig_spike_ac_unit_count),'black',1,5); %all trials
dofill(f2,sum(list_sig_spike_ac_count./list_sig_spike_ac_unit_count),'red',1,5); %all trials
hold off
xlim([100 450])
xlabel('Frequency (Hz)')
ylabel('% of Units')
legend('Sequence','List')
title('All Significant Frequencies')

