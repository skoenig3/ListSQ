clar

task = 'ListSQ';
min_blks = 2;

all_list_eye_FFT = NaN(347,2^15+1);
all_list_spike_FFT =  NaN(347,2^15+1);
all_seq_eye_FFT = NaN(347,2^15+1);
all_seq_spike_FFT =  NaN(347,2^15+1);
unit_ind = zeros(1,2);;
for monkey = 2:-1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        %         listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 155;%155.85 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 135;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %         listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for sess = 1:length(session_data)
        task_file = get_task_data(session_data{sess},'ListsQ');
        if isempty(task_file)
            warning('No file could be found for specificed task. Exiting function...')
            continue
        end
        disp(task_file(1:8))
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat']);
        
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
            = get_task_data(session_data{sess},'ListSQ');
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        clear unit_names
        
        if num_units == 0; %if no units exit function
            disp([task_file(1:8) ': no units could be found. Exiting function...'])
            continue;
        end
        
        valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
        if all(isnan(valid_trials(1,:)))
            disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
            continue %since no valid units skip analysis
        end
        
        try
        load([data_dir task_file(1:8) '-Freq_Mod_Analysis.mat'],...
            'list_eye_ac_observed_FFTs','seq_eye_ac_observed_FFTs',...
            'list_spike_ac_observed_FFTs','seq_spike_ac_observed_FFTs');
        catch
            continue
        end
        
        for unit = 1:length(list_eye_ac_observed_FFTs)
            if ~isempty(list_eye_ac_observed_FFTs{unit})
                unit_ind(1) = unit_ind(1)+1;
                all_list_eye_FFT(unit_ind(1),:) = list_eye_ac_observed_FFTs{unit};
                all_list_spike_FFT(unit_ind(1),:) = list_spike_ac_observed_FFTs{unit};
            end
            if ~isempty(seq_eye_ac_observed_FFTs{unit})
                unit_ind(2) = unit_ind(2)+1;
                all_seq_eye_FFT(unit_ind(2),:)= seq_eye_ac_observed_FFTs{unit};
                all_seq_spike_FFT(unit_ind(2),:)= seq_spike_ac_observed_FFTs{unit};
            end
        end
    end
end

all_list_eye_FFT(unit_ind(1)+1:end,:) = [];
all_list_spike_FFT(unit_ind(1)+1:end,:) = [];
all_seq_eye_FFT(unit_ind(2)+1:end,:) = [];
all_seq_spike_FFT(unit_ind(2)+1:end,:) = [];
%%
Fs = 1000;
NFFT = 2^16;%number of frequency points from 0-Fs/2
f2 = Fs/2*linspace(0,1,NFFT/2+1);%frequency values
%%
avg_seq_eye = mean(all_seq_eye_FFT(:));
avg_eye = mean(all_seq_eye_FFT)/avg_seq_eye;
avg_seq_eye = avg_seq_eye-mean(avg_seq_eye);
%%
norm1 = mean(all_seq_spike_FFT)./avg_eye;
norm1 = norm1/max(norm1);
 
%%
norm2 = mean(all_seq_spike_FFT);
norm2 = norm2/max(norm2);
%%
norm3=mean(all_list_spike_FFT);
norm3 = norm3/max(norm3);
%%
norm4 = mean(all_list_spike_FFT)./avg_eye;
norm4 = norm4/max(norm4);
%%
whitening_filter = avg_seq_eye;
save('AutoCorrWhitenFilter','whitening_filter');