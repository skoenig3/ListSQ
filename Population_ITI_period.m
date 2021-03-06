% Population Visual Response
% written Seth Konig 9/1/16

clar
task = 'ListSQ';
Fs = 1000;%Hz
smval = 60;
imageX = 800;
imageY = 600;


ITI_firing = [];
cell_count = zeros(2,2);
sig_count = zeros(1,2);
all_unit_names = {};
FWHM = [];
spatialness = [];
for monkey = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for session = 1:length(session_data)
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{session},task);
        if isempty(task_file)
            continue
        end
        
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials','cfg',...
            'hdr','data','fixationstats','whole_session_mean_firing_rate');
        
           
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        
        if exist([data_dir task_file(1:8) '-ListSQ-ITI_results.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-ListSQ-ITI_results.mat']);
        else
            continue
        end
        
        %load spatial analysis data
        disp(task_file(1:8))
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'spatial_info')
        
        
        disp(['Importing ' task_file(1:8)])
        
        for unit = 1:num_units
            if ~isempty(time_locked_firing{unit})
                if (ITI_temporal_info.rate_prctile(unit) > 95) && (ITI_temporal_info.temporalstability_prctile(1,unit) > 95)
                    firing_rate = time_locked_firing{unit};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    
                     if nanmean(firing_rate(1:twin)) > nanmean(firing_rate(twin:end))
                         firing_rate = firing_rate-nanmean(firing_rate(twin:end));%maintains that first part is higher
                     else
                          firing_rate = firing_rate-nanmean(firing_rate(1:twin));
                     end
                     firing_rate = firing_rate/max(firing_rate);

%                       firing_rate = firing_rate-nanmean(firing_rate(1:twin));
% %                     firing_rate = firing_rate-whole_session_mean_firing_rate(unit);%remove average firing rate sort of like baseline
%                     if max(firing_rate) > abs(min(firing_rate)) %normalize to max
%                         firing_rate = firing_rate/max(abs(firing_rate));
%                     else%normalize to min %could have some neurons that only show supression
%                         firing_rate = firing_rate/min(firing_rate);
%                     end
                    ITI_firing = [ITI_firing ; firing_rate];
                    
                                        
                    %is unit spatially modulated
                    if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                            && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                        spatialness = [spatialness 1]; %place cell
                    else
                        spatialness = [spatialness 0]; %non-place cell
                    end

                    
                    [mx,max_tuning_ind] = max(firing_rate);
                    large_ind = find(firing_rate > 0.5);
                    gaps = findgaps(large_ind);
                    window = [];
                    for g =1:size(gaps,1)
                        gp = gaps(g,:);
                        gp(gp == 0) = [];
                        if any(gp == max_tuning_ind)
                            window = gp;
                            break
                        end
                    end
                    
                    FWHM = [FWHM [median(window); length(window)]];
                    
                    cell_count(1,monkey) = cell_count(monkey)+1;
                    sig_count(monkey) =sig_count(1,monkey)+1;
                    all_unit_names = [all_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}];
                else
                    cell_count(2,monkey) = cell_count(2,monkey)+1;
                end
            end
        end
    end
end

%% All Cells Aligned to Max
vals = ITI_firing(:);
vals = vals(vals < 0);

figure

[m,i] = max(ITI_firing,[],2);
[mm,ii] = sort(i);
imagesc([-twin:ITI_dur],[1:size(ITI_firing,1)],ITI_firing(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(ITI_firing,1)],'w--');
hold off
xlabel('Time from ITI start  (ms)')
ylabel('Neuron #')
title('ITI Period')
caxis([-std(vals) 1])

names = all_unit_names(ii);

figure
plot(-twin:999,nanmean(ITI_firing));
xlabel('Time From ITI Start (ms)')
ylabel('Normalized Firing Rate')
title('Normalized Population Average Response')