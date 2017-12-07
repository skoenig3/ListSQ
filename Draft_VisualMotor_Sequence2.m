% Population Saccadic Modulation
% written Seth Konig 8/12/16
%
% code imports data from List_Saccade_AnalysisV2 results and looks at how
% modulated the whole population is by saccades. Analysis includes...
% 1) Average saccade/fixation-locked firing rate
% 2) AP Differences
%   a) % of neurons significantly modulated by eye movements
%   b) degree of modulation e.g. % change in firing rate?
% 3) Relative timing of modulation i.e. lag
%   a) on average
%   b) distribution
clar
task = 'ListSQ';
Fs = 1000;%Hz
smval = 60;
imageX = 800;
imageY = 600;
numshuffs = 1000;

vizmotor_status = [];
t = -500:499;
avg_fixation_firing_rates = cell(1,2);
peak_fixation_firing_rates = cell(1,2);

spatialness = [];
visual_motor_visuomotor = [];
unit_names = [];
sig_event = [];
sig_fix = [];
bin_size = 100; %trials per bin
set(0,'DefaultFigureVisible','OFF');
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Population Figures\Place Cell Eye Movements\';
for monkey = 2:-1:1
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
        
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials','cfg','hdr','data','fixationstats');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        
        if exist([data_dir task_file(1:8) '-Fixation_locked_Sequence_results.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-Fixation_locked_Sequence_results.mat']);
        else
            continue
        end
        
        disp(['Importing ' task_file(1:8)])
        
        %load spatial analysis data
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'spatial_info')
        
        for unit = 1:num_units
            if ~isempty(fixation_locked_firing{unit})
                %response to all items
                if all_items_info.temporalstability_prctile(unit) > 95 &&  all_items_info.rate_prctile(unit) > 95
                    item_appearance = 1;
                else
                    item_appearance = 0;
                end
                if  fixation_info.rate_prctile(unit) > 95 && fixation_info.temporalstability_prctile(unit) > 95
                    fixation_on_item = 1;
                else
                    fixation_on_item = 0;
                end
                
                sig_event = [sig_event item_appearance];
                sig_fix = [sig_fix fixation_on_item];
                
                %is unit spatially modulated
                if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                        && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                    spatialness = [spatialness 1]; %place cell
                else
                    spatialness = [spatialness 0]; %non-place cell
                end
                
                if item_appearance == 1 || fixation_on_item == 1
                    %%
                    fixation_firing = fixation_locked_firing{unit}(:,:);
                    event_firing = event_locked_firing{unit}(:,:);
                    
                    rts = reaction_time{unit};
                    %                     [~,sorted_rts] = sort(rts);
                    %                     sorted_fixation_aligned = fixation_firing(sorted_rts,:);
                    %                     sorted_item_aligned = event_firing(sorted_rts,:);
                    [~,all_fix_curves] =  nandens(fixation_firing,smval,'gauss',Fs,'nanflt');
                    [~,all_event_curves] = nandens(event_firing,smval,'gauss',Fs,'nanflt');
                    which_seq = sequence_nums{unit};
                    which_item = item_nums{unit};
                    
                    %%
                    figure
                    fix_gain_modulation = NaN(2,4);
                    item_gain_modulation = NaN(2,4);
                    for seq = 1:2
                        for item = 1:4
                            subplot(1,2,1)
                            hold on
                            mean_firing = nanmean(all_fix_curves(which_seq == seq & which_item == item,:));
                            plot(mean_firing)
                            fix_gain_modulation(seq,item) = max(mean_firing)-min(mean_firing);
                            
                            subplot(1,2,2)
                            hold on
                            mean_firing = nanmean(all_event_curves(which_seq == seq & which_item == item,:));
                            plot(mean_firing)
                            item_gain_modulation(seq,item) = max(mean_firing)-min(mean_firing);
                        end
                    end
                    
                    fix_gain_modulation = fix_gain_modulation/sum(fix_gain_modulation(:));
                    item_gain_modulation = item_gain_modulation/sum(item_gain_modulation(:));
                    %%
                    
                    observed_fix_corr = NaN(2,4);
                    observed_event_corr = NaN(2,4);
                    for seq = 1:2
                        for item = 1:4
                            indeces = (which_seq == seq & which_item == item);
                            count_pct = round(sum(indeces)*.33);
                            these_rt = rts(indeces);
                            [~,these_sorted_rts] = sort(these_rt);
                            
                            these_fix = all_fix_curves(indeces,:);
                            these_fix = these_fix(these_sorted_rts,:);
                            these_event = all_event_curves(indeces,:);
                            these_event = these_event(these_sorted_rts,:);
                            
                            fix_c1 = nanmean(these_fix(1:count_pct,:));
                            fix_c2 = nanmean(these_fix(size(these_fix,1)-count_pct+1:end,:));
                            observed_fix_corr(seq,item) = corr(fix_c1(:),fix_c2(:),'row','pairwise','type','Spearman');
                            
                            event_c1 = nanmean(these_event(1:count_pct,:));
                            event_c2 = nanmean(these_event(size(these_event,1)-count_pct+1:end,:));
                            observed_event_corr(seq,item) = corr(event_c1(:),event_c2(:),'row','pairwise','type','Spearman');
                        end
                    end
                    
                    observed_fix_corr = sum(sum(observed_fix_corr.*fix_gain_modulation));
                    observed_event_corr = sum(sum(observed_event_corr.*item_gain_modulation));
                    %                     %%
                    %                     count_40pct = round(size(all_fix_curves,1)*.33);
                    %                     fix_c1 = nanmean(all_fix_curves(1:count_40pct,:));
                    %                     fix_c2 = nanmean(all_fix_curves(size(all_fix_curves,1)-count_40pct+1:end,:));
                    %                     observed_fix_corr = corr(fix_c1(:),fix_c2(:),'row','pairwise','type','Spearman');
                    %
                    %                     event_c1 = nanmean(all_event_curves(1:count_40pct,:));
                    %                     event_c2 = nanmean(all_event_curves(size(all_event_curves,1)-count_40pct+1:end,:));
                    %                     observed_event_corr = corr(event_c1(:),event_c2(:),'row','pairwise','type','Spearman');
                    %
                    %%
                    observed_corr_diff = observed_fix_corr-observed_event_corr;
                    
                    if isnan(observed_corr_diff)
                        titel_str = 'Unknown';
                        vizmotor_status = [vizmotor_status NaN];
                    else
                        
                        shuffled_diff = NaN(1,numshuffs);
                        
                        %                     parfor shuff = 1:numshuffs;
                        %                         ind = randperm(trial_count);
                        %                         shuff_fix = all_fix_curves(ind,:);
                        %                         shuff_event = all_event_curves(ind,:);
                        %
                        %                         event_c1 = nanmean(shuff_event(1:count_40pct,:));
                        %                         event_c2 = nanmean(shuff_event(size(all_event_curves,1)-count_40pct+1:end,:));
                        %                         shuff_event_corr = corr(event_c1(:),event_c2(:),'row','pairwise','type','Spearman');
                        %
                        %                         fix_c1 = nanmean(shuff_fix(1:count_40pct,:));
                        %                         fix_c2 = nanmean(shuff_fix(size(all_fix_curves,1)-count_40pct+1:end,:));
                        %                         shuff_fix_corr = corr(fix_c1(:),fix_c2(:),'row','pairwise','type','Spearman');
                        %
                        %
                        %                         shuffled_diff(shuff) = shuff_fix_corr-shuff_event_corr;
                        %                     end
                        
                        parfor shuff = 1:numshuffs;
                            fix_corr = NaN(2,4);
                            event_corr = NaN(2,4);
                            for seq = 1:2
                                for item = 1:4
                                    indeces = (which_seq == seq & which_item == item);
                                    count_pct = round(sum(indeces)*.33);
                                    these_rt = rts(indeces);
                                    %[~,these_sorted_rts] = sort(these_rt);
                                    these_sorted_rts =  randperm(sum(indeces));
                                    
                                    these_fix = all_fix_curves(indeces,:);
                                    these_fix = these_fix(these_sorted_rts,:);
                                    these_event = all_event_curves(indeces,:);
                                    these_event = these_event(these_sorted_rts,:);
                                    
                                    fix_c1 = nanmean(these_fix(1:count_pct,:));
                                    fix_c2 = nanmean(these_fix(size(these_fix,1)-count_pct+1:end,:));
                                    fix_corr(seq,item) = corr(fix_c1(:),fix_c2(:),'row','pairwise','type','Spearman');
                                    
                                    event_c1 = nanmean(these_event(1:count_pct,:));
                                    event_c2 = nanmean(these_event(size(these_event,1)-count_pct+1:end,:));
                                    event_corr(seq,item) = corr(event_c1(:),event_c2(:),'row','pairwise','type','Spearman');
                                end
                            end
                            
                            fix_corr = sum(sum(fix_corr.*fix_gain_modulation));
                            event_corr = sum(sum(event_corr.*item_gain_modulation));
                            shuffled_diff(shuff) = fix_corr-event_corr;
                        end
                        
                        if 100*sum(observed_corr_diff > shuffled_diff)/numshuffs > 97.5 %fixation so motor
                            title_str = 'Motor';
                            vizmotor_status = [vizmotor_status 1];
                        elseif 100*sum(observed_corr_diff > shuffled_diff)/numshuffs < 2.5 %event so visual
                            title_str = 'Visual';
                            vizmotor_status = [vizmotor_status 2];
                        else
                            title_str =  'VisuoMotor';
                            vizmotor_status = [vizmotor_status 3];
                        end
                    end
                    %%
                    figure
                    twin1 = 500;
                    subplot(2,2,1)
                    [trial,time] = find(event_firing > 0);
                    if ~isempty(trial)
                        plot(time-twin1,trial,'.k')
                        xlim([-500 500])
                        ylim([0 max(trial)+1]);
                    end
                    ylabel('Trial Ranked by RT')
                    xlabel('Time from Cross On (ms)')
                    box off
                    title(['\rho_{rt} = ' num2str(observed_event_corr,2)])
                    
                    subplot(2,2,2)
                    [trial,time] = find(fixation_firing > 0);
                    if ~isempty(trial)
                        plot(time-twin1,trial,'.k')
                        xlim([-500 500])
                        ylim([0 max(trial)+1]);
                    end
                    ylabel('Trial Ranked by RT')
                    xlabel('Time from Fixation Start (ms)')
                    box off
                    title(['\rho_{rt} = ' num2str(observed_fix_corr,2)])
                    
                    subplot(2,2,3)
                    hold on
                    dofill2([1:size(event_firing,2)]-500,event_firing,'black',1,smval); %smoothed cross on aligned curve
                    dofill2([1:size(fixation_firing,2)]-500,fixation_firing,'green',1,smval); %smoothed cross on aligned curve
                    hold off
                    xlabel('Event Onset (ms)')
                    ylabel('Firing Rate (Hz)')
                    xlim([-500 500])
                    legend('Item On','Fixation Start')
                    
                    if ~isnan(observed_corr_diff)
                        subplot(2,2,4)
                        hist(shuffled_diff,25)
                        hold on
                        yl = ylim;
                        plot([observed_corr_diff observed_corr_diff],[0 yl(2)],'r--')
                        hold off
                        box off
                        title(['\Delta \rho_{rt} = ' num2str(observed_fix_corr-observed_event_corr,2)...
                            ' (' num2str(100*sum(observed_corr_diff > shuffled_diff)/numshuffs,3) '%)'])
                    end
                    subtitle(title_str)
                    
                    %%
                    save_and_close_fig('C:\Users\seth.koenig\Desktop\Test_VizMotor\',[task_file(1:8) '_' unit_stats{1,unit}])
                    
                end
            end
        end
    end
end
set(0,'DefaultFigureVisible','ON');