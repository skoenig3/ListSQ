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

t = -500:499;
avg_fixation_firing_rates = cell(1,2);
peak_fixation_firing_rates = cell(1,2);

spatialness = [];
visual_motor_visuomotor = [];
unit_names = [];
sig_event = [];
sig_fix = [];
bin_size = 100; %trials per bin

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
                
                fixation_firing = fixation_locked_firing{unit};
                event_firing = event_locked_firing{unit};
                rts = reaction_time{unit};
                [~,sorted_rts] = sort(rts);
                fixation_firing = fixation_firing(sorted_rts,:);
                event_firing = event_firing(sorted_rts,:);
                
                if item_appearance == 1 || fixation_on_item == 1
                    
                    bins = floor(size(fixation_firing,1)/bin_size);
                    rem = size(fixation_firing,1)-bins*bin_size;
                    %skip the fastest since may be hardest to tell
                    %%
                    N = length(rts);
                    n_trials = floor(N/3);
                    
                    event_curves1 = nandens(event_firing(1:n_trials,:),smval,'gauss',Fs,'nanflt');
                    event_curves2 = nandens(event_firing(N-n_trials:end,:),smval,'gauss',Fs,'nanflt');
                    
                    observed_event_corr = corr(event_curves1(:),event_curves2(:),'row','pairwise','type','Spearman');
                    
                    fix_curves1 = nandens(fixation_firing(1:n_trials,:),smval,'gauss',Fs,'nanflt');
                    fix_curves2 = nandens(fixation_firing(N-n_trials:end,:),smval,'gauss',Fs,'nanflt');
                    
                    observed_fix_corr = corr(fix_curves1(:),fix_curves2(:),'row','pairwise','type','Spearman');
                    
                    if observed_fix_corr > 0 && observed_event_corr < -0.1
                        visuomotor = 1; %motor
                    elseif observed_fix_corr <-0.1 && observed_event_corr > 0
                        visuomotor = 2; %visual
                    else
                        visuomotor = 3;%visual motor
                    end
                    visual_motor_visuomotor = [visual_motor_visuomotor visuomotor];
                    unit_names = [unit_names [{task_file(1:8) '_3_' unit_stats{1,unit}}]];
                    
                    %                      [~,all_fix_curves] =  nandens(event_firing,smval,'gauss',Fs,'nanflt');
                    %                      [~,all_event_curves] = nandens(fixation_firing,smval,'gauss',Fs,'nanflt');
                    %%
                    %                      shuff_event_corrs = NaN(1,1000);
                    %                      shuff_fix_corrs = NaN(1,1000);
                    %                      N2 = size(all_fix_curves,1);
                    %                      parfor shuffs = 1:1000
                    %                          shuffled_ind = randperm(N2);
                    %                          shuffled_fix = all_fix_curves(shuffled_ind,:);
                    %                          shuffled_fix1 = mean(shuffled_fix(1:n_trials,:));
                    %                          shuffled_fix2 = mean(shuffled_fix(N-n_trials:end,:));
                    %                          shuffled_event = all_event_curves(shuffled_ind,:);
                    %                          shuffled_event1 = mean(shuffled_event(1:n_trials,:));
                    %                          shuffled_event2 = mean(shuffled_event(N-n_trials:end,:));
                    %
                    %                          shuff_fix_corrs(shuffs) = corr(shuffled_fix1(:),shuffled_fix2(:),'row','pairwise','type','Spearman');
                    %                          shuff_event_corrs(shuffs) = corr(shuffled_event1(:),shuffled_event2(:),'row','pairwise','type','Spearman');
                    %                      end
                    %                      fix_prctile = sum(observed_fix_corr > shuff_fix_corrs);
                    %                      event_prctile = sum(observed_event_corr > shuff_event_corrs);
                    %%
                    %                     event_curves = NaN(bins,size(event_firing,2));
                    %                     for b = 1:bins;
                    %                         ind = ((b-1)*bin_size:bin_size*b)+rem;
                    %                         event_curves(b,:) = nandens(event_firing(ind,:),smval,'gauss',Fs,'nanflt');%fixation aligned firing rate
                    %                     end
                    %
                    %                     fix_curves = NaN(bins,size(fixation_firing,2));
                    %                     for b = 1:bins;
                    %                         ind = ((b-1)*bin_size:bin_size*b)+rem;
                    %                         fix_curves(b,:) = nandens(fixation_firing(ind,:),smval,'gauss',Fs,'nanflt');%fixation aligned firing rate
                    %                     end
                    
                    
                    %
                    %                     [U,S,V] = pca([event_curves; fix_curves],2);
                    %                     PC_event = U(1:bins,:);
                    %                     PC_fix =  U(bins+1:end,:);
                    %
                    %                     %find the component with the greatest covariont between
                    %                     %the 2 types of aligend acitivity. For most neurons this
                    %                     %will be firing rate
                    %                     cv = NaN(2,2);
                    %                     for p = 1:2
                    %                         [r,pval] = corrcoef(PC_event(:,p),PC_fix(:,p));
                    %                         cv(1,p) = r(2);
                    %                         cv(2,p) = pval(2);
                    %                     end
                    %
                    %                     event_rrt = NaN(2,2);
                    %                     fix_rrt = NaN(2,2);
                    %                     for p = 1:2
                    %                         [r,pval] = corrcoef(1:bins,PC_event(:,p));
                    %                         event_rrt(1,p) = r(2);
                    %                         event_rrt(2,p) = pval(2);
                    %
                    %                         [r,pval] = corrcoef(1:bins,PC_fix(:,p));
                    %                         fix_rrt(1,p) = r(2);
                    %                         fix_rrt(2,p) = pval(2);
                    %                     end
                    %
                    %                     [~,largest_cov] = max(cv(1,:));
                    %                     if cv(2,largest_cov) > 0.05
                    %                        largest_cov = [];
                    %                     end
                    %
                    %                     event_rrt(:,largest_cov) = [];
                    %                     fix_rrt(:,largest_cov) = [];
                    %
                    %                     if (event_rrt(2) < 0.05) && (fix_rrt(2) > 0.05)
                    %                         visuomotor = 1; %motor
                    %                     elseif (event_rrt(2) > 0.05) && (fix_rrt(2) < 0.05)
                    %                          visuomotor = 2; %visual
                    %                     else
                    %                         visuomotor = 3;%visual motor
                    %                     end
                    %
                    
                    %                     if  fix_prctile > 5 &&  event_prctile < 5
                    %                           visomotor = 1; %motor
                    %                     elseif fix_prctile < 5 &&  event_prctile > 5
                    %                         visomotor = 2;%visual
                    %                     else
                    %                          visomotor = 3; %visual motor
                    %                     end
                    %                         fix_prctile = sum(observed_fix_corr > shuff_fix_corrs);
                    %                      event_prctile = sum(observed_event_corr > shuff_event_corrs);
                    
                    %
                    figure
                    subplot(2,2,1)
                    [trial,time] = find(fixation_firing == 1);
                    plot(time-twin,(trial),'.k')
                    ylabel('Ranked by RT')
                    xlabel('Time from fixation Start (ms)')
                    title('Fixation Aligned')
                    xlim([-twin twin])
                    if ~isempty(trial)
                        ylim([0 max(trial)+1])
                    end
                    box off
                    axis square
                    
                    subplot(2,2,3)
                    plot(fix_curves2,'r')
                    hold on
                    plot(fix_curves1,'k')
                    hold off
                    ylabel('Firing Rate (Hz)')
                    xlabel('Time from Fixation Start (ms)')
                    title(num2str(observed_fix_corr,2))
                    box off
                    axis square
                    
                    subplot(2,2,2)
                    [trial,time] = find(event_firing == 1);
                    plot(time-twin,(trial),'.k')
                    ylabel('Ranked by RT')
                    xlabel('Time from Item On(ms)')
                    title('Item Appearance Aligned')
                    xlim([-twin twin])
                    if ~isempty(trial)
                        ylim([0 max(trial)+1])
                    end
                    box off
                    axis square
                    
                    
                    subplot(2,2,4)
                    hold on
                    plot(event_curves2,'r')
                    plot(event_curves1,'k')
                    hold off
                    ylabel('Firing Rate (Hz)')
                    xlabel('Time from Fixation Start (ms)')
                    title(num2str(observed_event_corr,2))
                    box off
                    axis square
                    
                    if observed_fix_corr > 0 && observed_event_corr < -0.1
                        visname = 'Motor';
                    elseif observed_fix_corr <-0.1 && observed_event_corr > 0
                        visname = 'Visual';
                    else
                        visname = 'Visual-Motor';
                    end
                    
                    
                    
                    subtitle([task_file(1:8) ' ' unit_stats{1,unit} ' ' visname])
                    close
                else
                    if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                            && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                        spatialness = [spatialness 1]; %place cell
                    else
                        spatialness = [spatialness 0]; %non-place cell
                    end
                    visual_motor_visuomotor = [visual_motor_visuomotor NaN];
                    unit_names = [unit_names [{task_file(1:8) '_3_' unit_stats{1,unit}}]];
                end
            elseif ~isnan(spatial_info.shuffled_rate_prctile(unit))
                if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                        && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                    spatialness = [spatialness 1]; %place cell
                else
                    spatialness = [spatialness 0]; %non-place cell
                end
                visual_motor_visuomotor = [visual_motor_visuomotor NaN];
                unit_names = [unit_names [{task_file(1:8) '_3_' unit_stats{1,unit}}]];
            end
        end
    end
end