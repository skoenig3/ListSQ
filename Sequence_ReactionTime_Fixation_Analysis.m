function Sequence_ReactionTime_Fixation_Analysis(data_dir,figure_dir,session_data,predict_rt)
% written 2/14/2017 by Seth Konig
% code anlayzes fixation aligned firing rate curves for overlapping
% sequence locations and asks wether neurons differentiate between slow and
% fastst reaction times.
% loads data pre-processed by Sequence_Fixation_AnalysisV2.

% Inputs:
%   1) data_dir: direction where preprocessed data is located
%   2) figure_dir: location of where to put save figures
%   3) session_data: data containing relevant session data
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-sequence_contextual_fixation_results'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set Default Parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure_dir = [figure_dir 'Sequence Fixation Analysis\'];
numshuffs = 10000; %recommend this is between 100 & 1000, for bootstrapping to
Fs = 1000;
min_blks = 2;
predict_thresh = 10;%10% needed to do anlaysis/say did remember
upper_lower = 40;%40 perctile for anlaysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task = 'ListSQ';
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end

%load trial data
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials',...
    'hdr','fixationstats');
num_trials = length(cfg.trl);
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
if num_units == 0
    return
end

valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,'ListSQ');
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end

load([data_dir task_file(1:8) '-Fixation_locked_Sequence_results.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Permutation/Resampling Analysis for Significance---%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reaction_time_sig_times_fixation = cell(4,2,num_units);
reaction_time_sig_times_event = cell(4,2,num_units);
for unit = 1:num_units
    if ~isempty(fixation_locked_firing{unit})
        for c  = 1:4
            for s = 1:2
                rts = reaction_time{unit}(item_nums{unit} == c & sequence_nums{unit} == s);
                
                firing_rate = fixation_locked_firing{unit}(item_nums{unit} == c & sequence_nums{unit} == s,:);
                [~,firing_rate_curve] = nandens(firing_rate,smval,'gauss',Fs,'nanflt'); %trial-by-trial smoothed
                
                event_rate = event_locked_firing{unit}(item_nums{unit} == c & sequence_nums{unit} == s,:);
                [~,event_rate_curve] = nandens(event_rate,smval,'gauss',Fs,'nanflt'); %trial-by-trial smoothed
                
                lower_prctile = prctile(rts,upper_lower); %lower limit
                upper_prctile = prctile(rts,100-upper_lower); %upper limit
                
                all_curves = NaN(numshuffs,twin*2);
                all_curves2 = NaN(numshuffs,twin*2);
                parfor shuff = 1:numshuffs;
                    ind = randperm(length(rts));
                    shuff_rts = rts(ind);%randomly assign fixations as in or out
                    
                    %---Fixation Aligned Firing---%
                    shuff_curve1 = nanmean(firing_rate_curve(shuff_rts <= lower_prctile,:)); %sequence 1
                    shuff_curve2 = nanmean(firing_rate_curve(shuff_rts >= upper_prctile,:));%sequence 2
                    all_curves(shuff,:) = shuff_curve1-shuff_curve2; %calculate shuffled difference
                    
                    %---Event Aligned Firing---%
                    shuff_curve1 = nanmean(event_rate_curve(shuff_rts <= lower_prctile,:)); %sequence 1
                    shuff_curve2 = nanmean(event_rate_curve(shuff_rts >= upper_prctile,:));%sequence 2
                    all_curves2(shuff,:) = shuff_curve1-shuff_curve2; %calculate shuffled difference
                end
                
                %for fixations
                curve1 = nanmean(firing_rate_curve(rts <= lower_prctile,:));  %firing rate curve for fixations in the field
                curve2 = nanmean(firing_rate_curve(rts >= upper_prctile,:)); %firing rate fo fixations out of the filed
                observed_diff = curve1-curve2; %observed difference in firing rate
                [~,reaction_time_sig_times_fixation{c,s,unit}] = cluster_level_statistic(observed_diff,all_curves,2,smval); %multiple comparision corrected significant indeces
                
                %for events
                curve1 = nanmean(event_rate_curve(rts <= lower_prctile,:));  %firing rate curve for fixations in the field
                curve2 = nanmean(event_rate_curve(rts >= upper_prctile,:)); %firing rate fo fixations out of the filed
                observed_diff = curve1-curve2; %observed difference in firing rate
                [~,reaction_time_sig_times_event{c,s,unit}] = cluster_level_statistic(observed_diff,all_curves,2,smval); %multiple comparision corrected significant indeces
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%
%---Plot Results---%
%%%%%%%%%%%%%%%%%%%%
t = -twin+1:twin;
for unit = 1:num_units
    if ~isempty(fixation_locked_firing{unit})
        
        %---For Fixations---%
        figure
        for c  = 1:4
            for s = 1:2
                rts = reaction_time{unit}(item_nums{unit} == c & sequence_nums{unit} == s);
                firing_rate = fixation_locked_firing{unit}(item_nums{unit} == c & sequence_nums{unit} == s,:);
                
                lower_prctile = prctile(rts,upper_lower); %lower limit
                upper_prctile = prctile(rts,100-upper_lower); %upper limit
                
                subplot(2,4,c+4*(s-1))
                hold on
                dofill(t,firing_rate(rts <= lower_prctile,:),'green',1,smval);
                dofill(t,firing_rate(rts >= upper_prctile,:),'red',1,smval);
                yl = ylim;
                if yl(1) < 0
                    yl(1) = 0;
                    ylim([0 yl(2)]);
                end
                plot([0 0],[yl(1) yl(2)],'k--')
                gaps = findgaps(find(reaction_time_sig_times_fixation{c,s,unit}));
                if ~isempty(gaps)
                    for g = 1:size(gaps,1)
                        gp = gaps(g,:);
                        gp(gp == 0) = [];
                        h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin,...
                            [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
                        uistack(h,'down')
                        set(h,'facealpha',.25,'EdgeColor','None')
                    end
                end
                hold off
                xlabel('Time from Fixation Start (ms)')
                ylabel('Firing Rate (Hz)')
                xlim([-twin twin])
                title(['Item#' num2str(c) ', Seq #' num2str(s) ' ']);
            end
        end
        n_str = [' n = ' num2str(round(sum(sequence_nums{unit} == s)/4))]; %average number could be off a few
        if multiunit(unit)
            subtitle(['Fixation-Locked Reaction Time: Multiunit ' task_file(1:8) ' ' unit_stats{1,unit} n_str]);
        else
            subtitle(['Fixation-Locked Reaction Time: ' task_file(1:8) ' ' unit_stats{1,unit} n_str]);
        end
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Sequence_Fixation_Reaction_Time']);
        
        %---For events---%
        figure
        for c  = 1:4
            for s = 1:2
                rts = reaction_time{unit}(item_nums{unit} == c & sequence_nums{unit} == s);
                firing_rate = event_locked_firing{unit}(item_nums{unit} == c & sequence_nums{unit} == s,:);
                
                lower_prctile = prctile(rts,upper_lower); %lower limit
                upper_prctile = prctile(rts,100-upper_lower); %upper limit
                
                subplot(2,4,c+4*(s-1))
                hold on
                dofill(t,firing_rate(rts <= lower_prctile,:),'green',1,smval);
                dofill(t,firing_rate(rts >= upper_prctile,:),'red',1,smval);
                yl = ylim;
                if yl(1) < 0
                    yl(1) = 0;
                    ylim([0 yl(2)]);
                end
                plot([0 0],[yl(1) yl(2)],'k--')
                gaps = findgaps(find(reaction_time_sig_times_event{c,s,unit}));
                if ~isempty(gaps)
                    for g = 1:size(gaps,1)
                        gp = gaps(g,:);
                        gp(gp == 0) = [];
                        h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin,...
                            [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
                        uistack(h,'down')
                        set(h,'facealpha',.25,'EdgeColor','None')
                    end
                end
                hold off
                xlabel('Time from event Start (ms)')
                ylabel('Firing Rate (Hz)')
                xlim([-twin twin])
                xlim([-twin twin])
                title(['Item#' num2str(c) ', Seq #' num2str(s) ' ']);
            end
        end
        n_str = [' n = ' num2str(round(sum(sequence_nums{unit} == s)/4))]; %average number could be off a few
        if multiunit(unit)
            subtitle(['event-Locked Reaction Time: Multiunit ' task_file(1:8) ' ' unit_stats{1,unit} n_str]);
        else
            subtitle(['event-Locked Reaction Time: ' task_file(1:8) ' ' unit_stats{1,unit} n_str]);
        end
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Sequence_Event_Reaction_Time']);
        
        for unit = 1:num_units
            
            any_predictive = 0;
            sequences = [];
            figure
            for seq = 1:2
                for c = 1:4
                    if  (sequence_info.rate_prctile(c,s,unit) > 95) && (sequence_info.temporalstability_prctile(c,s,unit) > 95)
                        rts = reaction_time{unit}(item_nums{unit} == c & sequence_nums{unit} == s);
                        firing_rate = fixation_locked_firing{unit}(item_nums{unit} == c & sequence_nums{unit} == s,:);
                        predict_percent = 100*sum(rts < predict_rt)/length(rts);
                        
                        if predict_percent > predict_thresh
                            any_predictive = 1;
                            sequences = [sequences seq];
                            
                            predictive = find(rts < predict_rt);
                            if predict_percent < 33
                                upper_limit = prctile(rts,100-predict_percent);
                                average = find(rts > predict_rt & rts < upper_limit);
                                slowest = find(rts > upper_limit);
                            else
                                average = [];
                                slowest = find(rts > predict_rt);
                            end
                            
                            subplot(2,4,c+4*(s-1))
                            hold on
                            dofill(t,firing_rate(predictive,:),'green',1,smval);
                            dofill(t,firing_rate(slowest,:),'red',1,smval);
                            if predict_percent < 33
                                dofill(t,firing_rate(average,:),'black',1,smval);
                                legend('Predictive','Reactive','Slow Reactive')
                            else
                                legend('Predictive','Reactive')
                            end
                            yl = ylim;
                            if yl(1) < 0
                                yl(1) = 0;
                                ylim([0 yl(2)]);
                            end
                            plot([0 0],[yl(1) yl(2)],'k--')
                            hold off
                            xlabel('Time from Fixation Start (ms)')
                            ylabel('Firing Rate (Hz)')
                            xlim([-twin twin])
                            
                            title_str = ['Item#' num2str(c) ', Seq#' num2str(s) ', n_{predict} = ' num2str(predict_percent,2) '%%'];
                            title(sprintf(title_str));
                        end
                    end
                end
            end
            n_str = [' n = ' num2str(round(sum(sequence_nums{unit} == s)/4))]; %average number could be off a few
            if multiunit(unit)
                subtitle(['event-Locked Reaction Time: Multiunit ' task_file(1:8) ' ' unit_stats{1,unit} n_str]);
            else
                subtitle(['event-Locked Reaction Time: ' task_file(1:8) ' ' unit_stats{1,unit} n_str]);
            end
            
            if any_predictive == 0
                close
            else
                save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Sequence_Event_Predictive']);
            end
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%
%---Save Results---%
%%%%%%%%%%%%%%%%%%%%
save([data_dir task_file(1:8) '-sequence_fixation_Reaction_Time_results.mat'],...
    'numshuffs','reaction_time_sig_times_fixation','reaction_time_sig_times_event')