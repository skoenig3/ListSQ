function Sequence_Contextual_Fixation_Analysis(data_dir,figure_dir,session_data)
% written 2/14/2017 by Seth Konig
% code anlayzes fixation aligned firing rate curves for overlapping
% sequence locations and asks wether neurons differentiate between contexts
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
contextual_sig_times = cell(4,num_units);
for unit = 1:num_units
    if ~isempty(fixation_locked_firing{unit})
        for o = 1:length(overlap)
            c = overlap(o);
            seq_num = sequence_nums{unit}(item_nums{unit} == c);
            firing_rate = fixation_locked_firing{unit}(item_nums{unit} == c,:);
            [~,firing_rate_curve] = nandens(firing_rate,smval,'gauss',Fs,'nanflt'); %trial-by-trial smoothed
            all_curves = NaN(numshuffs,twin*2);
            
            parfor shuff = 1:numshuffs;
                ind = randperm(length(seq_num));
                shuff_nums = seq_num(ind);%randomly assign fixations as in or out
                shuff_curve1 = nanmean(firing_rate_curve(shuff_nums == 1,:)); %sequence 1
                shuff_curve2 = nanmean(firing_rate_curve(shuff_nums == 2,:));%sequence 2
                all_curves(shuff,:) = shuff_curve1-shuff_curve2; %calculate shuffled difference
            end
            
            %for fixations out->in vs out->out
            curve1 = nanmean(firing_rate_curve(seq_num == 1,:));  %firing rate curve for fixations in the field
            curve2 = nanmean(firing_rate_curve(seq_num == 2,:)); %firing rate fo fixations out of the filed
            observed_diff = curve1-curve2; %observed difference in firing rate
            [~,contextual_sig_times{c,unit}] = cluster_level_statistic(observed_diff,all_curves,2,smval); %multiple comparision corrected significant indeces
        end
    end
end

%%%%%%%%%%%%%%%%%%%%
%---Plot Results---%
%%%%%%%%%%%%%%%%%%%%
sbs = [1 2 4 5];
t = -twin+1:twin;
for unit = 1:num_units
    if ~isempty(fixation_locked_firing{unit})
        
        figure
        
        for o = 1:length(overlap)
            c = overlap(o);
            
            subplot(2,3,sbs(c))
            hold on
            dofill(t,fixation_locked_firing{unit}((sequence_nums{unit} == 1) & (item_nums{unit} == c),:),'green',1,smval);
            dofill(t,fixation_locked_firing{unit}((sequence_nums{unit} == 2) & (item_nums{unit} == c),:),'magenta',1,smval);
            yl = ylim;
            if yl(1) < 0
                yl(1) = 0;
                ylim([0 yl(2)]);
            end
            plot([0 0],[yl(1) yl(2)],'k--')
            gaps = findgaps(find(contextual_sig_times{c,unit}));
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
            title_str = ['Item #' num2str(c) ':'];
            if sequence_info.rate_prctile(c,1,unit) > 95
                title_str = [title_str 'fix_{1} = ' num2str(sequence_info.rate_prctile(c,1,unit),3) '%%'];
            end
            if sequence_info.temporalstability_prctile(c,1,unit) > 95
                title_str = [title_str ', \\rho_{1} = ' num2str(sequence_info.temporalstability_prctile(c,1,unit),3) '%%'];
            end
            if sequence_info.rate_prctile(c,2,unit) > 95
                title_str = [title_str 'fix_{2} = ' num2str(sequence_info.rate_prctile(c,2,unit),3) '%%'];
            end
            if sequence_info.temporalstability_prctile(c,2,unit) > 95
                title_str = [title_str ', \\rho_{2} = ' num2str(sequence_info.temporalstability_prctile(c,2,unit),3) '%%'];
            end
            if ~isempty(title_str)
                title(sprintf(title_str))
            end
            title(sprintf(title_str))
        end
        
        subplot(2,3,3)
        hold on
        dofill(t,fixation_locked_firing{unit}((sequence_nums{unit} == 1),:),'green',1,smval);
        dofill(t,fixation_locked_firing{unit}((sequence_nums{unit} == 2),:),'magenta',1,smval);
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
        title('All Sequence 1 vs All Sequence 2')

        n_str = [' n = ' num2str(round(sum(item_nums{unit} == c)/2))]; %average number could be off by 1 or 2
        if multiunit(unit)
            subtitle(['Fixation-Locked Multiunit ' task_file(1:8) ' ' unit_stats{1,unit} n_str]);
        else
            subtitle(['Fixation-Locked ' task_file(1:8) ' ' unit_stats{1,unit} n_str]);
        end
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Sequence_Fixation_Context']);
    end
end


%%%%%%%%%%%%%%%%%%%%
%---Save Results---%
%%%%%%%%%%%%%%%%%%%%
save([data_dir task_file(1:8) '-sequence_contextual_fixation_results.mat'],...
    'numshuffs','contextual_sig_times')