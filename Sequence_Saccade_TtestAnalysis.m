
function  Sequence_Saccade_TtestAnalysis(data_dir,preprocessed_data_file,figure_dir,task)
%written by Seth Konig on 5/4/15
%code runs a sliding window t-test to determine if context 1 is encoded
%different than context 2
 
%load previously collected data
load([data_dir preprocessed_data_file],'item_set','cfg','num_units','multiunit');
load([data_dir preprocessed_data_file(1:8) '-Eyemovement_Locked_Sequence_results.mat']);

[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
overlap = find((sequence_locations{1}(1,:) ==  sequence_locations{2}(1,:)) & ...
    (sequence_locations{1}(2,:) ==  sequence_locations{2}(2,:)));

%preallocate space and parallel structure of cfg
successful_sequence_trials = NaN(1,length(cfg.trl));
which_sequence = NaN(1,length(cfg.trl));
for t = 1:length(cfg.trl);
    if sum(cfg.trl(t).allval == 3) >= 6; %in which sequence trials were rewarded
        which_sequence(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
        successful_sequence_trials(t) = t;
    end
end
successful_sequence_trials = laundry(successful_sequence_trials);
which_sequence = laundry(which_sequence);

saccade_t_test_pvals = cell(4,num_units);
fixation_t_test_pvals = cell(4,num_units);

for epic = 1:4
    for unit = 1:num_units
        for period = 1:91
            temp_center_index = (period+4)*10;
            t_window = temp_center_index-49:temp_center_index+49;
            sac1 = saccade_locked_firing{epic,unit}(which_sequence == 1,t_window);% spikes for sequence 1 within temporal window
            sac2 = saccade_locked_firing{epic,unit}(which_sequence == 2,t_window);% spikes for sequence 2 within temporal window
            sac1 = sum(sac1'); %total spikes within window per trial equiv to firing rate
            sac2 = sum(sac2'); %total spikes within window per trial equiv to firing rate
            [~,saccade_t_test_pvals{epic,unit}(period)] = ttest2(sac1,sac2);
            
            fix1 = fixation_locked_firing{epic,unit}(which_sequence == 1,t_window);% spikes for sequence 1 within temporal window
            fix2 = fixation_locked_firing{epic,unit}(which_sequence == 2,t_window);% spikes for sequence 2 within temporal window
            fix1 = sum(fix1'); %total spikes within window per trial equiv to firing rate
            fix2 = sum(fix2'); %total spikes within window per trial equiv to firing rate
            [~,fixation_t_test_pvals{epic,unit}(period)] = ttest2(fix1,fix2);
        end
    end
end

t = -twin:twin-1;
unit_names = cfg.channel;
smval = 120;
num_trials = length(which_sequence);
for unit = 1:num_units
    figure
    ylims = NaN(1,4);
    for c = 1:4
        subplot(2,2,c)
        hold on
        if ~all(all(isnan(saccade_locked_firing{c,unit}(which_sequence == 1,:))))
            dofill(t,saccade_locked_firing{c,unit}(which_sequence == 1,:),'red',1,smval); %reactive sequence 1
        end
        if ~all(all(isnan(saccade_locked_firing{c,unit}(which_sequence == 2,:))))
            dofill(t,saccade_locked_firing{c,unit}(which_sequence == 2,:),'blue',1,smval);%reactive sequence 2
        end
        hold off
        xlabel('Time from Eye Movement (ms)')
        ylabel('Firing Rate (Hz)')
        yl = ylim;
        ylims(c) = yl(2);
        
        bitstr = [];
        if saccade_info(c,unit) > saccade_info_95(c,unit)
            bitstr = [bitstr 'sac_{95} = ' num2str(saccade_info(c,unit))];
        elseif saccade_info(c,unit) > saccade_info_90(c,unit)
            bitstr = [bitstr 'sac_{90} = ' num2str(saccade_info(c,unit))];
        end
        
        if any(overlap == c)
            title(['Overlapping Item # ' num2str(c) ' ' bitstr])
        else
            title(['Item # ' num2str(c) ' ' bitstr])
        end
        
    end
    maxy = max(ylims);
    for c = 1:4;
        subplot(2,2,c)
        sigp = find(saccade_t_test_pvals{c,unit} < 0.05);
        if length(sigp) > 1
            new_pthresh = 0.05/length(sigp);%simple multicomparison correction
        else
            new_pthresh = 0.05;
        end
        sigp = find(saccade_t_test_pvals{c,unit} < new_pthresh);
        if ~isempty(sigp)
            [broken_ind]=findgaps(sigp);
            if isempty(broken_ind)
                broken_ind = sigp;
            end
            hold on
            for b = 1:size(broken_ind,1)
                ind = broken_ind(b,:);
                ind(ind == 0) = [];
                temp_center_index = (ind+4)*10;
                t_window = [temp_center_index(1)-9 temp_center_index(end)+9]-twin;
                h = fill([t_window(1) t_window(2) t_window(2) t_window(1) t_window(1)],...
                    [0 0 maxy maxy 0],'k');
                uistack(h,'bottom')
                set(h,'facealpha',.25,'EdgeColor','None')
            end
            hold off
        end
        ylim([0 maxy]);
    end
    
    if multiunit(unit)
        subtitle(['Saccade-Locked Multiunit ' unit_names{unit} ' n =' num2str(num_trials)]);
    else
        subtitle(['Saccade-Locked ' unit_names{unit} ' n =' num2str(num_trials)]);
    end
    save_and_close_fig(figure_dir,[preprocessed_data_file(1:end-11) '_' unit_names{unit} '_Sequence_Saccade_Locked_Ttest_analysis']);
    
    figure
    ylims = NaN(1,4);
    for c = 1:4
        subplot(2,2,c)
        hold on
        if ~all(all(isnan(fixation_locked_firing{c,unit}(which_sequence == 1,:))))
            dofill(t,fixation_locked_firing{c,unit}(which_sequence == 1,:),'red',1,smval); %reactive sequence 1
        end
        if ~all(all(isnan(fixation_locked_firing{c,unit}(which_sequence == 2,:))))
            dofill(t,fixation_locked_firing{c,unit}(which_sequence == 2,:),'blue',1,smval);%reactive sequence 2
        end
        hold off
        xlabel('Time from Eye Movement (ms)')
        ylabel('Firing Rate (Hz)')
        yl = ylim;
        ylims(c) = yl(2);
        
        bitstr = [];
        if fixation_info(c,unit) > fixation_info_95(c,unit)
            bitstr = [bitstr 'fix_{95} = ' num2str(fixation_info(c,unit))];
        elseif fixation_info(c,unit) > fixation_info_90(c,unit)
            bitstr = [bitstr 'fix_{90} = ' num2str(fixation_info(c,unit))];
        end
        
        if any(overlap == c)
            title(['Overlapping Item # ' num2str(c) ' ' bitstr])
        else
            title(['Item # ' num2str(c) ' ' bitstr])
        end
        
    end
    maxy = max(ylims);
    for c = 1:4;
        subplot(2,2,c)
        sigp = find(fixation_t_test_pvals{c,unit} < 0.05);
        if length(sigp) > 1
            new_pthresh = 0.05/length(sigp);%simple multicomparison correction
        else
            new_pthresh = 0.05;
        end
        sigp = find(fixation_t_test_pvals{c,unit} < new_pthresh);
        if ~isempty(sigp)
            [broken_ind]=findgaps(sigp);
            if isempty(broken_ind)
                broken_ind = sigp;
            end
            hold on
            for b = 1:size(broken_ind,1)
                ind = broken_ind(b,:);
                ind(ind == 0) = [];
                temp_center_index = (ind+4)*10;
                t_window = [temp_center_index(1)-9 temp_center_index(end)+9]-twin;
                h = fill([t_window(1) t_window(2) t_window(2) t_window(1) t_window(1)],...
                    [0 0 maxy maxy 0],'k');
                uistack(h,'bottom')
                set(h,'facealpha',.25,'EdgeColor','None')
            end
            hold off
        end
        ylim([0 maxy]);
    end
    
    if multiunit(unit)
        subtitle(['fixation-Locked Multiunit ' unit_names{unit} ' n =' num2str(num_trials)]);
    else
        subtitle(['fixation-Locked ' unit_names{unit} ' n =' num2str(num_trials)]);
    end
    save_and_close_fig(figure_dir,[preprocessed_data_file(1:end-11) '_' unit_names{unit} '_Sequence_fixation_Locked_Ttest_analysis']);
end

save([data_dir preprocessed_data_file(1:8) '-Eyemovement_Locked_Sequence_Ttest_results.mat'],...
    'saccade_t_test_pvals','fixation_t_test_pvals');
end