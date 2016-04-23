function  Sequence_VisualResponse_TtestAnalysis(data_dir,preprocessed_data_file,figure_dir,task)
%written by Seth Konig on 5/5/15
%code runs a sliding window t-test to determine if context 1 is encoded
%different than context 2 locked to item onset/visually responsive

%load previously collected data
load([data_dir preprocessed_data_file],'item_set','cfg','num_units','multiunit');
load([data_dir preprocessed_data_file(1:8) '-ListSQ-time_locked_results.mat']);

[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
overlap = find((sequence_locations{1}(1,:) ==  sequence_locations{2}(1,:)) & ...
    (sequence_locations{1}(2,:) ==  sequence_locations{2}(2,:)));

events = [2 5 8 11 13]; %item 1 on, item 2 on, item 3 on, item 4 on, reward period

t_test_pvals = cell(5,num_units);


%for items 1-4
for epic = 1:4
    for unit = 1:num_units
        for period = 1:91
            temp_center_index = (period+4)*10;
            t_window = temp_center_index-49:temp_center_index+49;
            seq1 = time_lock_firing{events(epic),unit}(which_sequence == 1,t_window);% spikes for sequence 1 within temporal window
            seq2 = time_lock_firing{events(epic),unit}(which_sequence == 2,t_window);% spikes for sequence 1 within temporal window
            seq1 = sum(seq1'); %total spikes within window per trial equiv to firing rate
            seq2 = sum(seq2'); %total spikes within window per trial equiv to firing rate
            [~,t_test_pvals{epic,unit}(period)] = ttest2(seq1,seq2);
        end
    end
end

%for reward period
epic = 5;
for unit = 1:num_units
    for period = 1:floor(size(time_lock_firing{13,1},2)/10-9)
        temp_center_index = (period+4)*10;
        t_window = temp_center_index-49:temp_center_index+49;
        seq1 = time_lock_firing{events(epic),unit}(which_sequence == 1,t_window);% spikes for sequence 1 within temporal window
        seq2 = time_lock_firing{events(epic),unit}(which_sequence == 2,t_window);% spikes for sequence 1 within temporal window
        seq1 = sum(seq1'); %total spikes within window per trial equiv to firing rate
        seq2 = sum(seq2'); %total spikes within window per trial equiv to firing rate
        [~,t_test_pvals{epic,unit}(period)] = ttest2(seq1,seq2);
    end
end

xlabels = {'Time from Item 1 On (ms)','Time from Item 2 On (ms)','Time from Item 3 On (ms)', ...
    'Time from Item 4 On (ms)','Time from Reward Period Start (ms)'};
plotnums = {1,2,3,4,5:6};
twin = 500;
timewindow = {-twin:twin-1,-twin:twin-1,-twin:twin-1,-twin:twin-1,...
    -twin:size(time_lock_firing{13,1},2)-twin-1};
unit_names = cfg.channel;
smval = 120;
num_trials = length(which_sequence);
for unit = 1:num_units
    figure
    ylims = NaN(1,4);
    for epic = 1:5
        subplot(2,3,plotnums{epic})
        hold on
        if ~all(all(isnan(time_lock_firing{events(epic),unit}(which_sequence == 1,:))))
            dofill(timewindow{epic},time_lock_firing{events(epic),unit}(which_sequence == 1,:),'blue',1,smval); %reactive sequence 1
        end
        if ~all(all(isnan(time_lock_firing{events(epic),unit}(which_sequence == 1,:))))
            dofill(timewindow{epic},time_lock_firing{events(epic),unit}(which_sequence == 2,:),'red',1,smval); %reactive sequence 1
        end
        hold off
        xlabel(xlabels{epic})
        ylabel('Firing Rate (Hz)')
        yl = ylim;
        ylims(epic) = yl(2);
        
        if any(overlap == epic)
            title(['Overlapping Item # ' num2str(epic)])
        end
    end
    
    maxy = max(ylims);
    for epic = 1:5;
        subplot(2,3,plotnums{epic})
        sigp = find(t_test_pvals{epic,unit} < 0.05);
        if length(sigp) > 1
            new_pthresh = 0.05/length(sigp);%simple multicomparison correction
        else
            new_pthresh = 0.05;
        end
        sigp = find(t_test_pvals{epic,unit} < new_pthresh);
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
        subtitle(['Multiunit ' unit_names{unit} ' n =' num2str(num_trials)]);
    else
        subtitle([unit_names{unit} ' n =' num2str(num_trials)]);
    end
    save_and_close_fig(figure_dir,[preprocessed_data_file(1:end-11) '_' unit_names{unit} '_Sequence_VisualResponsive_Locked_Ttest_analysis']);
    
    
end

save([data_dir preprocessed_data_file(1:8) '-VisualResponse_Locked_Sequence_Ttest_results.mat'],...
    't_test_pvals');
end