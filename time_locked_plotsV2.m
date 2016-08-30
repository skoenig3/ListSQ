function time_locked_plotsV2(figure_dir,task_file,time_lock_firing,epoch_data,...
    infodata,taskdata,unit_names,smval,fr_threshold)
% written by Seth Konig August, 2014. Updated to V2 during January 2016 to
% handle partial session data for when each unit was stable SDK
% generates and save plots for time locked spike analysis
%
% Inputs:
%   1) figure_dir: directory where to save figures
%   2) task_file: name of the preprocessed_file
%   2) time_lock_firing: matrix containing spike times locked to events by trial
%   3) epoch_data...
%       a) epoch_data.firing_rates: contains distribution of spikes per trial for each
%       b) epoch_data.dur: duration of epoch
%       c) epoch_data.num_trials: number of trials in the epoch
%       d) epoch_data.pvals: contains p-vals for whether distribution of
%       spikes varies between 2 or more conditions (i.e. sequence or
%       novelty) in any epoch
%       e) epoch_data.event_names: epoch name
%       f) epoch_data.event_codes: epoch event cortex code
%       g) epoch_data.event_durs: epoch duration
%       h) epoch_data.event_t0: epoch start index in time_lock_firing matrix
%   4) infodata...
%       a) infodata.rate: the observed information rate in bits/sec
%       b) infodata.shuffled_info_rate: information rate expected by chance
%       c) infodata.shuffled_95_percentile: 95-perecentile of boostrapped
%       information rate (chance/shuffled information rate)
%       d) infodata.shuffled_90_percentile: 90-perecentile of boostrapped
%   5) task_type: type of plot to create by task and/or subtask...
%       a) 'ListSQ_Sequence': plot for sequence trials in ListSQ task
%       b) 'ListSQ_List': plot for List images trials in ListSQ task
%       c) ''cvtnew': plot for cvtnew task
%   6) unit_names...
%       a) unit_names.name: name of unit e.g. sig001a
%       b) unit_names.multiunit: 1 for multiunit 0 for single unit
%   7) smval: smoothing parameters of varying lenghts depending on task
%
% Outputs:
%   1) saved figures into figure_dir

tooslow_firing_folder = [figure_dir '\FiringRateTooSlow_Temporal\'];

num_units = size(time_lock_firing,2);
task_type = taskdata.task_type;
peak_firing_rate = taskdata.peak_firing_rate;

switch task_type
    case 'ListSQ_Sequence'
        
        which_sequence = taskdata.which_sequence;
        for unit = 1:num_units
            if isempty(time_lock_firing{1,unit})
                continue
            end
            figure
            
            ylims = NaN(1,13);
            
            [s(1),ylims(1)]=make_subplot([4,4],1,time_lock_firing{1,unit},...
                epoch_data.event_t0(1),250,epoch_data.event_names{1}{1},....
                infodata.rate(1,unit,1),infodata.shuffled_95_percentile(1,unit,1),...
                infodata.shuffled_90_percentile(1,unit,1),smval);
            
            for event = 2:13
                time_seq_1 = time_lock_firing{event,unit}(which_sequence{unit} == 1,:);
                time_seq_2 = time_lock_firing{event,unit}(which_sequence{unit} == 2,:);
                [s(event),ylims(event)]= make_subplot_2CND([4,4],event,time_seq_1,time_seq_2,...
                    epoch_data.event_t0(event),250,epoch_data.event_names{event},....
                    infodata.rate(event,unit,1),infodata.shuffled_95_percentile(event,unit,1),...
                    infodata.shuffled_90_percentile(event,unit,1),smval);
            end
            
            means = NaN(2,14);
            stds = NaN(2,14);
            numpoints = NaN(2,14);
            total_spikes = cell(1,2);
            for seq = 1:2
                for event = 1:13;
                    if event == 1
                        dur = epoch_data.dur(event,unit,1);
                        numpoints(seq,event) = epoch_data.num_trials(event,unit,1);
                        means(seq,event) = nanmean(epoch_data.firing_rates{event,unit,1})/dur;
                        stds(seq,event) = nanstd(epoch_data.firing_rates{event,unit,1})/dur;
                    else
                        dur = epoch_data.dur(event,unit,1);
                        means(seq,event) = nanmean(epoch_data.firing_rates{event,unit,seq})/dur;
                        stds(seq,event) = nanstd(epoch_data.firing_rates{event,unit,seq})/dur;
                        numpoints(seq,event) = epoch_data.num_trials(event,unit,seq);
                    end
                end
            end
            
            if any(peak_firing_rate(:,unit) > fr_threshold)
                temp_data1 =[];
                temp_data2 =[];
                for event = 2:13;
                    temp_data1 = [temp_data1; epoch_data.firing_rates{event,unit,1}/epoch_data.dur(event,unit,1)];
                    temp_data2 = [temp_data2; epoch_data.firing_rates{event,unit,2}/epoch_data.dur(event,unit,1)];
                end
                means(1,end) = nanmean(temp_data1);
                means(2,end) = nanmean(temp_data2);
                stds(1,end) = nanstd(temp_data1);
                stds(2,end) = nanstd(temp_data2);
                numpoints(1,end) = sum(~isnan(temp_data1));
                numpoints(2,end) = sum(~isnan(temp_data2));
                [~,all_pval] =kstest2(temp_data1,temp_data2); %determine if across all epochs there a difference in firing rate
                
                subplot(4,4,[14 15 16])
                hold on
                bar(means')
                errorb(means',stds'./sqrt(numpoints'));
                ylabel('Firing Rate')
                xlabel('Epoch')
                title('Average Firing Rate by Epoch')
                set(gca,'Xtick',1:14)
                set(gca,'XtickLabel',[num2cell(1:13) {'all'}])
                xlim([0 15])
                yl = ylim;
                ylim([0 yl(2)]);
                
                % plot astriks (*) above epochs in which firing rate was
                % significantly different between sequences in an epoch
                if all_pval < 0.05
                    plot(14,yl(2)-yl(2)/10,'*k')
                end
                sig_sequence_effect = find(epoch_data.pvals(:,unit) < 0.05);
                for sse = 1:length(sig_sequence_effect)
                    plot(sig_sequence_effect(sse),yl(2)-yl(2)/10,'*k')
                end
                hold off
            end
            
            
            max_y = 5*round(max(ylims)/5);
            max_y(max_y < 1) = 1;
            for i = 1:13;
                set(s(i),'ylim',[0 max_y])
            end
            
            if unit_names.multiunit(unit)
                subtitle(['Multiunit ' unit_names.name{unit} ]);
            else
                subtitle(unit_names.name{unit});
            end
            
            if all(peak_firing_rate(:,unit) < fr_threshold) %silent neuron
                save_and_close_fig(tooslow_firing_folder,[task_file(1:end-11) '_' unit_names.name{unit} '_Sequence-time_locked_analysis']);
            else
                if unit_names.multiunit(unit)
                    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names.name{unit} '_List-time_locked_analysis']);
                else
                    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names.name{unit} '_Sequence-time_locked_analysis']);
                end
            end
        end
        
    case 'ListSQ_List'
        
        
        smval_novrep = smval(2);
        smval = smval(1);
        novel_vs_repeat = taskdata.novel_vs_repeat;
        for unit = 1:num_units
            if isempty(time_lock_firing{1,unit})
                continue
            end
            
            figure
            
            ylims = NaN(1,6);
            s = [];
            
            [s(1),ylims(1)]=make_subplot([2,4],1,time_lock_firing{14,unit},...
                epoch_data.event_t0(14),250,epoch_data.event_names{14},....
                infodata.rate(14,unit,1),infodata.shuffled_95_percentile(14,unit,1),...
                infodata.shuffled_90_percentile(14,unit,1),smval);
            
            [s(2),ylims(2)]=make_subplot([2,4],5,time_lock_firing{15,unit},...
                epoch_data.event_t0(15),250,epoch_data.event_names{15},....
                infodata.rate(15,unit,1),infodata.shuffled_95_percentile(15,unit,1),...
                infodata.shuffled_90_percentile(15,unit,1),smval);
            
            novel_time_lock_firing = time_lock_firing{16,unit}(novel_vs_repeat{unit} == 1,:);
            repeat_time_lock_firing = time_lock_firing{16,unit}(novel_vs_repeat{unit} == 2,:);
            [s(3),ylims(3)]=make_subplot_2CND([2,4],[2 3],novel_time_lock_firing,repeat_time_lock_firing,...
                epoch_data.event_t0(16),500,epoch_data.event_names{16},infodata.rate(16,unit,end-2:end),...
                infodata.shuffled_95_percentile(16,unit,end-2:end),...
                infodata.shuffled_90_percentile(16,unit,end-2:end),smval_novrep);
            
            novel_time_lock_firing = time_lock_firing{17,unit}(novel_vs_repeat{unit} == 1,:);
            repeat_time_lock_firing = time_lock_firing{17,unit}(novel_vs_repeat{unit} == 2,:);
            [s(4),ylims(4)]=make_subplot_2CND([2,4],[6 7],novel_time_lock_firing,repeat_time_lock_firing,...
                epoch_data.event_t0(17),500,epoch_data.event_names{17},infodata.rate(17,unit,end-2:end),...
                infodata.shuffled_95_percentile(17,unit,end-2:end),...
                infodata.shuffled_90_percentile(17,unit,end-2:end),smval_novrep);
            
            [s(5),ylims(5)]=make_subplot([2,4],4,time_lock_firing{18,unit},...
                epoch_data.event_t0(18),250,epoch_data.event_names{18},....
                infodata.rate(18,unit,1),infodata.shuffled_95_percentile(18,unit,1),...
                infodata.shuffled_90_percentile(18,unit,1),smval);
            
            [s(6),ylims(6)]=make_subplot([2,4],8,time_lock_firing{19,unit},...
                epoch_data.event_t0(19),250,epoch_data.event_names{19},....
                infodata.rate(19,unit,1),infodata.shuffled_95_percentile(19,unit,1),...
                infodata.shuffled_90_percentile(19,unit,1),smval);
            
            max_y = 5*round(max(ylims)/5);
            max_y(max_y < 1) = 1;
            for i = 1:6;
                set(s(i),'ylim',[0 max_y])
            end
            
            if unit_names.multiunit(unit)
                subtitle(['Multiunit ' unit_names.name{unit} ]);
            else
                subtitle(unit_names.name{unit});
            end
            
             if unit_names.multiunit(unit)
                 subtitle(['Multiunit ' unit_names.name{unit} ]);
             else
                 subtitle(unit_names.name{unit});
             end
             if all(peak_firing_rate(:,unit) < fr_threshold) %silent neuron
                 save_and_close_fig(tooslow_firing_folder,[task_file(1:end-11) '_' unit_names.name{unit} '_List-time_locked_analysis']);
             else
                 if unit_names.multiunit(unit)
                     save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names.name{unit} '_List-time_locked_analysis']);
                 else
                     save_and_close_fig([figure_dir '\MultiUnit\'],[task_file(1:end-11) '_' unit_names.name{unit} '_List_time_locked_analysis']);
                 end
             end
             
             
            
            %plot image data locked to crosshair appearance and fixation
            figure
            ylims = NaN(1,2);
            s = [];
            
            novel_time_lock_firing = time_lock_firing{20,unit}(novel_vs_repeat{unit} == 1,:);
            repeat_time_lock_firing = time_lock_firing{20,unit}(novel_vs_repeat{unit} == 2,:);
            [s(1),ylims(1)]=make_subplot_2CND([2,1],1,novel_time_lock_firing,repeat_time_lock_firing,...
                epoch_data.event_t0(20),500,epoch_data.event_names{20},infodata.rate(20,unit,end-2:end),...
                infodata.shuffled_95_percentile(20,unit,end-2:end),...
                infodata.shuffled_90_percentile(20,unit,end-2:end),smval_novrep);
            
            novel_time_lock_firing = time_lock_firing{21,unit}(novel_vs_repeat{unit} == 1,:);
            repeat_time_lock_firing = time_lock_firing{21,unit}(novel_vs_repeat{unit} == 2,:);
            [s(2),ylims(2)]=make_subplot_2CND([2,1],2,novel_time_lock_firing,repeat_time_lock_firing,...
                epoch_data.event_t0(21),500,epoch_data.event_names{21},infodata.rate(21,unit,end-2:end),...
                infodata.shuffled_95_percentile(21,unit,end-2:end),...
                infodata.shuffled_90_percentile(21,unit,end-2:end),smval_novrep);
            
            max_y = 5*round(max(ylims)/5);
            max_y(max_y < 1) = 1;
            for i = 1:2;
                set(s(i),'ylim',[0 max_y])
            end
            
            if unit_names.multiunit(unit)
                subtitle(['Multiunit ' unit_names.name{unit} ]);
            else
                subtitle(unit_names.name{unit});
            end
            if all(peak_firing_rate(:,unit) < fr_threshold) %silent neuron
                save_and_close_fig(tooslow_firing_folder,[task_file(1:end-11) '_' unit_names.name{unit} '_List2-time_locked_analysis']);
            else
                if unit_names.multiunit(unit)
                    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names.name{unit} '_List2-time_locked_analysis']);
                else
                    save_and_close_fig([figure_dir '\MultiUnit\'],[task_file(1:end-11) '_' unit_names.name{unit} '_List2_time_locked_analysis']);
                end
            end
        end
        
    case {'cvtnew','CVTNEW'}
        
        for unit = 1:num_units
            figure
            
            ylims = NaN(1,7);
            
            [s(1),ylims(1)]=make_subplot([3,3],1,time_lock_firing{1,unit},...
                epoch_data.event_t0(1),250,epoch_data.event_names{1}{1},....
                infodata.rate(1,unit,1),infodata.shuffled_95_percentile(1,unit,1),...
                infodata.shuffled_90_percentile(1,unit,1),smval);
            
            for event = [2:3 5:7]
                [s(event),ylims(event)] = make_subplot([3,3],event,time_lock_firing{event,unit},...
                    epoch_data.event_t0(event),250,epoch_data.event_names{event},....
                    infodata.rate(event,unit),infodata.shuffled_95_percentile(event,unit),...
                    infodata.shuffled_90_percentile(event,unit),smval);
            end
            
            %since event 4 is the dot on-dot color change the duration of
            %this event is uniformly random and there is no other locking
            %possible, it gets it's own subplot function
            event = 4;
            [s(event),ylims(event)] = make_subplot4([3,3],event,time_lock_firing{event,unit},...
                epoch_data.event_t0(event),250,epoch_data.event_names{event},....
                infodata.rate(event,unit),infodata.shuffled_95_percentile(event,unit),...
                infodata.shuffled_90_percentile(event,unit),smval,epoch_data.trial_duration{unit});
            
            
            means = NaN(1,8);
            stds = NaN(1,8);
            numpoints = NaN(1,8);
            total_spikes = cell(1,2);
            for event = 1:7;
                dur = epoch_data.dur(event,unit,1);
                numpoints(event) = epoch_data.num_trials(event,unit);
                means(event) = nanmean(epoch_data.firing_rates{event,unit})/dur;
                stds(event) = nanstd(epoch_data.firing_rates{event,unit})/dur;
            end
            
            temp_data =[];
            for event = 2:7;
                temp_data = [temp_data; epoch_data.firing_rates{event,unit}/epoch_data.dur(event,unit)];
            end
            means(8) = nanmean(temp_data);
            stds(8) = nanstd(temp_data);
            numpoints(8) = sum(~isnan(temp_data));
            
            subplot(3,3,[8 9])
            hold on
            bar(means)
            errorb(means,stds./sqrt(numpoints));
            hold off
            ylabel('Firing Rate')
            xlabel('Epoch')
            title('Average Firing Rate by Epoch')
            set(gca,'Xtick',1:8)
            set(gca,'XtickLabel',[num2cell(1:7) {'all'}])
            xlim([0 9])
            yl = ylim;
            ylim([0 yl(2)]);
            
            max_y = max(ylims);
            max_y(max_y < 1) = 1;
            for i = 1:7;
                set(s(i),'ylim',[0 max_y])
            end
            
            if unit_names.multiunit(unit)
                subtitle(['Multiunit ' unit_names.name{unit} ]);
            else
                subtitle(unit_names.name{unit});
            end
            save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names.name{unit} '_cvtnew-time_locked_analysis']);
        end
end
end

function [s,ylimit]=make_subplot(src,plot_nums,time_matrix,t0,timestep,epoch_name,inforate,info95,info90,smval)
% src: subplot row and column
% plot_nums: plot number or numbers
% inforate: observed inforamtion rate
% info95: bootstrapped 95 percentile

s = subplot(src(1),src(2),plot_nums);
t = 1:size(time_matrix,2);
if ~isempty(t)
    [~,~,~,y] =dofill(t,time_matrix,'black',1,smval);
    set(gca,'Xtick',0:timestep:size(time_matrix,2))
    set(gca,'XtickLabel',num2cell((0:timestep:size(time_matrix,2))-t0));
    ylabel('Firing Rate (Hz)')
    xlabel(['Time from ' epoch_name ' (ms)'])
    num_trial_str = ['n = ' num2str(size(time_matrix,1))];
    if any(inforate > info95)
        title([num_trial_str ', Bits_{95} = ' num2str(inforate(end))]);
    elseif any(inforate > info90)
        title([num_trial_str ', Bits_{90} = ' num2str(inforate(end))]);
    else
        title(num_trial_str);
    end
    ylimit = max(y)*1.1;
else
    ylimit = NaN;
end
end

function [s,ylimit]=make_subplot4(src,plot_nums,time_matrix,t0,timestep,epoch_name,inforate,info95,info90,smval,trial_duration,twin)
%subplotting function just for 4th event
% src: subplot row and column
% plot_nums: plot number or numbers
% inforate: observed inforamtion rate
% info95: bootstrapped 95 percentile

%CVTNEW variable trial lenght
%time is + 300 ms
short = [700 1133]; %short duration trials
medium = [1134 1567];
long = [1568 2000];%cap should be 2000 but cortex can have some lag

dur_thresh = [short; medium; long];


s = subplot(src(1),src(2),plot_nums);
t = 1:size(time_matrix,2);

clr = ['bgr'];
ylimit = 0;
if ~isempty(t) && ~isempty(time_matrix)
    for dur = 1:size(dur_thresh,1)
        these_length_trials = find(trial_duration >= dur_thresh(dur,1) & trial_duration <= dur_thresh(dur,2));
        
        these_time_matrix = time_matrix(these_length_trials,1:t0+dur_thresh(dur,2)); %plot window + 500 ms end
        these_time = t(1:t0+dur_thresh(dur,2));
        
        %find time points in which we have engouh data to plot as to not
        %create artifiacts from low data counts
        total_trials = sum(~isnan(these_time_matrix));
        enough_data = find(total_trials > 0.33*total_trials(1));
        these_time = these_time(enough_data);
        these_time_matrix = these_time_matrix(:,enough_data);
        
        [~,~,~,y] =dofill(these_time,these_time_matrix,clr(dur),1,smval);
        if 1.1*max(y) < 2.5*mean(y)
            ylimit = max(ylimit,1.1*max(y));
        else
            ylimit = max(ylimit,2.5*mean(y));
        end
    end
end

set(gca,'Xtick',0:timestep:size(time_matrix,2))
set(gca,'XtickLabel',num2cell((0:timestep:size(time_matrix,2))-t0));
ylabel('Firing Rate (Hz)')
xlabel(['Time from ' epoch_name ' (ms)'])
if any(inforate > info95)
    title(['Bits_{95} = ' num2str(inforate(end))]);
elseif any(inforate) > info90
    title(['Bits_{90} = ' num2str(inforate(end))]);
end

end

function [s,ylimit]= make_subplot_2CND(src,plot_nums,time_matrix1,time_matrix2,t0,timestep,epoch_name,inforate,info95,info90,smval)
% essentailly same as above just 2 time_matrices
s = subplot(src(1),src(2),plot_nums);
t = 1:size(time_matrix1,2);
if ~isempty(t)
    dofill(t,time_matrix1,'blue',1,smval);
    dofill(t,time_matrix2,'red',1,smval);
end
set(gca,'Xtick',0:timestep:size(time_matrix1,2))
set(gca,'XtickLabel',num2cell((0:timestep:size(time_matrix1,2))-t0));
ylabel('Firing Rate (Hz)')
xlabel(['Time from ' epoch_name ' (ms)'])

num_trial_str = ['n_1 = ' num2str(size(time_matrix1,1)) ' n_2 = ' num2str(size(time_matrix2,1))];
if any(inforate > info95)
    title([num_trial_str ', Bits_{95} = ' num2str(inforate(end))]);
elseif any(inforate > info90)
    title([num_trial_str ', Bits_{90} = ' num2str(inforate(end))]);
else
    title(num_trial_str)
end
ylimit = ylim;
ylimit = ylimit(2);
end