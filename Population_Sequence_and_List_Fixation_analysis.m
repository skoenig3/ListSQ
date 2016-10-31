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

list_fixation_firing_rates = [];
seq_fixation_firing_rates = [];
list_mean_fixation_firing_rates = [];
seq_mean_fixation_firing_rates = [];

figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Population Figures\Place Cell Eye Movements\';
for monkey = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        listsq_read_excel(data_dir,excel_file);
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
            seq_fix = fixation_locked_firing; %since same variable
            seq_info = fixation_info;
            seq_info2 = fixation_info2;
            load([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat']);
        else
            continue
        end
        
        disp(['Importing ' task_file(1:8)])
        
        for unit = 1:num_units
            if ~isempty(fixation_locked_firing{unit})
                if ((fixation_info.rate(unit) >  fixation_info_95.rate(unit)) &&  (saccade_info.rate(unit) >  saccade_info_95.rate(unit))) || ...
                        ((fixation_info2.rate(unit) >  fixation_info2_95.rate(unit)) &&  (saccade_info2.rate(unit) >  saccade_info2_95.rate(unit)))
                    if seq_info.rate_prctile(unit) > 95 ||  seq_info2.rate_prctile(unit) > 95
                        col = 1;
                    else
                        continue
                    end
                else
                    continue
                end
                
                
                fixation_firing = fixation_locked_firing{unit}(fixation_information{unit}(:,4) > 2*twin,:);                fixation_firing(nansum(fixation_firing,2) == 0,:) = [];%remove fixations without spikes
                [firing_rate,~]= nandens(fixation_firing,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                firing_rate = firing_rate-nanmean(firing_rate);
                list_mean_fixation_firing_rates = [list_mean_fixation_firing_rates; firing_rate];
                seq_firing_rate = firing_rate;
                fr = firing_rate;
%                 if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                    firing_rate = firing_rate/max(abs(firing_rate));
%                 else%normalize to min %could have some neurons that only show supression
%                     firing_rate = firing_rate/min(firing_rate);
%                 end
                list_fixation_firing_rates = [list_fixation_firing_rates; firing_rate];
                
                
                fixation_firing = seq_fix{unit};
                fixation_firing(nansum(fixation_firing,2) == 0,:) = [];%remove fixations without spikes
                [firing_rate,~]= nandens(fixation_firing,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                firing_rate = firing_rate-nanmean(firing_rate);
                seq_mean_fixation_firing_rates = [seq_mean_fixation_firing_rates; firing_rate];
                list_firing_rate = firing_rate;
%                 if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                    firing_rate = firing_rate/max(abs(firing_rate));
%                 else%normalize to min %could have some neurons that only show supression
%                     firing_rate = firing_rate/min(firing_rate);
%                 end
                seq_fixation_firing_rates = [seq_fixation_firing_rates; firing_rate];
                
          
                
            end
        end
    end
end
%% All Cells Aligned to Max
figure
subplot(1,3,1)
[m,i] = max(list_fixation_firing_rates(:,450:750),[],2);
[mm,ii] = sort(i);
ii3 = ii;
imagesc([-500:499],[1:size(list_fixation_firing_rates,1)],list_fixation_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(list_fixation_firing_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('List')

subplot(1,3,2)
[m,i] = max(seq_fixation_firing_rates(:,:),[],2);
[mm,ii] = sort(i);
imagesc([-500:499],[1:size(seq_fixation_firing_rates,1)],seq_fixation_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(seq_fixation_firing_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Sequence')

subplot(1,3,3)
imagesc([-500:499],[1:size(seq_fixation_firing_rates,1)],seq_fixation_firing_rates(ii3,:))
hold on
plot([0 0],[1 size(seq_fixation_firing_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('List sorted as  Sequence')
%%
for i = 50:60%size(seq_fixation_firing_rates,1)
    figure
    plot(t,list_mean_fixation_firing_rates(i,:))
    hold on
    plot(t,seq_mean_fixation_firing_rates(i,:))
    yl = ylim;
    plot([0 0],[yl(1) yl(2)],'k--')
    hold off
    legend('List','Seq')
    xlabel('Time from Fixation Start (ms)')
    ylabel('Zero-mean firing rate')
end
%%
%% All Cells Aligned to Min
figure
subplot(2,2,1)
[m,i] = min(list_fixation_firing_rates(:,450:850),[],2);
[mm,ii] = sort(i);
ii3 = ii;
imagesc([-500:499],[1:size(list_fixation_firing_rates,1)],list_fixation_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(list_fixation_firing_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('List')

subplot(2,2,2)
[m,i] = min(seq_fixation_firing_rates(:,:),[],2);
[mm,ii] = sort(i);
imagesc([-500:499],[1:size(seq_fixation_firing_rates,1)],seq_fixation_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(seq_fixation_firing_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Sequence')

subplot(2,2,3)
imagesc([-500:499],[1:size(list_fixation_firing_rates,1)],list_fixation_firing_rates(ii,:))
hold on
plot([0 0],[1 size(seq_fixation_firing_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('List sorted as  Sequence')

subplot(2,2,4)
imagesc([-500:499],[1:size(seq_fixation_firing_rates,1)],seq_fixation_firing_rates(ii3,:))
hold on
plot([0 0],[1 size(seq_fixation_firing_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Sequence sorted as List')




