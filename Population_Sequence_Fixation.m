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

sig_count = zeros(1,3);

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
        else
            continue
        end
        
        disp(['Importing ' task_file(1:8)])
        
        for unit = 1:num_units
            if ~isempty(fixation_locked_firing{unit})
                if fixation_info.rate_prctile(unit) > 95
                    col = 1;%significant modulation
                    sig_count(1) = sig_count(1)+1;
                elseif fixation_info2.rate_prctile(unit) > 95
                    col = 1;%significant modulation
                    sig_count(2) = sig_count(2)+1;
                else
                    col = 2;%not saccade modulated
                     sig_count(3) = sig_count(3)+1;
                     continue
                end
                
                
                
                fixation_firing = fixation_locked_firing{unit};
                fixation_firing(nansum(fixation_firing,2) == 0,:) = [];%remove fixations without spikes
                
                
                [firing_rate,~]= nandens(fixation_firing,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                firing_rate = firing_rate-nanmean(firing_rate);
                fr = firing_rate;
                if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                    firing_rate = firing_rate/max(firing_rate);
                else%normalize to min %could have some neurons that only show supression
                    firing_rate = firing_rate/min(firing_rate);
                end
                
                avg_fixation_firing_rates{col} =  [avg_fixation_firing_rates{col}; firing_rate];
                fr = fr/max(fr);
                peak_fixation_firing_rates{col} = [peak_fixation_firing_rates{col}; fr];
                
            end
        end
    end
end
%% All Cells
figure
[m,i] = max(avg_fixation_firing_rates{1}(:,:),[],2);
[mm,ii] = sort(i);
imagesc([-500:499],[1:size(avg_fixation_firing_rates{1},1)],avg_fixation_firing_rates{1}(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(avg_fixation_firing_rates{1},1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
%%
% figure
% [m,i] = max(peak_fixation_firing_rates{1},[],2);
% [mm,ii] = sort(i);
% imagesc([-500:499],[1:size(peak_fixation_firing_rates{1},1)],peak_fixation_firing_rates{1}(ii,:))
% colormap('jet')
% hold on
% plot([0 0],[1 size(peak_fixation_firing_rates{1},1)],'w--');
% hold off
% xlabel('Time from Fixation Start')
% ylabel('Neuron #')
%%
[U,S,V] = pca(avg_fixation_firing_rates{1});
T = kmeans(U,4);
t = -500:499;
figure
hold all
for i = 1:max(T)
    if sum(T == i) > 1
        plot(t,mean(avg_fixation_firing_rates{1}(T == i,:)))
    else
        plot(t,avg_fixation_firing_rates{1}(T == i,:));
    end
end
plot([0 0],[-1 1],'k--')
plot([-500 500],[0 0],'k')
xlabel('Time from Fixation Start')
ylabel('Normalized Firing')
%%

t = -500:499;
figure
plot(t,mean(avg_fixation_firing_rates{1}))
hold on
plot([0 0],[-1 1],'k--')
plot([-500 500],[0 0],'k')
hold off
xlabel('Time from Fixation Start')
ylabel('Normalized Firing')
title('Average across All Neurons')