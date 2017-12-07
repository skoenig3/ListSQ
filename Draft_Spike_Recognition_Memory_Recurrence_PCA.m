% Population Visual Response
% Code below creates population summary for Significnat Visually Reponsive Neurons
% Written by Seth Konig  written Seth Konig 9/1/16, updated 1/16/2017
% Code does the following
% 1) Summarizes visual responses to images for short and long windows
% 2) Determines if neurons may be sequentially organized in time
% 3) Determines whether place cells are also visually responsive
% 4) Determines if visually responsive neurons are also modulated by novel/repeat images
% 5) Tracks AP location, unit counts, and which monkey (not currently used)
% 6) Copies relevant figures to summary directory

%Code rechecked by SDK on 1/16/2017

clar %clear,clc
set(0,'DefaultFigureVisible','OFF');

task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 500; %vertical image size

% time_locked_firing since ITI period is defined by 2 events (15 & 16)
image_on_code = 23;
image_off_code = 24;
reward_code = 3;
trial_start_code = 15; %ITI start

minimum_firing_rate = 2;%Hz so can see up and downs in firnig rates based on memory
minimum_number_of_images = 24;%number of novel repeat pairs to do analysis

%---Misc. Parameters---%
monkey_all_unit_count = zeros(2,2);%row 1 place row 2 non-place, column by monkey
all_multi_unit_count = zeros(1,2); %all_multi_unit_countunits
all_unit_names = {}; %place cell unit names
all_monkeys = []; %1s and 2s
AP_location = []; %AP location of recorded place cell

%---Unit Test Significance---%
monkeys = {'Vivian','Tobii'};
figure_dir = {};
unit_count = 0;
image_count = [];
low_high_mem_modulation = NaN(2,347);
for monk =2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir{1} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir{2} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    end
    
    for sess = 1:length(session_data)
        %read in task file data
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        %load task file data
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials',...
            'stability_attribute','cfg','hdr','data','fixationstats');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        if num_units == 0
            continue
        end
        
        num_trials = length(cfg.trl); %number of trials
        valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
        if all(isnan(valid_trials(1,:)))
            continue
        end
        
        disp(task_file(1:8))
        if exist([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat']) %visual response analysis
            load([data_dir task_file(1:8) '-ListSQ-Visual_Response_Memory_results']) %memory visual response analysis
        else
            if num_units ~= 0
                error('Where is this file')
            end
            continue
        end
        
        
        for unit = 1:num_units
            if ~isempty(time_lock_firing{unit,1})
                
                %---for image on short 1 second window---%
                if (epoch_data.rate_prctile(unit,3) > 95) && (epoch_data.temporalstability_prctile(unit,3) > 95)
                    visual_response_stats_short = 1; %signficant
                else
                    visual_response_stats_short = 0; %not signficant
                end

                %---for long image on 5 second period---%
                if (epoch_data.rate_prctile(unit,5) > 95) && (epoch_data.temporalstability_prctile(unit,5) > 95)
                    visual_response_stats_long = 1;%signficant
                else
                    visual_response_stats_long = 0;%not signficant
                end
                
                if visual_response_stats_short == 0 && visual_response_stats_long == 0
                    continue
                end
                
                %---combine sig times together across short and long windows---%
                sig_times = zeros(1,twin1+twin4);%full image window period
                if visual_response_stats_short; %short window
                    sig_times(1:twin1+twin2) = sig_short{unit};
                end
                if visual_response_stats_long %long window
                    sig_times = sig_times+sig_long{unit};
                end
                sig_times(sig_times > 2) = 2; %ensure overlappping points only count once
                
                %remove any time points that include the time before the
                %image appears
                gaps = findgaps(find(sig_times));
                for g = 1:size(gaps,1)
                    gp = gaps(g,:);
                    gp(gp == 0) = [];
                    if any(gp < twin1)
                        if all(gp < twin1)
                            sig_times(gp) = 0; %remove these points
                        else
                            gp(gp < twin1) = [];
                            if length(gp) < smval %if too short then remove
                                sig_times(gp) = 0; %remove these points
                            else
                                sig_times(1:twin1) = 0;
                            end
                        end
                    end
                end
                
%                 time_window = [];
%                 if sum(sig_times) == 0
%                     continue %if not memory modulated then go to next neuron
%                 end
%       
%                 time_window = find(sig_times);    
%                 
%                 if length(time_window) < 50
%                     disp('now')
%                 end
%                     
%                 average_firing_rate = mean(sum(time_lock_firing{unit,5}(:,time_window),2))*1000/length(time_window);
%                 if average_firing_rate < minimum_firing_rate
%                    disp('Firing Rate Too low')
%                    continue
%                 end
                
                if sum(nvr{unit} == 1) < minimum_number_of_images
                   disp('Too few trials')
                   continue
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%---Get Fixation Duration Information---%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fixation_durs = NaN(1,length(trial_number{unit}));
                all_recurrence_measures = NaN(length(trial_number{unit}),5);
                
                for tr = 1:length(trial_number{unit})
                    t = trial_number{unit}(tr); %trial number within cfg struct
                    
                    %---Trial Event Codes/Times---%
                    trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code); %start of ITI
                    trial_end = cfg.trl(t).alltim(end)-trial_start; %end of trial after image off
                    imgon = cfg.trl(t).alltim(cfg.trl(t).allval == image_on_code)-trial_start; %when image turns on
                    imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == image_off_code)-trial_start; %when image turns off
                    
                    %---get fixation/saccaqde information---%
                    fixationtimes = fixationstats{t}.fixationtimes; %fixtaion start and end times
                    fixations = fixationstats{t}.fixations; %fixtaion start and end times
                    
                    %remove fixations started ending after image turned off
                    invalid= find(fixationtimes(2,:) > imgoff);
                    fixationtimes(:,invalid) = [];
                    fixations(:,invalid) = [];

                    %remove fixations that started before image turned on
                    invalid= find(fixationtimes(1,:) < imgon);
                    fixationtimes(:,invalid) = [];
                    
                    fixdurs = diff(fixationtimes)+1; %calculate fixation duration
                    fixation_durs(tr) = nanmean(fixdurs(4:end));
                    
                     [recurrence_rate,~,corm,laminarity,~,forward_trace,~,reverse_trace,~] = ...
                        calculate_auto_recurrence(fixations,40,48);
                    
                    all_recurrence_measures(tr,1) = recurrence_rate;
                    all_recurrence_measures(tr,2) = corm;
                    all_recurrence_measures(tr,3) = laminarity;
                    all_recurrence_measures(tr,4) = forward_trace;
                    all_recurrence_measures(tr,5) = reverse_trace;
 
                end
                
                %%
                %---Determine high versus low memory trials---%
                change_dur = NaN(96,6);
                for img = 1:96
                    ind = find(which_images{unit} == img);
                    if ~isempty(ind)
                        if length(ind) ~= 2
                            error('what')
                        end
                        
                        mean_fix_durs = fixation_durs(ind(1));
                        rep_recurrence_measures = all_recurrence_measures(ind(1),:);
%                         rep_recurrence_measures(:,1) = [];
%                         rep_recurrence_measures(:,3) = [];
%                         
                        nandurs = find(isnan(mean_fix_durs));
                        nanmeasures = find(isnan(sum(rep_recurrence_measures',1)));
                        nanvals = unique([nandurs,nanmeasures]);
                        change_dur(img,1) = mean_fix_durs; 
                        change_dur(img,2:6) = rep_recurrence_measures;
                        
                   
%                         change_dur(img) =  fixation_durs(ind(2))- fixation_durs(ind(1));
%                         rep_durs(img) = fixation_durs(ind(2));
%                         nov_durs(img) = fixation_durs(ind(1));
                    end
                end
                %%
                image_count = [image_count  sum(~isnan(unique(which_images{unit})))];
                
                non_nan_rows = find(sum(isnan(change_dur')) == 0);
                change_dur = change_dur(non_nan_rows,:);
                [PCAU,~,~] = pca(change_dur ,1);    
                
                change_thresh_high = prctile(PCAU,67);
                change_thresh_low = prctile(PCAU,33);
                
                low_change = find(PCAU <= change_thresh_low);
                high_change = find(PCAU >= change_thresh_high);
                
                %---for long image on 5 second period---%
                %firing_rate = time_lock_firing{unit,5};
                
                low_change_firing_nov = [];
                low_change_firing_rep = [];
                for l = 1:length(low_change);
                    ind = find(which_images{unit} == non_nan_rows(low_change(l)) & nvr{unit} == 1);
                    low_change_firing_nov = [low_change_firing_nov; time_lock_firing{unit,5}(ind,:)];
                    
                    ind = find(which_images{unit} == non_nan_rows(low_change(l)) & nvr{unit} == 2);
                    low_change_firing_rep = [low_change_firing_rep; time_lock_firing{unit,5}(ind,:)];
                end
                
                high_change_firing_nov = [];
                high_change_firing_rep = [];
                for l = 1:length(high_change);
                    ind = find(which_images{unit} == non_nan_rows(high_change(l)) & nvr{unit} == 1);
                    high_change_firing_nov = [high_change_firing_nov; time_lock_firing{unit,5}(ind,:)];
                    
                    ind = find(which_images{unit} == non_nan_rows(high_change(l)) & nvr{unit} == 2);
                    high_change_firing_rep = [high_change_firing_rep; time_lock_firing{unit,5}(ind,:)];
                end
             
      %%
                if visual_response_stats_long
                    smval_win = smval;
                    xl = [-twin1 twin2];
                else
                    smval_win = smval2;
                    xl = [-twin1 twin4];
                end
                
                ynov = nandens(time_lock_firing{unit,5}(nvr{unit} == 1,:),smval_win,'gauss',1000,'nanflt');
                yrep = nandens(time_lock_firing{unit,5}(nvr{unit} == 2,:),smval_win,'gauss',1000,'nanflt');
%                 if mean(ynov(time_window)) > mean(yrep(time_window))
%                     posneg = 1;
%                 else
%                     posneg = -1;
%                 end
                
                
                if visual_response_stats_long
                    smval_win = smval;
                    xl = [-twin1 twin2];
                else
                    smval_win = smval2;
                    xl = [-twin1 twin4];
                end
                
                %time_2_firing_rate = 1000/length(time_window); %spike count 2 firing rate based on window size             
%%
                unit_count = unit_count+1;
                ylims = NaN(4,2);
                
                figure
                subplot(2,2,1)
                hold on
                [~,~,~,y1] =dofill(-twin1:twin4-1,time_lock_firing{unit,5}(nvr{unit} == 1,:),'green',1,smval_win);
                [~,~,~,y2] =dofill(-twin1:twin4-1,time_lock_firing{unit,5}(nvr{unit} == 2,:),'black',1,smval_win);
                hold off
                xlim(xl)
                ylabel('Firing Rate (Hz)')
                xlabel('Time from Image On (ms)')
                title(['All Novel & All Repeat, n = ' num2str(image_count(end))])
                
                subplot(2,2,3)
                hold on
                [~,~,~,y3] = dofill(-twin1:twin4-1,low_change_firing_nov,'blue',1,smval_win); %all trials
                [~,~,~,y4] = dofill(-twin1:twin4-1,low_change_firing_rep,'red',1,smval_win); %all trials
                hold off
                ylabel('Firing Rate (Hz)')
                xlabel('Time from Image On (ms)')
                ylims(3,:) = ylim;
                xlim(xl)
                
                if visual_response_stats_long
                    low_high_mem_modulation(1,unit) = 100*sum(abs(y3(twin1:twin2+twin1)-y4(twin1:twin2+twin1)))...
                        /(sum(ynov(twin1:twin2+twin1)+yrep(twin1:twin2+twin1))/2);
                else
                    low_high_mem_modulation(1,unit) =  100*sum(abs(y3(twin1:end)-y4(twin1:end)))...
                        /(sum(ynov(twin1:end)+yrep(twin1:end))/2);
                end
                title(['Low: ' num2str(low_high_mem_modulation(1,unit),2)])

                subplot(2,2,4)
                hold on
                [~,~,~,y5] = dofill(-twin1:twin4-1,high_change_firing_nov,'blue',1,smval_win); %all trials
                [~,~,~,y6] = dofill(-twin1:twin4-1,high_change_firing_rep,'red',1,smval_win); %all trials
                hold off
                ylabel('Firing Rate (Hz)')
                xlabel('Time from Image On (ms)')
                ylims(4,:) = ylim;
                xlim(xl)
                ymin = min(ylims(:,1));
                ymin(ymin < 0) = 0;
                ymax = max(ylims(:,2));
                
                if visual_response_stats_long
                    low_high_mem_modulation(2,unit) = 100*sum(abs(y5(twin1:twin2+twin1)-y6(twin1:twin2+twin1)))...
                        /(sum(ynov(twin1:twin2+twin1)+yrep(twin1:twin2+twin1))/2);
                else
                    low_high_mem_modulation(2,unit) =  100*sum(abs(y5(twin1:end)-y6(twin1:end)))...
                        /(sum(ynov(twin1:end)+yrep(twin1:end))/2);
                end
                title(['High: ' num2str(low_high_mem_modulation(2,unit),2)])
                
                for sb = [1 3 4]
                    subplot(2,2,sb)
%                     hold on
% %                   h = fill([time_window(1) time_window(end) time_window(end) time_window(1) time_window(1)]-twin1,...
% %                         [ymin ymin ymax ymax ymin],'k');
%                     uistack(h,'bottom')
%                     set(h,'facealpha',.25,'EdgeColor','None')
%                     hold off
                    ylim([ymin ymax])
                end
                
                subtitle([task_file(1:8) ' ' unit_stats{1,unit}]);
                %%
                save_and_close_fig([ pwd 'Temp\'],[task_file(1:8) '_' unit_stats{1,unit}]);
                
            end
        end
    end
end
%%
p1 = signrank(high_vs_low_recognition_change_durs(1,:),high_vs_low_recognition_change_durs(2,:))
p2 = signrank(high_vs_low_recognition_repeat_durs(1,:),high_vs_low_recognition_repeat_durs(2,:))
p3 = signrank(high_vs_low_recognition_novel_durs(1,:),high_vs_low_recognition_novel_durs(2,:))
%%

figure
high_vs_low_recognition_change_durs(high_vs_low_recognition_change_durs < 0) = 0;%for ploting
high_vs_low_recognition_change_durs(high_vs_low_recognition_change_durs > 10) = 10;%for ploting

plot(high_vs_low_recognition_change_durs(1,:),high_vs_low_recognition_change_durs(2,:),'.')
hold on
plot([0 10],[0 10],'k--')
xlabel('Low Recognition: Change (Hz)')
ylabel('HIgh Recognition: Change (Hz)')
box off

%%

figure
high_vs_low_recognition_repeat_durs(high_vs_low_recognition_repeat_durs < 0) = 0;%for ploting
high_vs_low_recognition_repeat_durs(high_vs_low_recognition_repeat_durs > 10) = 10;%for ploting

plot(high_vs_low_recognition_repeat_durs(1,:),high_vs_low_recognition_repeat_durs(2,:),'.')
hold on
plot([0 10],[0 10],'k--')
xlabel('Low Recognition: Change (Hz)')
ylabel('HIgh Recognition: Change (Hz)')
box off
