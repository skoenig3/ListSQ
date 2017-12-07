%adapted from various code to determine if neurons are visual, motor, or
%visual-motor responses during cross hair period
set(0,'DefaultFigureVisible','OFF');
clar %clear,clc

min_rt = 75;
max_rt = 400;

numshuffs = 10000;

%where to store spatial analysis figure copies
summary_directory = 'C:\Users\seth.koenig\Desktop\Significant Units\Visual Response\';
if ~isdir(summary_directory)
    mkdir(summary_directory)
end

task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size

% time_locked_firing since ITI period is defined by 2 events (15 & 16)
event_codes = [35 8 23 24 23];
trial_start_code = 15;
task = 'ListSQ';
min_blks = 2;

%---Misc. Parameters---%
monkey_all_unit_count = zeros(2,2);%row 1 place row 2 non-place, column by monkey
all_multi_unit_count = zeros(1,2); %all_multi_unit_countunits
all_unit_names = {}; %place cell unit names
all_monkeys = []; %1s and 2s
AP_location = []; %AP location of recorded place cell

%---Unit Test Significance---%
visual_response_stats_short = []; %short 1 second window, 1 for sig, 0 for not sig
visual_response_stats_long = []; %long 5 second window, 1 for sig, 0 for not sig
visual_response_stats_short_sliding = []; %sliding window analysis for short 1 second window, 1 for sig, 0 for not sig
visual_response_stats_long_sliding = []; %sliding window analysis for long 5 second window, 1 for sig, 0 for not sig
spatialness = []; %1 for place cells, 0 for non place cells
fixation_on_cross_status = []; %1 for sig, 0 for not sig
image_off_status = [];%1 for sig, 0 for not sig
memory_short = [];%significant indeces within the first 1 second window
memory_long = [];%significant indeces for 5 second window

%---Firing Rate Curves---%
cross_fixation_firing_rates = [];%fixation on across
image_on_firing_rates = []; %image on 1 second window
long_image_on_firing_rates = [];%image on 5 second window
image_off_firing_rates = []; %image off

sig_event = [];
sig_fix = [];
vizmotor_status = [];

monkeys = {'Vivian','Tobii'};
figure_dir = {};
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
            'stability_attribute','cfg','hdr','data');
        
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
        else
            if num_units ~= 0
                error('Where is this file')
            end
            continue
        end
        
        %get trials in which spikes are valid i.e. trials for which neuron is
        %stable for
        valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
        if all(isnan(valid_trials(1,:)))
            disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
            return %since no valid units skip analysis
        end
        
        %get important task specific information
        [itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,23);
        
        
        for unit = 1:num_units
            if ~isempty(time_lock_firing{unit,1})
                
                if epoch_data.temporalstability_prctile(unit,1) > 95 && epoch_data.rate_prctile(unit,1) > 95
                    item_appearance = 1;
                else
                    item_appearance = 0;
                end
                if epoch_data.temporalstability_prctile(unit,2) > 95 && epoch_data.rate_prctile(unit,2) > 95
                    fixation_on_item = 1;
                else
                    fixation_on_item = 0;
                end
                
                sig_event = [sig_event item_appearance];
                sig_fix = [sig_fix fixation_on_item];
                
                if item_appearance == 1 || fixation_on_item == 1
                    
                    fixation_rt = NaN(1,192);
                    trial_number= NaN(1,192);
                    fixation_aligned = NaN(192,1150);
                    item_aligned = NaN(192,1150);
                    for t = 1:num_trials
                        if any(cfg.trl(t).allval == event_codes(3)); %in which image was displayed or 1st item in sequence was displayed
                            if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(end) %then sequence trial
                                continue % go to the next trial
                            end
                            
                            
                            img_index = find(cfg.trl(t).cnd == img_cnd); %image index
                            if any(isnan(which_img(img_index)))
                                continue
                            end
                            
                            trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code); %start of ITI
                            trial_end = cfg.trl(t).alltim(end)-trial_start; %end of trial after image off
                            crosson = cfg.trl(t).alltim(cfg.trl(t).allval == 35)-trial_start; %cross hair on
                            fix_cross = cfg.trl(t).alltim(cfg.trl(t).allval == 8)-trial_start; %fixation on cross hair according to cortex
                            fix_cross = fix_cross(1);%since fixation on image also counts as event 8
                            imgon = cfg.trl(t).alltim(cfg.trl(t).allval == 23)-trial_start; %when image turns on
                            imgon2 = imgon;
                            fixation_rt(img_index) = fix_cross-crosson; %reaction time to fixate on cross
                            trial_number(img_index) = t; %trial #
             
                            spikes = find(data(unit).values{t}); %spike trains for this trial
                            
                            twin1 = 400;
                            if imgon-fix_cross > 750;
                                imgon = fix_cross+750;
                            end
                            %---fix cross locked firing---%
                            %variable duration event
                            event_spikes = spikes(spikes > fix_cross-twin1 & spikes <= imgon)-fix_cross+twin1;
                            timevec = [fix_cross-twin1+1 imgon]-fix_cross+twin1;
                            tempvec = zeros(1,timevec(2)-timevec(1)+1);
                            tempvec(event_spikes) = 1;
                            fixation_aligned(img_index,timevec(1):timevec(2)) = tempvec;
                            
                            %---cross on locked firing---%
                            if imgon2-crosson > 950;
                                imgon = crosson+950;
                            else
                                imgon = imgon2;
                            end
                            twin1 = 200;
                            event_spikes = spikes(spikes > crosson-twin1 & spikes <= imgon)-crosson+twin1;
                            timevec = [crosson-twin1+1 imgon]-crosson+twin1;
                            tempvec = zeros(1,timevec(2)-timevec(1)+1);
                            tempvec(event_spikes) = 1;
                            item_aligned(img_index,timevec(1):timevec(2)) = tempvec;
                            
                        end
                    end
                    %%
                    %clean up ultra fast and slow trials
                    too_fast = find(fixation_rt < min_rt);
                    fixation_rt(too_fast) = [];
                    item_aligned(too_fast,:) = [];
                    fixation_aligned(too_fast,:) = [];
                    
                    too_slow = find(fixation_rt > max_rt);
                    fixation_rt(too_slow) = [];
                    item_aligned(too_slow,:) = [];
                    fixation_aligned(too_slow,:) = [];
                    
                    
                    fixation_rt = laundry(fixation_rt);
                    item_aligned = laundry(item_aligned);
                    fixation_aligned = laundry(fixation_aligned);

                    
                    %%
                    [~,ranked_rt] = sort(fixation_rt);
                    
                    sorted_fixation_aligned = fixation_aligned(ranked_rt,:);
                    enough_samples = find(sum(~isnan(sorted_fixation_aligned)) > 20); %variable event length so cut time points without enough samples
                    sorted_fixation_aligned = sorted_fixation_aligned(:,enough_samples);
                    
                    sorted_item_aligned = item_aligned(ranked_rt,:);
                    enough_samples = find(sum(~isnan(sorted_item_aligned)) > 20); %variable event length so cut time points without enough samples
                    sorted_item_aligned = sorted_item_aligned(:,enough_samples);
                    %%
                    [~,all_fix_curves] =  nandens(sorted_fixation_aligned,smval,'gauss',Fs,'nanflt');
                    [~,all_event_curves] = nandens(sorted_item_aligned,smval,'gauss',Fs,'nanflt');

                    count_40pct = round(size(all_fix_curves,1)*.4);
                    fix_c1 = nanmean(all_fix_curves(1:count_40pct,:));
                    fix_c2 = nanmean(all_fix_curves(size(all_fix_curves,1)-count_40pct+1:end,:));
                    observed_fix_corr = corr(fix_c1(:),fix_c2(:),'row','pairwise','type','Spearman');
             
                    event_c1 = nanmean(all_event_curves(1:count_40pct,:));
                    event_c2 = nanmean(all_event_curves(size(all_event_curves,1)-count_40pct+1:end,:));
                    observed_event_corr = corr(event_c1(:),event_c2(:),'row','pairwise','type','Spearman');
                    
                    %%
                    observed_corr_diff = observed_fix_corr-observed_event_corr;
                    shuffled_diff = NaN(1,numshuffs);
                    trial_count = size(all_fix_curves,1);
                    parfor shuff = 1:numshuffs;
                        ind = randperm(trial_count);
                        shuff_fix = all_fix_curves(ind,:);
                        shuff_event = all_event_curves(ind,:); 
                        
                        event_c1 = nanmean(shuff_event(1:count_40pct,:));
                        event_c2 = nanmean(shuff_event(size(all_event_curves,1)-count_40pct+1:end,:));
                        shuff_event_corr = corr(event_c1(:),event_c2(:),'row','pairwise','type','Spearman');
                        
                        fix_c1 = nanmean(shuff_fix(1:count_40pct,:));
                        fix_c2 = nanmean(shuff_fix(size(all_fix_curves,1)-count_40pct+1:end,:));
                        shuff_fix_corr = corr(fix_c1(:),fix_c2(:),'row','pairwise','type','Spearman');
                        
                        
                        shuffled_diff(shuff) = shuff_fix_corr-shuff_event_corr;
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

                    %%
                    figure
                    twin1 = 200;
                    subplot(2,2,1)
                    [trial,time] = find(sorted_item_aligned == 1);
                    if ~isempty(trial)
                        plot(time-twin1,trial,'.k')
                        xlim([-200 950])
                        ylim([0 max(trial)+1]);
                    end
                    set(gca,'Xtick',[-200 0 200 600])
                    ylabel('Trial #')
                    xlabel('Time from Cross On (ms)')
                    box off
                    title(['\rho_{rt} = ' num2str(observed_event_corr,2)])
                    
                    twin1 = 400;
                    subplot(2,2,2)
                    [trial,time] = find(sorted_fixation_aligned == 1);
                    if ~isempty(trial)
                        plot(time-twin1,trial,'.k')
                        xlim([-400 750])
                        ylim([0 max(trial)+1]);
                    end
                    set(gca,'Xtick',[-200 0 200 600])
                    ylabel('Trial #')
                    xlabel('Time from Fixation Start (ms)')
                    box off
                    title(['\rho_{rt} = ' num2str(observed_fix_corr,2)])

                    subplot(2,2,3)
                    hold on
                    dofill2([1:size(sorted_item_aligned,2)]-200,sorted_item_aligned,'black',1,smval); %smoothed cross on aligned curve
                    dofill2([1:size(sorted_fixation_aligned,2)]-200,sorted_fixation_aligned,'green',1,smval); %smoothed cross on aligned curve
                    hold off
                    xlabel('Event Onset (ms)')
                    ylabel('Firing Rate (Hz)')
                    xlim([-200 750])
                    legend('Item On','Fixation Start')
                    
                    subplot(2,2,4)
                    hist(shuffled_diff,25)
                    hold on
                    yl = ylim;
                    plot([observed_corr_diff observed_corr_diff],[0 yl(2)],'r--')
                    hold off
                    box off
                    title(['\Delta \rho_{rt} = ' num2str(observed_fix_corr-observed_event_corr,2)...
                        ' (' num2str(100*sum(observed_corr_diff > shuffled_diff)/numshuffs,3) '%)'])
                    
                    subtitle(title_str)
                    
                   save_and_close_fig('C:\Users\seth.koenig\Desktop\Test_VizMotor\',[task_file(1:8) '_' unit_stats{1,unit}])
                    %%
                end
            end
        end
    end
end
set(0,'DefaultFigureVisible','ON');