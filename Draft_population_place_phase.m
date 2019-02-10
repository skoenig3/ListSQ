% Code below creates population summary for Significnat Place Cells
% Written by Seth Konig
% Code does the following
% 1) Summarizes spatial correlations for place cell and non-place cells
% 2) Calculates fixation aligned firing rate curves for place cells
% 3) Calculates eye coverage and place field coverage for place cells
% 4) Tracks AP location, unit counts, and which monkey (not currently used)
% 5) Contextual differences between list and sequence task
% 6) Copies relevant figures for place cells to summary directory

%Code rechecked by SDK on 1/5/2017
%Small bug fixed on 2/6/18 on designating inside and outside sequence
%trials SDK

clar %clear,clc


min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect

task = 'ListSQ';
min_blks = 2; %only anaedlyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size

imageX = 800;
imageY = 600;
img_on_code = 23; %cortex code when image turns on
img_off_code = 24; %cortex code when image turns off
ITIstart_code = 15; %start of ITI/trial
numshuffs = 10000; %number of shuffles for resmampling

min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect

Fs = 1000; %Hz sampling frequency
fixwin = 5;%size of fixation window on each crosshair
smval = 30;%2*std of gaussian kernel so 15 ms standard deviation
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
twin1 = 200;%200 ms before fixation
twin2 = 400;%400 ms after start of fixation
image_on_twin = 500;%how much time to ignore eye movements for to remove strong visual response though some may last longer


%---Misc. Parameters (Mostly Place Cells)---%
monkey_all_unit_count = zeros(2,2);%row 1 place row 2 non-place, column by monkey

%first set for out2in only and 2nd set for all fixations
cell_ind = 1;
skagg_info_raw_warped = NaN(4,109); 

STA_LFP = cell(2,109);
LFP_phase_vals = cell(2,109);

monkeys = {'Vivian','Tobii'};
figure_dir = {};
for monk =2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    end
    
    for sess = 45%1:length(session_data)
        %read in task file data
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,...
            sorting_quality,waveform_count,lfp_quality,comments] = get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        %load task file data
        load([data_dir task_file(1:end-11) '-preprocessed.mat']);
        
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
                
        %load spatial analysis data
        disp(task_file(1:8))
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
            'spatial_info','spike_times','eyepos','binsize','filter_width')
        
        if exist([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat']) %visual response analysis
            load([data_dir task_file(1:8) '-ListSQ-Visual_Response_Memory_results']) %memory visual response analysis
        else
            if num_units ~= 0
                error('Where is this file')
            end
            continue
        end
        
        %load Place Cell Fixation Analysis data
        load([data_dir task_file(1:8) '-Place_Cell_Analysis.mat'])
        if numshuffs < 1000
            error('Should have 1000 shuffles')%make sure file has right number of shuffles
        end
        if smval ~= 30
            error('Smoothing Value (2xStd) does not match expectations!')
        end
        
        if any(spatial_info.shuffled_rate_prctile> 95  & spatial_info.spatialstability_halves_prctile > 95)
        else
            continue
        end
        
        %Save as new variableso can reload later...kind of dumb but that was how it was written
        absolute_fixationstats = fixationstats;
        absolute_cfg = cfg;
        
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
            = get_task_data(session_data{sess},task);
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        clear unit_names
        
        %remove units with too few trials
        %these are the absolute minimum data required to do data analysis
        %set the image duration
        if str2double(task_file(3:8)) < 140805
            imgdur = 7000;
        else
            imgdur = 5000;
        end
        
        %get important task specific information
        [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);
        
        %remove bad LFP channels
        LFPchannels = find_desired_channels(cfg,'LFP');
        bad_channels = [];
        for channel = 1:4
            if cell2mat(strfind(hdr.label,['AD0' num2str(channel)])) %make sure have recorded channel
                if  lfp_quality(channel) == 0; %if it is bad
                    bad_channels = [bad_channels channel];
                end
            end
        end
        LFPchannels(bad_channels) = NaN;
        if all(isnan(LFPchannels))
            continue %no analyses to do
        end
        
        for unit = 1:num_units
            if ~isempty(spatial_info.shuffled_info_rate{unit})&& length(spatial_info.shuffled_info_rate{unit}) < 1000
                error('Should have 1000 shuffles') %make sure file has right number of shuffles
            elseif isempty(spatial_info.shuffled_info_rate{unit})
                continue
            end
            
            if isnan(stats_across_tasks(1,unit))%no peak detected
                continue
            end
            
            if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                    && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                
                place_field_matrix = all_place_field_matrix{unit};
                %note ignores incomplete coverage

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%---Calculate Firing Rate Locked to Fixations---%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                fix_in_out = NaN(1,3000); %firing rate for fixations inside or outside  of place field
                %1) first fixation in: out-> in
                %2) fixation in field but not first: in -> in
                %3) first fixation out of field: in -> out
                %4) fixation out of field but not first: out-> out
                
                spike_count = zeros(1,2);
                STA_LFP{1,cell_ind} = zeros(1,401);
                STA_LFP{2,cell_ind} = zeros(1,401);
                LFP_phase_vals{1,cell_ind} = NaN(28,10000);
                LFP_phase_vals{2,cell_ind} = NaN(28,10000);
                
                fix_locked_firing = NaN(1000,(twin1+twin1)); %spike trains locked to fixations
                fix_locked_firing_warped = NaN(1000,(twin1+twin1)); %spike trains locked to fixations
                
                all_fix_dur = NaN(1,1000);
                all_fix_locked_firing = NaN(1000,(twin1+twin1)); %spike trains locked to fixations
                all_fix_locked_firing_warped = NaN(1000,(twin1+twin1)); %spike trains locked to fixations

                fix_ind = 1; %fixation # so can track in variables above 
                all_fix_ind = 1; %fixation # so can track in variables above 

                fixationstats = absolute_fixationstats; %reload because written over below
                cfg = absolute_cfg; %reload because written over below
                num_trials = length(cfg.trl);%number of trials
                for t = 1:num_trials
                    if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
                        if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                            trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                            imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %when image turns on
                            imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start; %when image turns off
                            
                            %---image info---%
                            img_index = find(cfg.trl(t).cnd == img_cnd); %image index
                            
                            % if monkey isn't paying attention and looked away image presentation
                            % is now longer than imgdur (because of cumulative looking time)
                            % so data isn't probably worth much plus have to cut off somewhere
                            if imgoff-imgon > 1.5*imgdur-1 %cut off trial at 1.5x length of desired looking time
                                imgoff = imgon+1.5*imgdur-1;
                            end
                            imgon = imgon+image_on_twin; %to avoid visual response and strong central bias
                            
                            fixationtimes = fixationstats{t}.fixationtimes; %fixation start and end times
                            saccadetimes = fixationstats{t}.saccadetimes; %saccade start and end times
                            fixations = round(fixationstats{t}.fixations); %mean fixation location
                            xy = fixationstats{t}.XY; %xy eye trace
                            
                            %find fiations and saccades that did not occur during the image period;
                            %should also take care of the 1st fixation on the crosshair
                            
                            %fixation started before image turned on
                            invalid= find(fixationtimes(1,:) < imgon);
                            fixationtimes(:,invalid) = [];
                            fixations(:,invalid) = [];
                            
                            %fixation started after the image turned off and/or firing rate could corrupted by image turning off
                            invalid= find(fixationtimes(1,:) > imgoff-twin2);
                            fixationtimes(:,invalid) = [];
                            fixations(:,invalid) = [];
                            
                            %saccade started before image turned on
                            invalid= find(saccadetimes(1,:) < imgon);
                            saccadetimes(:,invalid) = [];
                            
                            %saccade started after the image turned off and/or firing rate could corrupted by image turning off
                            invalid= find(saccadetimes(1,:) > imgoff-twin2);
                            saccadetimes(:,invalid) = [];
                            
                            %remove fixations that are too short to really look at for analysis
                            fixdur = fixationtimes(2,:)-fixationtimes(1,:)+1;
                            fixationtimes(:,fixdur < min_fix_dur) = [];
                            fixations(:,fixdur < min_fix_dur) = [];
                            
                            img_index = find(img_cnd == cfg.trl(t).cnd);
                            %which image monkey viewed if nan or empty then skip cuz
                            %bad trial
                            if isempty(img_index) || any(isnan(which_img(img_index)))
                                continue
                            end
                            
                            spikes = find(data(unit).values{t}); %spike trains for this trial
                            
                            if sum(~isnan(LFPchannels)) == 1 %chose the only available channel
                                LFPs = data(LFPchannels(~isnan(LFPchannels))).values{t};
                            else
                                unit_channel = str2double(unit_stats{1,unit}(6));
                                if ~isnan(LFPchannels(unit_channel))
                                    LFPs = data(LFPchannels(unit_channel)).values{t};
                                else
                                    first_non_nan = find(~isnan(LFPchannels));
                                    first_non_nan = first_non_nan(1);
                                    LFPs = data(LFPchannels(first_non_nan)).values{t};
                                end
                            end
                            [~,trialphase,wfq,~,~,~,~] = waveletanalysis(LFPs);
                            if length(wfq) ~= 28
                                error('what?');
                            end
                            
                            for f = 2:size(fixations,2)%ignore first fixation not sure where it was/possibly contaminated anyway
                                fixdur = fixationtimes(2,f)-fixationtimes(1,f);
                                prior_sac = find(saccadetimes(2,:) == fixationtimes(1,f)-1);%next fixation should start immediately after
                                if isempty(prior_sac) %no prior saccade so was proabbly looking off screen
                                    continue;
                                end
                                sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %saccade amplitude
                                if sacamp < min_sac_amp %prior saccade is too small so ignore
                                    continue
                                end
                                
                                
                                prior_fix_in_out = NaN;
                                %determine if prior fixation was in or out of place field
                                fixx = fixations(1,f-1);
                                fixx(fixx < 1) = 1;
                                fixx(fixx > imageX) = imageX;
                                fixy = imageY-fixations(2,f-1);
                                fixy(fixy < 1) = 1;
                                fixy(fixy > imageY) = imageY;
                                
                                if place_field_matrix(fixy,fixx) == 1 %then prior fixation was already inside
                                    prior_fix_in_out = 1;
                                else
                                    prior_fix_in_out = 0;
                                end
                                last_fixx = fixx;
                                last_fixy = fixy;
                                
                                %determine if this fixation was in our out of place field
                                fixx = fixations(1,f);
                                fixx(fixx < 1) = 1;
                                fixx(fixx > imageX) = imageX;
                                fixy = imageY-fixations(2,f);
                                fixy(fixy < 1) = 1;
                                fixy(fixy > imageY) = imageY;
                                if place_field_matrix(fixy,fixx) == 1 %then inside
                                    if prior_fix_in_out == 1%prior fixation was inside so in->in
                                        fix_in_out(fix_ind) = 2;
                                    else %out->in
                                        fix_in_out(fix_ind) = 1;
                                    end
                                elseif place_field_matrix(fixy,fixx) == 0 %then inside, NaNs are for locations not analyzed
                                    if prior_fix_in_out == 1%prior fixation was inside so in->out
                                        fix_in_out(fix_ind) = 3;
                                    else %out->out
                                        fix_in_out(fix_ind) = 4;
                                    end
                                else %not a valid fixation location too sparse of coverage to analyze
                                    continue
                                end
                                
                                if fix_in_out(fix_ind) == 1 %out2in
                                    %get firing rate locked to fixation
                                    fixt = fixationtimes(1,f);%start of fixation
                                    fix_spikes = spikes(spikes > fixt-twin1 & spikes <= fixt+twin1)-fixt+twin1;
                                    temp = zeros(1,twin1+twin1);
                                    temp(fix_spikes) = 1;
                                    fix_locked_firing(fix_ind,:) = temp;
                                    
                                    if fixdur == 200
                                        fix_locked_firing_warped(fix_ind,:) = temp;
                                    else
                                        fix_locked_firing_warped(fix_ind,1:twin1) = temp(1:twin1);
                                        fix_spikes2 = spikes(spikes > fixt & spikes <= fixt + fixdur)-fixt;
                                        temp2 = zeros(1,fixdur);
                                        temp2(fix_spikes2) = 1;
                                        resampled = round(linspace(1,fixdur,200));
                                        temp2 = temp2(resampled);
                                        fix_locked_firing_warped(fix_ind,twin1+1:end) = temp2;
                                    end
                                    fix_ind = fix_ind+1;
                                    
                                    fix_spikes = spikes(spikes > fixt-twin1 & spikes <= fixt+twin1);
                                    for fs = 1:length(fix_spikes)
                                       if fix_spikes(fs) < imgoff-200 %%already has to be image on
                                           spike_count(1) = spike_count(1)+1;
                                           STA_LFP{1,cell_ind} = STA_LFP{1,cell_ind} + LFPs(fix_spikes(fs)-200:fix_spikes(fs)+200);
                                           for wf = 1:length(wfq)
                                                 LFP_phase_vals{1,cell_ind}(wf,spike_count(1)) = trialphase(wf,fix_spikes(fs));
                                           end
                                       end
                                    end
                                end
                                    
                                all_fix_dur(all_fix_ind) = fixdur;
                                
                                fixt = fixationtimes(1,f);%start of fixation
                                fix_spikes = spikes(spikes > fixt-twin1 & spikes <= fixt+twin1)-fixt+twin1;
                                temp = zeros(1,twin1+twin1);
                                temp(fix_spikes) = 1;
                                all_fix_locked_firing(all_fix_ind,:) = temp;
                                
                                if fixdur == 200
                                    all_fix_locked_firing_warped(all_fix_ind,:) = temp;
                                else
                                    all_fix_locked_firing_warped(all_fix_ind,1:twin1) = temp(1:twin1);
                                    fix_spikes2 = spikes(spikes > fixt & spikes <= fixt + fixdur)-fixt;
                                    temp2 = zeros(1,fixdur);
                                    temp2(fix_spikes2) = 1;
                                    resampled = round(linspace(1,fixdur,200));
                                    temp2 = temp2(resampled);
                                    all_fix_locked_firing_warped(all_fix_ind,twin1+1:end) = temp2;
                                end
                                
                                all_fix_ind = all_fix_ind+1; %fixation # so can track in variables above
                                
                                
                                fix_spikes = spikes(spikes > fixt-twin1 & spikes <= fixt+twin1);
                                for fs = 1:length(fix_spikes)
                                    if fix_spikes(fs) < imgoff-200 %%already has to be image on
                                        spike_count(2) = spike_count(2)+1;
                                        STA_LFP{2,cell_ind} = STA_LFP{1,cell_ind} + LFPs(fix_spikes(fs)-200:fix_spikes(fs)+200);
                                        for wf = 1:length(wfq)
                                            LFP_phase_vals{2,cell_ind}(wf,spike_count(2)) = trialphase(wf,fix_spikes(fs));
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                STA_LFP{1,cell_ind} = STA_LFP{1,cell_ind}/spike_count(1);
                STA_LFP{2,cell_ind} = STA_LFP{2,cell_ind}/spike_count(2);
                
                %remove excess NaNs
                fix_locked_firing = laundry(fix_locked_firing);
                fix_locked_firing_warped = laundry(fix_locked_firing_warped);
                
                all_fix_locked_firing = laundry(all_fix_locked_firing);
                all_fix_locked_firing_warped = laundry(all_fix_locked_firing_warped);
                all_fix_dur = laundry(all_fix_dur);
               
                LFP_phase_vals{1,cell_ind} = laundry(LFP_phase_vals{1,cell_ind});
                LFP_phase_vals{2,cell_ind} = laundry(LFP_phase_vals{2,cell_ind});
                  
                %%
                figure
                
                subplot(2,3,1)
                [trial,time] = find(fix_locked_firing == 1);
                plot(time-twin1,trial,'.k')
                [trialw,timew] = find(fix_locked_firing_warped == 1);
                hold on
                plot(timew-twin1,trialw+max(trial),'.b')
                hold off
                ylim([0 max(trialw) + max(trial)])
                xlabel('Time from Fixation Start (ms)')
                ylabel('Fixation #')
                title('Rasters for Out2In Only')
                box off
                
                yraw = nandens(fix_locked_firing(:,twin1-50:end),smval,'gauss',Fs,'nanflt');
                ywarped = nandens(fix_locked_firing_warped(:,twin1-50:end),smval,'gauss',Fs,'nanflt');
                                 
                %total_time = 1:2*twin1; 
                total_time = 1:twin1+51; 
                p_x = total_time/sum(total_time);
                lambda_x = yraw; 
                lambda = nansum(nansum(lambda_x.*p_x));
                plogp = lambda_x.*log2(lambda_x/lambda);
                raw_skaggs = nansum(nansum(plogp.*p_x)); %bits/second
                
                 %total_time = 1:2*twin1; 
                total_time = 1:twin1+51; 
                p_x = total_time/sum(total_time);
                lambda_x = ywarped; 
                lambda = nansum(nansum(lambda_x.*p_x));
                plogp = lambda_x.*log2(lambda_x/lambda);
                warped_skaggs = nansum(nansum(plogp.*p_x)); %bits/second
                
                skagg_info_raw_warped(1,cell_ind) = raw_skaggs; 
                skagg_info_raw_warped(2,cell_ind) = warped_skaggs; 
                
                subplot(2,3,2)
                t11 = -twin1:twin1-1;
                hold on
                dofill(t11,fix_locked_firing,'black',1,smval);
                dofill(t11,fix_locked_firing_warped,'blue',1,smval);
                hold off
                ylabel('Firing Rate (Hz)')
                xlabel('Time from Fixation Start (ms)')
                xlim([-twin1 twin1])
                title(['Out2In Only: Raw_{skaggs} = ' num2str(raw_skaggs,3) ', Warped_{skaggs} = ' num2str(warped_skaggs,3)])
                legend('Raw','Warped')
               
                [s,si] = sort(all_fix_dur);
                all_fix_locked_firing = all_fix_locked_firing(si,:);
                [trial,time] = find(all_fix_locked_firing == 1);
                all_fix_locked_firing_warped = all_fix_locked_firing_warped(si,:);
                [trialw,timew] = find(all_fix_locked_firing_warped == 1);
            
                subplot(2,3,4)
                plot(time-twin1,trial,'.k')
                hold on
                plot(timew-twin1,trialw+max(trial),'.b')
                hold off
                ylim([0 max(trialw) + max(trial)])
                xlabel('Time from Fixation Start (ms)')
                ylabel('Fixation #')
                title('Rasters for All Fixations, sorted by fixation duration')
                box off
                
                yraw = nandens(all_fix_locked_firing(:,twin1-50:end),smval,'gauss',Fs,'nanflt');
                ywarped = nandens(all_fix_locked_firing_warped(:,twin1-50:end),smval,'gauss',Fs,'nanflt');
                                 
                %total_time = 1:2*twin1; 
                total_time = 1:twin1+51; 
                p_x = total_time/sum(total_time);
                lambda_x = yraw; 
                lambda = nansum(nansum(lambda_x.*p_x));
                plogp = lambda_x.*log2(lambda_x/lambda);
                raw_skaggs = nansum(nansum(plogp.*p_x)); %bits/second
                
                 %total_time = 1:2*twin1; 
                total_time = 1:twin1+51; 
                p_x = total_time/sum(total_time);
                lambda_x = ywarped; 
                lambda = nansum(nansum(lambda_x.*p_x));
                plogp = lambda_x.*log2(lambda_x/lambda);
                warped_skaggs = nansum(nansum(plogp.*p_x)); %bits/second
                
                skagg_info_raw_warped(3,cell_ind) = raw_skaggs; 
                skagg_info_raw_warped(4,cell_ind) = warped_skaggs; 
                
                subplot(2,3,5)
                t11 = -twin1:twin1-1;
                hold on
                dofill(t11,all_fix_locked_firing,'black',1,smval);
                dofill(t11,all_fix_locked_firing_warped,'blue',1,smval);
                hold off
                ylabel('Firing Rate (Hz)')
                xlabel('Time from Fixation Start (ms)')
                xlim([-twin1 twin1])
                title(['All fixations: Raw_{skaggs} = ' num2str(raw_skaggs,3) ', Warped_{skaggs} = ' num2str(warped_skaggs,3)])
                
                if sqrt(mean(STA_LFP{1,cell_ind}.^2))/2 > sqrt(mean(STA_LFP{2,cell_ind}.^2))
                    scale = sqrt(mean(STA_LFP{1,cell_ind}.^2))/ sqrt(mean(STA_LFP{2,cell_ind}.^2));
                else
                    scale = 1;
                end 
                subplot(2,3,3)
                plot(-200:200,STA_LFP{1,cell_ind},'r')
                hold on
                plot(-200:200,scale*STA_LFP{2,cell_ind},'b')
                hold off
                xlabel('Time from Spike (ms)')
                ylabel('scaled LFP (uV)')
                legend('Out2In','All')
                box off
                
                subplot(2,3,6)
                circ_var = NaN(2,length(wfq));
                for wf = 1:length(wfq)
                    circ_var(1,wf) = circ_r(LFP_phase_vals{1,cell_ind}(wf,:));
                    circ_var(2,wf) = circ_r(LFP_phase_vals{2,cell_ind}(wf,:));
                end
                plot(wfq,circ_var(1,:),'r');
                hold on
                plot(wfq,circ_var(2,:),'b');
                hold off
                xlabel('Frequency (Hz)')
                ylabel('Circular Variance')
                box off
               
                subtitle([task_file(1:end-11) '_' unit_stats{1,unit}])
                %%
                save_and_close_fig('C:\Users\seth.koenig\Desktop\New folder2\',[task_file(1:end-11) '_' unit_stats{1,unit} '_view_cell__time_warping']);
                
                %%
                cell_ind = cell_ind +1;
            else
                continue
            end
        end
    end
end
%%
[~,p12] = ttest(skagg_info_raw_warped(1,:),skagg_info_raw_warped(2,:));
[~,p34] = ttest(skagg_info_raw_warped(3,:),skagg_info_raw_warped(4,:));

change = 100*(skagg_info_raw_warped(1,:)-skagg_info_raw_warped(2,:))./skagg_info_raw_warped(1,:);
change(change < -500) = -500;
figure
subplot(1,2,1)
hist(change,25)
box off
xlabel('% Change in Skagg Info')
ylabel('View Cell Count')
title(['Out2in Only: mean change = ' num2str(100*nanmean((skagg_info_raw_warped(1,:)-skagg_info_raw_warped(2,:))./skagg_info_raw_warped(1,:)),3) ...
    ' (paired test, p = ' num2str(p12,3) ')'])

subplot(1,2,2)
hist(100*(skagg_info_raw_warped(3,:)-skagg_info_raw_warped(4,:))./skagg_info_raw_warped(3,:),25)
box off
xlabel('% Change in Skagg Info')
ylabel('View Cell Count')
title(['All Fixations: mean change = ' num2str(100*nanmean((skagg_info_raw_warped(3,:)-skagg_info_raw_warped(4,:))./skagg_info_raw_warped(3,:)),3) ...
    ' (paired test, p = ' num2str(p34,3) ')'])

%%
avg_STA_LFP = zeros(2,401);
for c = 1:cell_ind-1
    avg_STA_LFP(1,:) =  avg_STA_LFP(1,:) + STA_LFP{1,c};
    avg_STA_LFP(2,:) =  avg_STA_LFP(2,:) + STA_LFP{2,c};
end
avg_STA_LFP = avg_STA_LFP/(cell_ind-1);


    
if sqrt(mean(avg_STA_LFP(1,:).^2))/2 > sqrt(mean(avg_STA_LFP(2,:).^2))
    scale = sqrt(mean(avg_STA_LFP(1,:).^2))/ sqrt(mean(avg_STA_LFP(2,:).^2));
else
    scale = 1;
end

figure
subplot(1,2,1)
plot(avg_STA_LFP(1,:),'r')
hold on
plot(scale*avg_STA_LFP(2,:),'b')
hold off
xlabel('Time from Spike (ms)')
ylabel('scaled LFP (uV)')
legend('Out2In','All')
box off


all_circ_var1 = NaN(length(wfq),cell_ind-1);
all_circ_var2 = NaN(length(wfq),cell_ind-1);
for wf = 1:length(wfq)
    for c = 1:cell_ind-1
        all_circ_var1(wf,c) = circ_r(LFP_phase_vals{1,c}(wf,:));
        all_circ_var2(wf,c) = circ_r(LFP_phase_vals{2,c}(wf,:));
    end
end
%%
subplot(1,2,2)
plot(wfq,mean(all_circ_var1'),'r');
hold on
plot(wfq,mean(all_circ_var2'),'b');
hold off
xlabel('Frequency (Hz)')
ylabel('Circular Variance')
box off
