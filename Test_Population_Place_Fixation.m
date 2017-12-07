% Population Fixation Modulation
% written Seth Konig 8/12/16
%
% code imports data from List_fixation_AnalysisV results and looks at how
% modulated the whole population is by fixations.
% Code does the following
% 1) Summarizes fixation modulation for for place cell and non-place cells
% 2) Tracks AP location, unit counts, and which monkey (not currently used)
% 3) Asks how modulate whole population of all units is
% 6) Copies relevant figures for eye movment modulated cells to summary directory

clar %clear, clc
task = 'ListSQ';
Fs = 1000;%Hz
min_num_fix = 250; %at least 250 fixatoins with a certain duration to analyze for time limited to fixation duration
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)

imageX = 800;
imageY = 600;
img_on_code = 23; %cortex code when image turns on
img_off_code = 24; %cortex code when image turns off
ITIstart_code = 15; %start of ITI/trial
numshuffs = 10000; %number of shuffles for resmampling

min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect
min_sac_dur = 44;

Fs = 1000; %Hz sampling frequency
fixwin = 5;%size of fixation window on each crosshair
smval = 30;%2*std of gaussian kernel so 15 ms standard deviation
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
twin1 = 200;%200 ms before fixation
twin2 = 400;%400 ms after start of fixation
image_on_twin = 500;%how much time to ignore eye movements for to remove strong visual response though some may last longer

all_fields = [];
all_field_peaks = [];
all_cropped_fields = [];
all_cropped_peaks = [];
all_cropped_fields2 = [];
all_cropped_peaks2 = [];
all_warped_fields = [];
all_warped_peaks = [];

monkeys = {'Vivian','Tobii'};
figure_dir = {};
for monk = 2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for sess data---%%%
    %only need to run when somethings changed or sesss have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir{1} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_session_Unit_Data_Vivian.mat'])
        
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
        load([data_dir 'Across_session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    end
    
    for sess = 1:length(session_data)
        
        %read in task file data
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import task and unit data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(task_file)
            warning('No file could be found for specificed task. Exiting function...')
            continue
        end
        try
            load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat']);
        catch
            disp('Spatial Analysis File not found continueing')
            continue
        end
        
        %load task file data
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials',...
            'stability_attribute','cfg','hdr','data','fixationstats');
        unit_names = unit_names.name;
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        clear unit_names
        
        if num_units == 0; %if no units exit function
            disp([task_file(1:8) ': no units could be found. Exiting function...'])
            continue;
        end
        
        valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
        if all(isnan(valid_trials(1,:)))
            disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
            continue %since no valid units skip analysis
        end
        
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
        [which_img,~,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);
        
        %filter parameters
        H = define_spatial_filter(filter_width);
        
        for unit = 1:num_units
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%---Determine if Unit Passes 95% for both Skaggs and Stability--%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isnan(spatial_info.shuffled_rate_prctile(unit))
                continue %spatial analysis wasn't run on this unit
            elseif (spatial_info.shuffled_rate_prctile(unit) < 95) || (spatial_info.spatialstability_halves_prctile(unit) < 95)
                continue %unit not spatial
            end
            
            disp(num2str(spatial_info.shuffled_rate_prctile(unit)))
            disp(spatial_info.spatialstability_halves_prctile(unit))
            if isempty(eyepos{unit})
                continue %no data for this neuron
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Determine List Firing rate Locked to Fixations Inside vs Outside Field---%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate Firing Rate Map---%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            firing_rate_map = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all'); %get firing rate map
            place_field_matrix = define_place_field(firing_rate_map,imageX,imageY); %place field location
            all_place_field_matrix{unit} = place_field_matrix; %save place field location
            area(unit) = 100*nansum(nansum(place_field_matrix))/(imageX*imageY); %calculate area as % of screen size
            %note ignores incomplete coverage
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Calculate Firing Rate Locked to Fixations---%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fix_in_out = NaN(1,3000); %firing rate for fixations inside or outside  of place field
            %1) first fixation in: out-> in
            %2) fixation in field but not first: in -> in
            %3) first fixation out of field: in -> out
            %4) fixation out of field but not first: out-> out
            fix_locked_firing = NaN(3000,(twin1+twin2)); %spike trains locked to fixations
            
            fix_info = NaN(3000,2);
            %1) prior fix start (i.e. sacdur+fixdur)
            %2) current fix dur
            
            fix_ind = 1; %fixation # so can track in variables above
            num_trials = length(cfg.trl);%number of trials
            for t = 1:num_trials
                if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
                    if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                        imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %when image turns on
                        imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start; %when image turns off
                        
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
                        for f = 2:size(fixations,2)%ignore first fixation not sure where it was/possibly contaminated anyway
                            prior_sac = find(saccadetimes(2,:) == fixationtimes(1,f)-1);%next fixation should start immediately after
                            if isempty(prior_sac) %no prior saccade so was proabbly looking off screen
                                continue;
                            end
                            sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %saccade amplitude
                            if sacamp < min_sac_amp %prior saccade is too small so ignore
                                continue
                            end
                            
                            prior_fix = find(saccadetimes(1,prior_sac)-1 == fixationtimes(2,:));
                            if isempty(prior_fix)
                                continue
                            end
                            
                            prior_fix_in_out = [];
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
                            
                            %get firing rate locked to fixation
                            fixt = fixationtimes(1,f);%start of fixation
                            fix_spikes = spikes(spikes > fixt-twin1 & spikes <= fixt+twin2)-fixt+twin1;
                            temp = zeros(1,twin1+twin2);
                            temp(fix_spikes) = 1;
                            fix_locked_firing(fix_ind,:) = temp;
                            
                            fix_info(fix_ind,1) = fixationtimes(1,f)-fixationtimes(1,prior_fix);
                            fix_info(fix_ind,2) = diff(fixationtimes(:,f))+1;
                            
                            fix_ind = fix_ind+1;
                        end
                    end
                end
            end
            
            %%
            out_2_in = find(fix_in_out == 1);
            fixdurs = fix_info(out_2_in,1);
            these_spikes = fix_locked_firing(out_2_in,:);
            
            median_dur = median(fixdurs);
            
            warped_spikes = NaN(size(these_spikes,1),twin1+twin2);
            for t = 1:size(these_spikes,1)
                this_dur = fixdurs(t);
                if this_dur == median_dur
                    warped_spikes(t,:) = these_spikes(t,:);
                elseif this_dur > median_dur
                    ind = round(linspace(1,this_dur,median_dur));
                    warped_spikes(t,1:twin1) = these_spikes(t,1:twin1);
                    if this_dur < 400
                        warped_spikes(t,twin1:twin1+length(ind)-1) = these_spikes(t,twin1+ind);
                    else
                        ind(ind > 400) = [];
                        warped_spikes(t,twin1:twin1+length(ind)-1) = these_spikes(t,twin1+ind);
                    end
                    if this_dur < 400
                        warped_spikes(t,(twin1+median_dur:end)) = these_spikes(t,median_dur:ind(end)+twin2-ind(end));
                    end
                else% this_dur < median_dur
                    ind = round(linspace(1,this_dur,median_dur));
                    ind2 = round(linspace(median_dur+1,twin2,twin2-median_dur+1));
                    
                    warped_spikes(t,1:twin1) = these_spikes(t,1:twin1);
                    warped_spikes(t,twin1:twin1+length(ind)-1) = these_spikes(t,twin1+ind);
                    warped_spikes(t,(twin1+median_dur):end) = these_spikes(t,twin1+ind2);
                end
            end
            %%
            out_2_in = find(fix_in_out == 1);
            median_back_dur = round(nanmedian(fix_info(out_2_in,1)));%how far back can we go in time
            if median_back_dur < min_sac_dur+min_fix_dur;
                median_back_dur =  min_sac_dur+min_fix_dur;
            elseif median_back_dur > twin1
                median_back_dur = twin1;
            end
            
            median_forward_dur = round(nanmedian(fix_info(out_2_in,2)));%how far back can we go in time
            if median_forward_dur < min_fix_dur
                median_forward_dur = min_fix_dur;
            elseif median_forward_dur > twin1;
                median_forward_dur = twin1;
            end
            
            these_spikes = fix_locked_firing(out_2_in,:);
            these_fix_info = fix_info(out_2_in,:);
            
            
            start_window = twin1-median_back_dur;%how early before saccade can you look for direction tuning
            end_window = twin2+median_forward_dur;%how late after saccade can you look for direction tuning,
            
            cropped_spikes = these_spikes;
            cropped_spikes(these_fix_info(:,1) < median_back_dur,:) = NaN;
            cropped_spikes(these_fix_info(:,2) < median_forward_dur,:) = NaN;
            %%
            cropped_spikes2 = these_spikes;
            for f = 1:size(cropped_spikes2)
                if these_fix_info(f,1) < median_back_dur
                    cropped_spikes2(f,1:twin1-these_fix_info(f,1)-1) = NaN;
                end
                if these_fix_info(f,2) < median_forward_dur
                    cropped_spikes2(f,twin1+these_fix_info(f,2)+1:end) = NaN;
                end
            end
            
%             t = -twin1:twin2-1;
%             figure
%             subplot(2,2,1)
%             hold on
%             dofill(t,these_spikes,'red',1,smval);%out-> in
%             dofill(t,cropped_spikes,'blue',1,smval);%out-> in cropped for longer fixation durations
%             dofill2(t,cropped_spikes2,'green',1,smval);%out-> in cropped for longer fixation durations
%             hold off
%             xlabel('Time from Fixation Start (ms)')
%             ylabel('Firing Rate (Hz)')
%             legend('Raw','Cropped to median of all','Cropped to individual')
%             
%             subplot(2,2,3)
%             [trial,time] = find(these_spikes == 1);
%             plot(time-twin1,trial,'r.')
%             box off
%             xlabel('Time from Fixation Start (ms)')
%             ylabel('Occurence')
%             title('Cropped to Median of All')
%             
%             subplot(2,2,2)
%             hold on
%             dofill(t,these_spikes,'red',1,smval);%out-> in
%             dofill(t,cropped_spikes,'blue',1,smval);%out-> in cropped for longer fixation durations
%             dofill2(t,warped_spikes,'black',1,smval);%out-> in cropped for longer fixation durations
%             hold off
%             xlabel('Time from Fixation Start (ms)')
%             ylabel('Firing Rate (Hz)')
%             legend('Raw','Cropped to median of all','Warped to Median')
%             
%             subplot(2,2,4)
%             [trial,time] = find(warped_spikes == 1);
%             plot(time-twin1,trial,'k.')
%             box off
%             xlabel('Time from Fixation Start (ms)')
%             ylabel('Occurence')
%             title('Cropped to Median of All')
            
            
            
            %%
            avg_these = nandens(these_spikes,smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
            avg_these = avg_these-nanmean(avg_these(1:twin1)); %remove base line
            avg_these = avg_these/max(avg_these);
            [PKS,LOCS]= findpeaks(avg_these,'MinPeakWidth',smval); %find peaks > 2*std in width
            if ~isempty(LOCS)
                %remove peaks less than 2/3 the max
                LOCS(PKS < 0.66) = [];
            end
            if ~isempty(LOCS)
                all_fields = [all_fields; avg_these./avg_these(LOCS(1))];
                all_field_peaks = [all_field_peaks LOCS(1)];
            end
            
            
            %%
            avg_cropped = nandens(cropped_spikes,smval,'gauss',Fs,'nanflt');
            avg_cropped = avg_cropped-nanmean(avg_cropped(1:twin1));
            avg_cropped = avg_cropped/max(avg_cropped);
            [PKS,LOCS]= findpeaks(avg_cropped,'MinPeakWidth',smval); %find peaks > 2*std in width
            if ~isempty(LOCS)
                %remove peaks less than 2/3 the max
                LOCS(PKS < 0.66) = [];
            end
            if ~isempty(LOCS)
                all_cropped_fields = [all_cropped_fields; avg_cropped./avg_cropped(LOCS(1))];
                all_cropped_peaks = [all_cropped_peaks LOCS(1)];
            end
            
                        avg_warped = nandens(warped_spikes,smval,'gauss',Fs,'nanflt');
            avg_warped = avg_warped-nanmean(avg_warped(1:twin1));
            avg_warped = avg_warped/max(avg_warped);
            [PKS,LOCS]= findpeaks(avg_warped,'MinPeakWidth',smval); %find peaks > 2*std in width
            if ~isempty(LOCS)
                %remove peaks less than 2/3 the max
                LOCS(PKS < 0.66) = [];
            end
            if ~isempty(LOCS)
                all_warped_fields = [all_warped_fields; avg_warped./avg_warped(LOCS(1))];
                all_warped_peaks = [all_warped_peaks LOCS(1)];
            end
            
            avg_cropped2 = nandens3(cropped_spikes2,smval,Fs);
            avg_cropped2 = avg_cropped2-nanmean(avg_cropped2(1:twin1));
            avg_cropped2 = avg_cropped2/max(avg_cropped2);
            [PKS,LOCS]= findpeaks(avg_cropped2,'MinPeakWidth',smval); %find peaks > 2*std in width
            if ~isempty(LOCS)
                %remove peaks less than 2/3 the max
                LOCS(PKS < 0.66) = [];
            end
            if ~isempty(LOCS)
                all_cropped_fields2 = [all_cropped_fields2; avg_cropped2./avg_cropped2(LOCS(1))];
                all_cropped_peaks2 = [all_cropped_peaks2 LOCS(1)];
            end
            
        end
    end
end
%%
t = -twin1:twin2-1;
figure
hold on
plot(t,nanmean(all_fields))
plot(t,nanmean(all_cropped_fields))
plot(t,nanmean(all_cropped_fields2))
plot(t,nanmean(all_warped_fields))
plot([-200 400],[0 0],'k--')
plot([0 0],[-0.1 0.6],'k--')
plot([-44 -44],[-0.1 0.6],'k--')
hold off
xlabel('Time from Fixation Onset (ms)')
ylabel('Normalized Firing Rate')
legend('Raw','Cropped to median of all','Cropped to individual','Warped')

%%


figure
subplot(2,2,1)
vals = all_fields(:,1:twin1); %"baseline" out of field firing rate
vals(isnan(vals)) = [];
%[mx,mxi] = max(all_fields');
[~,place_order] = sort(all_field_peaks); %sort order by peak firing time
imagesc([-twin1:twin2-1],[1:size(all_fields,1)],all_fields(place_order,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_fields,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Raw')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar

subplot(2,2,2)
vals = all_cropped_fields(:,1:twin1); %"baseline" out of field firing rate
vals(isnan(vals)) = [];
% [mx,mxi] = max(all_cropped_fields');
[~,place_order] = sort(all_cropped_peaks); %sort order by peak firing time
imagesc([-twin1:twin2-1],[1:size(all_cropped_fields,1)],all_cropped_fields(place_order,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_cropped_fields,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Cropped to median of all')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar

subplot(2,2,3)
vals = all_cropped_fields2(:,1:twin1); %"baseline" out of field firing rate
vals(isnan(vals)) = [];
% [mx,mxi] = max(all_cropped_fields2');
[~,place_order] = sort(all_cropped_peaks2); %sort order by peak firing time
imagesc([-twin1:twin2-1],[1:size(all_cropped_fields2,1)],all_cropped_fields2(place_order,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_cropped_fields2,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Cropped to individual')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar

subplot(2,2,4)
vals = all_warped_fields(:,1:twin1); %"baseline" out of field firing rate
vals(isnan(vals)) = [];
% [mx,mxi] = max(all_warped_fields');
[~,place_order] = sort(all_warped_peaks); %sort order by peak firing time
imagesc([-twin1:twin2-1],[1:size(all_warped_fields,1)],all_warped_fields(place_order,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_warped_fields,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Cropped to individual')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar