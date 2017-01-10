function Visual_Response_AnalysisV2(data_dir,figure_dir,session_data)
% written by Seth Konig July 11, 2016 updated to V2 August 31, 2016 to focus on
% visual response and time delay so can look at lag for eye movements and
% events.

figure_dir = [figure_dir 'Visual Response\'];
twin1 = 200;% how much time to take before event cross appears and how much to ignore during fixation
twin2 = 1000;%how much time to look at after stimulus onset
twin3 = 500;%for image off
twin4 = 5000; %for long window on image on
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
numshuffs = 1000; %number of shuffles to do for bootstrapping
smval =60; %temporal 1/2 width of gaussian smoothing filter
smval2 = 300;%for long window around imageonset


% time_locked_firing since ITI period is defined by 2 events (15 & 16)
event_codes = [35 8 23 24 23];
trial_start_code = 15;
task = 'ListSQ';
min_blks = 2;

sliding_window = twin1;
sliding_step = 25;  %step size in ms
p_thresh = 0.01;%significance level

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials','hdr');
%grab unit data
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

if num_units == 0
    return
end
num_trials = length(cfg.trl);

valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end

%get important task specific information
[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,23);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import & reformat data so that spikes are locked to events---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Aligning spike times to trial events')

time_lock_firing = cell(num_units,length(event_codes));
which_images = cell(1,num_units); %image #
trial_number = cell(1,num_units); %# in cfg so can import fixstats at later time
nvr = cell(1,num_units);%novel vs repeat
fix_rt = cell(1,num_units); %how long does it take on average to fixate cross once on
%preallocate space and parallel structure of cfg
for unit = 1:num_units
    time_lock_firing{unit,1} = NaN(192,twin1+twin3);%cross hair on
    time_lock_firing{unit,2} = NaN(192,twin1+twin2);%fix crosshair
    time_lock_firing{unit,3} = NaN(192,twin1+twin2);%imgon
    time_lock_firing{unit,4} = NaN(192,2*twin3);%imgoff.
    time_lock_firing{unit,5} = NaN(192,twin1+twin4);%img on long
    fix_rt{unit}=  NaN(1,192);
    nvr{unit} = NaN(1,length(cfg.trl));
    which_images{unit} = NaN(1,length(cfg.trl));
    trial_number{unit} = NaN(1,length(cfg.trl));
end

fix_rt = cell(1,num_units);
for t = 1:num_trials
    if any(cfg.trl(t).allval == event_codes(3)); %in which image was displayed or 1st item in sequence was displayed
        if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(end) %then sequence trial
            continue % go to the next trial
        end
        
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
        trial_end = cfg.trl(t).alltim(end)-trial_start;
        crosson = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(1))-trial_start;
        fix_cross = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(2))-trial_start;
        fix_cross = fix_cross(1);%since fixation on image also counts as event 8
        imgon = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(3))-trial_start;
        
        if imgon-fix_cross > 750%cortex sometimes has minimal lag upto ~20 ms
            imgon2 =  fix_cross+750;
        else
            imgon2 = imgon;
        end
        imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(4))-trial_start;
        for unit = 1:num_units
            if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only put data in for valid trials
                
                %---image info---%
                img_index = find(cfg.trl(t).cnd == img_cnd);
                if any(isnan(which_img(img_index)))
                    continue
                end
                nvr{unit}(img_index) = novel_vs_repeat(img_index);
                which_images{unit}(img_index) = which_img(img_index);
                fix_rt{unit}(img_index) = fix_cross-crosson;
                trial_number{unit}(img_index) = t;
                
                spikes = find(data(unit).values{t});
                
                %---crosson locked firing---%
                event_spikes = spikes(spikes > crosson-twin1 & spikes <= crosson+twin3)-crosson+twin1;
                timevec = [crosson-twin1+1 crosson+twin3]-crosson+twin1;
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                time_lock_firing{unit,1}(img_index,timevec(1):timevec(2)) = tempvec;
                
                %---fix cross locked firing---%
                event_spikes = spikes(spikes > fix_cross-twin1 & spikes <= imgon2)-fix_cross+twin1;
                timevec = [fix_cross-twin1+1 imgon2]-fix_cross+twin1;
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                time_lock_firing{unit,2}(img_index,timevec(1):timevec(2)) = tempvec;
                
                %---image on locked firing---%
                event_spikes = spikes(spikes > imgon-twin1 & spikes <= imgon+twin2)-imgon+twin1;
                timevec = [imgon-twin1+1 imgon+twin2]-imgon+twin1;
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                time_lock_firing{unit,3}(img_index,timevec(1):timevec(2)) = tempvec;
                
                %---image off locked firing---%
                if t ~= length(cfg.trl)
                    %will need to look in next trial to get rest of data
                    %since usually only get ~200 ms of data
                    rest_time = twin3-(trial_end-imgoff);
                    spikes2 = find(data(unit).values{t+1});
                    spikes2(spikes2 > rest_time) = [];
                    spikes2 = spikes2+(twin3-rest_time)+twin3;
                    
                    event_spikes = spikes(spikes > imgoff-twin3 & spikes <= imgoff+twin3)-imgoff+twin3;
                    event_spikes = [event_spikes spikes2];
                    timevec = [imgoff-twin3+1 imgoff+twin3]-imgoff+twin3;
                    tempvec = zeros(1,timevec(2)-timevec(1)+1);
                    tempvec(event_spikes) = 1;
                    time_lock_firing{unit,4}(img_index,timevec(1):timevec(2)) = tempvec;
                end
                
                
                %---long image on locked firing---%
                event_spikes = spikes(spikes > imgon-twin1 & spikes <= imgon+twin4)-imgon+twin1;
                timevec = [imgon-twin1+1 imgon+twin4]-imgon+twin1;
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                time_lock_firing{unit,5}(img_index,timevec(1):timevec(2)) = tempvec;
                
            end
        end
    end
end

%remove excess NaNs associated with error trials
time_lock_firing = laundry(time_lock_firing);
nvr = laundry(nvr);
which_images = laundry(which_images);
fix_rt = cellfun(@nanmean,fix_rt);
trial_number = laundry(trial_number);

%for comparing novel to repeat only want the trials in which both
%images were shown
for unit = 1:num_units
    rmv = [];
    for img = 1:96
        ind = find(which_images{unit} == img);
        if length(ind) == 1 %so either novel or repeat but not both
            rmv = [rmv ind];
        end
    end
    which_images{unit}(rmv) = NaN;
    nvr{unit}(rmv) = NaN;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Perform Temporal analysis and signifance testing---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Determining if spikes are locked to trial events')
Fs = data(1).fsample; %should be 1000
info_type = 'temporal';

% perform mutual information theory analysis to determine
% if neurons encode info about a particular time within an epoch
epoch_data = [];
epoch_data.rate = NaN(num_units,length(event_codes));
epoch_data.temporalstability = NaN(num_units,length(event_codes));
epoch_data.shuffled_rate = cell(num_units,length(event_codes));
epoch_data.shuffled_temporalstability = cell(num_units,length(event_codes));
epoch_data.rate_prctile = NaN(num_units,length(event_codes));
epoch_data.temporalstability_prctile = NaN(num_units,length(event_codes));
%
for unit = 1:num_units
    for event = 1:length(event_codes);
        if isempty(time_lock_firing{unit,event})
            continue
        elseif nansum(nansum(time_lock_firing{unit,event})) == 0; %no spikes
            epoch_data.rate(unit,event) = 0;
            epoch_data.rate_prctile(unit,event) = 0;
            epoch_data.temporalstability_prctile(unit,event) = 0;
        else
            if event == 5 %long image period
                [temporal_info,shuffled_info_rate] = ...
                    estimated_mutual_information(time_lock_firing{unit,event},numshuffs,info_type,smval2,Fs);
            elseif event == 2 %fixation on cross variable duration until image is on so only take first part
                if nansum(nansum(time_lock_firing{unit,event}(:,1:twin3+twin1))) > 0
                    [temporal_info,shuffled_info_rate] = ...
                        estimated_mutual_information(time_lock_firing{unit,event}(:,1:twin3+twin1),numshuffs,info_type,smval,Fs);
                else
                    epoch_data.rate(unit,event) = 0;
                    epoch_data.rate_prctile(unit,event) = 0;
                    epoch_data.temporalstability_prctile(unit,event) = 0;
                end
            else
                [temporal_info,shuffled_info_rate] = ...
                    estimated_mutual_information(time_lock_firing{unit,event},numshuffs,info_type,smval,Fs);
            end
            
            epoch_data.rate(unit,event) = temporal_info.skaggs;
            epoch_data.temporalstability(unit,event) = temporal_info.temporalstability(1);
            epoch_data.shuffled_rate{unit,event} = shuffled_info_rate.skaggs;
            epoch_data.shuffled_temporalstability{unit,event} = shuffled_info_rate.temporalstability(1,:);
            
            epoch_data.rate_prctile(unit,event) = 100*sum(epoch_data.rate(unit,event)>...
                epoch_data.shuffled_rate{unit,event})/numshuffs;
            epoch_data.temporalstability_prctile(unit,event) = 100*sum(epoch_data.temporalstability(unit,event)>...
                epoch_data.shuffled_temporalstability{unit,event})/numshuffs;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Determine if Neuron is Visually Responsive---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sig_nov_rep_times =  zeros(num_units,size(time_lock_firing{unit,3},2)); %for short window
% sig_nov_rep_times2 =  zeros(num_units,size(time_lock_firing{unit,3},2)); %for long window
% p_nov_rep = NaN(num_units,50); %for short window
% p_nov_rep2 = NaN(num_units,50); %for long window
sig_visual_response = zeros(num_units,size(time_lock_firing{unit,3},2)); %for short window
sig_visual_response2 = zeros(num_units,size(time_lock_firing{unit,3},2)); %for long window
baseline_firing_rate = NaN(2,num_units); %mean and std
p_visual_response = NaN(num_units,50);
for unit = 1:num_units
    if ~isempty(time_lock_firing{unit,3})
        
        %define baseline firing rate of neuron during central fixation ignoring
        %first twin1 ~200 ms in case there a fixation response
        
        %twin1 get to event start then want to go another twin1 to ignore 200 ms
        avg_firing_rate = 1000*nansum(time_lock_firing{unit,3}(:,1:twin1),2)./...
            sum(~isnan(time_lock_firing{unit,3}(:,1:twin1)),2);
        baseline_firing_rate(1,unit) = mean(avg_firing_rate);
        baseline_firing_rate(2,unit) = std(avg_firing_rate);
        fr_distribution = sum(time_lock_firing{unit,3}(:,1:twin1),2);
        
        
        %---Analysis of Visual Response for Short Window---%
        total_time = size(time_lock_firing{unit,3},2);
        down = 1000/sliding_window;%scaling of average firing rate to amount of time in sliding window
        for step = 1:(total_time/sliding_step)-sliding_window/sliding_step+1
            ind = sliding_step*(step-1)+1:sliding_window+sliding_step*(step-1);
            [~,p_visual_response(unit,step)] = kstest2(sum(time_lock_firing{unit,3}(:,ind)'),...
                fr_distribution);
            if p_visual_response(unit,step) < p_thresh
                sig_visual_response(unit,ind) = 2;
            end
        end
        
        %---Analysis for Long Window---%
        total_time = size(time_lock_firing{unit,5},2);
        for step = 1:(total_time/sliding_step)-sliding_window/sliding_step+1
            ind = sliding_step*(step-1)+1:sliding_window+sliding_step*(step-1);
            [~,p_visual_response2(unit,step)] = kstest2(sum(time_lock_firing{unit,5}(:,ind)'),...
                fr_distribution);
            if p_visual_response2(unit,step) < p_thresh
                sig_visual_response2(unit,ind) = 2;
            end
        end
    end
end

% p_nov_rep = cell(2,num_units);% for short window
% p_nov_rep2 = cell(2,num_units); %for long window
% for unit = 1:num_units
%     if ~isempty(time_lock_firing{unit,3})
%
%         %for short window
%         all_curves = NaN(10*numshuffs,twin1+twin2);
%         nov_rep_ind = find(~isnan(nvr{unit})); %can have NaNs which are meaningful here since unpaired novel/repeat images
%         parfor shuff = 1:10*numshuffs
%             ind = randperm(length(nov_rep_ind ));
%             shuffled_ind = nov_rep_ind(ind);
%             nov_curve = nandens(time_lock_firing{unit,3}(nvr{unit}(shuffled_ind) == 1,:),smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
%             rep_curve = nandens(time_lock_firing{unit,3}(nvr{unit}(shuffled_ind) == 2,:),smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
%             all_curves(shuff,:) = nov_curve-rep_curve;
%         end
%         p_nov_rep{1,unit} = prctile(all_curves,97.5,1);
%         p_nov_rep{2,unit} = prctile(all_curves,2.5,1);
%
%         %for short window
%         all_curves = NaN(10*numshuffs,twin1+twin4);
%         nov_rep_ind = find(~isnan(nvr{unit})); %can have NaNs which are meaningful here since unpaired novel/repeat images
%         parfor shuff = 1:10*numshuffs
%             ind = randperm(length(nov_rep_ind ));
%             shuffled_ind = nov_rep_ind(ind);
%             nov_curve = nandens(time_lock_firing{unit,5}(nvr{unit}(shuffled_ind) == 1,:),smval2,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
%             rep_curve = nandens(time_lock_firing{unit,5}(nvr{unit}(shuffled_ind) == 2,:),smval2,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
%             all_curves(shuff,:) = nov_curve-rep_curve;
%         end
%         p_nov_rep2{1,unit} = prctile(all_curves,99,1);
%         p_nov_rep2{2,unit} = prctile(all_curves,1,1);
%
%     end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot and Save Figures of Results---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t12 = -twin1:twin2-1;
t13 = -twin1:twin3-1;
t133 = -twin1:twin3-1+250;
t3 = -twin3:twin3-1;
t14 = -twin1:twin4-1;
unit_names = unit_stats(1,:);
for unit = 1:num_units
    if ~isempty(time_lock_firing{unit,1})
        figure
        
        %---Firing rate locked to Cross and Fixation on Cross---%
        subplot(2,2,1)
        hold on
        [~,~,~,y] =dofill(t13-round(fix_rt(unit)),time_lock_firing{unit,1},'black',1,smval);
        [~,~,~,y1] =dofill2(t133,time_lock_firing{unit,2},'green',1,smval);
        plot([-twin3 twin3+250],[baseline_firing_rate(1,unit) baseline_firing_rate(1,unit)],'k--')
        hold off
        xl = xlim;
        xlim([xl(1) twin3+250])
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Fixation on Cross (ms)')
        legend('Cross On','Fix on Cross','"Baseline"')
        if epoch_data.rate_prctile(unit,2) > 90 || epoch_data.temporalstability_prctile(unit,2) > 90
            title(['bit ' num2str(epoch_data.rate_prctile(unit,2),3) '% ' ...
                '\rho = ' num2str(epoch_data.temporalstability(unit,2),2) ' ' ...
                num2str(epoch_data.temporalstability_prctile(unit,2),3)  '%' ])
        end
        
        %---Raster aligned to fixation on cross---%
        subplot(2,2,3)
        [trial,time] = find(time_lock_firing{unit,2}(:,1:twin3+twin1) == 1);
        if ~isempty(trial)
            plot(time-twin1,trial,'.k')
            xlim([xl(1) twin3])
            ylim([0 max(trial)+1]);
        end
        ylabel('Trial #')
        xlabel('Time from Fixation on Cross (ms)')
        
        
        %---Firing rate locked to Image Off---%
        subplot(2,2,2)
        [~,~,~,y2] =dofill(t3,time_lock_firing{unit,4},'black',1,smval);
        xlim([-twin3 twin3])
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Image Off (ms)')
        if epoch_data.rate_prctile(unit,4) > 90 || epoch_data.temporalstability_prctile(unit,4) > 90
            title(['bit ' num2str(epoch_data.rate_prctile(unit,4),3) '% ' ...
                '\rho = ' num2str(epoch_data.temporalstability(unit,4),2) ' ' ...
                num2str(epoch_data.temporalstability_prctile(unit,4),3)  '%' ])
        end
        
        %---Raster aligned to Image off---%
        subplot(2,2,4)
        [trial,time] = find(time_lock_firing{unit,4} == 1);
        if ~isempty(trial)
            plot(time-twin3,trial,'.k')
            xlim([-twin3 twin3])
            ylim([0 max(trial)+1]);
        end
        ylabel('Trial #')
        xlabel('Time from Image Off (ms)')
        
        
        %---Scale plots to be the same---%
        ymin = 0.85*min([y y1 y2]);
        if ymin < 0.1
            ymin = 0;
        end
        ymax = 1.2*max([y y1 y2]);
        if ymin ~= ymax
            subplot(2,2,1)
            ylim([ymin ymax])
            subplot(2,2,2)
            ylim([ymin ymax])
        end
        
        n_str = ['Cross and Image Off, n_ =' num2str(size(time_lock_firing{unit,1},1))];
        if multiunit(unit)
            subtitle(['Multiunit ' unit_names{unit} n_str]);
        else
            subtitle(['' unit_names{unit} n_str]);
        end
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_Visual Response_CrossHair']);
        
        
        figure
        
        %---Firing rate locked to Image On---%
        subplot(2,2,1)
        hold on
        [~,~,~,y3] =dofill(t12,time_lock_firing{unit,3},'black',1,smval);
        plot([-twin1 twin2],[baseline_firing_rate(1,unit) baseline_firing_rate(1,unit)],'k--')
        plot([-twin1 twin2],[baseline_firing_rate(1,unit) baseline_firing_rate(1,unit)],'k--')
        yl = ylim;
        gaps = findgaps(find(( sig_visual_response(unit,:))));
        if ~isempty(gaps)
            for g = 1:size(gaps,1)
                gp = gaps(g,:);
                gp(gp == 0) = [];
                h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
                    [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
                uistack(h,'down')
                set(h,'facealpha',.25,'EdgeColor','None')
            end
        end
        hold off
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Image On (ms)')
        xlim([-twin1 twin2])
        if epoch_data.rate_prctile(unit,3) > 90 || epoch_data.temporalstability_prctile(unit,3) > 90
            title(['bit ' num2str(epoch_data.rate_prctile(unit,3),3) '% ' ...
                '\rho = ' num2str(epoch_data.temporalstability(unit,3),2) ' ' ...
                num2str(epoch_data.temporalstability_prctile(unit,3),3)  '%' ])
        end
        
        %---raster for image on---%
        %in case we want to look at latency or something
        subplot(2,2,3)
        [trial,time] = find(time_lock_firing{unit,3} == 1); %novel
        if ~isempty(trial)
            plot(time-twin1,trial,'.k')
            xlim([-twin1 twin2])
            ylim([0 max(trial)+1]);
        end
        ylabel('Trial #')
        xlabel('Time from Image On (ms)')
        
     
        %---Firing rate locked to Long Image On---%
        subplot(2,2,2)
        [~,~,~,y4] =dofill(t14,time_lock_firing{unit,5},'black',1,smval2);
        hold on
        plot([-twin1 twin4],[baseline_firing_rate(1,unit) baseline_firing_rate(1,unit)],'k--')
          yl = ylim;
        gaps = findgaps(find((sig_visual_response2(unit,:))));
        if ~isempty(gaps)
            for g = 1:size(gaps,1)
                gp = gaps(g,:);
                gp(gp == 0) = [];
                h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
                    [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
                uistack(h,'down')
                set(h,'facealpha',.25,'EdgeColor','None')
            end
        end
        hold off
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Image On (ms)')
        xlim([-twin1 twin4])
        if epoch_data.rate_prctile(unit,5) > 90 || epoch_data.temporalstability_prctile(unit,5) > 90
            title(['bit ' num2str(epoch_data.rate_prctile(unit,5),3) '% ' ...
                '\rho = ' num2str(epoch_data.temporalstability(unit,5),2) ' ' ...
                num2str(epoch_data.temporalstability_prctile(unit,5),3)  '%' ])
        end
        
        
     
        %---Raster for Long Image On---%
        subplot(2,2,4)
        [trial,time] = find(time_lock_firing{unit,5} == 1);
        if ~isempty(trial)
            plot(time-twin1,trial,'.k')
            xlim([-twin1 twin4])
            ylim([0 max(trial)+1]);
        end
        ylabel('Trial #')
        xlabel('Time from Image On (ms)')

        %---Scale plots to be the same---%
        ymin = 0.85*min([y3 y4]);
        if ymin < 0.1
            ymin = 0;
        end
        ymax = 1.2*max([y3 y4]);
        if ymin ~= ymax
            subplot(2,2,1)
            ylim([ymin ymax])
            subplot(2,2,2)
            ylim([ymin ymax])
        end
        
        
        n_str = ['Image Onset, n_ =' num2str(size(time_lock_firing{unit,1},1))];
        if multiunit(unit)
            subtitle(['Multiunit ' unit_names{unit} n_str]);
        else
            subtitle(['' unit_names{unit} n_str]);
        end
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_Image_Visual Response']);
        
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Finally save all the data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat'],'time_lock_firing',...
    'smval','smval2','epoch_data','unit_stats','twin1','twin2','twin3','twin4',...
    'p_thresh','sliding_window','sliding_step','baseline_firing_rate',...
    'sig_visual_response','p_visual_response','sig_visual_response2','p_visual_response2',...
    'unit_names','nvr','which_images','fix_rt','trial_number')
%   'p_nov_rep2','sig_nov_rep_times2','p_nov_rep','sig_nov_rep_times'
disp(['Time Locked Data Analyis for ' task_file(1:8) ' saved']);
end