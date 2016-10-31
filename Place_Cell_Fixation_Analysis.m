function  Place_Cell_Fixation_Analysis(data_dir,figure_dir,session_data)
%written by Seth Konig 4/25/16 updated August 30, 2016
%code only analyzes units with significant (95+ percentile) skaggs information
%score and Miriam's spatial stability score. Code determines/verifies
%that firing rate within the place field is great than firing rate outside
%of the place field. A common theme off
%view cells (based on an initial pass and a logical assumption one would make)
%is that they show strong perisaccadic modulation usually firing more after
%a few fixation has started (hence based on viewing location). This code also
%analyzes firing rates during sequence trials for these neurons to
%determine if there is a correlation in firing rate of the view cell during
% images and during sequence trials.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set important task parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
place_cell_dir = [figure_dir 'Place Cells Fixation Analysis\']; %mother dir
best_place_cell_dir = [place_cell_dir 'Best Place Cells Fixation Analysis\']; %super conservative: skagg+corr 1/2+corr e/o
skagg_cell_dir99 = [place_cell_dir 'Skagg 99\'];%want to see what significant skagg scores > 99%-tile
skagg_cell_dir = [place_cell_dir 'Skagg only\'];%want to see what significant skagg scores shows
spatial_corr_dir = [place_cell_dir 'Spatial Correlation Only\']; %want to see what spatially consistent score shows
task = 'ListSQ';

colors =['rgkm'];
shapes = ['xo'];

imageX = 800;
imageY = 600;
img_on_code = 23;
img_off_code = 24;
ITIstart_code = 15;

min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect

Fs = 1000;
sliding_window = 200; %width in ms
sliding_step = 50;  %step size in ms
p_thresh = 0.05;%significance level
fixwin = 5;%size of fixation window on each crosshair
event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
smval = 60;%0.5*smval is gaussian std for smoothing
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
twin1 = 200;%200 ms before fixation
twin2 = 400;%400 ms after start of fixation
image_on_twin = 500;%how much time to ignore eye movements for to remove strong visual response though some may last longer

statistical_window = twin1+1:twin1+twin2;

task_file = get_task_data(session_data,task);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
try
    load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat']);
catch
    disp('Spatial Analysis File not found continueing')
    return
end
load([data_dir task_file(1:end-11) '-preprocessed.mat']);
absolute_fixationstats = fixationstats;
absolute_cfg = cfg;

[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

if num_units == 0; %if no units exit function
    disp([task_file(1:8) ': no units could be found. Exiting function...'])
    return;
end

valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
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
%get important task specific information
[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[which_img,~,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);


%---determine number of blocks unit was stable for
%remove units with too few trials
%these are the absolute minimum data required to do data analysis may want
%to be more strigent later but not worth doing analysis (especially
%shuffling) on such few trials for these neurons


%filter parameters
H = define_spatial_filter(filter_width);

%---Pre-allocate space for Fixations Inside vs Outside Field Analysis---%
list_fixation_locked_firing = cell(1,num_units); %firing rate locked to fixations
in_out = cell(1,num_units); %firing rate for fixations inside (1) or outside (0) of place field
area = NaN(1,num_units);
% p_in_out = NaN(num_units,(twin1+twin2)/sliding_step); %is firing rate higher inside or outside the firing field
p_in_out = NaN(num_units,2); %is firing rate higher inside or outside the firing field


%---Pre-allocate space for Fixations during List vs Sequence Trials Analysis---%
sequence_fixation_locked_firing = cell(4,num_units);
trial_nums = cell(1,num_units);
corr_coeff_seqeuence_list = NaN(1,num_units);
p_corr_seqeuence_list = NaN(1,num_units);
slope_seqeuence_list = NaN(1,num_units);
sequences_inside = NaN(2,4);
% seq_p_in_out = NaN(num_units,(twin1+twin2)/sliding_step); %is firing rate higher inside or outside the firing field
seq_p_in_out = NaN(1,num_units); %is firing rate higher inside or outside the firing field
all_place_field_matrix = cell(1,num_units);

contextual_gain = NaN(5,num_units);
%1) list peak location
%2) list peak
%3) sequence peak location
%4) sequence peak
%5) gain: peak list/peak sequence

for unit = 1:num_units
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Unit Passes 95% for both Skaggs and Stability--%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isnan(spatial_info.shuffled_rate_prctile(unit))
        continue %analysis wasn't run on this unit
    elseif (spatial_info.shuffled_rate_prctile(unit) < 95) && (spatial_info.spatialstability_halves_prctile(unit) < 95) ...
            && (spatial_info.spatialstability_even_odd_prctile(unit) < 95);
        continue %unit not spatial or correlated in any sense
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine Firing rate Locked to Fixations Inside vs Outside Field---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate Firing Rate Map---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    firing_rate_map = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all');
    place_field_matrix = define_place_field(firing_rate_map,imageX,imageY); %place field location
    all_place_field_matrix{unit} = place_field_matrix;
    area(unit) = 100*nansum(nansum(place_field_matrix))/(imageX*imageY);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate Firing Rate Locked to Fixations---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fix_in_out = NaN(1,3000); %firing rate for fixations inside or outside  of place field
    %1) first fixation in: out-> in
    %2) fixation in field but not first: in -> in
    %3) first fixation out of field: in -> out
    %4) fixation out of field but not first: out-> out
    fix_locked_firing = NaN(3000,(twin1+twin2)); %firing rate locked to fixations
    
    fix_ind = 1; %fixation # so can track in variables above
    fixationstats = absolute_fixationstats; %reload because written over below
    cfg = absolute_cfg;
    num_trials = length(cfg.trl);
    %easier to start over rather than figure out format from imported data
    for t = 1:num_trials
        if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
            if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %want to avoid visual response
                imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start;
                
                % if monkey isn't paying attention and looked away image presentation
                % is now longer than imgdur (because of cumulative looking time)
                % so data isn't probably worth much plus have to cut off somewhere
                if imgoff-imgon > 1.5*imgdur-1 %cut off trial at 1.5x length of desired looking time
                    imgoff = imgon+1.5*imgdur-1;
                end
                imgon = imgon+image_on_twin;
                
                fixationtimes = fixationstats{t}.fixationtimes;
                saccadetimes = fixationstats{t}.saccadetimes;
                fixations = fixationstats{t}.fixations;
                xy = fixationstats{t}.XY;
                
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
                if isempty(img_index) || any(isnan(which_img(img_index)))
                    continue
                end
                
                spikes = find(data(unit).values{t});
                
                for f = 2:size(fixations,2)%ignore first fixation not sure where it was/possibly contaminated anyway
                    prior_sac = find(saccadetimes(2,:) == fixationtimes(1,f)-1);%next fixation should start immediately after
                    if isempty(prior_sac) %no prior saccade so was proabbly looking off screen
                        continue;
                    end
                    sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %saccade amplitude
                    if sacamp < min_sac_amp %saccade is too small so ignore
                        continue
                    end
                    
                    
                    prior_fix_in_out = [];
                    %determine if prior fixation was in or out of place field
                    fixx = round(fixations(1,f-1));
                    fixx(fixx < 1) = 1;
                    fixx(fixx > imageX) = imageX;
                    fixy = imageY-round(fixations(2,f-1));
                    fixy(fixy < 1) = 1;
                    fixy(fixy > imageY) = imageY;
                    
                    if place_field_matrix(fixy,fixx) == 1 %then prior fixation was already inside
                        prior_fix_in_out = 1;
                    else
                        prior_fix_in_out = 0;
                    end
                    
                    
                    %determine if fixation was in our out of place field
                    fixx = round(fixations(1,f));
                    fixx(fixx < 1) = 1;
                    fixx(fixx > imageX) = imageX;
                    fixy = imageY-round(fixations(2,f));
                    fixy(fixy < 1) = 1;
                    fixy(fixy > imageY) = imageY;
                    if place_field_matrix(fixy,fixx) == 1 %then inside
                        if prior_fix_in_out == 1%so prior fixation was inside so in->in
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
                    
                    fix_ind = fix_ind+1;
                end
            end
        end
    end
    
    %remove excess NaNs associated with error trials
    fix_in_out = laundry(fix_in_out);
    fix_locked_firing = laundry(fix_locked_firing);
    
    
    %store variables across units for later access
    list_fixation_locked_firing{unit} = fix_locked_firing;
    in_out{unit} = fix_in_out;
    
    %     for step = 1:((twin1+twin2)/sliding_step)-sliding_window/sliding_step+1
    %         ind = twin1+1:twin1+twin2;
    ind = statistical_window;
    %for out->in vs out->out
    [~,p_in_out(unit,1)] = kstest2(nansum(fix_locked_firing(fix_in_out == 1,ind)'),nansum(fix_locked_firing(fix_in_out == 4,ind)'),'tail','smaller');
    if sum(fix_in_out  == 1) > 0  || sum(fix_in_out == 4) > 0
        %for ?->in vs ?->out
        [~,p_in_out(unit,2)] = kstest2(nansum(fix_locked_firing(fix_in_out == 1 | fix_in_out == 2,ind)'),...
            nansum(fix_locked_firing(fix_in_out == 4 | fix_in_out == 3,ind)'),'tail','smaller');
    end
    
    %     sig_time_in_out = zeros(1,twin1+twin2);
    %     sig_time_nov_rep = zeros(1,twin1+twin2);
    %     for step = 1:((twin1+twin2)/sliding_step)-sliding_window/sliding_step+1
    %         ind = sliding_step*(step-1)+1:sliding_window+sliding_step*(step-1);
    %         if p_in_out(unit,step) < p_thresh
    %             sig_time_in_out(ind) = 2;
    %         end
    %         if p_nov_rep(unit,step) < p_thresh
    %             sig_time_nov_rep(ind) = 2;
    %         end
    %     end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%
    %%%---Make Plot---%%%
    %%%%%%%%%%%%%%%%%%%%%
    t = -twin1:twin2-1;
    maxfr = prctile(firing_rate_map(:),95);
    
    figure
    
    %---plot firing rate map for all images---%
    subplot(2,3,1)
    h = imagesc(firing_rate_map);
    set(h,'alphadata',~isnan(firing_rate_map));
    title('All images')
    axis off
    axis equal
    colorbar
    colormap('jet')
    clim = caxis;
    caxis([clim(1) maxfr])
    
    title_str = ['All images, peak rate = ' num2str(max(firing_rate_map(:)),3) ' Hz\n'];
    if spatial_info.shuffled_rate_prctile(unit) >= 95;
        title_str = [title_str 'Bits = ' num2str(spatial_info.rate(unit),3) ...
            '(' num2str(spatial_info.shuffled_rate_prctile(unit),3) '%%) '];
    end
    if (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
        title_str = [title_str '\n \\rho_{1/2} = ' num2str(spatial_info.spatialstability_halves(1,unit),2) ...
            '(' num2str(spatial_info.spatialstability_halves_prctile(1,unit),3) '%%)'];
    end
    if (spatial_info.spatialstability_even_odd_prctile(1,unit) > 95)
        title_str = [title_str ' \\rho_{e/o} = ' num2str(spatial_info.spatialstability_even_odd(1,unit),2) ...
            '(' num2str(spatial_info.spatialstability_even_odd_prctile(1,unit),3) '%%)'];
        
    end
    title(sprintf(title_str));
    
    
    %---draw place field---%
    subplot(2,3,4)
    imagesc(place_field_matrix)
    axis off
    axis equal
    title(sprintf(['place field Location \n Area = ' num2str(area(unit),2)]));
    
    %---Fixations in->out vs out->out---%
    subplot(2,3,2)
    [trial,time] = find(fix_locked_firing(fix_in_out == 4,:) == 1);
    plot(time-twin1,(trial),'.r')
    hold on
    if ~isempty(trial)
        b4 = max(trial);
    else
        b4 = 0;
    end
    [trial,time] = find(fix_locked_firing(fix_in_out == 1,:) == 1);
    trial = trial+b4;
    plot(time-twin1,(trial),'.b')
    if ~isempty(trial)
        ylim([0 max(trial)])
    else
        ylim([0 b4])
    end
    ylabel('Occurence #')
    set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
    xlabel('Time from Fixatoin Start')
    title('Fixation Aligned Rasters')
    
    subplot(2,3,5)
    hold on
    [~,~,~,y_list,~] = dofill(t,fix_locked_firing(fix_in_out == 1,:),'blue',1,smval);%out-> in
    dofill(t,fix_locked_firing(fix_in_out == 4,:),'red',1,smval);%out->out
    yl = ylim;
    if yl(1) < 0
        yl(1) = 0;
        ylim(yl);
    end
    plot([0 0],[yl(1) yl(2)],'k--')
    %     gaps = findgaps(find(sig_time_in_out));
    %     if ~isempty(gaps)
    %         for g = 1:size(gaps,1)
    %             gp = gaps(g,:);
    %             gp(gp == 0) = [];
    %             h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
    %                 [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
    %             uistack(h,'down')
    %             set(h,'facealpha',.25,'EdgeColor','None')
    %         end
    %     end
    if p_in_out(unit,1) < p_thresh
        h = fill([twin1 twin1+twin2  twin1+twin2 twin1 twin1]-twin1,...
            [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
        set(h,'facealpha',.25,'EdgeColor','None')
    end
    xlim([-twin1 twin2]);
    hold off
    legend('out->in','out->out','Location','NorthWest')
    set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
    title(sprintf(['n_{out->in} = ' num2str(sum(fix_in_out == 1))...
        ', n_{out->out} = ' num2str(sum(fix_in_out == 4))]));
    
    %find peaks for list task
    [pks,locs] = findpeaks(y_list,'MinPeakWidth',50);
    pks(locs < twin1) = [];
    locs(locs < twin1) = [];
    list_pk = locs(pks == max(pks));
    list_max = max(pks);
    
    %---Fixations ?->in vs ?-> out---%
    subplot(2,3,6)
    hold on
    dofill(t,fix_locked_firing(fix_in_out == 1 | fix_in_out == 2,:),'black',1,smval);%?->in
    dofill(t,fix_locked_firing(fix_in_out == 3 | fix_in_out == 4,:),'green',1,smval);%?->out
    yl = ylim;
    if yl(1) < 0
        yl(1) = 0;
        ylim(yl);
    end
    plot([0 0],[yl(1) yl(2)],'k--')
    if p_in_out(unit,2) < p_thresh
        h = fill([twin1 twin1+twin2  twin1+twin2 twin1 twin1]-twin1,...
            [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
        set(h,'facealpha',.25,'EdgeColor','None')
    end
    hold off
    xlabel('Time from Fixation Onset (ms)')
    ylabel('Firing Rate (Hz)')
    xlim([-twin1 twin2])
    set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
    legend('?->in','?->out','Location','NorthWest')
    title(sprintf(['n_{?->in} = ' num2str(sum(fix_in_out == 1 | fix_in_out == 2))...
        ', n_{?->out} = ' num2str(sum(fix_in_out == 3 | fix_in_out == 4))]));
    
    subplot(2,3,3)
    [trial,time] = find(fix_locked_firing(fix_in_out == 3 | fix_in_out == 4,:) == 1);
    plot(time-twin1,(trial),'.g')
    hold on
    if ~isempty(trial)
        b4 = max(trial);
    else
        b4 = 0;
    end
    [trial,time] = find(fix_locked_firing(fix_in_out == 1 | fix_in_out == 2,:) == 1);
    trial = trial+b4;
    plot(time-twin1,(trial),'.k')
    if ~isempty(trial)
        ylim([0 max(trial)])
    else
        ylim([0 b4])
    end
    ylabel('Occurence #')
    set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
    xlabel('Time from Fixatoin Start')
    title('Fixation Aligned Rasters')
    xlim([-twin1 twin2])
    hold off
    
    num_trials_str = [' n_{img} = ' num2str(size(spike_times{unit},1))];
    if multiunit(unit)
        multi_str = 'Multiunit ';
    else
        multi_str = ' ';
    end
    subtitle(['Spatial Plots' num_trials_str multi_str ' ' task_file(1:8) ' ' unit_stats{1,unit}]);
    
    
    if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
            && (spatial_info.spatialstability_even_odd_prctile(unit) > 95) ... %spatial consistency
            && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
        if p_in_out(unit,2) > 0.05 && p_in_out(unit,1) > 0.05
            disp('why')
        end
    end
    
    if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95) ...
            && (spatial_info.spatialstability_even_odd_prctile(1,unit) > 95);
        save_and_close_fig(best_place_cell_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_fixation_analysis']);
    elseif (spatial_info.shuffled_rate_prctile(unit) > 95) && ((spatial_info.spatialstability_halves_prctile(1,unit) > 95) ...
            || (spatial_info.spatialstability_even_odd_prctile(1,unit) > 95));
        save_and_close_fig(place_cell_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_fixation_analysis']);
    elseif (spatial_info.shuffled_rate_prctile(unit) > 99)
        save_and_close_fig( skagg_cell_dir99,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_fixation_analysis']);
    elseif (spatial_info.shuffled_rate_prctile(unit) > 95)
        save_and_close_fig(skagg_cell_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_fixation_analysis']);
    else
        save_and_close_fig(spatial_corr_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_fixation_analysis']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Correlation during List and Sequence trials---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %preallocate space and parallel structure of cfg
    successful_sequence_trials = NaN(1,length(cfg.trl));
    which_sequence = NaN(1,length(cfg.trl));
    for t = 1:length(cfg.trl);
        if sum(cfg.trl(t).allval == 3) == 6; %in which sequence trials were rewarded
            which_sequence(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
            successful_sequence_trials(t) = t;
        end
    end
    successful_sequence_trials = laundry(successful_sequence_trials);
    which_sequence = laundry(which_sequence);
    
    num_trials = length(successful_sequence_trials);
    event_times = NaN(length(successful_sequence_trials),8);
    event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
    for t = 1:num_trials
        trial_start = cfg.trl(successful_sequence_trials(t)).alltim(cfg.trl(successful_sequence_trials(t)).allval == ITIstart_code);
        for event = 1:length(event_codes);
            event_times(t,event) = cfg.trl(successful_sequence_trials(t)).alltim(cfg.trl(successful_sequence_trials(t)).allval == event_codes(event))-trial_start;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---process eye data locked to trial events---%%%
    fixationstats = absolute_fixationstats; %reload because will be written over below
    cfg = absolute_cfg; %reload because will be written over below
    fixationstats = fixationstats(successful_sequence_trials);
    cfg.trl = cfg.trl(successful_sequence_trials);
    
    saccade_start_time = NaN(length(fixationstats),4);%when did saccade to item start
    fixation_start_time = NaN(length(fixationstats),4);%when did fixation on item start
    reaction_time = NaN(length(fixationstats),4); %how long after a stimulus disappears do they saccade
    time_to_fixation = NaN(length(fixationstats),4); %how fast did they get to it the first time around
    fixation_accuracy = NaN(length(fixationstats),4); %how far off
    fixation_duration = NaN(length(fixationstats),4); %fixation duration
    extrafixations = NaN(length(fixationstats),4); %how many noncross hair fixations do they mak
    for trial = 1:num_trials
        locs = sequence_locations{which_sequence(trial)};
        
        %convert to DVA for this analysis
        locs(1,:) = (locs(1,:)-400)/24;
        locs(2,:) = (locs(2,:)-300)/24;
        fixationstats{trial}.fixations(1,:) = (fixationstats{trial}.fixations(1,:)-400)/24;
        fixationstats{trial}.fixations(2,:) = (fixationstats{trial}.fixations(2,:)-300)/24;
        
        event_codes = cfg.trl(trial).allval;
        event_codes(event_codes == 100)= 0;
        event_codes(1) = 100;%eye data starts for recording right away
        event_times = cfg.trl(trial).alltim;
        
        trialdata = analyze_sequence_trial(fixationstats{trial},locs,fixwin,...
            event_codes,event_times);
        
        time_to_fixation(trial,:) = trialdata.t2f;
        fixation_duration(trial,:) = trialdata.fixation_duration;
        reaction_time(trial,:)  = trialdata.time_to_leave;
        extrafixations(trial,:) = trialdata.extrafixations;
        fixation_accuracy(trial,:) =  trialdata.accuracy;
        
        fixation_numbers = trialdata.fixationnums; %fixation number for each item
        fixationtimes = fixationstats{trial}.fixationtimes;
        saccadetimes = fixationstats{trial}.saccadetimes;
        for item = 1:4
            if ~isnan(fixation_numbers(item))
                fixation_start_time(trial,item) = fixationtimes(1,fixation_numbers(item));
                saccadeind = find(saccadetimes(2,:)+1 ==  fixation_start_time(trial,item));
                if ~isempty(saccadeind)
                    saccade_start_time(trial,item) = saccadetimes(1,saccadeind);
                end
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate Firing Rate Locked to Eye Data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for c = 1:4
        sequence_fixation_locked_firing{c,unit} = NaN(num_trials,twin1+twin2);
        trial_nums{c,unit} = NaN(1,num_trials);
    end
    
    for trial = 1:num_trials
        if successful_sequence_trials(trial) >= valid_trials(1,unit) && successful_sequence_trials(trial) <= valid_trials(2,unit)
            spikes = find(data(unit).values{successful_sequence_trials(trial)});
            for c = 1:4;
                fixt = fixation_start_time(trial,c);
                sact = saccade_start_time(trial,c);
                
                %if no saccade detected on item want to check if previous
                %trial has saccade/begining of fixation
                if isnan(sact) && isnan(fixt)
                    continue
                elseif ~isnan(sact) && isnan(fixt)
                    error(['Could not identify fixation but could locate saccade! Trial#:' num2str(trial)])
                elseif ~isnan(fixt)
                    fix_spikes = spikes(spikes > fixt-twin1 & spikes <= fixt+twin2)-fixt+twin1;
                    temp = zeros(1,twin1+twin2);
                    temp(fix_spikes) = 1;
                    sequence_fixation_locked_firing{c,unit}(trial,:) = temp;
                    trial_nums{c,unit}(trial) = trial;
                end
            end
        end
    end
    sequence_fixation_locked_firing(:,unit) = laundry(sequence_fixation_locked_firing(:,unit));
    trial_nums(:,unit) = laundry(trial_nums(:,unit));
    
    %---Calculate Observed correlation/slope and shuffled corrleation/slope---%
    %using trial by trial permutation
    temp_list = []; %need to define for parfor even if don't use,
    %for sequence locations in list pixels wihtout enough coverage
    
    %observed data
    imgy = size(firing_rate_map,1); %vertical size of firing rate matrix
    fr_seq = NaN(2,4);
    fr_list = NaN(2,4);
    all_sequence_fixation_locked_firing = [];
    all_sequences = [];
    for c = 1:4
        for seq = 1:2
            [y,~] = nandens(sequence_fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == seq ,twin1:end),smval,'gauss',Fs,'nanflt');%nandens made by nathan killian
            binned_location = round([sequence_locations{seq}(1,c)/binsize,imgy-sequence_locations{seq}(2,c)/binsize]);
            fr_seq(seq,c) =  max(y);
            fr_list(seq,c) = firing_rate_map(binned_location(2),binned_location(1));
            if isnan(fr_list(seq,c)) %check surrounding bins since binning and rounding
                for i = -1:1
                    for ii = -1:1
                        if ~isnan(firing_rate_map(binned_location(2)+i,binned_location(1)+ii))
                            fr_list(seq,c) = firing_rate_map(binned_location(2)+i,binned_location(1)+ii);
                            break
                        end
                    end
                end
            end
            
            all_sequence_fixation_locked_firing = [all_sequence_fixation_locked_firing; sequence_fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}),twin1:end)];
            all_sequences = [all_sequences ((seq-1)*4+c)*ones(1,length((which_sequence(trial_nums{c,unit}))))];
        end
    end
    if any(isnan(fr_list(:)))
        temp_list = fr_list(:);
        temp_list(isnan(temp_list)) = [];
        temp_seq = fr_seq(:);
        temp_seq(isnan(fr_list(:))) = [];
        
        [r,rp] =corrcoef(temp_seq(:),temp_list(:));
        corr_coeff_seqeuence_list(unit) =r(2);
        p_corr_seqeuence_list(unit)=rp(2);
        P_fit = polyfit(temp_seq(:),temp_list(:),1);
        slope_seqeuence_list(unit) = P_fit(1);
        
    else
        [r,rp] =corrcoef(fr_seq(:),fr_list(:));
        corr_coeff_seqeuence_list(unit) =r(2);
        p_corr_seqeuence_list(unit)=rp(2);
        P_fit = polyfit(fr_seq(:),fr_list(:),1);
        slope_seqeuence_list(unit) = P_fit(1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Plot Corrrelation Between List and Sequence Firing Rates---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    imgy = size(place_field_matrix,1);
    %Determine if any of the items are inside the matrix or not?
    sequence_inside = NaN(2,4);
    for c = 1:4
        for seq = 1:2
            if place_field_matrix(imgy-sequence_locations{seq}(2,c),sequence_locations{seq}(1,c)) == 1
                sequence_inside(seq,c) =1;
            end
        end
    end
    
    if any(sequence_inside(:) == 1) && ~all(sequence_inside(:) == 1)
        sequence_insides = sequence_inside;
        
        fixation_locked_inside = [];
        fixation_locked_outside = [];
        for c = 1:4
            for seq = 1:2
                if  sequence_inside(seq,c) == 1;
                    fixation_locked_inside =[fixation_locked_inside; sequence_fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == seq ,:)];
                else
                    fixation_locked_outside =[fixation_locked_outside; sequence_fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == seq ,:)];
                end
            end
        end
        
        %         for step = 1:((twin1+twin2)/sliding_step)-sliding_window/sliding_step+1
        %             ind = sliding_step*(step-1)+1:sliding_window+sliding_step*(step-1);
        %             [~,seq_p_in_out(unit,step)] = kstest2(nanmean(fixation_locked_inside(:,ind)'),nanmean(fixation_locked_outside(:,ind)'));
        %         end
        ind = statistical_window;
        [~,seq_p_in_out(unit)] = kstest2(nanmean(fixation_locked_inside(:,ind)'),nanmean(fixation_locked_outside(:,ind)'));
        
        %         sig_time_in_out = zeros(1,twin1+twin2);
        %         for step = 1:((twin1+twin2)/sliding_step)-sliding_window/sliding_step+1
        %             ind = sliding_step*(step-1)+1:sliding_window+sliding_step*(step-1);
        %             if seq_p_in_out(unit,step) < p_thresh
        %                 sig_time_in_out(ind) = 2;
        %             end
        %         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Compare Sequence Firing Rates for Items Inside vs Outside Field---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        colors = ['rgbm'];
        t = -twin1:twin2-1;
        
        figure
        a = subplot(2,3,1);
        h = imagesc(firing_rate_map);
        set(h,'alphadata',~isnan(firing_rate_map));
        axis off
        axis equal
        colormap(a,'jet')
        colorbar
        caxis([clim(1) maxfr])
        title_str = ['All images, peak rate = ' num2str(max(firing_rate_map(:)),3) ' Hz\n'];
        if spatial_info.shuffled_rate_prctile(unit) >= 95;
            title_str = [title_str 'Bits = ' num2str(spatial_info.rate(unit),3) ...
                '(' num2str(spatial_info.shuffled_rate_prctile(unit),3) '%%) '];
        end
        if (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
            title_str = [title_str '\n \\rho_{1/2} = ' num2str(spatial_info.spatialstability_halves(1,unit),2) ...
                '(' num2str(spatial_info.spatialstability_halves_prctile(1,unit),3) '%%)'];
        end
        if (spatial_info.spatialstability_even_odd_prctile(1,unit) > 95)
            title_str = [title_str ' \\rho_{e/o} = ' num2str(spatial_info.spatialstability_even_odd(1,unit),2) ...
                '(' num2str(spatial_info.spatialstability_even_odd_prctile(1,unit),3) '%%)'];
            
        end
        title(sprintf(title_str));
        
        
        b = subplot(2,3,2);
        imagesc(place_field_matrix);
        hold on
        for c = 1:4
            for seq = 1:2
                plot(sequence_locations{seq}(1,c),imgy-sequence_locations{seq}(2,c),[colors(c) shapes(seq)],'markersize',12)
            end
        end
        hold off
        xlim([0 800])
        ylim([0 600])
        axis equal
        axis off
        colormap(b,'gray')
        title(sprintf(['place field Location \n Area = ' num2str(area(unit),2)]));
        
        subplot(2,3,3)
        
        hold on
        for seq = 1:2
            for c = 1:4
                [Dmean,~] = nandens (sequence_fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == seq ,:), smval, 'gauss',Fs);
                plot(t(1:10:end),Dmean(1:10:end),[colors(c) shapes(seq) '-']);
            end
        end
        xlabel('Time from Fixation (ms)')
        ylabel('Firing Rate (Hz)')
        yl = ylim;
        if yl(1) < 0
            yl(1) = 0;
            ylim(yl);
        end
        plot([0 0],[yl(1) yl(2)],'k--')
        hold off
        title('Firing Rate Curves for Each Item')
        
        subplot(2,3,4)
        hold on
        [~,~,~,y_seq,~] = dofill(t,fixation_locked_inside,'blue',1,smval);
        dofill(t,fixation_locked_outside,'red',1,smval);
        yl = ylim;
        if yl(1) < 0
            yl(1) = 0;
            ylim(yl);
        end
        if seq_p_in_out(unit) < p_thresh
            h = fill([twin1 twin1+twin2  twin1+twin2 twin1 twin1]-twin1,...
                [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
            set(h,'facealpha',.25,'EdgeColor','None')
        end
        plot([0 0],[yl(1) yl(2)],'k--')
        xlabel('Time from Fixation (ms)')
        ylabel('Firing Rate (Hz)')
        hold off
        legend('Items inside','Items Outside','Location','NorthWest')
        
        %find peaks for list task
        [pks,locs] = findpeaks(y_seq,'MinPeakWidth',50);
        pks(locs < twin1) = [];
        locs(locs < twin1) = [];
        seq_pk = locs(pks == max(pks));
        seq_max = max(pks);
        title(['Sequence Trials: peak of ' num2str(seq_max,3) 'Hz @ ' num2str(seq_pk-twin1)])
        
        
        subplot(2,3,5)
        hold on
        [~,~,~,y_list,~]=dofill(t,list_fixation_locked_firing{unit}(in_out{unit} == 1,:),'blue',1,smval);%out-> in
        dofill(t,list_fixation_locked_firing{unit}(in_out{unit} == 4,:),'red',1,smval);%out->out
        yl = ylim;
        if yl(1) < 0
            yl(1) = 0;
            ylim(yl);
        end
        plot([0 0],[yl(1) yl(2)],'k--')
        if p_in_out(unit,2) < p_thresh
            h = fill([twin1 twin1+twin2  twin1+twin2 twin1 twin1]-twin1,...
                [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
            set(h,'facealpha',.25,'EdgeColor','None')
        end
        hold off
        xlabel('Time from Fixation Onset (ms)')
        ylabel('Firing Rate (Hz)')
        legend('out->in','in->in','out->out','Location','NorthWest')
        
        title(['List Image Trials: peak of ' num2str(list_max,3) 'Hz @ ' num2str(list_pk-twin1)])
        
        subplot(2,3,6)
        hold on
        dofill(t,list_fixation_locked_firing{unit}(in_out{unit} == 1,:),'black',1,smval);%list out-> in
        dofill(t,fixation_locked_inside,'green',1,smval);%sequence
        yl = ylim;
        if yl(1) < 0
            yl(1) = 0;
            ylim(yl);
        end
        plot([0 0],[yl(1) yl(2)],'k--')
        hold off
        xlabel('Time from Fixation Onset (ms)')
        ylabel('Firing Rate (Hz)')
        legend('List','Sequence','Location','NorthWest')
        
        if ~isempty(list_pk)
            list_pk = list_pk(1);
            list_max = list_max(1);
            contextual_gain(1,unit) = list_pk;
            contextual_gain(2,unit) = list_max;
        end
        if ~isempty(seq_pk)
            seq_pk = seq_pk(1);
            seq_max = seq_max(1);
            contextual_gain(3,unit) = seq_pk;
            contextual_gain(4,unit) = seq_max;
        end
        if ~isempty(list_pk) && ~isempty(seq_pk)
            contextual_gain(5,unit) = list_max/seq_max;
            title(['Contextual Gain: ' num2str(100*(list_max-seq_max)/list_max,3)])
        end
        
        subtitle([task_file(1:8) ' ' unit_stats{1,unit}]);
        if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95) ...
                && (spatial_info.spatialstability_even_odd_prctile(1,unit) > 95);
            save_and_close_fig(best_place_cell_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_Sequence_InSideOutside']);
        elseif (spatial_info.shuffled_rate_prctile(unit) > 95) && ((spatial_info.spatialstability_halves_prctile(1,unit) > 95) ...
                || (spatial_info.spatialstability_even_odd_prctile(1,unit) > 95));
            save_and_close_fig(place_cell_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_Sequence_InSideOutside']);
        elseif (spatial_info.shuffled_rate_prctile(unit) > 99)
            save_and_close_fig( skagg_cell_dir99,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_Sequence_InSideOutside']);
        elseif (spatial_info.shuffled_rate_prctile(unit) > 95)
            save_and_close_fig(skagg_cell_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_Sequence_InSideOutside']);
        else
            save_and_close_fig(spatial_corr_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_Sequence_InSideOutside']);
        end
    end
end

save([data_dir task_file(1:8) '-Place_Cell_Analysis.mat'],...
    'twin1','twin2','smval','list_fixation_locked_firing','in_out','p_in_out',...
    'sliding_step','sliding_window','p_thresh','task_file','area','contextual_gain',...
    'unit_stats','sequence_fixation_locked_firing',...
    'all_place_field_matrix','seq_p_in_out','statistical_window')

end