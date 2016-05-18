function  View_Cell_Fixation_Analysis(data_dir,figure_dir,session_data,task)
%written by Seth Konig 4/25/16
%code only analyzes units with significant (95+ percentile) skaggs information
%score and Miriam's spatial stability score. Code determines/verifies
%that firing rate within the view field is great than firing rate outside
%of the view field. Additionally, the code determines whether there is a
%differeince in firing rate for novel vs repeat images. A common theme of
%view cells (based on an initial pass and a logical assumption one would make)
%is that they show strong perisaccadic modulation usually firing more after
%a few fixation has started (hence based on viewing location). This code
%compared in-and-out of field firing and novel-vs-repeat firing based on
%firing rates locked to fixations (could also try saccades). This code also
%analyzes firing rates during sequence trials for these neurons to
%determine if there is a correlation in firing rate of the view cell during
% images and during sequence trials.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set important task parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
twin = 500;% here how much after image onset to ignore data
view_cell_dir = [figure_dir 'View Cells\'];
imageX = 800;
imageY = 600;
img_on_code = 23;
img_off_code = 24;
ITIstart_code = 15;
minimum_fix_duration = 100;%miniumum fixation duration to look at data for
smval = 60;
Fs = 1000;
numshuffs = 0;
sliding_window = 100; %width in ms
sliding_step = 10;  %step size in ms
p_thresh = 0.05;%significance level
min_bin_dur = 0.100; %minimum of 500 ms in each bin to use so no outlier
fixwin = 5;%size of fixation window on each crosshair
event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
smval = 60;%gaussian 1/2 width for smoothing

task_file = get_task_data(session_data,task);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat']);
load([data_dir task_file(1:end-11) '-preprocessed.mat']);
absolute_fixationstats = fixationstats;
absolute_cfg = cfg;

[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
[multiunit,~,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);

%set the image duration
if str2double(task_file(3:8)) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end

%get important task specific information
%get important task specific information
[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
overlap = find((sequence_locations{1}(1,:) ==  sequence_locations{2}(1,:)) & ...
    (sequence_locations{1}(2,:) ==  sequence_locations{2}(2,:)));
[which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);

%NaNs are for start and end trials otherwise cut
num_trials = length(cfg.trl);
valid_trials(1,isnan(valid_trials(1,:))) = 1;
valid_trials(2,isnan(valid_trials(2,:))) = num_trials;

%these are the absolute minimum data required to do data analysis may want
%to be more strigent later but not worth doing analysis (especially
%shuffling) on such few trials for these neurons
if str2double(task_file(3:8)) < 140805 %7 second viewing with fewere sequence trials
    minimum_trials_1 = 101; %for session start minimum number of trials: includes fam block and 16 novel/repeat images + sequence trials
    minimum_trials_2 =  80;%for rest of session: includes 16 novel/repeat images + sequence trials
else
    minimum_trials_1 = 117; %for session start minimum number of trials: includes fam block and 16 novel/repeat images + sequence trials
    minimum_trials_2 =  96;%for rest of session: includes 16 novel/repeat images + sequence trials
end


for unit = 1:size(valid_trials,2)
    %---determine number of blocks unit was stable for
    start_end = valid_trials(:,unit);
    if isnan(start_end(1))
        start_end(1) = 1;
    end
    if isnan(start_end(2))
        start_end(2) = length(cfg.trl);
    end
    start_end(start_end == 0) = 1;
    min_trial = cfg.trl(start_end(1)).cnd-1000; %get the first condition number
    max_trial = cfg.trl(start_end(2)).cnd-1000; %get the last condition number
    
    if min_trial < 22 %includes fam block
        min_trial = 22; %remove then count from there
    end
    num_blks = floor((max_trial-min_trial+1)/minimum_trials_2);
    
    if num_blks < 2
        valid_trials(:,unit) = 0;
    end
end


%filter parameters
filter_size = filter_width*10;
H = fspecial('gaussian',filter_size,filter_width);

%---Pre-allocate space for Fixations Inside vs Outside Field Analysis---%
list_fixation_locked_firing = cell(1,num_units); %firing rate locked to fixations
in_out = cell(1,num_units); %firing rate for fixations inside (1) or outside (0) of view filed
view_cell_nvr = cell(1,num_units); %firing rate for fixations during novel (1) and repeat (2) images only inside
area = NaN(1,num_units);
p_in_out = NaN(num_units,2*twin/sliding_step); %is firing rate higher inside or outside the firing field
p_nov_rep = NaN(num_units,2*twin/sliding_step); %is firing rate inside the field different for novel vs repeat
n_in_not_first = 0; %number of fixations inside view field but prior fixations was also inside

%---Pre-allocate space for Fixations during List vs Sequence Trials Analysis---%
fixation_locked_firing = cell(4,num_units);
trial_nums = cell(1,num_units);
corr_coeff_seqeuence_list = NaN(1,num_units);
p_corr_seqeuence_list = NaN(1,num_units);
slope_seqeuence_list = NaN(1,num_units);
shuffled_slope = cell(1,num_units);
shuffled_corr_coeff = cell(1,num_units);
sequences_inside = NaN(2,4);
seq_p_in_out = NaN(num_units,2*twin/sliding_step); %is firing rate higher inside or outside the firing field

go_no_go = 0; %did we even do an analysis for this file
for unit = 1:num_units
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Unit Passes 95% for both Skaggs and Stability--%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if all(isnan(peak_firing_rate(:,unit))) || all(peak_firing_rate(:,unit) < 1)
        continue %unit doesn't fire enough go to next unit
    elseif all((spatial_info.shuffled_rate_prctile(:,unit) < 95) & (spatial_info.shuffled_spatialstability_prctile(:,unit) < 95))
        continue %unit likely not spatial
    elseif  ~any((spatial_info.shuffled_rate_prctile(:,unit) > 95) & (spatial_info.shuffled_spatialstability_prctile(:,unit) > 95))
        continue %only process cells that pass both criterion under at least 1 condition
    end
    
    go_no_go = 1; %yes we are doing analyses
    %copy jittered spike location plots and heat maps to view cell folder
    %     try
    %         if multiunit(unit)
    %             copyfile([figure_dir '\MultiUnit\' task_file(1:end-11) '_' unit_names{unit} '_List_spatial_analysis.bmp'],...
    %                 [view_cell_dir task_file(1:end-11) '_' unit_names{unit} '_List_spatial_analysis.bmp']);
    %             copyfile([figure_dir '\MultiUnit\' task_file(1:end-11) '_' unit_names{unit} '_List_spatial_analysis.fig'],...
    %                 [view_cell_dir task_file(1:end-11) '_' unit_names{unit} '_List_spatial_analysis.fig']);
    %         else
    %             copyfile([figure_dir task_file(1:end-11) '_' unit_names{unit} '_List_spatial_analysis.bmp'],...
    %                 [view_cell_dir task_file(1:end-11) '_' unit_names{unit} '_List_spatial_analysis.bmp']);
    %             copyfile([figure_dir task_file(1:end-11) '_' unit_names{unit} '_List_spatial_analysis.fig'],...
    %                 [view_cell_dir task_file(1:end-11) '_' unit_names{unit} '_List_spatial_analysis.fig']);
    %         end
    %     end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine Firing rate Locked to Fixations Inside vs Outside Field---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate View Field Location---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %draw spikes in locations above threshold (currently > 20% of max)
    which_condition_skaggs = find(spatial_info.shuffled_spatialstability_prctile(:,unit) >= 95);
    which_condition_stable = find(spatial_info.shuffled_spatialstability_prctile(:,unit) >= 95);
    which_condition_both = which_condition_skaggs & which_condition_stable;
    if sum(which_condition_both > 1) %if more than 1 condition produces significant results select all images
        condition = 3;
    elseif sum(which_condition_both) == 1;
        %Calculated smoothed firing rate for all images
        condition = find(which_condition_both);
    else %both scores are not significant
        if ~isempty(which_condition_skaggs)
            if length(which_condition_skaggs) > 1
                condition = 3;
            else
                condition = which_condition_skaggs;
            end
        else
            if length(which_condition_stable) > 1
                condition = 3;
            else
                condition = which_condition_stable;
            end
        end
    end
    
    if condition == 3 %if condition is 3 then all images
        filtered_time = filter_time(eyepos{unit},imageX,imageY,Fs,binsize,H);
        filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
        filtered_space = filter_space(eyepos{unit},spike_times{unit},imageX,imageY,binsize,H);
    elseif condition == 1 %if condition is 1 then for novel images
        filtered_time = filter_time(select_eyepos(eyepos{unit},nvr{unit}==1),imageX,imageY,Fs,binsize,H);
        filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
        filtered_space = filter_space(select_eyepos(eyepos{unit},nvr{unit}==1),spike_times{unit}(nvr{unit} ==1,:),imageX,imageY,binsize,H);
    else %condition == 2
        filtered_time = filter_time(select_eyepos(eyepos{unit},nvr{unit}==2),imageX,imageY,Fs,binsize,H);
        filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
        filtered_space = filter_space(select_eyepos(eyepos{unit},nvr{unit}==2),spike_times{unit}(nvr{unit} ==2,:),imageX,imageY,binsize,H);
    end
    
    firing_rate = filtered_space./filtered_time;
    fr = sort(firing_rate(1:end));
    fr(isnan(fr)) = [];
    maxfr = fr(round(0.95*length(fr)));% the ~95%-tile
    
    median_firing_rate = median(fr)+0.5*std(fr);%leaving nmae
    [r,c] = find(firing_rate > median_firing_rate);
    firing_ind = sub2ind(size(firing_rate),r,c);
    threshold_matrix = zeros(size(firing_rate));
    threshold_matrix(firing_ind) = 1;
    threshold_matrix(isnan(firing_rate)) = NaN; %remove locations with too little data from further analysis
    if sum(sum(~isnan(threshold_matrix))) < 20
        continue %way too little area to process
    end
    threshold_matrix = imresize(threshold_matrix,[imageY,imageX],'method','nearest');%upsample to images size in pixel coordinates
    area(unit) = 100*nansum(nansum(threshold_matrix))/sum(sum(~isnan(threshold_matrix)));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate Firing Rate Locked to Fixations---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fix_in_out = NaN(1,1000); %firing rate for fixations inside (1) or outside (0) of view filed
    fix_nvr = NaN(1,1000);  %firing rate for fixations during novel (1) and repeat (2) images only inside
    fix_locked_firing = NaN(1000,2*twin); %firing rate locked to fixations
    
    fix_ind = 1; %fixation # so can track in variables above
    fixationstats = absolute_fixationstats; %reload because written over below
    cfg = absolute_cfg;
    num_trials = length(cfg.trl);
    %easier to start over rather than figure out format from imported data
    for t = 1:num_trials
        if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
            if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start+twin; %want to avoid visual response
                imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start;
                
                % if monkey isn't paying attention data isn't probably
                % worth much plus have to cut off somewhere
                if imgoff-imgon > 1.5*imgdur-1
                    imgoff = imgon+1.5*imgdur-1;
                end
                
                %---Get eye data and remove fixations before/after image period---%
                fixationtimes = fixationstats{t}.fixationtimes;
                fixations = fixationstats{t}.fixations;
                %fixation started before image turned on
                invalid= find(fixationtimes(1,:) < imgon);
                fixationtimes(:,invalid) = [];
                fixations(:,invalid) = [];
                %fixation started after the image turned off and/or firing rate could corrupted by image turning off
                invalid= find(fixationtimes(1,:) > imgoff);
                fixationtimes(:,invalid) = [];
                fixations(:,invalid) = [];
                
                %remove fixations that are too short to really look at for analysis
                fixdur = fixationtimes(2,:)-fixationtimes(1,:)+1;
                fixationtimes(:,fixdur < minimum_fix_duration) = [];
                fixations(:,fixdur < minimum_fix_duration) = [];
                
                img_index = find(img_cnd == cfg.trl(t).cnd);
                if length(img_index) > 1
                    %                     emailme(['Spatial Analysis Importing Data found 2 image presentations. Img condition ' num2str(img_cnd(img_index(1)))...
                    %                         ' ' task_file])
                    %remove that image from analysis. For TO set 17 seems
                    %to be a cortex error...
                    imgnum = which_img(img_index(1));
                    img_index = find(which_img == imgnum);
                    which_img(img_index) = NaN;
                    img_cnd(img_index) = NaN;
                    novel_vs_repeat(img_index) = NaN;
                    img_index = find(img_cnd == cfg.trl(t).cnd);
                end
                if isempty(img_index)
                    continue
                end
                
                
                novel_vs_repeat(img_index);
                spikes = find(data(unit).values{t});
                
                for f = 1:size(fixations,2)
                    
                    if f > 1
                        %determine if prior fixation was in our out of view field
                        fixx = round(fixations(1,f-1));
                        fixx(fixx < 1) = 1;
                        fixx(fixx > imageX) = imageX;
                        fixy = imageY-round(fixations(2,f-1));
                        fixy(fixy < 1) = 1;
                        fixy(fixy > imageY) = imageY;
                        
                        if threshold_matrix(fixy,fixx) == 1 %then prior fixation was already inside
                            n_in_not_first = n_in_not_first+1;
                            continue %move to the next one
                        end
                    end
                    
                    %determine if fixation was in our out of view field
                    fixx = round(fixations(1,f));
                    fixx(fixx < 1) = 1;
                    fixx(fixx > imageX) = imageX;
                    fixy = imageY-round(fixations(2,f));
                    fixy(fixy < 1) = 1;
                    fixy(fixy > imageY) = imageY;
                    if threshold_matrix(fixy,fixx) == 1 %then inside
                        fix_in_out(fix_ind) = 1;
                    elseif threshold_matrix(fixy,fixx) == 0 %then inside, NaNs are for locations not analyzed
                        fix_in_out(fix_ind) = 0;
                    else %not a valid fixation location too sparse of coverage to analyze
                        continue
                    end
                    
                    
                    %get firing rate locked to fixation
                    fixt = fixationtimes(1,f);%start of fixation
                    fix_spikes = spikes(spikes > fixt-twin & spikes <= fixt+twin)-fixt+twin;
                    temp = zeros(1,twin*2);
                    temp(fix_spikes) = 1;
                    fix_locked_firing(fix_ind,:) = temp;
                    
                    fix_nvr(fix_ind) = novel_vs_repeat(img_index); %novel or repeat image
                    fix_ind = fix_ind+1;
                end
            end
        end
    end
    
    %remove excess NaNs associated with error trials
    fix_in_out = laundry(fix_in_out);
    fix_nvr = laundry(fix_nvr);
    fix_locked_firing = laundry(fix_locked_firing);
    
    %store variables across units for later access
    list_fixation_locked_firing{unit} = fix_locked_firing;
    in_out{unit} = fix_in_out;
    view_cell_nvr{unit} = fix_nvr;
    
    for step = 1:(2*twin/sliding_step)-sliding_window/sliding_step+1
        ind = sliding_step*(step-1)+1:sliding_window+sliding_step*(step-1);
        [~,p_in_out(unit,step)] = ttest2(nanmean(fix_locked_firing(fix_in_out == 1,ind)'),nanmean(fix_locked_firing(fix_in_out == 0,ind)'));
        [~,p_nov_rep(unit,step)] = ttest2(nanmean(fix_locked_firing(fix_nvr == 1,ind)'),nanmean(fix_locked_firing(fix_nvr == 2,ind)'));
    end
    
    sig_time_in_out = zeros(1,2*twin);
    sig_time_nov_rep = zeros(1,2*twin);
    for step = 1:(2*twin/sliding_step)-sliding_window/sliding_step+1
        ind = sliding_step*(step-1)+1:sliding_window+sliding_step*(step-1);
        if p_in_out(unit,step) < p_thresh
            sig_time_in_out(ind) = 2;
        end
        if p_nov_rep(unit,step) < p_thresh
            sig_time_nov_rep(ind) = 2;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    %%%---Make Plot---%%%
    %%%%%%%%%%%%%%%%%%%%%
    figure
    clims = NaN(2,3);
    
    %---plot firing rate map for all images---%
    filtered_time = filter_time(eyepos{unit},imageX,imageY,Fs,binsize,H);
    filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
    filtered_space = filter_space(eyepos{unit},spike_times{unit},imageX,imageY,binsize,H);
    firing_rate = filtered_space./filtered_time;
    
    subplot(2,3,1)
    h = imagesc(filtered_space./filtered_time);
    set(h,'alphadata',~isnan(filtered_time));
    title('All images')
    axis off
    axis equal
    colorbar
    
    fr = sort(firing_rate(1:end));
    fr(isnan(fr)) = [];
    clims(:,1) = caxis;
    if length(fr) > 20
        clims(2,1) = fr(round(0.90*length(fr)));% the ~95%-tile
    end
    
    title_str = ['All images, peak rate = ' num2str(max(fr),3) ' Hz\n'];
    if spatial_info.shuffled_rate_prctile(3,unit) >= 95;
        title_str = [title_str 'Bits = ' num2str(spatial_info.rate(3,unit),3) ...
            '(' num2str(spatial_info.shuffled_rate_prctile(3,unit),3) '%%) '];
    end
    if spatial_info.shuffled_spatialstability_prctile(3,unit) >= 95
        title_str = [title_str 'r = ' num2str(spatial_info.spatialstability(3,unit),2) ...
            '(' num2str(spatial_info.shuffled_spatialstability_prctile(3,unit),3) '%%)'];
    end
    title(sprintf(title_str));
    
    
    %---plot firing rate map for novel images---%
    filtered_time = filter_time(select_eyepos(eyepos{unit},nvr{unit}==1),imageX,imageY,Fs,binsize,H);
    filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
    filtered_space = filter_space(select_eyepos(eyepos{unit},nvr{unit}==1),spike_times{unit}(nvr{unit} ==1,:),imageX,imageY,binsize,H);
    subplot(2,3,2)
    h = imagesc(filtered_space./filtered_time);
    set(h,'alphadata',~isnan(filtered_time));
    title('All images')
    axis off
    axis equal
    
    fr = filtered_space./filtered_time;
    fr = sort(fr(1:end));
    fr(isnan(fr)) = [];
    clims(:,2) = caxis;
    if length(fr) > 20
        clims(2,2) = fr(round(0.90*length(fr)));% the ~95%-tile
    end
    
    title_str = ['Novel images, peak rate = ' num2str(max(fr),3) ' Hz\n'];
    if spatial_info.shuffled_rate_prctile(1,unit) >= 95;
        title_str = [title_str 'Bits = ' num2str(spatial_info.rate(1,unit),3) ...
            '(' num2str(spatial_info.shuffled_rate_prctile(1,unit),3) '%%) '];
    end
    if spatial_info.shuffled_spatialstability_prctile(1,unit) >= 95
        title_str = [title_str 'r = ' num2str(spatial_info.spatialstability(1,unit),2) ...
            '(' num2str(spatial_info.shuffled_spatialstability_prctile(1,unit),3) '%%)'];
    end
    title(sprintf(title_str));
    
    
    %---plot firing rate map for repeat images---%
    filtered_time = filter_time(select_eyepos(eyepos{unit},nvr{unit}==2),imageX,imageY,Fs,binsize,H);
    filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
    filtered_space = filter_space(select_eyepos(eyepos{unit},nvr{unit}==2),spike_times{unit}(nvr{unit}==2,:),imageX,imageY,binsize,H);
    subplot(2,3,3)
    h = imagesc(filtered_space./filtered_time);
    set(h,'alphadata',~isnan(filtered_time));
    title('All images')
    axis off
    axis equal
    
    fr = filtered_space./filtered_time;
    fr = sort(fr(1:end));
    fr(isnan(fr)) = [];
    clims(:,3) = caxis;
    if length(fr) > 20
        clims(3,2) = fr(round(0.90*length(fr)));% the ~95%-tile
    end
    
    title_str = ['Repeat images, peak rate = ' num2str(max(fr),3) ' Hz\n'];
    if spatial_info.shuffled_rate_prctile(2,unit) >= 95;
        title_str = [title_str 'Bits = ' num2str(spatial_info.rate(2,unit),3) ...
            '(' num2str(spatial_info.shuffled_rate_prctile(2,unit),3) '%%) '];
    end
    if spatial_info.shuffled_spatialstability_prctile(2,unit) >= 95
        title_str = [title_str 'r = ' num2str(spatial_info.spatialstability(2,unit),2) ...
            '(' num2str(spatial_info.shuffled_spatialstability_prctile(2,unit),3) '%%)'];
    end
    title(sprintf(title_str));
    
    minc = nanmin(clims(1,:));
    maxc = nanmax(clims(2,:));
    diffc = maxc-minc;
    minc = minc - 0.1*diffc;
    for sp = 4:6
        subplot(2,3,sp)
        colormap('jet')
        caxis([minc maxc])
        c = colormap;
        c(1,:) = [1 1 1];%turn nans in to white pixels
        colormap(c);
    end
    
    %---draw view field---%
    subplot(2,3,4)
    imagesc(threshold_matrix)
    axis off
    axis equal
    title(sprintf(['View Field Location \n Area = ' num2str(area(unit),2)]));
    
    t = -twin:twin-1;
    subplot(2,3,5)
    hold on
    dofill(t,fix_locked_firing(fix_in_out == 1,:),'blue',1,smval);%inside
    dofill(t,fix_locked_firing(fix_in_out == 0,:),'red',1,smval);%outside
    xlabel('Time from Fixation Onset (ms)')
    ylabel('Firing Rate (Hz)')
    legend('Fix Inside','Fix Outside','Location','NorthEast')
    title(sprintf(['n_{in} = ' num2str(sum(fix_in_out == 1)) ', n_{in but not first} = ' num2str(n_in_not_first) ', n_{out} = ' num2str(sum(fix_in_out == 0))]))
    if sum( sig_time_in_out) > 1
        yl = ylim;
        gaps = findgaps(find(sig_time_in_out));
        for g = 1:size(gaps,1)
            gp = gaps(g,:);
            gp(gp == 0) = [];
            h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin,...
                [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
            uistack(h,'down')
            set(h,'facealpha',.25,'EdgeColor','None')
        end
    end
    xlim([-twin twin]);
    hold off
    
    subplot(2,3,6)
    hold on
    dofill(t,fix_locked_firing(fix_nvr == 1 & fix_in_out == 1,:),'blue',1,smval);%inside
    dofill(t,fix_locked_firing(fix_nvr == 2 & fix_in_out == 1,:),'red',1,smval);%outside
    xlabel('Time from Fixation Onset (ms)')
    ylabel('Firing Rate (Hz)')
    legend('Novel','Repeat','Location','NorthEast')
    title(sprintf(['n_{nov} = ' num2str(sum(fix_nvr == 1 & fix_in_out == 1)) ' n_{rep} = ' num2str(sum(fix_nvr == 2 & fix_in_out == 1))]))
    if sum(sig_time_nov_rep) > 1
        yl = ylim;
        gaps = findgaps(find(sig_time_nov_rep));
        for g = 1:size(gaps,1)
            gp = gaps(g,:);
            gp(gp == 0) = [];
            h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin,...
                [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
            uistack(h,'bottom')
            set(h,'facealpha',.25,'EdgeColor','None')
        end
    end
    xlim([-twin twin])
    hold off
    
    num_trials_str = [' n_ = ' num2str(size(spike_times{unit},1))];
    if multiunit(unit)
        multi_str = 'Multiunit ';
    else
        multi_str = ' ';
    end
    
    subtitle(['Spatial Plots' num_trials_str multi_str unit_names{unit}]);
    save_and_close_fig(view_cell_dir,[task_file(1:end-11) '_' unit_names{unit} '_View_Cell_fixation_analysis']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Correlation during List and Sequence trials---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %preallocate space and parallel structure of cfg
    successful_sequence_trials = NaN(1,length(cfg.trl));
    which_sequence = NaN(1,length(cfg.trl));
    for t = 1:length(cfg.trl);
        if sum(cfg.trl(t).allval == 3) >= 6; %in which sequence trials were rewarded
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
    fixstats = fixationstats;
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
        fixation_locked_firing{c,unit} = NaN(num_trials,twin*2);
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
                elseif isnan(sact) && ~isnan(fixt)
                    if fixt == 1 %only fixation is already started when the trial starts
                        %disp('Fixation already on item when trial started. Looking back to previous trial')
                        if successful_sequence_trials(trial) > 1 %can't look back if it's on the first trial there's "no" eye data
                            if ~isnan(fixstats{successful_sequence_trials(trial)-1}.XY(1,end)) %is the monkey looking at the screen at the end of the trial
                                fixations1 = fixstats{successful_sequence_trials(trial)-1}.fixations(:,end);%last fixation previous trial
                                fixations2 =  fixstats{successful_sequence_trials(trial)}.fixations(:,1);%first fixation this trial
                                fixationtimes1 = fixstats{successful_sequence_trials(trial)-1}.fixationtimes(:,end); %fixation times of last fixation previous trial
                                saccadetimes = fixstats{successful_sequence_trials(trial)-1}.saccadetimes(:,end); %saccade times of last fixation previous trial
                                last_location = fixstats{successful_sequence_trials(trial)-1}.XY(:,end); %last eye location from previous trial
                                if fixationtimes(2,end) > saccadetimes(2,end) %Yes, then did the trial end on fixation and already there
                                    %if it did are the eyes in the same location
                                    if sqrt(sum((fixations1-fixations2).^2)) < 36 %within 1.5 dva so approx yes
                                        spikes2 = find(data(unit).values{successful_sequence_trials(trial)-1});
                                        fixt = fixationtimes1(1);
                                        fix_spikes = spikes2(spikes2 > fixt-twin & spikes2 <= fixt+twin)-fixt+twin;
                                        fixt = fixationtimes1(1)-length(data(unit).values{successful_sequence_trials(trial)-1});%should be negative
                                        fix_spikes = [fix_spikes find(spikes < fixt+twin)-fixt+twin];
                                        temp = zeros(1,twin*2);
                                        temp(fix_spikes) = 1;
                                        fixation_locked_firing{c,unit}(trial,:) = temp;
                                        continue %so we don't overwrite this trial in the next section
                                    else %if not possibly blinked at end of trial no good estimate of eye movements possible
                                        continue
                                    end
                                end
                            else %possibly blinked at end of trial no good estimate of eye movements possible
                                continue
                            end
                        else
                            continue %just ignore this trial then for now
                        end
                    else
                        %warning(['Could not identify saccade but could locate fixation! Trial#:' num2str(trial)])
                    end
                elseif ~isnan(sact) && isnan(fixt)
                    error(['Could not identify fixation but could locate saccade! Trial#:' num2str(trial)])
                end
                
                if ~isnan(fixt)
                    fix_spikes = spikes(spikes > fixt-twin & spikes <= fixt+twin)-fixt+twin;
                    temp = zeros(1,twin*2);
                    temp(fix_spikes) = 1;
                    fixation_locked_firing{c,unit}(trial,:) = temp;
                end
                trial_nums{c,unit}(trial) = trial;
            end
        end
    end
    fixation_locked_firing(:,unit) = laundry(fixation_locked_firing(:,unit));
    trial_nums(:,unit) = laundry(trial_nums(:,unit));
    
    %---Calculate Observed correlation/slope and shuffled corrleation/slope---%
    %using trial by trial permutation
    temp_list = []; %need to define for parfor even if don't use,
    %for sequence locations in list pixels wihtout enough coverage
    
    %observed data
    imgy = size(firing_rate,1); %vertical size of firing rate matrix
    fr_seq = NaN(2,4);
    fr_list = NaN(2,4);
    all_fixation_locked_firing = [];
    all_sequences = [];
    for c = 1:4
        for seq = 1:2
            [y,~] = nandens(fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == seq ,twin:end),smval,'gauss',Fs,'nanflt');%nandens made by nathan killian
            binned_location = round([sequence_locations{seq}(1,c)/binsize,imgy-sequence_locations{seq}(2,c)/binsize]);
            fr_seq(seq,c) =  mean(y(75:300));%alternatively the max(y);
            fr_list(seq,c) = firing_rate(binned_location(2),binned_location(1));
            if isnan(fr_list(seq,c)) %check surrounding bins since binning and rounding
                for i = -1:1
                    for ii = -1:1
                        if ~isnan(firing_rate(binned_location(2)+i,binned_location(1)+ii))
                            fr_list(seq,c) = firing_rate(binned_location(2)+i,binned_location(1)+ii);
                            break
                        end
                    end
                end
            end
            
            all_fixation_locked_firing = [all_fixation_locked_firing; fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}),twin:end)];
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
    
    shl = NaN(1,numshuffs); %shuffled slope values
    shc = NaN(1,numshuffs); %shuffled correlation values
    %permutation analysis
    
    
    num_seq = length(all_sequences);
    parfor shuff = 1:numshuffs
        randseq = all_sequences(randperm(num_seq));
        
        rand_fr = NaN(8,1);
        for s = 1:8
            [y,~] = nandens(all_fixation_locked_firing(randseq==s,:),smval,'gauss',Fs,'nanflt');%nandens made by nathan killian
            rand_fr(s) = max(y);
        end
        
        if any(isnan(fr_list(:)))
            rand_fr(isnan(fr_list(:))) = [];
            sr = corrcoef(rand_fr,temp_list);
            P = polyfit(rand_fr,temp_list,1);
        else
            sr = corrcoef(rand_fr,fr_list(:));
            P = polyfit(rand_fr,fr_list(:),1);
        end
        shc(shuff) = sr(2);
        shl(shuff) = P(1);
        
    end
    shuffled_slope{unit} = shl;
    shuffled_corr_coeff{unit} = shc;
    
    
    
    colors =['rgkm'];
    shapes = ['xo'];
    figure
    
    subplot(2,2,1)
    filtered_time = filter_time(eyepos{unit},imageX,imageY,Fs,binsize,H);
    filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
    filtered_space = filter_space(eyepos{unit},spike_times{unit},imageX,imageY,binsize,H);
    firing_rate = filtered_space./filtered_time;
    
    h = imagesc(firing_rate);
    set(h,'alphadata',~isnan(filtered_time));
    title('All images')
    axis off
    axis equal
    hold on
    for c = 1:4
        for seq = 1:2
            plot(sequence_locations{seq}(1,c)/binsize,imgy-sequence_locations{seq}(2,c)/binsize,[colors(c) shapes(seq)],'markersize',12)
        end
    end
    hold off
    colormap('jet')
    colorbar
    
    ylims = NaN(2,2);
    subplot(2,2,2)
    t = -twin:twin-1;
    hold on
    seq = 1;
    for c = 1:4
        dofill(t,fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == seq ,:),colors(c),1,smval); %reactive sequence 1
    end
    xlabel('Time from Fixation (ms)')
    ylabel('Firing Rate (Hz)')
    ylims(1,:) = ylim;
    hold off
    title(['Sequence 1 (' shapes(1) ')'])
    
    subplot(2,2,4)
    t = -twin:twin-1;
    hold on
    seq = 2;
    for c = 1:4
        dofill(t,fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == seq ,:),colors(c),1,smval); %reactive sequence 1
    end
    xlabel('Time from Fixation (ms)')
    ylabel('Firing Rate (Hz)')
    ylims(1,:) = ylim;
    hold off
    title(['Sequence 2 (' shapes(2) ')'])
    
    ymin = min(ylims(:,1));
    ymin(ymin < 0) = 0;
    ymax = 1.2*max(ylims(:,2));
    
    subplot(2,2,2)
    ylim([ymin ymax])
    subplot(2,2,4)
    ylim([ymin ymax])
    
    subplot(2,2,3)
    hold on
    for c = 1:4
        for seq = 1:2
            plot(fr_seq(seq,c),fr_list(seq,c),[colors(c) shapes(seq)],'markersize',12)
        end
    end
    plot([min(fr_seq) max(fr_seq)],polyval(P_fit,[min(fr_seq) max(fr_seq)]),'k')
    hold off
    xlabel('Sequence Peak Firing rate (Hz)')
    ylabel('List Firing Rate (after fixation onset) (Hz)')
    title_str = ['r = ' num2str(corr_coeff_seqeuence_list(unit),2)];
    if p_corr_seqeuence_list(unit) < 0.05
        title_str = [title_str '(p = ' num2str(p_corr_seqeuence_list(unit),2) ')'];
    end
    if corr_coeff_seqeuence_list(unit) > prctile(shuffled_corr_coeff{unit},95)
        percentile = 100*sum(corr_coeff_seqeuence_list(unit) > shuffled_corr_coeff{unit})/numshuffs;
        title_str = [title_str '(' num2str(percentile,3) '%%)'];
    end
    title_str = [title_str ' m = ' num2str(slope_seqeuence_list(unit),3)];
    if slope_seqeuence_list(unit) > prctile(shuffled_slope{unit},95)
        percentile = 100*sum(slope_seqeuence_list(unit) > shuffled_slope{unit})/numshuffs;
        title_str = [title_str '(' num2str(percentile,3) '%%)'];
    end
    title(sprintf(title_str))
    
    save_and_close_fig(view_cell_dir,[task_file(1:end-11) '_' unit_names{unit} '_View_Cell_Sequence_fixation_analysis'])
    
    colors =['rgcm'];
    imgy = size(threshold_matrix,1);
    
    %Determine if any of the items are inside the matrix or not?
    sequence_inside = NaN(2,4);
    for c = 1:4
        for seq = 1:2
            if threshold_matrix(imgy-sequence_locations{seq}(2,c),sequence_locations{seq}(1,c)) == 1
                sequence_inside(seq,c) =1;
            end
        end
    end
    
    if any(sequence_inside(:) == 1)
        sequence_insides = sequence_inside;
        
        fixation_locked_inside = [];
        fixation_locked_outside = [];
        for c = 1:4
            for seq = 1:2
                if  sequence_inside(seq,c) == 1;
                    fixation_locked_inside =[fixation_locked_inside; fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == seq ,:)];
                else
                    fixation_locked_outside =[fixation_locked_outside; fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == seq ,:)];
                end
            end
        end
        
        for step = 1:(2*twin/sliding_step)-sliding_window/sliding_step+1
            ind = sliding_step*(step-1)+1:sliding_window+sliding_step*(step-1);
            [~,seq_p_in_out(unit,step)] = ttest2(nanmean(fixation_locked_inside(:,ind)'),nanmean(fixation_locked_outside(:,ind)''),'tail','right');
        end
        
        sig_time_in_out = zeros(1,2*twin);
        for step = 1:(2*twin/sliding_step)-sliding_window/sliding_step+1
            ind = sliding_step*(step-1)+1:sliding_window+sliding_step*(step-1);
            if seq_p_in_out(unit,step) < p_thresh
                sig_time_in_out(ind) = 2;
            end
        end

        figure
        subplot(1,2,1)
        imagesc(threshold_matrix)
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
        colormap('gray')
        title(sprintf(['View Field Location \n Area = ' num2str(area(unit),2)]));
        
        subplot(1,2,2)
        hold on
        dofill(t,fixation_locked_inside,'blue',1,smval);
        dofill(t,fixation_locked_outside,'red',1,smval);
        if sum( sig_time_in_out) > 1
            yl = ylim;
            gaps = findgaps(find(sig_time_in_out));
            for g = 1:size(gaps,1)
                gp = gaps(g,:);
                gp(gp == 0) = [];
                h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin,...
                    [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
                uistack(h,'down')
                set(h,'facealpha',.25,'EdgeColor','None')
            end
        end
        xlabel('Time from Fixation (ms)')
        ylabel('Firing Rate (Hz)')
        hold off
        legend('Items inside','Items Outside')
        axis square
        
        save_and_close_fig(view_cell_dir,[task_file(1:end-11) '_' unit_names{unit} '_View_Cell_Sequence_InSideOutside'])

    end
end

if go_no_go %then we did some analysis so save the data, otherwise don't
    save([data_dir task_file(1:8) '-View_Cell_Analysis.mat'],...
        'twin','smval','fixation_locked_firing','in_out','view_cell_nvr',...
        'p_nov_rep','p_in_out','sliding_step','sliding_window','p_thresh',...
        'n_in_not_first','corr_coeff_seqeuence_list','p_corr_seqeuence_list',...
        'slope_seqeuence_list','shuffled_corr_coeff','shuffled_slope',...
        'list_fixation_locked_firing','multiunit','sequences_inside',...
        'seq_p_in_out','task_file')
end
end