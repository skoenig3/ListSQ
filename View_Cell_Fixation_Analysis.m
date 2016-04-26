function  View_Cell_Fixation_Analysis(data_dir,figure_dir,session_data,task)
%written by Seth Konig 4/25/16
%code anlayzes units with significant (95+ percentile) skaggs information
%score, Miriam's spatial stability score, or both. Code determines/verifies
%that firing rate within the view field is great than firing rate outside
%of the view field. Additionally, the code determines whether there is a
%differeince in firing rate for novel vs repeat images. A common theme of
%view cells (based on an initial pass and a logical assumption one would make)
%is that they show strong perisaccadic modulation usually firing more after
%a few fixation has started (hence based on viewing location). This code
%compared in-and-out of field firing and novel-vs-repeat firing based on
%firing rates locked to fixations (could also try saccades).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set important task parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
twin = 750;% here how much after image onset to ignore data
view_cell_dir = [figure_dir 'View Cells\'];
imageX = 800;
imageY = 600;
img_on_code = 23;
img_off_code = 24;
ITIstart_code = 15;
minimum_fix_duration = 100;%miniumum fixation duration to look at data for
smval = 60;
Fs = 1000;
sliding_window = 100; %width in ms
sliding_step = 10;  %step size in ms
p_thresh = 0.001;%significance level

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat']);
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'fixationstats','cfg','data','hdr','valid_trials');

[multiunit,~,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);

%set the image duration
if str2double(task_file(3:8)) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end

%get important task specific information
[itmlist,sequence_items] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);

%NaNs are for start and end trials otherwise cut
num_trials = length(cfg.trl);
valid_trials(1,isnan(valid_trials(1,:))) = 1;
valid_trials(2,isnan(valid_trials(2,:))) = num_trials;

%filter parameters
filter_size = filter_width*10;
H = fspecial('gaussian',filter_size,filter_width);

fixation_locked_firing = cell(1,num_units); %firing rate locked to fixations
in_out = cell(1,num_units); %firing rate for fixations inside (1) or outside (0) of view filed
view_cell_nvr = cell(1,num_units); %firing rate for fixations during novel (1) and repeat (2) images only inside
area = NaN(1,num_units);
p_in_out = NaN(num_units,2*twin/sliding_step); %is firing rate higher inside or outside the firing field
p_nov_rep = NaN(num_units,2*twin/sliding_step); %is firing rate inside the field different for novel vs repeat
for unit = 1:num_units
    if all(isnan(peak_firing_rate(:,unit))) || all(peak_firing_rate(:,unit) < 1)
        continue %unit doesn't fire enough go to next unit
    elseif all(spatial_info.shuffled_rate_prctile(:,unit) < 95) && all(spatial_info.shuffled_spatialstability_prctile(:,unit) < 95)
        continue %unit likely not spatial
    end
    
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
        filtered_time(filtered_time < 0.025) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
        filtered_space = filter_space(eyepos{unit},spike_times{unit},imageX,imageY,binsize,H);
    elseif condition == 1 %if condition is 1 then for novel images
        filtered_time = filter_time(select_eyepos(eyepos{unit},nvr{unit}==1),imageX,imageY,Fs,binsize,H);
        filtered_time(filtered_time < 0.025) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
        filtered_space = filter_space(select_eyepos(eyepos{unit},nvr{unit}==1),spike_times{unit}(nvr{unit} ==1,:),imageX,imageY,binsize,H);
    else %condition == 2
        filtered_time = filter_time(select_eyepos(eyepos{unit},nvr{unit}==2),imageX,imageY,Fs,binsize,H);
        filtered_time(filtered_time < 0.025) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
        filtered_space = filter_space(select_eyepos(eyepos{unit},nvr{unit}==2),spike_times{unit}(nvr{unit} ==2,:),imageX,imageY,binsize,H);
    end
    
    firing_rate = filtered_space./filtered_time;
    fr = sort(firing_rate(1:end));
    fr(isnan(fr)) = [];
    maxfr = fr(round(0.95*length(fr)));% the ~95%-tile
    
    median_firing_rate = median(fr)+std(fr);%leaving nmae
    [r,c] = find(firing_rate > median_firing_rate);
    firing_ind = sub2ind(size(firing_rate),r,c);
    threshold_matrix = NaN(size(firing_rate));
    threshold_matrix(firing_ind) = 1;
    threshold_matrix = imresize(threshold_matrix,[imageY,imageX],'method','nearest');%upsample to images size in pixel coordinates
    area(unit) = 100*nansum(nansum(threshold_matrix))/(imageX*imageY);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate Firing Rate Locked to Fixations---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fix_in_out = NaN(1,1000); %firing rate for fixations inside (1) or outside (0) of view filed
    fix_nvr = NaN(1,1000);  %firing rate for fixations during novel (1) and repeat (2) images only inside
    fix_locked_firing = NaN(1000,2*twin); %firing rate locked to fixations
    
    fix_ind = 1; %fixation # so can track in variables above
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
                    emailme(['Spatial Analysis Importing Data found 2 image presentations. Img condition ' num2str(img_cnd(img_index(1)))...
                        ' ' task_file])
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
                    %get firing rate locked to fixation
                    fixt = fixationtimes(1,f);%start of fixation
                    fix_spikes = spikes(spikes > fixt-twin & spikes <= fixt+twin)-fixt+twin;
                    temp = zeros(1,twin*2);
                    temp(fix_spikes) = 1;
                    fix_locked_firing(fix_ind,:) = temp;
                    
                    fix_nvr(fix_ind) = novel_vs_repeat(img_index); %novel or repeat image
                    
                    %determine if fixation was in our out of view field
                    fixx = round(fixations(1,f));
                    fixx(fixx < 1) = 1;
                    fixx(fixx > imageX) = imageX;
                    fixy = imageY-round(fixations(2,f));
                    fixy(fixy < 1) = 1;
                    fixy(fixy > imageY) = imageY;
                    if threshold_matrix(fixy,fixx) == 1 %then inside
                        fix_in_out(fix_ind) = 1;
                    else
                        fix_in_out(fix_ind) = 0;
                    end
                    
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
    fixation_locked_firing{unit} = fix_locked_firing;
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
    filtered_time(filtered_time < 0.025) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
    filtered_space = filter_space(eyepos{unit},spike_times{unit},imageX,imageY,binsize,H);
    firing_rate = filtered_space./filtered_time;
    
    subplot(2,3,1)
    h = imagesc(filtered_space./filtered_time);
    set(h,'alphadata',~isnan(filtered_time));
    title('All images')
    axis off
    axis equal
    
    fr = sort(firing_rate(1:end));
    fr(isnan(fr)) = [];
    clims(:,1) = caxis;
    if length(fr) > 20
        clims(2,1) = fr(round(0.95*length(fr)));% the ~95%-tile
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
    filtered_time(filtered_time < 0.025) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
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
        clims(2,2) = fr(round(0.95*length(fr)));% the ~95%-tile
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
    filtered_time(filtered_time < 0.025) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
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
        clims(3,2) = fr(round(0.95*length(fr)));% the ~95%-tile
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
    title(sprintf(['n_{in} = ' num2str(sum(fix_in_out == 1)) ' n_{out} = ' num2str(sum(fix_in_out == 0))]))
    if sum( sig_time_in_out) > 1
        yl = ylim;
        gaps = findgaps(find(sig_time_in_out))-twin;
        for g = 1:size(gaps,1)
            gp = gaps(g,:);
            gp(gp == 0) = [];
            h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)],...
                [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
            uistack(h,'down')
            set(h,'facealpha',.25,'EdgeColor','None')
        end
    end
    xlim([-twin twin]);
    hold off
    
    subplot(2,3,6)
    hold on
    dofill(t,fix_locked_firing(fix_nvr == 1,:),'blue',1,smval);%inside
    dofill(t,fix_locked_firing(fix_nvr == 2,:),'red',1,smval);%outside
    xlabel('Time from Fixation Onset (ms)')
    ylabel('Firing Rate (Hz)')
    legend('Novel','Repeat','Location','NorthEast')
    title(sprintf(['n_{nov} = ' num2str(sum(fix_nvr == 1)) ' n_{rep} = ' num2str(sum(fix_nvr == 2))]))
    if sum(sig_time_nov_rep) > 1
        yl = ylim;
        gaps = findgaps(find(sig_time_nov_rep))-twin;
        for g = 1:size(gaps,1)
            gp = gaps(g,:);
            gp(gp == 0) = [];
            h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)],...
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
    save_and_close_fig(view_cell_dir,[task_file(1:end-11) '_' unit_names{unit} '_View_Cell_analysis']);
end

save([data_dir task_file(1:8) '-View_Cell_Analysis.mat'],...
    'twin','smval','fixation_locked_firing','in_out','view_cell_nvr',...
    'p_nov_rep','p_in_out','sliding_step','sliding_window','p_thresh')
end

function [filtered_time] = filter_time(eyepos,imageX,imageY,Fs,binsize,H)
%calculate the total time spent at any locaitons in binned pixels
spatial_time = time_per_pixel(eyepos,imageX,imageY,Fs);
filtered_time = bin2(spatial_time,binsize,binsize);
filtered_time = imfilter(filtered_time,H);
filtered_time(filtered_time == 0) = NaN; %can cause aribitrarily high firing rates
filtered_time = filtered_time(end:-1:1,:); %flip so rightside up
end

function [filtered_space] = filter_space(eyepos,spike_times,imageX,imageY,binsize,H)
%caluclate total spikes over space
[firing_location] = pixel_spike_location(eyepos,spike_times,imageX,imageY);
filtered_space = bin2(firing_location,binsize,binsize);
filtered_space = imfilter(filtered_space,H);
% filtered_space(filtered_space == 0) = NaN;
filtered_space = filtered_space(end:-1:1,:); %flip so rightside up
end

function [selected] = select_eyepos(eyepos,select_rows)
%x eye eyepos is in odd rows and y ye eyepos is in even rows
%select rows should be a logical index
x = eyepos(1:2:end,:);
y = eyepos(2:2:end,:);
x = x(select_rows,:);
y = y(select_rows,:);
selected = NaN(sum(select_rows)*2,size(x,2));
selected(1:2:end,:) = x;
selected(2:2:end,:) = y;
end
