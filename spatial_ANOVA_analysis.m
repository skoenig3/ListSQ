function spatial_ANOVA_analysis(data_dir,preprocessed_data_file,figure_dir)
%written 5/4/15 by Seth Konig
% code analyzes whether fixations in the locations of the items during the
% sequence task are modulated by locaiton using an ANOVA analysis.

tolerance = 2.5*24;%fixwin in dva*pixels/dva
twin = 250; %temporal window

load([data_dir preprocessed_data_file],'data','cfg','multiunit','num_units','item_set','fixationstats');
disp('Collecting Spike Locations')

%set/get some general important info
[eyechans] = find_desired_channels(data(1).hdr,data,'eye');
num_trials = length(cfg.trl);
Fs = data(1).fsample; %should be 1000
trial_start_code = 15;
img_on_code= 23;
img_off_code = 24;
imageX = 800; %horizontal size of the image
imageY = 600; %horizontal size of the image
smval = 60; %for temporal filtering
binsize = 12; %pixels per bin spatial bin in either dimension ~1/2 dva
filter_width = 3; %std of 2D guassian filter ~ 2 dva


%set the image duration
if str2double(preprocessed_data_file(3:8)) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end

%get important task specific information
[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
[which_img,novel_vs_repeat] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);

overlapping_location = find((sequence_locations{1}(1,:) == sequence_locations{2}(1,:)) & ...
    (sequence_locations{1}(2,:) == sequence_locations{2}(2,:)));
all_locations = [sequence_locations{1} sequence_locations{2}];
all_locations(:,overlapping_location) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preallocate space and parallel structure of cfg
num_trials = length(cfg.trl);
image_trials = zeros(length(which_img));
for t = 1:num_trials
    if any(cfg.trl(t).allval == 23) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end);
        image_trials(t) = 1;
    end
end

%remove excess data associated with non image trials
which_img = laundry(which_img);
image_trials = find(image_trials);
num_trials = length(image_trials);


fixation_locked_data = cell(1,25);
saccade_locked_data = cell(1,25);
for fixation = 1:25
    fixation_locked_data{fixation} = NaN(num_trials,1000);
    saccade_locked_data{fixation} = NaN(num_trials,1000);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---process eye data locked to Fixations---%%%
fixationstats = fixationstats(image_trials);

fixation_start_time = cell(1,length(all_locations)); %1 cell per location
for t = 1:num_trials
    fixationtimes = fixationstats{t}.fixationtimes;
    fixations = fixationstats{t}.fixations;
    
    trial_start = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == trial_start_code);
    imgon = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == img_on_code)-trial_start;
    imgoff = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == img_off_code)-trial_start;
    
    %find fiations and saccades that did not occur during the image period;
    %should also take care of the 1st fixation on the crosshair
    invalid= find(fixationtimes(1,:) < imgon);
    fixationtimes(:,invalid) = [];
    fixations(:,invalid) = [];
    invalid= find(fixationtimes(1,:) > imgoff);
    fixationtimes(:,invalid) = [];
    fixations(:,invalid) = [];
    
    for f = 1:size(fixations,2);
        fixinfixwin = find((fixations(1,f) < all_locations(1,:) + tolerance & ...
            fixations(1,f) > all_locations(1,:)-tolerance) & ...
            (fixations(2,f) < all_locations(2,:) + tolerance & ...
            fixations(2,f) > all_locations(2,:)-tolerance));%find fixations within fixation window
        if isempty(fixinfixwin)
            continue
        elseif length(fixinfixwin) > 1
            disp('too many potential locations...skipping this fixation')
            continue
        else
            fixation_start_time{fixinfixwin} = [fixation_start_time{fixinfixwin}; ...
                [t fixationtimes(1,f)]];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Firing Rate Locked to Fixations---%%%
fixation_locked_firing = cell(length(all_locations),num_units);
for c = 1:length(all_locations)
    for t = 1:size(fixation_start_time{c},1);
        trial = fixation_start_time{c}(t,1);
        fixt = fixation_start_time{c}(t,2);
        for unit = 1:num_units;
            spikes = find(data(unit).values{image_trials(trial)});
            fix_spikes = spikes(spikes > fixt-twin & spikes <= fixt+twin)-fixt+twin;
            temp = zeros(1,twin*2);
            temp(fix_spikes) = 1;
            fixation_locked_firing{c,unit} = [fixation_locked_firing{c,unit}; temp];
        end
    end
end

Hz = 1000/(twin*2); %to get firing rate by trial from number of spikes
means = NaN(length(all_locations),num_units);
stds = NaN(length(all_locations),num_units);
numpts = NaN(length(all_locations),num_units);
all_trials = cell(1,num_units);
all_items =  cell(1,num_units);
for c =  1:length(all_locations)
    for unit = 1:num_units
        numpts(c,unit) = size(fixation_locked_firing{c,unit},1);
        mn  = 1000*mean(fixation_locked_firing{c,unit}')*Hz;
        means(c,unit) = mean(mn);
        stds(c,unit)  = std(mn);
        all_trials{unit} = [all_trials{unit}; mn'];
        all_items{unit} = [all_items{unit}; c*ones(length(mn),1)];
    end
end

%get pvals for each unit to determine if fixation locked eye data is
%significnatly modulated by  location
for unit = 1:num_units
    p_anova(unit) = anovan(all_trials{unit},all_items{unit},'display','off');
end

%load spatial analysis data
unit_name = cfg.channel;
load([data_dir preprocessed_data_file(1:10) '-spatial_analysis_results.mat']);
H = fspecial('gaussian',[filter_width*10+5 filter_width*10+5],filter_width);
trial_type = novel_vs_repeat;
item_style = 'XO';

for unit = 1:num_units
    clims = NaN(2,3);
    
    figure
    
    make_spike_jittered_plot(eyepos,spike_times{unit},[3,4],1)
    hold on
    for sequence = 1:2
        for item = 1:4
            plot(sequence_locations{sequence}(1,item),...
                sequence_locations{sequence}(2,item),['k' item_style(sequence)],...
                'MarkerSize',10)
        end
    end
    hold off
    xlim([-50 850])
    ylim([-50 650])
    title('All images')
    set(gca,'Xcolor','w')
    set(gca,'Ycolor','w')
    
    subplot(3,4,2)
    filtered_time = filter_time(eyepos,imageX,imageY,Fs,binsize,H);
    filtered_space = filter_space(eyepos,spike_times{unit},imageX,imageY,binsize,H);
    imagesc(filtered_space./filtered_time);
    hold on
    for sequence = 1:2
        for item = 1:4
            plot(sequence_locations{sequence}(1,item),...
                600-sequence_locations{sequence}(2,item),['k' item_style(sequence)],...
                'MarkerSize',5)
        end
    end
    hold off
    axis off
    if spatial_info.rate(1,unit) > spatial_info.shuffled_95_percentile(1,unit);
        title(['All images, max firing rate = ' num2str(max(max(filtered_space./filtered_time))) ...
            ' Bits = ' num2str(spatial_info.rate(1,unit))])
    else
        title(['All images, max firing rate = ' num2str(max(max(filtered_space./filtered_time)))])
    end
    clims(:,1) = caxis;
    
    subplot(3,4,3)
    filtered_time = filter_time(select_eyepos(eyepos,trial_type == 1),imageX,imageY,Fs,binsize,H);
    filtered_space = filter_space(select_eyepos(eyepos,trial_type == 1),spike_times{unit}(trial_type == 1,:),imageX,imageY,binsize,H);
    imagesc(filtered_space./filtered_time);
    hold on
    for sequence = 1:2
        for item = 1:4
            plot(sequence_locations{sequence}(1,item),...
                600-sequence_locations{sequence}(2,item),['k' item_style(sequence)],...
                'MarkerSize',5)
        end
    end
    hold off
    axis off
    if spatial_info.rate(2,unit) > spatial_info.shuffled_95_percentile(2,unit);
        title(['Novel images, max firing rate = ' num2str(max(max(filtered_space./filtered_time))) ...
            ' Bits = ' num2str(spatial_info.rate(2,unit))])
    else
        title(['Novel images, max firing rate = ' num2str(max(max(filtered_space./filtered_time)))])
    end
    clims(:,2) = caxis;
    
    subplot(3,4,4)
    filtered_time = filter_time(select_eyepos(eyepos,trial_type == 2),imageX,imageY,Fs,binsize,H);
    filtered_space = filter_space(select_eyepos(eyepos,trial_type == 2),spike_times{unit}(trial_type == 2,:),imageX,imageY,binsize,H);
    imagesc(filtered_space./filtered_time);
    hold on
    for sequence = 1:2
        for item = 1:4
            plot(sequence_locations{sequence}(1,item),...
                600-sequence_locations{sequence}(2,item),['k' item_style(sequence)],...
                'MarkerSize',5)
        end
    end
    hold off
    axis off
    if spatial_info.rate(3,unit) > spatial_info.shuffled_95_percentile(3,unit);
        title(['Repeat images, max firing rate = ' num2str(max(max(filtered_space./filtered_time))) ...
            ' Bits = ' num2str(spatial_info.rate(3,unit))])
    else
        title(['Repeat images, max firing rate = ' num2str(max(max(filtered_space./filtered_time)))])
    end
    clims(:,3) = caxis;
    
    for sp = 2:4
        subplot(3,4,sp)
        caxis([min(clims(1,:)) max(clims(2,:))])
        c = colormap;
        c(1,:) = [1 1 1];%turn nans in to white pixels
        colormap(c);
    end
    
    subplot(3,4,5)
    errorbar(means(:,unit),stds(:,unit)./sqrt(numpts(:,unit)));
    xlabel('Item Location')
    ylabel('Firing Rate (Hz)')
    if p_anova(unit) < 0.05
        title([' anova_p = ' num2str(p_anova(unit))]);
    end
    
   t = -twin:twin-1;
   ylimits = [];
    for c = 1:length(all_locations)
        subplot(3,4,5+c)
        dofill(t,fixation_locked_firing{c,unit},'blue',1,smval);
        xlim([-250 250])
        xlabel('Time from Eye Movement (ms)')
        ylabel('Firing Rate (Hz)')
        title(['Item Location #' num2str(c)])
        yl = ylim;
        ylimits(c) = yl(2);
    end
    
    for c = 1:length(all_locations)
        subplot(3,4,5+c)
        ylim([0 max(ylimits)])
    end
    
    save_and_close_fig(figure_dir,[preprocessed_data_file(1:end-11) '_' unit_name{unit} '_List_Spatial_ANOVA_analysis']);
end

save([data_dir preprocessed_data_file(1:8) '-List_Spatial_ANOVA_analysis.mat'],...
    'twin','smval','p_anova','fixation_locked_firing')

end

function make_spike_jittered_plot(eyepos,spike_times,src,subnum)
%src: subplot row and subplot col
%subnum: subplot number
jitter = 6; %1/4 dva

x = eyepos(1:2:end,:);
y = eyepos(2:2:end,:);

[trial,time] = find(spike_times);
spikeind = sub2ind(size(x),trial,time);

xs = x(spikeind);
ys = y(spikeind);
[xs,ys] = remove_nans(xs,ys);

if size(xs,1) == 1; %should only happen if theres 1 trial
    xs = xs';
    ys = ys';
end

xs = xs+randi(jitter,length(xs),1);
ys = ys+randi(jitter,length(ys),1);

subplot(src(1),src(2),subnum)
hold on
plot(x',y','color',[0.8 0.8 0.8])
plot(xs,ys,'.r')
hold off
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
set(gca,'Xcolor',[1 1 1]);
set(gca,'Ycolor',[1 1 1]);

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