%% Visually Check trial by trial spike variability aligned to fixation to check if neuron was stable
for unit = 1:num_units
    allspikes = zeros(length(cfg.trl),7000);
    for t = 1:length(cfg.trl);
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == 15);
        fixation = cfg.trl(t).alltim(cfg.trl(t).allval == 8)-trial_start;
        spikes = find(data(unit).values{t});
        spikes = spikes-fixation;
        spikes(spikes < 1) = [];
        allspikes(t,spikes) = 1;
    end
    
    figure
    [trial,time] = find(allspikes(type == i,:));
    plot(time,trial,'.k')
    xlable('Time from Stimulus Onset (ms)')
    ylabel('Trial #')
    
    if multiunit(unit)
        subtitle(['Multiunit ' cfg.channel{unit}]);
    else
        subtitle(cfg.channel{unit});
    end
end
%%
multiunit = [1 0 0 1 0 1];

imageY = 600;
imageX = 800;
H = fspecial('gaussian',[25 25],2);
twin = 500; %how much time to take before or after event starts or ends respectively
t0 = twin ; %which index to define when an event occured

ITIstart_code = 15;
ITIend_code =16;
fixspot_on_code = 35;
fixation_occured_code = 8;
dot_on_code = 25;
dot_clrchng_code = 27;
bar_code_response = 4; %monkey made a move
reward_code = 3;

timed_firing = cell(7,num_units); %spikes related to particular events
spatial_firing = cell(1,num_units); %spikes during moving dot period by spatial location
direction= cell(1,num_units); %direction of the path the item took in the last 250 ms before firing
dot_location = NaN(length(cfg.trl)*2,3000); %location of the dot on the screen for the trial
for unit =1:num_units
    for event = 1:size(timed_firing,1)
        timed_firing{event,unit} = NaN(length(cfg.trl),4000);
    end
    spatial_firing{unit} = NaN(length(cfg.trl),3000);
end

total_spikes = zeros(1,6);
trial_dur = NaN(1,length(cfg.trl));

for t = 1:length(cfg.trl);
    if any(cfg.trl(t).allval == reward_code); %only take correct trials
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
        ITIperiod = [cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code) cfg.trl(t).alltim(cfg.trl(t).allval == ITIend_code)]-trial_start;
        fixspot_on = cfg.trl(t).alltim(cfg.trl(t).allval == fixspot_on_code)-trial_start;
        fixation = cfg.trl(t).alltim(cfg.trl(t).allval == fixation_occured_code)-trial_start;
        pathon =  cfg.trl(t).alltim(cfg.trl(t).allval == dot_on_code)-trial_start; %meta data starts at 1st 666 which appears before dot actually turns on in event array
        dot_clrchng = cfg.trl(t).alltim(cfg.trl(t).allval == dot_clrchng_code)-trial_start;
        responded = cfg.trl(t).alltim(cfg.trl(t).allval == bar_code_response)-trial_start;
        reward = cfg.trl(t).alltim(cfg.trl(t).allval == reward_code)-trial_start;
        reward = [reward(1) reward(end)];
        
        xn = interp1(meta(t).sample,meta(t).x,meta(t).sample(1):meta(t).sample(end),'cubic');
        yn = interp1(meta(t).sample,meta(t).y,meta(t).sample(1):meta(t).sample(end),'cubic');
        
        dot_location(2*t-1,1:length(xn)) = round(xn);
        dot_location(2*t,1:length(xn)) = round(yn);
        
        for unit = 1:num_units
            
            spikes = find(data(unit).values{t});
            total_spikes(unit) = total_spikes(unit)+length(spikes);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %---spatial information---%
            
            spikeind = spikes(spikes >= pathon & spikes <= responded)-pathon;
            spikeind(spikeind < 1) = []; %should only happen when spikes occur when path turns on
            
            temp = zeros(1,length(xn));
            temp(spikeind) = 1;
            spatial_firing{unit}(t,1:length(temp)) = temp;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %---temporal information---%
            trial_dur(t) = length(xn);
            
            %event 1: ITI period
            if t == 1;
                event_1_spikes = spikes(spikes <= ITIperiod(2) +twin)+t0;
                timevec1 = [1 ITIperiod+twin]+t0;
            else %look back in time to end of previous trial
                prespikes = find(data(unit).values{t-1});
                prespikes = prespikes(prespikes > length(data(unit).values{t-1})-twin) - length(data(unit).values{t-1});
                event_1_spikes = [prespikes spikes(spikes <= ITIperiod(2) +twin)]+t0;
                timevec1 = [1-twin ITIperiod(2)+twin]+t0;
            end
            tempvec = zeros(1,timevec1(2)-timevec1(1)+1);
            tempvec(event_1_spikes) = 1;
            timed_firing{1,unit}(t,timevec1(1):timevec1(2)) = tempvec;
            
            %event 2: fixation spot turns on
            event_2_spikes = spikes(spikes > fixspot_on-twin & spikes <= fixspot_on+twin)-fixspot_on+t0;
            timevec2 = [1-twin twin]+t0;
            tempvec = zeros(1,timevec2(2)-timevec2(1)+1);
            tempvec(event_2_spikes) = 1;
            timed_firing{2,unit}(t,timevec2(1):timevec2(2)) = tempvec;
            
            %event 3: fixatoin occurs
            event_3_spikes = spikes(spikes > fixation-twin & spikes <= fixation+twin)-fixation+t0;
            timevec3 = [1-twin twin]+t0;
            tempvec = zeros(1,timevec3(2)-timevec3(1)+1);
            tempvec(event_3_spikes) = 1;
            timed_firing{3,unit}(t,timevec3(1):timevec3(2)) = tempvec;
            
            %event 4: dot turns on
            event_4_spikes = spikes(spikes > pathon-twin & spikes <= pathon+twin)-pathon+t0;
            timevec4 = [1-twin twin]+t0;
            tempvec = zeros(1,timevec4(2)-timevec4(1)+1);
            tempvec(event_4_spikes) = 1;
            timed_firing{4,unit}(t,timevec4(1):timevec4(2)) = tempvec;
            
            %event 5: dot changes color.
            event_5_spikes = spikes(spikes > dot_clrchng-twin & spikes <= dot_clrchng+twin)-dot_clrchng+t0;
            timevec5 = [1-twin twin]+t0;
            tempvec = zeros(1,timevec5(2)-timevec5(1)+1);
            tempvec(event_5_spikes) = 1;
            timed_firing{5,unit}(t,timevec5(1):timevec5(2)) = tempvec;
            
            %event 6: responded to color change
            event_6_spikes = spikes(spikes > responded-twin & spikes <= responded+twin)-responded+t0;
            timevec6 = [1-twin twin]+t0;
            tempvec = zeros(1,timevec6(2)-timevec6(1)+1);
            tempvec(event_6_spikes) = 1;
            timed_firing{6,unit}(t,timevec6(1):timevec6(2)) = tempvec;
            
            %event 7: reward period
            if t == length(cfg.trl);
                if length(data(unit).values{t}) > reward(2)+twin;
                    rewardstop = reward(2)+twin;
                else
                    rewardstop = length(data(unit).values{t});
                end
                event_7_spikes = spikes(spikes > reward(1)-twin & spikes <= rewardstop)-reward(1)+t0;
                timevec7 = [1-twin rewardstop-reward(1)]+t0;
            else %look to the next trial
                if length(data(unit).values{t}) > reward(2)+twin;
                    rewardstop = reward(2)+twin;
                    event_7_spikes = spikes(spikes > reward(1)-twin & spikes <= rewardstop)-reward(1)+t0;
                    timevec7 = [1-twin rewardstop-reward(1)]+t0;
                else
                    rewardstop = twin-(length(data(unit).values{t})-reward(2))+1;
                    postspikes = find(data(unit).values{t+1});
                    postspikes = postspikes(postspikes <= rewardstop);
                    event_7_spikes = [spikes(spikes > reward(1)-twin)-reward(1)...
                        postspikes+(length(data(unit).values{t})-reward(2))+(reward(2)-reward(1))]+t0;
                    timevec7 = [1-twin reward(2)-reward(1)+twin]+t0;
                end
            end
            tempvec = zeros(1,timevec7(2)-timevec7(1)+1);
            tempvec(event_7_spikes) = 1;
            timed_firing{7,unit}(t,timevec7(1):timevec7(2)) = tempvec;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finally trim down the execess data to save processing time
for unit =1:num_units
    for event = 1:size(timed_firing,1)
        if event == 1;
            timed_firing{event,unit} = timed_firing{event,unit}(:,1:2000);
        elseif event == 7;
            timed_firing{event,unit} = timed_firing{event,unit}(:,1:1800);
        end
        nancols = find(nansum(isnan(timed_firing{event,unit})) == size(timed_firing{event,unit},1));
        nancols(nancols < twin) = [];
        timed_firing{event,unit}(:,nancols) = [];
        nanrows = find(nansum(isnan(timed_firing{event,unit}),2) == size(timed_firing{event,unit},2));
        timed_firing{event,unit}(nanrows,:) = [];
    end
end

%% Calculate Temporal Information
smval = 60;
fs = 1000;
time_observed_info = NaN(size(timed_firing,1),num_units);
time_shuffled_info = cell(size(timed_firing,1),num_units);
numshuffs = 100;
for unit = 1:num_units
    for event = 1:size(timed_firing,1);
        total_time = sum(~isnan(timed_firing{event,unit}));
        prob_time = total_time/sum(total_time);
        lambda = nansum(nansum(timed_firing{event,unit}))/sum(sum(~isnan(timed_firing{event,unit})))*1000;
        [lambda_x,~]= nandens(timed_firing{event,unit},smval,'gauss',fs,'nanflt');
        plogp = lambda_x.*log2(lambda_x/lambda);
        plogp(isnan(plogp)) = 0;
        time_observed_info(event,unit) = sum(lambda_x.*plogp.*prob_time);
        
        time_shuffled_info{event,unit} = NaN(1,numshuffs);
        for shuffled = 1:numshuffs;
            shuffed_firing = circshift_row(timed_firing{event,unit});
            [lambda_x,~]= nandens(shuffed_firing,smval,'gauss',fs,'nanflt');
            plogp = lambda_x.*log2(lambda_x/lambda);
            plogp(isnan(plogp)) = 0;
            time_shuffled_info{event,unit}(shuffled) = sum(lambda_x.*plogp.*prob_time);
        end
    end
end

time_shuffled_95_percentile = NaN(size(timed_firing,1),num_units);

for unit = 1:num_units
    for event = 1:size(timed_firing,1);
        time_shuffled_95_percentile(event,unit) = prctile(time_shuffled_info{event,unit},95);
    end
end
%% 
smval = 60;
for unit = 1:num_units;
    yls = NaN(2,7);
    
    figure
    
    s = subplot(2,4,1);
    t = 1:size(timed_firing{1,unit},2);
    dofill(t,timed_firing{1,unit},'red',1,smval);
    if shuffled_95_percentile(1,unit) < time_observed_info(1,unit)
        title(['Spikes = ' num2str(nansum(nansum(timed_firing{1,unit})))  'Bits: ' num2str(time_observed_info(1,unit))])
    else
        title(['Spikes = ' num2str(nansum(nansum(timed_firing{1,unit})))])
    end
    xlabel('Time from ITI start (ms)')
    ylabel('Firing rate (Hz)')
    xlim([0 2000])
    set(s,'Xtick',[0:500:2000])
    set(s,'XtickLabel',num2cell([0:500:2000]-500));
    yls(:,1) = ylim;
    
    s = subplot(2,4,2);
    t = 1:size(timed_firing{2,unit},2);
    dofill(t,timed_firing{2,unit},'red',1,smval);
    if shuffled_95_percentile(2,unit) < time_observed_info(2,unit)
        title(['Spikes = ' num2str(nansum(nansum(timed_firing{2,unit})))  'Bits: ' num2str(time_observed_info(2,unit))])
    else
        title(['Spikes = ' num2str(nansum(nansum(timed_firing{2,unit})))])
    end
    xlabel('Time from Fixation Spot On (ms)')
    ylabel('Firing rate (Hz)')
    xlim([0 1000])
    set(s,'Xtick',[0:250:1000])
    set(s,'XtickLabel',num2cell([0:250:1000]-500));
    yls(:,2) = ylim;
    
    s = subplot(2,4,3);
    t = 1:size(timed_firing{3,unit},2);
    dofill(t,timed_firing{3,unit},'red',1,smval);
    if shuffled_95_percentile(3,unit) < time_observed_info(3,unit)
        title(['Spikes = ' num2str(nansum(nansum(timed_firing{3,unit})))  'Bits: ' num2str(time_observed_info(3,unit))])
    else
        title(['Spikes = ' num2str(nansum(nansum(timed_firing{3,unit})))])
    end
    xlabel('Time from Fixation on Cross Hair (ms)')
    ylabel('Firing rate (Hz)')
    xlim([0 1000])
    set(s,'Xtick',[0:250:1000])
    set(s,'XtickLabel',num2cell([0:250:1000]-500));
    yls(:,4) = ylim;
    
    s = subplot(2,4,4);
    t = 1:size(timed_firing{4,unit},2);
    dofill(t,timed_firing{4,unit},'red',1,smval);
    if shuffled_95_percentile(4,unit) < time_observed_info(4,unit)
        title(['Spikes = ' num2str(nansum(nansum(timed_firing{4,unit})))  'Bits: ' num2str(time_observed_info(4,unit))])
    else
        title(['Spikes = ' num2str(nansum(nansum(timed_firing{4,unit})))])
    end
    xlabel('Time from Dot Appearing (ms)')
    ylabel('Firing rate (Hz)')
    xlim([0 1000])
    set(s,'Xtick',[0:250:1000])
    set(s,'XtickLabel',num2cell([0:250:1000]-500));
    yls(:,4) = ylim;
    
    s = subplot(2,4,5);
    t = 1:size(timed_firing{5,unit},2);
    dofill(t,timed_firing{5,unit},'red',1,smval);
    if shuffled_95_percentile(5,unit) < time_observed_info(5,unit)
        title(['Spikes = ' num2str(nansum(nansum(timed_firing{5,unit})))  'Bits: ' num2str(time_observed_info(5,unit))])
    else
        title(['Spikes = ' num2str(nansum(nansum(timed_firing{5,unit})))])
    end
    xlabel('Time from Dot Color Change (ms)')
    ylabel('Firing rate (Hz)')
    xlim([0 1000])
    set(s,'Xtick',[0:250:1000])
    set(s,'XtickLabel',num2cell([0:250:1000]-500));
    yls(:,5) = ylim;
    
    s = subplot(2,4,6);
    t = 1:size(timed_firing{6,unit},2);
    dofill(t,timed_firing{6,unit},'red',1,smval);
    if shuffled_95_percentile(6,unit) < time_observed_info(6,unit)
        title(['Spikes = ' num2str(nansum(nansum(timed_firing{6,unit})))  'Bits: ' num2str(time_observed_info(6,unit))])
    else
        title(['Spikes = ' num2str(nansum(nansum(timed_firing{6,unit})))])
    end
    xlabel('Time from Response to Color Change (ms)')
    ylabel('Firing rate (Hz)')
    xlim([0 1000])
    set(s,'Xtick',[0:250:1000])
    set(s,'XtickLabel',num2cell([0:250:1000]-500));
    yls(:,16) = ylim;
    
    s = subplot(2,4,7);
    t = 1:size(timed_firing{7,unit},2);
    dofill(t,timed_firing{7,unit},'red',1,smval);
    if shuffled_95_percentile(7,unit) < time_observed_info(7,unit)
        title(['Spikes = ' num2str(nansum(nansum(timed_firing{7,unit})))  'Bits: ' num2str(time_observed_info(7,unit))])
    else
        title(['Spikes = ' num2str(nansum(nansum(timed_firing{7,unit})))])
    end
    xlabel('Time from Reward onset(ms)')
    ylabel('Firing rate (Hz)')
    xlim([0 2000])
    set(s,'Xtick',[0:500:2000])
    set(s,'XtickLabel',num2cell([0:500:2000]-500));
    yls(:,7) = ylim;
    
    %by trial
    means = NaN(1,7);
    stds = NaN(1,7);
    numpoints = NaN(1,7);
    for event = 1:7;
        trial_spikes = nansum(timed_firing{event,unit}(:,twin:end),2); %ignore the 1st 500 ms this is from the previous event
        dur = (size(timed_firing{event,unit},2)-500)/1000; %duration in seconds
        means(event) = nanmean(trial_spikes)/dur;
        stds(event) = nanstd(trial_spikes)/dur;
        numpoints(event) = sum(~isnan(trial_spikes));
    end

    subplot(2,4,8)
    hold on
    bar(1:7,means)
    errorb(1:7,means,stds./sqrt(numpoints));
    hold off
    ylabel('Firing Rate')
    xlabel('Epoch')
    title('Average Firing Rate by Epoch')
    
    
    max_y = max(yls(2,:));
    for i = 1:7;
        subplot(2,4,i);
        ylim([0 max_y])
    end
    
    if multiunit(unit)
        subtitle(['Multiunit ' cfg.channel{unit} ', total spikes: ' num2str(total_spikes(unit))]);
    else
        subtitle([cfg.channel{unit} ', total spikes: ' num2str(total_spikes(unit))]);
    end
    
end
%% Calculate Spatial information
filter_width = 4;
binsize = 4;

H = fspecial('gaussian',[filter_width*10+5 filter_width*10+5],filter_width);

remove_x = floor(135/binsize);
remove_y = floor(35/binsize);

%calculate the total time spent at any locaitons in binned pixels
spatial_time = zeros(600,800);
dot_x = dot_location(1:2:end,:);
dot_x = dot_x(1:end);
dot_x(isnan(dot_x)) = [];
dot_y = dot_location(2:2:end,:);
dot_y = dot_y(1:end);
dot_y(isnan(dot_y)) = [];
xy_ind = sub2ind(size(spatial_time),round(dot_y),round(dot_x));
for i = 1:length(xy_ind);
    spatial_time(xy_ind(i)) =  spatial_time(xy_ind(i)) +1;
end
spatial_time = spatial_time/fs; %convert from ms to sec

tt = bin2(spatial_time,binsize,binsize);
tt = imfilter(tt,H);
tt(tt == 0) = NaN;

%because boundaries don't got that far out make sure filtering isn't reporting fake data
tt(:,1:remove_x) = NaN;
tt(:,end-remove_x+1:end) = NaN;
tt(1:remove_y,:) = NaN;
tt(end-remove_y+1:end,:) = NaN;
tt = tt(end:-1:1,:); %flip to reflex flip in firing rate to make rightside up images


firing_location = cell(1,num_units);
spatial_observed_info = NaN(1,num_units);
spatial_shuffled_info = cell(1,num_units);
spatial_shuffled_95_percentile = NaN(1,num_units);

direction = cell(1,num_units);

directional_info = NaN(1,num_units);
shuffled_directional_info = cell(1,num_units);
shuffled_directional_95_percentile = NaN(1,num_units);

rayleigh_vector = NaN(1,num_units);
shuffled_rayleigh_vector = cell(1,num_units);
shuffled_rayleigh_vector_95_percentile = NaN(1,num_units);

numshuffs = 100;
for unit = 1:num_units
    
    %%%---to compute spatial firing rate---%%%
    firing_location{unit} = zeros(600,800);
    [trial,time] = find(spatial_firing{unit} == 1);
    for i = 1:length(trial);
        xs = round(dot_location(2*trial(i)-1,time(i)));
        ys = round(dot_location(2*trial(i),time(i)));
        firing_location{unit}(ys,xs) =  firing_location{unit}(ys,xs) +1;
    end
    
    sf = bin2(firing_location{unit},binsize,binsize);
    sf = imfilter(sf,H);
    sf(sf == 0) = NaN;
    %because boundaries don't got that far out make sure filtering isn't reporting fake data
    sf(:,1:remove_x) = NaN;
    sf(:,end-remove_x+1:end) = NaN;
    sf(1:remove_y,:) = NaN;
    sf(end-remove_y+1:end,:) = NaN;
    sf = sf(end:-1:1,:);
    
    lambda_x = sf./tt; %observed firing rate over space
    prob_time = tt/nansum(nansum(tt));
    lambda = nansum(nansum(lambda_x))/sum(sum(~isnan(lambda_x)));
    plogp = lambda_x.*log2(lambda_x/lambda);
    plogp(isnan(plogp)) = 0;
    spatial_observed_info(unit) = nansum(nansum(lambda_x.*plogp.*prob_time));
    
    %%%---to compute directional selectivity---%%%
    % going to ignore spikes that shouldn't be related to the
    % direcitonal movement of the dot i.e. within 150 ms of the stimulus
    % turning on
    trial(time <= 150) = [];
    time(time <= 150) = [];
    direction{unit} = NaN(1,length(time));
    for t = 1:length(trial);
        dy = dot_location(2*trial(t),time(t))-dot_location(2*trial(t),time(t)-150);
        dx = dot_location(2*trial(t)-1,time(t))-dot_location(2*trial(t)-1,time(t)-150);
        direction{unit}(t) = atan2(dy,dx);
    end
    
    [r_theta_j,theta] = rose(direction{unit});
    r_theta_j = r_theta_j(1:2:end); %for plotting purposes every other is a 0
    theta = theta(1:2:end); %for plotting purposes every other is a 0
    n = length(theta);
    rayleigh_vector(unit) = abs(pi/(n*sin(pi/n))*sum(r_theta_j.*exp(-j*theta))/sum(r_theta_j));
    
    plogp = r_theta_j.*log2(r_theta_j/sum(r_theta_j));
    plogp(isnan(plogp)) = 0;
    directional_info(unit) = sum(r_theta_j.*plogp*1/n);%i think this is right .... idk could be scaled entropy otherwise
    
    spatial_shuffled_info{unit} = NaN(1,numshuffs);
    for shuffled = 1:numshuffs;
        shuffed_firing = circshift_row(spatial_firing{unit}); %rotate spike times within a trial randomly
        
        %%%---compute shuflfed spatial information---%%%
        shuffled_firing_location = zeros(600,800);
        [trial,time] = find(shuffed_firing == 1);
        for i = 1:length(trial);
            xs = round(dot_location(2*trial(i)-1,time(i)));
            ys = round(dot_location(2*trial(i),time(i)));
            shuffled_firing_location(ys,xs) = shuffled_firing_location(ys,xs) +1;
        end
        
        sf = bin2(shuffled_firing_location,binsize,binsize);
        sf = imfilter(sf,H);
        sf(sf == 0) = NaN;
        %because boundaries don't got that far out make sure filtering isn't reporting fake data
        sf(:,1:remove_x) = NaN;
        sf(:,end-remove_x+1:end) = NaN;
        sf(1:remove_y,:) = NaN;
        sf(end-remove_y+1:end,:) = NaN;
        sf = sf(end:-1:1,:);
        
        lambda_x = sf./tt;
        plogp = lambda_x.*log2(lambda_x/lambda);
        plogp(isnan(plogp)) = 0;
        spatial_shuffled_info{unit}(shuffled) = nansum(nansum(lambda_x.*plogp.*prob_time));
        
        %%%---compute shuffled directional information---%%%
        trial(time <= 150) = [];
        time(time <= 150) = [];
        shuffled_directions = NaN(1,length(time));
        for t = 1:length(trial);
            dy = dot_location(2*trial(t),time(t))-dot_location(2*trial(t),time(t)-150);
            dx = dot_location(2*trial(t)-1,time(t))-dot_location(2*trial(t)-1,time(t)-150);
            shuffled_directions(t) = atan2(dy,dx);
        end
        
        [r_theta_j,theta] = rose(shuffled_directions);
        r_theta_j = r_theta_j(1:2:end); %for plotting purposes every other is a 0
        theta = theta(1:2:end); %for plotting purposes every other is a 0
        n = length(theta);
        shuffled_rayleigh_vector{unit}(shuffled) = abs(pi/(n*sin(pi/n))*sum(r_theta_j.*exp(-j*theta))/sum(r_theta_j));
        
        plogp = r_theta_j.*log2(r_theta_j/sum(r_theta_j));
        plogp(isnan(plogp)) = 0;
        shuffled_directional_info{unit}(shuffled) = sum(r_theta_j.*plogp*1/n);%i think this is right .... idk could be scaled entropy otherwise
    end
end

for unit = 1:num_units
    spatial_shuffled_95_percentile(unit) = prctile(spatial_shuffled_info{unit},95);
    shuffled_directional_95_percentile(unit) = prctile(shuffled_directional_info{unit},95);
    shuffled_rayleigh_vector_95_percentile(unit) = prctile(shuffled_rayleigh_vector{unit},95);
end

for unit = 1:num_units
    figure
    subplot(1,2,1)
    sf = bin2(firing_location{unit},binsize,binsize);
    sf = imfilter(sf,H);
    sf(sf == 0) = NaN;
    
    %because boundaries don't got that far out make sure filtering isn't reporting fake data
    sf(:,1:remove_x) = NaN;
    sf(:,end-remove_x+1:end) = NaN;
    sf(1:remove_y,:) = NaN;
    sf(end-remove_y+1:end,:) = NaN;
    
    sf = sf(end:-1:1,:);
    r = sf./tt;
    pcolor(r),shading interp
    axis square
    axis off
    if spatial_shuffled_95_percentile(unit) < spatial_observed_info(unit)
         title(['Max Firing Rate ' num2str(max(max(r))) ' Bits: ' num2str(spatial_observed_info(unit))])
    else
        title(['Max Firing Rate ' num2str(max(max(r)))])
    end

    
    subplot(1,2,2)
    rose(direction{unit})
    titlevec= [];
    if shuffled_directional_95_percentile(unit) < directional_info(unit)
        titlevec = [titlevec ' Bits: ' num2str(directional_info(unit))];
    end
    if shuffled_rayleigh_vector_95_percentile(unit) < rayleigh_vector(unit)
        titlevec = [titlevec ' RV: ' num2str(rayleigh_vector(unit))];
    end
    title(titlevec)

    if multiunit(unit)
        subtitle(['Multiunit ' cfg.channel{unit}])
    else
        subtitle(cfg.channel{unit})
    end
end