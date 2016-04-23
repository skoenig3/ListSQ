% %% Visually Check trial by trial spike variability aligned to trial event to check if neuron was stable
% for unit = 1:num_units
%     type = NaN(1,length(cfg.trl));
%     allspikes = zeros(length(cfg.trl),7000);
%     for t = 1:length(cfg.trl);
%         trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == 15);
%         event = cfg.trl(t).alltim(cfg.trl(t).allval == 23)-trial_start;
%         
%         if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(1) % for the sequence trials
%             if  any(cfg.trl(t).allval == reward_code)
%                 type(t) = 1;
%             else
%                 continue;
%             end
%         elseif  itmlist(cfg.trl(t).cnd-1000) <= sequence_items(2)
%              if  any(cfg.trl(t).allval == reward_code)
%                  type(t) = 2;
%              else
%                  continue;
%              end
%         else
%             if  any(cfg.trl(t).allval == 23) %if image was shown
%                 type(t) = 3;
%             else
%                 continue;
%             end
%         end
%        
%         spikes = find(data(unit).values{t});
%         spikes = spikes-event;
%         spikes(spikes < 1) = [];
%         spikes(spikes > 7000) = [];
%         allspikes(t,spikes) = 1;
%     end
%     
%     allspikes(isnan(type),:) = [];
%     type(isnan(type)) = [];
%     
%     figure
%     for i = 1:3;
%         subplot(1,3,i)
%         [trial,time] = find(allspikes(type == i,:));
%         plot(time,trial,'.k')
%     end
%     xlabel('Time from Stimulus Onset (ms)')
%     ylabel('Trial #')
%     
%     if multiunit(unit)
%         subtitle(['Multiunit ' cfg.channel{unit}]);
%     else
%         subtitle(cfg.channel{unit});
%     end
% end
%%
multiunit = [zeros(1,9) 1];
ITMFile = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Item and Conditions Files\ListSQ08.itm';
CNDFile = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Item and Conditions Files\ListSQ08.cnd';

itmfil=[];
h =fopen(ITMFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(itmfil)
        if length(tline)>size(itmfil,2)
            tline=tline(1:size(itmfil,2));
        end
    end
    tline = [tline ones(1,(size(itmfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        itmfil=[itmfil; tline];
    else
        break
    end
end
fclose(h);

cndfil=[];
h=fopen(CNDFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(cndfil)
        if length(tline)>size(cndfil,2)
            tline=tline(1:size(cndfil,2));
        end
    end
    tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        cndfil=[cndfil; tline];
    else
        break
    end
end
fclose(h);

itmlist = zeros(size(cndfil,1)-1,1);
for i = 2:size(cndfil,1);
    str = textscan(cndfil(i,:),'%d');
    itmlist(i-1) = str{1}(end);
end

sequence_items = [11 19]; %what item starts the sequence

twin = 500; %ms before or after event is over in which to keep the data


timed_firing = cell(19,num_units); % unit by epoch
% epoch 1 ITI period
% epoch 2 item 1 in sequence is displayed
% epoch 3 item 1 in sequence is fixated
% epoch 4 item 1 turned off/IIT (inter item interval)
% epoch 5 item 2 in sequence is displayed
% epoch 6 item 2 in sequence is fixated
% epoch 7 item 2 turned off/IIT (inter item interval)
% epoch 8 item 3 in sequence is displayed
% epoch 9 item 3 in sequence is fixated
% epoch 10 item 3 turned off/IIT (inter item interval)
% epoch 11 item 4 in sequence is displayed
% epoch 12 item 4 in sequence is fixated
% epoch 13 reward period during sequence trials
% epoch 14 fixation cross is displayed
% epoch 15 fixation on crosshair to initialize image trials
% epoch 16 image viewing period start
% epoch 17 image viewing period end
% epoch 18 broke fixation while viewing image
% epoch 19 started viewing image again after breaking fixation

for i = 1:size(timed_firing,1);
    for ii = 1:size(timed_firing,2)
        if i == 1 || i == 13; %ITI and reward period
            timed_firing{i,ii} = NaN(size(cfg.trl,2),1000+2*twin);
        elseif i == 16 || i == 17 %image viewing period
            timed_firing{i,ii} = NaN(size(cfg.trl,2),7000+2*twin);
        else
            timed_firing{i,ii} = NaN(size(cfg.trl,2),2*twin);
        end
    end
end

eye_pos = NaN(2*length(cfg.trl),10000);
spatial_firing = cell(1,num_units);
for unit = 1:num_units
    spatial_firing{unit} = NaN(length(cfg.trl),10000);
end

ITIstart_code = 15;
ITIend_code = 16;
crosses_on_off = [23 25 27 29;
    24 26 28 30];
fixation_occured_code = 8;
fixspot_on_code = 35;
img_on_code= 23;
img_off_code = 24;
reward_code = 3;
break_fixation = 209; %on image trials
trial_type = NaN(1,length(cfg.trl)); %sequences = 1, images = 2
which_sequence = NaN(1,length(cfg.trl));
total_spikes = zeros(1,num_units);
broke_row = zeros(1,num_units);
for t = 1:length(cfg.trl);
    if any(cfg.trl(t).allval == img_on_code); %in which image was displayed or 1st item in sequence was displayed
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
        ITIperiod = [cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code) cfg.trl(t).alltim(cfg.trl(t).allval == ITIend_code)]-trial_start;
        ITIperiod(ITIperiod > 1000) = 1000; %some inconsitency in duration so keep at 1000
        reward = cfg.trl(t).alltim(cfg.trl(t).allval == reward_code)-trial_start;
        if ~isempty(reward) %therefore successful sequence trial
            reward = [reward(1) reward(end)];
            reward(reward > 800) = 800;
            crosses_on = NaN(1,length(crosses_on_off));
            crosses_off = NaN(1,length(crosses_on_off));
            fix_cross = NaN(1,length(crosses_on_off));
            for c = 1:length(crosses_on_off)
                crosses_on_ind = cfg.trl(t).allval == crosses_on_off(1,c);
                crosses_on(c) = cfg.trl(t).alltim(crosses_on_ind)-trial_start;
                crosses_off(c) = cfg.trl(t).alltim(cfg.trl(t).allval == crosses_on_off(2,c))-trial_start;
                possfix = find(cfg.trl(t).allval == fixation_occured_code);
                if ~isempty(possfix)
                    possfix(possfix < find(crosses_on_ind)) = [];
                    if ~isempty(possfix)
                        fix_cross(c) = cfg.trl(t).alltim(possfix(1))-trial_start;
                    end
                end
            end
            
            if itmlist(cfg.trl(t).cnd-1000) == sequence_items(1)
                which_sequence(t) = 1;
            elseif  itmlist(cfg.trl(t).cnd-1000) == sequence_items(2)
                which_sequence(t) = 2;
            end
            trial_type(t)=1;
        elseif itmlist(cfg.trl(t).cnd-1000) > sequence_items(2) %image trial
            fixspot_on = cfg.trl(t).alltim(cfg.trl(t).allval == fixspot_on_code)-trial_start;
            imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start;
            imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start;
            imgfixation = cfg.trl(t).alltim(cfg.trl(t).allval == fixation_occured_code)-trial_start;
            imgfixation = imgfixation(1); %this is fixation on 1st item
            broke = cfg.trl(t).alltim(cfg.trl(t).allval == break_fixation)-trial_start;
            broke(broke < imgon) = [];% just in case
            refixated = cfg.trl(t).alltim(cfg.trl(t).allval == fixation_occured_code)-trial_start;
            refixated(refixated < imgon+50) = [];% remove 1st fixation on crosshair sometimes fixation code gets displayed for image trial too?
            duration_outside = refixated-broke;
            %only take trials in which the monkey looked outside for long
            %durations setting to 500 since this is clearly more than 1 fixation
            refixated(duration_outside < 500) = [];
            broke(duration_outside < 500) = [];
            %remove break fixations that occurred when the monkey wasn't
            %paying attention and the image was up for a really long time >
            %10 secs
            refixated(broke > 10000) = [];
            broke(broke > 10000) = [];
            trial_type(t)=2;
        else
            continue
        end
        
        if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(2) %so sequence trials
            tstart = crosses_on(1,1);
            tend = reward(1);
        else
            tstart = imgon;
            if imgoff-imgon <= 10000
                tend = imgoff;
            else
                tend = imgon+9999;
            end
        end
        xn = data(end-1).values{t}(tstart:tend);
        yn = data(end).values{t}(tstart:tend);
        
        eye_pos(2*t-1,1:length(xn)) = round(xn);
        eye_pos(2*t,1:length(xn)) = round(yn);
        
        for unit = 1:num_units
            
            spikes = find(data(unit).values{t});
            total_spikes(unit) = total_spikes(unit)+length(spikes);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %---spatial information---%
            
            spikeind = spikes(spikes >= tstart & spikes <= tend)-tstart;
            spikeind(spikeind < 1) = []; %should only happen when spikes occur at the same time as tstart
            
            temp = zeros(1,length(xn));
            temp(spikeind) = 1;
            spatial_firing{unit}(t,1:length(temp)) = temp;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %---temporal information---%
            
            %event 1: ITI period
            if t == 1;
                event_spikes = spikes(spikes <= ITIperiod(2)+twin)+twin;
                timevec = [1 ITIperiod(2)+twin]+twin;
                tempvec = [NaN(1,twin) zeros(1,timevec(2)-timevec(1)+1)];
                tempvec(event_spikes) = 1;
                tempvec(1:twin) = [];
            else %look back in time to end of previous trial
                prespikes = find(data(unit).values{t-1});
                prespikes = prespikes(prespikes > length(data(unit).values{t-1})-twin) - length(data(unit).values{t-1});
                event_spikes = [prespikes spikes(spikes <= ITIperiod(2)+twin)]+twin;
                timevec = [1-twin ITIperiod(2)+twin]+twin;
                tempvec =  zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
            end
            timed_firing{1,unit}(t,timevec(1):timevec(2)) = tempvec;
            
            if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(2) %so sequence trials
                %event 2-9: evens item in sequence displayed, item in
                %sequence fixated, and item in sequence turned off
                for c = 1:length(crosses_on)
                    event_spikes = spikes(spikes > crosses_on(c)-twin & spikes <= crosses_on(c)+twin)-crosses_on(c)+twin;
                    timevec = [1-twin twin]+twin;
                    tempvec = zeros(1,timevec(2)-timevec(1)+1);
                    tempvec(event_spikes) = 1;
                    timed_firing{3*c-1,unit}(t,timevec(1):timevec(2)) = tempvec;
                    
                    event_spikes = spikes(spikes > fix_cross(c)-twin & spikes <= fix_cross(c)+twin)-fix_cross(c)+twin;
                    timevec = [1-twin twin]+twin;
                    tempvec = zeros(1,timevec(2)-timevec(1)+1);
                    tempvec(event_spikes) = 1;
                    timed_firing{3*c,unit}(t,timevec(1):timevec(2)) = tempvec;
                    
                    if c ~= 4
                        event_spikes = spikes(spikes > crosses_off(c)-twin & spikes <= crosses_off(c)+twin)-crosses_off(c)+twin;
                        timevec = [1-twin twin]+twin;
                        tempvec = zeros(1,timevec(2)-timevec(1)+1);
                        tempvec(event_spikes) = 1;
                        timed_firing{3*c+1,unit}(t,timevec(1):timevec(2)) = tempvec;
                    end
                        
                end
                
                %event 13: reward period for sequence trials
                if ~isempty(reward)
                    if t == length(cfg.trl);
                        if length(data(unit).values{t}) > reward(2)+twin;
                            rewardstop = reward(2)+twin;
                        else
                            rewardstop = length(data(unit).values{t});
                        end
                        event_spikes = spikes(spikes > reward(1)-twin & spikes <= rewardstop)-reward(1)+twin;
                        timevec = [1-twin rewardstop-reward(1)]+twin;
                    else %look to the next trial
                        if length(data(unit).values{t}) > reward(2)+twin;
                            rewardstop = reward(2)+twin;
                            event_spikes = spikes(spikes > reward(1)-twin & spikes <= rewardstop)-reward(1)+twin;
                            timevec = [1-twin rewardstop-reward(1)]+twin;
                        else
                            rewardstop = twin-(length(data(unit).values{t})-reward(2))+1;
                            postspikes = find(data(unit).values{t+1});
                            postspikes = postspikes(postspikes <= rewardstop);
                            event_spikes = [spikes(spikes > reward(1)-twin)-reward(1)...
                                postspikes+(length(data(unit).values{t})-reward(2))+(reward(2)-reward(1))]+twin;
                            timevec = [1-twin reward(2)-reward(1)+twin]+twin;
                        end
                    end
                    tempvec = zeros(1,timevec(2)-timevec(1)+1);
                    tempvec(event_spikes) = 1;
                    timed_firing{13,unit}(t,timevec(1):timevec(2)) = tempvec;
                end
                
            else %img trial
                %event 14: fixation spot turns on
                event_spikes = spikes(spikes > fixspot_on-twin & spikes <=fixspot_on+twin)-fixspot_on+twin;
                timevec = [1-twin twin]+twin;
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                timed_firing{14,unit}(t,timevec(1):timevec(2)) = tempvec;
                
                %event 15: fixation on fixspot occurs
                event_spikes = spikes(spikes > imgfixation-twin & spikes <= imgfixation+twin)-imgfixation+twin;
                timevec = [1-twin twin]+twin;
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                timed_firing{15,unit}(t,timevec(1):timevec(2)) = tempvec;
                
                %event 16: locked to start of image viewing period
                event_spikes = spikes(spikes > tstart-twin & spikes <= tend+twin)-tstart+twin;
                timevec = [1 8000];
                %limit spikes in case of image being displayed for longer than 7000 secs
                event_spikes(event_spikes > size(timed_firing{16,unit},2)) = [];
                timevec(timevec > size(timed_firing{16,unit},2)) = size(timed_firing{16,unit},2);
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                timed_firing{16,unit}(t,timevec(1):timevec(2)) = tempvec;
                
                %event 17: locked to end of image viewing period since
                %monkeys can look away for variable amounts of time
                event_spikes = spikes(spikes > imgoff-7000-twin & spikes <= imgoff+twin)-imgoff+8000-twin;
                event_spikes(event_spikes < 1) = [];   %limit spikes in case of image being displayed for longer than 7000 secs
                timevec = [1 8000];
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                timed_firing{17,unit}(t,timevec(1):timevec(2)) = tempvec;
                
                %event 18-19: break fixation and refixation on image period
                if ~isempty(broke);
                    for b = 1:length(broke)
                        broke_row(unit) = broke_row(unit)+1;
                        event_spikes = spikes(spikes > broke(b)-twin & spikes <= broke(b)+twin)-broke(b)+twin;
                        timevec = [1-twin twin]+twin;
                        tempvec = zeros(1,timevec(2)-timevec(1)+1);
                        tempvec(event_spikes) = 1;
                        timed_firing{18,unit}(broke_row(unit),timevec(1):timevec(2)) = tempvec;
                        
                        event_spikes = spikes(spikes > refixated(b)-twin & spikes <= refixated(b)+twin)- refixated(b)+twin;
                        timevec = [1-twin twin]+twin;
                        tempvec = zeros(1,timevec(2)-timevec(1)+1);
                        tempvec(event_spikes) = 1;
                        timed_firing{19,unit}(broke_row(unit),timevec(1):timevec(2)) = tempvec;
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finally trim down the execess data to save processing time
for unit =1:num_units
    for event = 1:size(timed_firing,1)
        nancols = find(nansum(isnan(timed_firing{event,unit})) == size(timed_firing{event,unit},1));
        nancols(nancols < twin) = [];
        timed_firing{event,unit}(:,nancols) = [];
        nanrows = find(nansum(isnan(timed_firing{event,unit}),2) == size(timed_firing{event,unit},2));
        timed_firing{event,unit}(nanrows,:) = [];
    end
    nancols = find(nansum(isnan(spatial_firing{unit})) == size(spatial_firing{unit},1));
    nancols(nancols < twin) = [];
    spatial_firing{unit}(:,nancols) = [];
    nanrows = find(nansum(isnan(spatial_firing{unit}),2) == size(spatial_firing{unit},2));
    spatial_firing{unit}(nanrows,:) = [];
end
eye_pos(sum(isnan(eye_pos),2) == size(eye_pos,2),:) = [];
which_sequence(isnan(which_sequence)) = [];
trial_type(isnan(trial_type)) = [];

which_img = NaN(1,96*2);
img_count = 1;
for t = 1:length(cfg.trl);
    if itmlist(cfg.trl(t).cnd-1000) > sequence_items(2)
        if any(cfg.trl(t).allval == img_on_code);
            which_img(img_count) = itmlist(cfg.trl(t).cnd-1000);
            img_count = img_count+1;
        end
    end
end

novel_vs_repeat = NaN(1,96*2);
for img = 1:max(which_img)
    imgind = find(which_img == img);
    dimgind = find(diff(imgind) > 1);
    if ~isempty(imgind);
        novel_vs_repeat(imgind(1:dimgind)) = 1;
        novel_vs_repeat(imgind(dimgind+1:end)) = 2;
    end
end
%% Plot Firing rates for image trials
smval_novrep = 250;
smval = 60;
for unit = 1:num_units
    figure
    
    yls = NaN(2,6);
    
    s(1) = subplot(2,4,1);
    t = 1:size(timed_firing{11,unit},2);
    dofill(t,timed_firing{14,unit},'black',1,smval);
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Onset of Fixation Cross(ms)')
    yls(:,1) = ylim;
    
    s(2) = subplot(2,4,5);
    t = 1:size(timed_firing{12,unit},2);
    dofill(t,timed_firing{15,unit},'black',1,smval);
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Fixation on Cross Onset(ms)')
    yls(:,2) = ylim;
    
    s(3) = subplot(2,4,[2 3]);
    novel_timed_firing = timed_firing{16,unit}(novel_vs_repeat == 1,:);
    repeat_timed_firing = timed_firing{16,unit}(novel_vs_repeat == 2,:);
    t = 1:size(novel_timed_firing,2);
    hold on
    dofill(t,novel_timed_firing,'blue',1,smval_novrep);
    dofill(t,repeat_timed_firing,'red',1,smval_novrep);
    set(gca,'Xtick',[0:500:8000])
    set(gca,'XtickLabel',num2cell([0:500:8000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Image Onset(ms)')
    legend('Novel','Repeat')
    xlim([0 4000]) %plot the 1st half of the image
    yls(:,3) = ylim;
    
    s(4) = subplot(2,4,[6 7]);
    novel_timed_firing = timed_firing{17,unit}(novel_vs_repeat == 1,:);
    repeat_timed_firing = timed_firing{17,unit}(novel_vs_repeat == 2,:);
    t = 1:size(novel_timed_firing,2);
    hold on
    dofill(t,novel_timed_firing,'blue',1,smval_novrep);
    dofill(t,repeat_timed_firing,'red',1,smval_novrep);
    set(gca,'Xtick',[0:500:8000])
    set(gca,'XtickLabel',num2cell([0:500:8000]-7500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Image Offset(ms)')
    legend('Novel','Repeat')
    xlim([4000 8000]) %plot the 1st half of the image
    yls(:,4) = ylim;
    
    s(5) = subplot(2,4,4);
    t = 1:size(timed_firing{18,unit},2);
    dofill(t,timed_firing{18,unit},'black',1,smval);
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time fom break fixation (ms)')
    yls(:,5) = ylim;
    
    s(6) = subplot(2,4,8);
    t = 1:size(timed_firing{19,unit},2);
    dofill(t,timed_firing{19,unit},'black',1,smval);
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time fom refixation (ms)')
    yls(:,6) = ylim;
    
    max_y = max(yls(2,:));
    for i = 1:6;
        set(s(i),'ylim',[0 max_y])
    end
    
    if multiunit(unit)
        subtitle(['Multiunit ' cfg.channel{unit} ]);
    else
        subtitle(cfg.channel{unit});
    end
end
%% Plot Temporal information for all Sequence Trials

for unit = 1:num_units
    figure
    
    yls = NaN(2,10);
    
    s(1) = subplot(4,4,1);
    t = 1:size(timed_firing{1,unit},2);
    hold on
    dofill(t,timed_firing{1,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{1,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:2000])
    set(gca,'XtickLabel',num2cell([0:250:2000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from ITI onset(ms)')
    yls(:,1) = ylim;
    
    s(2) = subplot(4,4,2);
    t = 1:size(timed_firing{2,unit},2);
    hold on
    dofill(t,timed_firing{2,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{2,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Item 1 Displayed(ms)')
    yls(:,2) = ylim;
    
    s(3) = subplot(4,4,3);
    t = 1:size(timed_firing{3,unit},2);
    hold on
    dofill(t,timed_firing{3,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{3,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Fixated Item 1 (ms)')
    yls(:,3) = ylim;
    
    s(4) = subplot(4,4,4);
    t = 1:size(timed_firing{4,unit},2);
    hold on
    dofill(t,timed_firing{4,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{4,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from IIT 1 start (ms)')
    yls(:,4) = ylim;
    
    s(5) = subplot(4,4,5);
    t = 1:size(timed_firing{5,unit},2);
    hold on
    dofill(t,timed_firing{5,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{5,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Item 2 Displayed(ms)')
    yls(:,5) = ylim;
    
    s(6) = subplot(4,4,6);
    t = 1:size(timed_firing{5,unit},2);
    hold on
    dofill(t,timed_firing{6,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{6,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Fixated Item 2 (ms)')
    yls(:,6) = ylim;
    
    s(7) = subplot(4,4,7);
    t = 1:size(timed_firing{7,unit},2);
    hold on
    dofill(t,timed_firing{7,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{7,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from IIT 2 start (ms)')
    yls(:,7) = ylim;
    
    s(8) = subplot(4,4,8);
    t = 1:size(timed_firing{8,unit},2);
    hold on
    dofill(t,timed_firing{8,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{8,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Item 3 Displayed(ms)')
    yls(:,8) = ylim;
    
    s(9) = subplot(4,4,9);
    t = 1:size(timed_firing{9,unit},2);
    hold on
    dofill(t,timed_firing{9,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{9,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Fixated Item 3 (ms)')
    yls(:,9) = ylim;
    
    s(10) = subplot(4,4,10);
    t = 1:size(timed_firing{10,unit},2);
    hold on
    dofill(t,timed_firing{10,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{10,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from IIT 3 Start (ms)')
    yls(:,10) = ylim;
    
    s(11) = subplot(4,4,11);
    t = 1:size(timed_firing{11,unit},2);
    hold on
    dofill(t,timed_firing{11,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{11,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Item 4 Displayed(ms)')
    yls(:,11) = ylim;
    
    s(12) = subplot(4,4,12);
    t = 1:size(timed_firing{12,unit},2);
    hold on
    dofill(t,timed_firing{12,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{12,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Fixated Item 4 (ms)')
    yls(:,12) = ylim;
    
    s(13) = subplot(4,4,13);
    t = 1:size(timed_firing{13,unit},2);
    hold on
    dofill(t,timed_firing{13,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{13,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:2000])
    set(gca,'XtickLabel',num2cell([0:250:2000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Reward Onset (ms)')
    yls(:,13) = ylim;
    
    %by trial
    means = NaN(2,14);
    stds = NaN(2,14);
    numpoints = NaN(2,14);
    total_spikes = cell(1,2);
    for seq = 1:2
        for event = 1:13;
            trial_spikes = nansum(timed_firing{event,unit}(which_sequence == seq,twin:end),2); %ignore the 1st 500 ms this is from the previous event
            dur = (size(timed_firing{event,unit},2)-500)/1000; %duration in seconds
            means(seq,event) = nanmean(trial_spikes)/dur;
            stds(seq,event) = nanstd(trial_spikes)/dur;
            numpoints(seq,event) = sum(~isnan(trial_spikes));
            if event ~= 1 && event ~= 13
                total_spikes{seq} = [total_spikes{seq};trial_spikes];
            end
        end
    end
    dur = 0.5; %for all epochs but the ITI and the reward
    for seq = 1:2
        means(seq,end) = nanmean(total_spikes{seq})/dur;
        stds(seq,end) = nanstd(total_spikes{seq})/dur;
        numpoints(seq,end) =  sum(~isnan(total_spikes{seq}));
    end
    subplot(4,4,[14 15 16])
    hold on
    bar(means')
    errorb(means',stds'./sqrt(numpoints'));
    hold off
    ylabel('Firing Rate')
    xlabel('Epoch')
    title('Average Firing Rate by Epoch')
    set(gca,'Xtick',1:14)
    set(gca,'XtickLabel',[num2cell(1:13) {'all'}])
    xlim([0 15])
    yl = ylim;
    ylim([0 yl(2)]);
    
    max_y = max(yls(2,:));
    for i = 1:13;
        set(s(i),'ylim',[0 max_y])
    end
    
    if multiunit(unit)
        subtitle(['Multiunit ' cfg.channel{unit} ]);
    else
        subtitle(cfg.channel{unit});
    end
end
%% Plot temporal info by sequence
for unit = 9%1:num_units
    figure
    
    yls = NaN(2,10);
    
    s(1) = subplot(4,4,1);
    t = 1:size(timed_firing{1,unit},2);
    dofill(t,timed_firing{1,unit},'black',1,smval);
    set(gca,'Xtick',[0:250:2000])
    set(gca,'XtickLabel',num2cell([0:250:2000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from ITI onset(ms)')
    yls(:,1) = ylim;
    
    s(2) = subplot(4,4,2);
    t = 1:size(timed_firing{2,unit},2);
    hold on
    dofill(t,timed_firing{2,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{2,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Item 1 Displayed(ms)')
    yls(:,2) = ylim;
    
    s(3) = subplot(4,4,3);
    t = 1:size(timed_firing{2,unit},2);
    hold on
    dofill(t,timed_firing{3,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{3,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Fixated Item 1 (ms)')
    yls(:,3) = ylim;
    
    s(4) = subplot(4,4,4);
    t = 1:size(timed_firing{4,unit},2);
    hold on
    dofill(t,timed_firing{4,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{4,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Item 2 Displayed(ms)')
    yls(:,4) = ylim;
    
    s(5) = subplot(4,4,5);
    t = 1:size(timed_firing{5,unit},2);
    hold on
    dofill(t,timed_firing{5,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{5,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Fixated Item 2 (ms)')
    yls(:,5) = ylim;
    
    s(6) = subplot(4,4,6);
    t = 1:size(timed_firing{6,unit},2);
    hold on
    dofill(t,timed_firing{6,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{6,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Item 3 Displayed(ms)')
    yls(:,6) = ylim;
    
    s(7) = subplot(4,4,7);
    t = 1:size(timed_firing{7,unit},2);
    hold on
    dofill(t,timed_firing{7,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{7,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Fixated Item  3 (ms)')
    yls(:,7) = ylim;
    
    s(8) = subplot(4,4,8);
    t = 1:size(timed_firing{8,unit},2);
    hold on
    dofill(t,timed_firing{8,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{8,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Item 4 Displayed(ms)')
    yls(:,8) = ylim;
    
    s(9) = subplot(4,4,9);
    t = 1:size(timed_firing{9,unit},2);
    hold on
    dofill(t,timed_firing{9,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{9,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:1000])
    set(gca,'XtickLabel',num2cell([0:250:1000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Fixated Item 4 (ms)')
    yls(:,9) = ylim;
    
    s(10) = subplot(4,4,10);
    t = 1:size(timed_firing{10,unit},2);
    hold on
    dofill(t,timed_firing{10,unit}(which_sequence == 1,:),'blue',1,smval);
    dofill(t,timed_firing{10,unit}(which_sequence == 2,:),'red',1,smval);
    hold off
    set(gca,'Xtick',[0:250:2000])
    set(gca,'XtickLabel',num2cell([0:250:2000]-500));
    ylabel('Firing Rate (Hz)')
    xlabel('Time from Reward Onset (ms)')
    yls(:,10) = ylim;
    
    %by trial
    means = NaN(2,10);
    stds = NaN(2,10);
    numpoints = NaN(2,10);
    for event = 1:10;
        for seq = 1:2
            if event ~= 1 
                trial_spikes = nansum(timed_firing{event,unit}(which_sequence == seq,twin:end),2); %ignore the 1st 500 ms this is from the previous event
            else %ITI
                trial_spikes = nansum(timed_firing{event,unit}(:,twin:end),2); %ignore the 1st 500 ms this is from the previous event
            end
            dur = (size(timed_firing{event,unit},2)-twin)/1000; %duration in seconds
            means(seq,event) = nanmean(trial_spikes)/dur;
            stds(seq,event) = nanstd(trial_spikes)/dur;
            numpoints(seq,event) = sum(~isnan(trial_spikes));
        end
    end
    
    subplot(4,4,[11 12])
    hold on
    bar(1:10,means')
    errorb(means',stds'./sqrt(numpoints'));
    hold off
    ylabel('Firing Rate')
    xlabel('Epoch')
    title('Average Firing Rate by Epoch')
    xlim([0 11])

    if unit == 1; %for some reason the code is being stupid on the first fig
            max_y = max(yls(2,:))/2;
    else
            max_y = max(yls(2,:));
    end
    for i = 1:10;
        set(s(i),'ylim',[0 max_y])
    end
    
    if multiunit(unit)
        subtitle(['Multiunit ' cfg.channel{unit} ]);
    else
        subtitle(cfg.channel{unit});
    end
end

%%
filter_width = 4;
binsize = 6;

H = fspecial('gaussian',[filter_width*10+5 filter_width*10+5],filter_width);
fs = 1000;

%calculate the total time spent at any locaitons in binned pixels
spatial_time = zeros(600,800);
dot_x = eye_pos(1:2:end,:);
dot_x = dot_x(trial_type == 2,:);

dot_y = eye_pos(2:2:end,:);
dot_y = dot_y(trial_type == 2,:);
% 
% figure
% for t = 1:size(dot_x,1);
%     hold on
%     plot(dot_x(t,:),dot_y(t,:),'color',[0.8 0.8 0.8])
% end

dot_x = dot_x(1:end);
dot_x(isnan(dot_x)) = [];
dot_y = dot_y(1:end);
dot_y(isnan(dot_y)) = [];
dot_x(dot_x == 0) = 1;
dot_y(dot_y == 0) = 1;
xy_ind = sub2ind(size(spatial_time),round(dot_y),round(dot_x));
for i = 1:length(xy_ind);
    spatial_time(xy_ind(i)) =  spatial_time(xy_ind(i)) +1;
end
spatial_time = spatial_time/fs; %convert from ms to sec

tt = bin2(spatial_time,binsize,binsize);
tt = imfilter(tt,H);
tt(tt == 0) = NaN;
tt = tt(end:-1:1,:);
%%

dot_x = eye_pos(1:2:end,:);
dot_x = dot_x(trial_type == 2,:);
dot_y = eye_pos(2:2:end,:);
dot_y = dot_y(trial_type == 2,:);

for unit = 1:num_units
    firing_location{unit} = zeros(600,800);

    temp_sf = spatial_firing{unit}(trial_type == 2,:);
    [trial,time] = find(temp_sf == 1);
    spikeind = sub2ind(size(dot_x),trial,time);
    xs = dot_x(spikeind);
    ys = dot_y(spikeind);
    
    [xs,ys] = remove_nans(xs,ys);
    for i = 1:length(xs);
        firing_location{unit}(ys(i),xs(i)) =  firing_location{unit}(ys(i),xs(i)) +1;
    end
    
    sf = bin2(firing_location{unit},binsize,binsize);
    sf = imfilter(sf,H);
    sf(sf == 0) = NaN;
    sf = sf(end:-1:1,:);
    figure
    imagesc(sf./tt)
end