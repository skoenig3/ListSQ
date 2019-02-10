clar

cortex_file_names = {'TO151204.3','TO151208.3','TO160120.3','TO160202.3','TO160322.3',...
    'PW140806.3','PW141007.3'};

samprate = 5;
imageX = 800;
imageY = 600;

Fs = 1000;
[bl,al] = butter(8,60/(Fs/2),'low');


for file = 5%1:length(cortex_file_names)
    
    if strcmpi(cortex_file_names{file}(1:2),'PW')
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
    elseif strcmpi(cortex_file_names{file}(1:2),'TO')
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
    else
        error('Wrong File name')
    end
    
    %load neural data
    load([data_dir cortex_file_names{file}(1:end-2) '_' cortex_file_names{file}(end) '-preprocessed.mat']);
    
    %load cortex data
    [time_arr,event_arr,eog_arr,~,~,~]  = get_ALLdata([data_dir cortex_file_names{file}]);
    
    [itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file);
    
    numrpt = size(event_arr,2);
    new_eog_arr=[];
    valrptcnt = 0;
    clear per clrchgind
    for rptlop = 1:numrpt
        if event_arr((find(event_arr(:,rptlop)>500,1,'last')),rptlop)-500 > 1 %1st block is color change always
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
            cnd = event_arr(cndnumind,rptlop)-1000;
            priorcndnumind = find(event_arr(:,rptlop-1) >= 1000 & event_arr(:,rptlop-1) <=2000);
            priorcnd =  event_arr(priorcndnumind,rptlop-1)-1000;
            
            perbegind = find(event_arr(:,rptlop) == 100); % eye data starts at 100; might have predictive looking
            perendind = find(event_arr(:,rptlop) == 101); % eye data stops collecting after rewards so can stop here
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
            begtimdum = time_arr(perbegind,rptlop);
            endtimdum = time_arr(perendind(end),rptlop);
            
            if endtimdum > begtimdum
                valrptcnt = valrptcnt + 1;
                per(valrptcnt).begsmpind = begtimdum;
                per(valrptcnt).endsmpind = endtimdum(1);
                per(valrptcnt).begpos = 1;
                per(valrptcnt).cnd = event_arr(cndnumind,rptlop)-1000;
                per(valrptcnt).blk = event_arr(blknumind,rptlop)-500;
                per(valrptcnt).allval = event_arr(:,rptlop);
                per(valrptcnt).alltim = time_arr(:,rptlop);
                per(valrptcnt).event = rptlop;
                new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
            end
        end
    end
    
    %---get eye data for only when fixation cross or picture is displayed---%
    eyedat = cell(1,length(per));
    cnd=[];
    teststart = [];
    for trlop=1:size(per,2)
        trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
        horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
        vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
        picstart=1*samprate;
        picend=per(trlop).endsmpind-per(trlop).begsmpind;
        picend(picend > samprate*length(horeog)) = length(horeog)*samprate;
        
        eyedat{trlop}(1,:) = horeog(round(picstart/samprate):floor(picend/samprate));
        eyedat{trlop}(2,:) = vrteog(round(picstart/samprate):floor(picend/samprate));
        cnd(trlop)=per(trlop).cnd;
    end
    %%
    %---Recalibrate and automatically scale eye data---%
    outside_xy = cell(1,length(eyedat)); %values outside image
    for t = 1:length(eyedat)
        x = eyedat{t}(1,:);
        y = eyedat{t}(2,:);
        
        outside_xy{t} = NaN(2,length(x));
        
        [x,y] = tformfwd(tform,x,y);
        
        x = 24*x; %convert from cortex dva to pixels
        y = 24*y; %convert from cortex dva to pixels
        x = x+imageX/2;
        y = y+imageY/2;
        
        %get x-values outside image border
        outside_xy{t}(1,x < 1) = x(x < 1);
        outside_xy{t}(2,x < 1) = y(x < 1);
        outside_xy{t}(1,x > imageX) = x(x > imageX);
        outside_xy{t}(2,x > imageX) = y(x > imageX);
        
        y(x < 1) = NaN;
        x(x < 1) = NaN;
        y(x > imageX) = NaN;
        x(x > imageX)= NaN;
        
        %get y-values outsie image border
        outside_xy{t}(1,y < 1) = x(y < 1);
        outside_xy{t}(2,y < 1) = y(y < 1);
        outside_xy{t}(1,y > imageY) = x(y > imageY);
        outside_xy{t}(2,y > imageY) = y(y > imageY);
        
        x(y < 1) = NaN;
        y(y < 1) = NaN;
        x(y > imageY) = NaN;
        y(y > imageY) = NaN;
        
        eyedat{t} = [x;y];
    end
    
    %%
    if length(per) ~= length(cfg.trl)
        if length(cfg.trl) < length(per)-1%%cutoff has occured before
            error('Number of trials in cortex data ~= number of trials in recording data');
        else%confirm that all cnds are the same
            for t= 1:min(length(cfg.trl),length(per))
                if cfg.trl(t).cnd-1000~=per(t).cnd
                    error('CNDs dont match')
                end
            end
            num_trials = min(length(cfg.trl),length(per));
        end
    else
        num_trials = length(per);
    end
    xxcorrs = NaN(length(per),501);
    yxcorrs = NaN(length(per),501);
    %%
    for t = 1:num_trials
        
        if per(t).alltim(1) > 1
            disp('now')
        end
        
        cortex_eye_start = per(t).alltim(per(t).allval == 100);
        cortex_stim_on = per(t).alltim(per(t).allval == 23)-cortex_eye_start;
        if isempty(cortex_stim_on) %no stimuli were turned on proabbly won't have good eye data then
            continue
        end
        
        if any(per(t).allval == 3)
            cortex_stim_off = per(t).alltim(per(t).allval == 3);
            cortex_stim_off = cortex_stim_off(end-1);
        else
            cortex_stim_off = per(t).alltim(per(t).allval == 24);
        end
        cortex_stim_off = cortex_stim_off-cortex_eye_start;
        
        %cortex_eye_dat = eyedat{t};
        
        %parsed_eyedat = preparse(eyedat);
        x = eyedat{t}(1,:);
        y = eyedat{t}(2,:);
        
        x = resample(x,5,1);%up sample to 1000 Hz
        y = resample(y,5,1);%up sample to 1000 Hz
        
        x = x(cortex_stim_on:cortex_stim_off);
        y = y(cortex_stim_on:cortex_stim_off);
        
        
        %%
        record_stim_on = per(t).alltim(per(t).allval == 23);
        if isempty(record_stim_on) %no stimuli were turned on proabbly won't have good eye data then
            %continue
        end
        
        if any(per(t).allval == 3)
            record_stim_off = per(t).alltim(per(t).allval == 3);
            record_stim_off = record_stim_off(end-1);
        else
            record_stim_off = per(t).alltim(per(t).allval == 24);
        end
        
        xx = fixationstats{t}.XY(1,record_stim_on:record_stim_off);
        yy = fixationstats{t}.XY(2,record_stim_on:record_stim_off);
        
        if length(x) < 1000 || length(xx) < 1000
            continue
        end
        %% Filter so can calculate velocity
        fx = x;
        fxx = xx;
        fy = y;
        fyy = yy;
        %         fx  = filtfilt(bl,al,x);
        %         fxx = filtfilt(bl,al,xx);
        %
        %         fy  = filtfilt(bl,al,y);
        %         fyy = filtfilt(bl,al,yy);
        %% Computer Velocity since cross corr deals with this better than xy data
        velx = diff(fx);
        velxx = diff(fxx);
        
        vely = diff(fy);
        velyy = diff(fyy);
        %% Find longest string of valid data between both recording and cortex data...should be similar timing within 100 ms
        valid = ~isnan(velx) & ~isnan(velxx);
        if sum(valid) == 0
            continue
        end
        gaps = findgaps(find(valid == 1));
        if size(gaps,2) < 1000 %% not enough eye data to match
            continue
        end
        if size(gaps,1) == 1;
            [rx,~]=xcorr(velx,velxx,250,'coeff');
            [ry,~]=xcorr(vely,velyy,250,'coeff');
        else
            sumzeros = sum(gaps' ~= 0);
            largest = find(sumzeros == max(sumzeros));
            velx = velx(gaps(largest));
            velxx = velxx(gaps(largest));
            vely = vely(gaps(largest));
            velyy = velyy(gaps(largest));
            
            [rx,~]=xcorr(velx,velxx,250,'coeff');
            [ry,~]=xcorr(vely,velyy,250,'coeff');
        end
        
        xxcorrs(t,:) = rx;
        yxcorrs(t,:) = ry;
    end
    %%
    xxcorrs = laundry(xxcorrs);
    yxcorrs = laundry(yxcorrs);
    
    figure
    subplot(2,2,1)
    imagesc(-250:250,1:size(xxcorrs,1),xxcorrs);
    hold on
    plot([0 0],[0 size(xxcorrs,1)],'k--')
    hold off
    xlabel('Lag (ms)')
    ylabel('Trial #')
    colormap('viridis')
    box off
    title('All Correlation Coefficient: Recording vs Cortex X Eye Data')
    
    subplot(2,2,2)
    imagesc(-250:250,1:size(yxcorrs,1),yxcorrs);
    hold on
    plot([0 0],[0 size(yxcorrs,1)],'k--')
    hold off
    xlabel('Lag (ms)')
    ylabel('Trial #')
    colormap('viridis')
    box off
    title('All Correlation Coefficient: Recording vs Cortex Y Eye Data')
    
    [PKS,LOCS] = findpeaks(mean(xxcorrs));
    LOCS = LOCS(PKS == max(PKS));
    subplot(2,2,3)
    plot(-250:250,mean(xxcorrs));
    hold on
    yl = ylim;
    plot([0 0],[yl(1) yl(2)],'k--')
    hold off
    xlim([-250 250])
    xlabel('Lag (ms)')
    ylabel('Correlation Coefficient')
    colormap('viridis')
    box off
    title(sprintf(['Average Correlation Coefficient: Recording vs Cortex X Eye Data \n maximum correlation @ ' num2str(LOCS+-250) ' ms']))
    
    [PKS,LOCS] = findpeaks(mean(yxcorrs));
    LOCS = LOCS(PKS == max(PKS));
    subplot(2,2,4)
    plot(-250:250,mean(yxcorrs));
    hold on
    yl = ylim;
    plot([0 0],[yl(1) yl(2)],'k--')
    hold off
    xlim([-250 250])
    xlabel('Lag (ms)')
    ylabel('Correlation Coefficient')
    colormap('viridis')
    box off
    title(sprintf(['Average Correlation Coefficient: Recording vs Cortex Y Eye Data \n maximum correlation @ ' num2str(LOCS+-250) ' ms']))
    
    subtitle(cortex_file_names{file})
    %%
end
%%
figure
subplot(1,2,1)
plot(1:length(fx),fx,'b')
hold on
plot(1:length(fxx),fxx,'r')
plot([1:length(fxx)]-80,fxx,'g')
hold off
xlabel('Time from Stim On (ms)')
ylabel('Eye Position (px)')
title('X')
box off
axis square

subplot(1,2,2)
plot(1:length(fy),fy,'b')
hold on
plot(1:length(fyy),fyy,'r')
plot([1:length(fyy)]-80,fyy,'g')
hold off
xlabel('Time from Stim On (ms)')
ylabel('Eye Position (px)')
title('Y')
box off
axis square

legend('Cortex','Recording','Recording Shifted')
subtitle(cortex_file_names{file})