%written by Seth Konig September 14, 2016
%code extacts fixation durations, saccade rate, saccade_durations, and saccade amplitudes
%
% code rechecked 1/9/2017 SDK

clar %clear, clc
task = 'ListSQ';
ITIstart_code = 15; %start of ITI
img_on_code= 23; %start of image presentation
img_off_code = 24; %end of image presentation

fltord = 60;
lowpasfrq = 30;
nyqfrq = 1000 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]); %30 Hz low pass filter
buffer = 100;

fltord = 60;
lowpasfrq = 100;
nyqfrq = 1000 ./ 2;
flt2 = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]); %100 Hz low pass filter

set = 0;%set index #

duration_range = [100 250;500 800; 1000 2000];
sequence_fixation_durations = cell(1,85);
sequence_saccade_durations = cell(1,85);
max_dispersion = cell(85,3); %max_xy-min_xy during fixation
smoothed_max_dispersion = cell(85,3); %low passed xy then max_xy-min_xy during fixation
over_smoothed_max_dispersion = cell(85,3);%much more low passed xy then max_xy-min_xy during fixation
total_drift = cell(1,85); %average xy at start of fixation minus average xy at end of fixation, averaged over 25 ms period
max_drift = cell(1,85); %average xy at start of fixation to maximum distance of xy from start of fixation, smoothed with 25 ms moving average

RMS_noise = cell(2,85); %variability in fixation position as a liberal estimate of eye tracking noise
%row 1 is x, row 2 is y
for monkey =1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 20]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for session = 1:length(session_data)
        task_file=get_task_data(session_data{session},task);
        if isempty(task_file)
            continue
        end
        
        %%%%%%%%%%%%%%%%%%%
        %---Import Data---%
        %%%%%%%%%%%%%%%%%%%
        set = set +1; %set index+1
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'cfg','fixationstats','item_file','cnd_file');
        
        %get important task specific information
        [itmlist,sequence_items] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);
                
        
        sequence_fixation_durations{set} = NaN(404,15);
        sequence_saccade_durations{set} = NaN(404,15);
        total_drift{set} = NaN(2,5000);
        max_drift{set} = NaN(2,5000);
        for durs = 1:size(duration_range,1)
            max_dispersion{set,durs}= NaN(1,500);
            smoothed_max_dispersion{set,durs}= NaN(1,500);
            over_smoothed_max_dispersion{set,durs}= NaN(1,500);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---Get Behavioral Stats---%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        num_trials = length(cfg.trl); %number of trials
        
        seq_ind = 1;
        dur_seq_ind = ones(1,size(duration_range,1));
        drift_seq_ind = 1;
        for t = 1:num_trials
            if sum(cfg.trl(t).allval == 3) == 6 %rewarded sequence trial
                
                %---get trial information---%
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code); %start time of trial
                item1on =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %when first item appeared
                reward = cfg.trl(t).alltim(cfg.trl(t).allval == 3)-trial_start; %when rewarded started
                reward = reward(1);
                
                %---get fixation/saccaqde information---%
                fixationtimes = fixationstats{t}.fixationtimes; %fixtaion start and end times
                saccadetimes = fixationstats{t}.saccadetimes; %saccade start and end times
                xy = fixationstats{t}.XY;
                
                %remove fixations that started before item 1 turned on
                invalid= find(fixationtimes(1,:) < item1on);
                fixationtimes(:,invalid) = [];
                invalid= find(saccadetimes(1,:) < item1on);
                saccadetimes(:,invalid) = [];
                
                %remove fixations started ending after reward started
                invalid= find(fixationtimes(2,:) > reward);
                fixationtimes(:,invalid) = [];
                invalid= find(saccadetimes(2,:) > reward);
                saccadetimes(:,invalid) = [];
                
                
                %---Smooth Eye Data---%
                parsed_eyedat = preparse(xy);
                smoothed_xy = [];
                over_smoothed_xy = [];
                for p = 1:length(parsed_eyedat)
                    if any(~isnan(parsed_eyedat{p}(1,:)))
                        
                        %raw velocity
                        x = parsed_eyedat{p}(1,:);
                        y = parsed_eyedat{p}(2,:);
                        
                        x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
                        y = [y(buffer:-1:1) y y(end:-1:end-buffer)];   %add buffer for filtering
                        
                        %low pass filtered velocity 30 Hz cutoff
                        xss = filtfilt(flt,1,x);
                        yss = filtfilt(flt,1,y);
                        xss = xss(101:end-101); %remove buffer after filtering
                        yss = yss(101:end-101); %remove buffer after filtering
                        
                        %low pass filtered velocity 100 Hz cutoff
                        xs = filtfilt(flt2,1,x);
                        ys = filtfilt(flt2,1,y);
                        xs = xs(101:end-101); %remove buffer after filtering
                        ys = ys(101:end-101); %remove buffer after filtering
                        
                        over_smoothed_xy =[over_smoothed_xy [xss;yss]];
                        smoothed_xy =[smoothed_xy [xs; ys]];
                    else
                        smoothed_xy =[smoothed_xy parsed_eyedat{p}];
                        over_smoothed_xy =[over_smoothed_xy parsed_eyedat{p}];
                    end
                end
                
                fixdurs = diff(fixationtimes)+1; %calculate fixation duration
                
                for f = 1:length(fixdurs)
                    if fixdurs(f) < duration_range(1,1)
                        continue
                    end
                    
                    xy_start = mean(xy(:,fixationtimes(1,f):fixationtimes(1,f)+24),2);
                    xy_end = mean(xy(:,fixationtimes(2,f)-24:fixationtimes(2,f)),2);                    
                    total_drift{set}(:,drift_seq_ind) = [fixdurs(f); sqrt(sum((xy_end-xy_start).^2))];
                    
                    distance_x = xy(1,fixationtimes(1,f):fixationtimes(2,f))-xy_start(1);
                    distance_y = xy(2,fixationtimes(1,f):fixationtimes(2,f))-xy_start(2);
                    distance_xy = sqrt(distance_x.^2+distance_y.^2);
                    smoothed_distance_xy = filtfilt(1/25*ones(1,25),1,distance_xy);
                    max_drift{set}(:,drift_seq_ind) = [fixdurs(f); max(smoothed_distance_xy(26:end))];
                    
                    drift_seq_ind = drift_seq_ind+1;
                    
                    range = find(fixdurs(f) > duration_range(:,1) & fixdurs(f) < duration_range(:,2));
                    if isempty(range)
                        continue
                    end
                    
                    fix_xy = xy(:,fixationtimes(1,f):fixationtimes(2,f));
                    disp_xy = max(fix_xy')-min(fix_xy');
                    max_dispersion{set,range}(dur_seq_ind(range)) = sqrt(sum(disp_xy.^2));
                    
                    fix_xy = smoothed_xy(:,fixationtimes(1,f):fixationtimes(2,f));
                    disp_xy = max(fix_xy')-min(fix_xy');
                    smoothed_max_dispersion{set,range}(dur_seq_ind(range)) = sqrt(sum(disp_xy.^2));
                    
                    fix_xy = over_smoothed_xy(:,fixationtimes(1,f):fixationtimes(2,f));
                    disp_xy = max(fix_xy')-min(fix_xy');
                    over_smoothed_max_dispersion{set,range}(dur_seq_ind(range)) = sqrt(sum(disp_xy.^2));
                    
                                    
                    
%                     if range == 3
%                         disp('really long "fixation"')
%                     end
                    
                    dur_seq_ind(range) = dur_seq_ind(range)+1;
                    
                    
                end
                
                
                if length(fixdurs) > 15
                    fixdurs = fixdurs(1:15);
                end
                sequence_fixation_durations{set}(seq_ind,1:length(fixdurs)) = fixdurs;
                
                sacdurs = diff(saccadetimes)+1; %calculate sacccade duration
                if length(sacdurs) > 15
                    sacdurs = sacdurs(1:15);
                end
                sequence_saccade_durations{set}(seq_ind,1:1:length(sacdurs)) = sacdurs;
                
                seq_ind = seq_ind+1;
            end
            
        end
        which_monkey(set) = monkey;
        trial_count(set) = num_trials;
    end
end

%% Fixation and Saaccade Durations

all_seq_fix_durs = [];
all_seq_sac_durs = [];

for set = 1:length(sequence_fixation_durations);
    durs = sequence_fixation_durations{set};
    durs(isnan(durs)) = [];
    all_seq_fix_durs = [all_seq_fix_durs durs];
    
    durs =  sequence_saccade_durations{set};
    durs(isnan(durs)) = [];
    all_seq_sac_durs = [all_seq_sac_durs durs];
end

all_seq_sac_durs(all_seq_sac_durs > 100) = [];
all_seq_fix_durs(all_seq_fix_durs > 1400) = [];

figure
hist(all_seq_sac_durs,25)
xlabel('Saccade Duration (ms)')
ylabel('Saccade Count')
box off

figure
hist(all_seq_fix_durs,50)
xlabel('Fixation Duration (ms)')
ylabel('Fixation Count')
box off

%% Dispersion of Eye during Various Lengthed Fixations
bin_vector = 0:0.1:5;
for range = 1:3;
    all_disp = [];
    all_smoothed_disp =[];
    all_over_smoothed_disp = [];
    for set = 1:length(max_dispersion);
        
        all_disp = [all_disp  max_dispersion{set,range}];
        all_smoothed_disp =[all_smoothed_disp smoothed_max_dispersion{set,range}];
        all_over_smoothed_disp =[all_over_smoothed_disp over_smoothed_max_dispersion{set,range}];
    end
        
    all_disp(all_disp > 5*24) = 5*24;
    [N1,X1] = hist(all_disp/24,bin_vector);

    all_smoothed_disp(all_smoothed_disp > 5*24) = 5*24;
    [N2,X2] = hist(all_smoothed_disp/24,bin_vector);
    
    all_over_smoothed_disp(all_over_smoothed_disp > 5*24) = 5*24;
     [N3,X3] = hist(all_over_smoothed_disp/24,bin_vector);
     
     subplot(2,2,1)
     hold on
     plot(X1,N1/sum(N1))
     hold off
     
     subplot(2,2,2)
     hold on
     plot(X2,N2/sum(N2))
     hold off
     
     subplot(2,2,3)
     hold on
     plot(X3,N3/sum(N3))
     hold off
     
end

subplot(2,2,1)
xlabel('Maximum Disperion (dva)')
ylabel('Probablity')
title('Raw Eye trace')
legend('100-250 ms','500-800 ms','1000+ ms')

subplot(2,2,2)
xlabel('Maximum Disperion (dva)')
ylabel('Probablity')
title('100 Hz Smoothed Eye trace')


subplot(2,2,3)
xlabel('Maximum Disperion (dva)')
ylabel('Probablity')
title('30 Hz Over Smoothed Eye trace')


subtitle('Maximum Disperion/Drift During Fixation of Various Durations')

%% Distance of Drift from start to end of Saccade
all_drift = [];
for set = 1:length(total_drift)
    temp = total_drift{set};
    temp(:,isnan(temp(1,:)))  = [];
    all_drift = [all_drift temp];
end
clear set

bin_durs = [100:100:1300];
binned_drift = cell(1,length(bin_durs)-1);
pct_5 = NaN(1,length(bin_durs)-1);
pct_95 = NaN(1,length(bin_durs)-1);
pct_25 = NaN(1,length(bin_durs)-1);
pct_75 = NaN(1,length(bin_durs)-1);
pct_90 = NaN(1,length(bin_durs)-1);

pct_greater_than_2_5 = NaN(1,length(bin_durs)-1);
pct_greater_than_3_5 = NaN(1,length(bin_durs)-1); %~2.5*sqrt(2)
for bin = 1:length(bin_durs)-1
    these_durs = find( all_drift(1,:) >= bin_durs(bin) & all_drift(1,:) < bin_durs(bin+1));
    binned_drift{bin} = all_drift(2,these_durs);
    
    pct_5(bin) = prctile(all_drift(2,these_durs),5);
    pct_95(bin) = prctile(all_drift(2,these_durs),95);
    pct_90(bin) = prctile(all_drift(2,these_durs),90);
    pct_25(bin) = prctile(all_drift(2,these_durs),25);
    pct_75(bin) = prctile(all_drift(2,these_durs),75);
    
    pct_greater_than_2_5(bin) = sum(all_drift(2,these_durs) > 2.5*24)/length(these_durs);
    pct_greater_than_3_5(bin) = sum(all_drift(2,these_durs) > 3.5*24)/length(these_durs);
end


figure
subplot(2,2,1)
plot(all_drift(1,:),all_drift(2,:)/24,'k.')
xlabel('Fixation Duration')
ylabel('Change in Position From Start-End (dva)')
xlim([100 1600])
ylim([0 5])
box off
title('Scatter of Drift vs Fixation Duration')


subplot(2,2,2)
hold on
plot(bin_durs(1:end-1)+50,cellfun(@mean,binned_drift)/24)
plot(bin_durs(1:end-1)+50,cellfun(@median,binned_drift)/24)
hold off
set(gca,'XMinorTick','on')
xlabel('Fixation Duration (ms)')
ylabel('"Drift" Amplitude (dva)')
legend('Mean','Median','Location','NorthWest')
xlim([100 1300])
title('Average Drift By Fixation Duration')

subplot(2,2,3)
hold on
plot(bin_durs(1:end-1)+50,pct_5/24)
plot(bin_durs(1:end-1)+50,pct_25/24)
plot(bin_durs(1:end-1)+50,cellfun(@median,binned_drift)/24)
plot(bin_durs(1:end-1)+50,pct_75/24)
plot(bin_durs(1:end-1)+50,pct_90/24)
plot(bin_durs(1:end-1)+50,pct_95/24)
hold off
set(gca,'XMinorTick','on')
xlabel('Fixation Duration (ms)')
ylabel('"Drift" Amplitude (dva)')
xlim([100 1300])
legend('5%','25%','50%','75%','90%','95%','Location','NorthWest')
title('Drift Percentile By Fixation Duration')


subplot(2,2,4)
hold on
plot(bin_durs(1:end-1)+50,100*pct_greater_than_2_5)
plot(bin_durs(1:end-1)+50,100*pct_greater_than_3_5)
hold off
set(gca,'XMinorTick','on')
xlabel('Fixation Duration (ms)')
ylabel('% of Fixations')
legend('2.5 dva','2.5*sqrt(2)','Location','NorthWest')
xlim([100 1300])
title('% of Fixations that Could have Left a Fixation Window')

subtitle('Drift From Start to End')
%% Maximum Drift/Distance from Start of Fixation
all_drift = [];
for set = 1:length(max_drift)
    temp = max_drift{set};
    temp(:,isnan(temp(1,:)))  = [];
    all_drift = [all_drift temp];
end
clear set

bin_durs = [100:100:1300];
binned_drift = cell(1,length(bin_durs)-1);
pct_5 = NaN(1,length(bin_durs)-1);
pct_95 = NaN(1,length(bin_durs)-1);
pct_25 = NaN(1,length(bin_durs)-1);
pct_75 = NaN(1,length(bin_durs)-1);
pct_90 = NaN(1,length(bin_durs)-1);

pct_greater_than_2_5 = NaN(1,length(bin_durs)-1);
pct_greater_than_3_5 = NaN(1,length(bin_durs)-1); %~2.5*sqrt(2)
for bin = 1:length(bin_durs)-1
    these_durs = find( all_drift(1,:) >= bin_durs(bin) & all_drift(1,:) < bin_durs(bin+1));
    binned_drift{bin} = all_drift(2,these_durs);
    
    pct_5(bin) = prctile(all_drift(2,these_durs),5);
    pct_95(bin) = prctile(all_drift(2,these_durs),95);
    pct_90(bin) = prctile(all_drift(2,these_durs),90);
    pct_25(bin) = prctile(all_drift(2,these_durs),25);
    pct_75(bin) = prctile(all_drift(2,these_durs),75);
    
    pct_greater_than_2_5(bin) = sum(all_drift(2,these_durs) > 2.5*24)/length(these_durs);
    pct_greater_than_3_5(bin) = sum(all_drift(2,these_durs) > 3.5*24)/length(these_durs);
end



figure
subplot(2,2,1)
plot(all_drift(1,:),all_drift(2,:)/24,'k.')
xlabel('Fixation Duration')
ylabel('Change in Position From Start-End (dva)')
xlim([100 1600])
ylim([0 5])
box off
title('Scatter of Drift vs Fixation Duration')


subplot(2,2,2)
hold on
plot(bin_durs(1:end-1)+50,cellfun(@mean,binned_drift)/24)
plot(bin_durs(1:end-1)+50,cellfun(@median,binned_drift)/24)
hold off
set(gca,'XMinorTick','on')
xlabel('Fixation Duration (ms)')
ylabel('"Drift" Amplitude (dva)')
legend('Mean','Median','Location','NorthWest')
xlim([100 1300])
title('Average Drift By Fixation Duration')

subplot(2,2,3)
hold on
plot(bin_durs(1:end-1)+50,pct_5/24)
plot(bin_durs(1:end-1)+50,pct_25/24)
plot(bin_durs(1:end-1)+50,cellfun(@median,binned_drift)/24)
plot(bin_durs(1:end-1)+50,pct_75/24)
plot(bin_durs(1:end-1)+50,pct_90/24)
plot(bin_durs(1:end-1)+50,pct_95/24)
hold off
set(gca,'XMinorTick','on')
xlabel('Fixation Duration (ms)')
ylabel('"Drift" Amplitude (dva)')
xlim([100 1300])
legend('5%','25%','50%','75%','90%','95%','Location','NorthWest')
title('Drift Percentile By Fixation Duration')


subplot(2,2,4)
hold on
plot(bin_durs(1:end-1)+50,100*pct_greater_than_2_5)
plot(bin_durs(1:end-1)+50,100*pct_greater_than_3_5)
hold off
set(gca,'XMinorTick','on')
xlabel('Fixation Duration (ms)')
ylabel('% of Fixations')
legend('2.5 dva','2.5*sqrt(2)','Location','NorthWest')
xlim([100 1300])
title('% of Fixations that Could have Left a Fixation Window')

subtitle('Maximum Drift from Start during Whole Fixation')