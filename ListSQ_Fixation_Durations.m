%written by Seth Konig September 14, 2016
%code extacts fixation durations, saccade rate, and saccade amplitudes
%
% code rechecked 1/9/2017 SDK

clar
task = 'ListSQ';
ITIstart_code = 15; %start of ITI
img_on_code= 23; %start of image presentation
img_off_code = 24; %end of image presentation

set = 0;%set index #
fixation_durations = cell(2,85);
saccade_amplitudes = cell(2,85);

for monkey = 1:2
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
        
        predict_rt = 206;%206 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML[
        
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
        
        
        %set the image duration
        if str2double(task_file(3:8)) < 140805
            imgdur = 7000;
        else
            imgdur = 5000;
        end
        
        %---Prealoccate Memory---%
        fixation_durations{1,set} = NaN(96,40);
        fixation_durations{2,set} = NaN(96,40);
        saccade_amplitudes{1,set} = NaN(96,40);
        saccade_amplitudes{2,set} = NaN(96,40);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---Get Behavioral Stats---%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        num_trials = length(cfg.trl); %number of trials
        for t = 1:num_trials
            if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                
                %---get trial information---%
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code); %start time of trial
                imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %when image turned on
                imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start; %when image turned off
                
                % if monkey isn't paying attention data isn't probably
                % worth much plus have to cut off somewhere
                if imgoff-imgon > 1.5*imgdur-1
                    imgoff = imgon+1.5*imgdur-1;
                end
                
                %get novel/repeat information
                img_index = find(img_cnd == cfg.trl(t).cnd); %get image index
                if any(isnan(which_img(img_index)))%presentation error so skip trial, see get_image_numbers.m
                    continue
                end
                which_image = which_img(img_index); %image #
                nvr = novel_vs_repeat(img_index); %whether image was novel or repeat
                
                
                %---get fixation/saccaqde information---%
                fixationtimes = fixationstats{t}.fixationtimes; %fixtaion start and end times
                fixations = fixationstats{t}.fixations; %fixation locations
                saccadetimes = fixationstats{t}.saccadetimes; %saccade start and end times
                xy = fixationstats{t}.XY;
                
                %remove fixations that started before image turned on
                invalid= find(fixationtimes(1,:) < imgon);
                fixationtimes(:,invalid) = [];
                fixations(:,invalid) = [];
                invalid= find(saccadetimes(1,:) < imgon);
                saccadetimes(:,invalid) = [];
                
                %remove fixations started ending after image turned off
                invalid= find(fixationtimes(2,:) > imgoff);
                fixationtimes(:,invalid) = [];
                fixations(:,invalid) = [];
                invalid= find(saccadetimes(2,:) > imgoff);
                saccadetimes(:,invalid) = [];
                
                fixdurs = diff(fixationtimes)+1; %calculate fixation duration
                fixation_durations{nvr,set}(which_image,1:size(fixationtimes,2)) = fixdurs;
                
                for f = 2:size(fixationtimes,2)%ignore first fixation not sure where it was/possibly contaminated anyway
                    prior_sac = find(saccadetimes(2,:) == fixationtimes(1,f)-1);%next fixation should start immediately after
                    if isempty(prior_sac) %no prior saccade so was proabbly looking off screen
                        continue;
                    end
                    sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %saccade amplitude
                    saccade_amplitudes{nvr,set}(which_image,f) = sacamp;
                end
            end
        end
    end
end
%% Fixuation Durations by Ordinal Fixation #
all_durs = cell(1,2); %all fixation durations
nov = NaN(size(fixation_durations,2),15);
rep = NaN(size(fixation_durations,2),15);
ordinal =  NaN(size(fixation_durations,2),15);
for set = 1:size(fixation_durations,2);
    all_durs{1} = [all_durs{1}; fixation_durations{1,set}(:,1:40)];
    all_durs{2} = [all_durs{2}; fixation_durations{2,set}(:,1:40)];
    nov(set,:) = nanmedian(fixation_durations{1,set}(:,1:15));
    rep(set,:) = nanmedian(fixation_durations{2,set}(:,1:15));
    ordinal(set,:) = 1:15;
end
% all_durs{1}(all_durs{1} > 500) = NaN;
% all_durs{2}(all_durs{2} > 500) = NaN;
%
figure
% subplot(1,2,1)
hold on
plot(nanmean(nov))
errorb(1:20,nanmean(nov),nanstd(nov)./sqrt(sum(~isnan(nov))),'color','b')
plot(nanmean(rep),'r')
errorb(1:20,nanmean(rep),nanstd(rep)./sqrt(sum(~isnan(rep))),'color','r')
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
xlim([0 21])
legend('Novel','Repeat')
% title('PW and TO')
title('Vivian')

durs = all_durs{1}(1:end);
durs(durs > 500) = [];
figure
hist(durs,475)
hold on
plot([round(nanmean(durs)) round(nanmean(durs))],[0 1000],'k--')
hold off
xlabel('Fixation Duration (ms)')
ylabel('Count')
title(sprintf(['Distribution of Fixation Durations \n'...
    'mean = ' num2str(round(nanmean(durs))) 'ms median = ' num2str(round(nanmedian(durs))) ' ms']))

%%
% group = [ones(numel(nov),1); 2*ones(numel(rep),1)];
% subject = [ordinal(:); ordinal(:)];
% groupsubj = [group subject];
% y = [nov(:); rep(:)];
%
% [P,T,STATS,TERMS] = anovan(y,groupsubj,'random',[2]);
% multcompare(STATS)

p_vals = [];
for f = 1:size(nov,2)
    [~,p_vals(f)] = ttest(nov(:,f),rep(:,f))
end
signrank(median(nov),median(rep))
friedman test
%%
all_rates = [];
for set = 1:size(fixation_durations,2);
    all_rates = [all_rates; saccade_rate{1,set}(:);saccade_rate{2,set}(:)];
end
all_rates(isnan(all_rates)) = [];
all_rates = 1000./all_rates;
all_rates(all_rates > 12) = [];
all_rates(all_rates < 1) = [];

subplot(1,2,2)
hist(all_rates,[1:0.5:12])
xlabel('Saccade Rate (Hz)')
ylabel('Count')
title(['Tobii Mean: ' num2str(mean(all_rates),3) ', Median: ' num2str(nanmedian(all_rates),3)])

%%
nov_durs = all_durs{1}(1:end);
nov_durs(nov_durs > 500) = [];
rep_durs = all_durs{2}(1:end);
rep_durs(rep_durs > 500) = [];

figure
N_nov = hist(nov_durs,[25:5:500]);
N_rep = hist(rep_durs,[25:5:500]);
hold on
plot([25:5:500],N_nov/sum(N_nov),'b')
plot([nanmedian(nov_durs) nanmedian(nov_durs)],[0 0.035],'b--')
plot([25:5:500],N_rep/sum(N_rep),'r')
plot([nanmedian(rep_durs) nanmedian(rep_durs)],[0 0.0325],'r--')
hold off
xlabel('Fixation Duration (ms)')
ylabel('Probability')
legend('Novel','Repeat')
title('PW and TO')



%% Saccade Durations (parallels increase in saccade amplitdue effect is small)
all_durs = cell(1,2);
for set = 1:size(saccade_durations,2);
    all_durs{1} = [all_durs{1}; saccade_durations{1,set}(:,1:40)];
    all_durs{2} = [all_durs{2}; saccade_durations{2,set}(:,1:40)];
end

figure
hold on
plot(nanmean(all_durs{1}(:,1:20)))
errorb(1:20,nanmean(all_durs{1}(:,1:20)),nanstd(all_durs{1}(:,1:20))./sqrt(sum(~isnan(all_durs{1}(:,1:20)))),'color','b')
plot(nanmean(all_durs{2}(:,1:20)),'r')
errorb(1:20,nanmean(all_durs{2}(:,1:20)),nanstd(all_durs{2}(:,1:20))./sqrt(sum(~isnan(all_durs{2}(:,1:20)))),'color','r')
hold off
xlabel('Ordinal saccade #')
ylabel('saccade Duration (ms)')
xlim([0 21])
legend('Novel','Repeat')
title('PW and TO')
%%
durs = all_durs{1}(1:end);
durs(durs > 80) = [];
figure
hist(durs,70)
hold on
plot([round(nanmean(durs)) round(nanmean(durs))],[0 8000],'k--')
hold off
xlabel('saccade Duration (ms)')
ylabel('Count')
title(sprintf(['Distribution of Sacccade Durations \n'...
    'mean = ' num2str(round(nanmean(durs))) 'ms median = ' num2str(round(nanmedian(durs))) ' ms']))

