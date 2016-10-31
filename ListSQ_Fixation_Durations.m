%written by Seth Konig September 14, 2016
%code extacts fixation durations

task = 'ListSQ';
ITIstart_code = 15;
img_on_code= 23;
img_off_code = 24;


set = 0;
fixation_durations = cell(2,85);
saccade_durations = cell(2,85);
saccade_rate = cell(2,85);
for monkey = 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML[
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for session = 1:length(session_data)
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{session},task);
        if isempty(task_file)
            continue
        end
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'cfg','fixationstats');
        
        %set the image duration
        if str2double(task_file(3:8)) < 140805
            imgdur = 7000;
        else
            imgdur = 5000;
        end
        
        
        
        set = set +1;
        fixation_durations{1,set} = NaN(96,40);
        fixation_durations{2,set} = NaN(96,40);
        
        saccade_durations{1,set} = NaN(96,40);
        saccade_durations{2,set} = NaN(96,40);
                
        saccade_rate{1,set} = NaN(96,40);
        saccade_rate{2,set} = NaN(96,40);
        
        %get important task specific information
        [itmlist,sequence_items] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);
        
        num_trials = length(cfg.trl);
        for t = 1:num_trials
            if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %want to avoid visual response
                imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start;
                
                % if monkey isn't paying attention data isn't probably
                % worth much plus have to cut off somewhere
                if imgoff-imgon > 1.5*imgdur-1
                    imgoff = imgon+1.5*imgdur-1;
                end
                
                img_index = find(img_cnd == cfg.trl(t).cnd);
                if isempty(img_index) || any(isnan(img_index))
                    continue
                end
                which_images = which_img(img_index);
                nvr = novel_vs_repeat(img_index);
                
                if any(isnan(which_img(img_index)))
                    continue
                end
                
                
                fixationtimes = fixationstats{t}.fixationtimes;
                saccadetimes = fixationstats{t}.saccadetimes;
                
                
                %find fiations and saccades that did not occur during the image period;
                %should also take care of the 1st fixation on the crosshair
                
                %fixation started before image turned on
                invalid= find(fixationtimes(1,:) < imgon);
                fixationtimes(:,invalid) = [];
                invalid= find(saccadetimes(1,:) < imgon);
                saccadetimes(:,invalid) = [];
                
                %fixation started after the image turned off and/or firing rate could corrupted by image turning off
                invalid= find(fixationtimes(1,:) > imgoff);
                fixationtimes(:,invalid) = [];
                invalid= find(saccadetimes(1,:) > imgoff);
                saccadetimes(:,invalid) = [];
                
                fixdurs = diff(fixationtimes)+1;
                fixation_durations{nvr,set}(which_images,1:size(fixationtimes,2)) = fixdurs;
                
                sacdurs = diff(saccadetimes)+1;
                saccade_durations{nvr,set}(which_images,1:size(saccadetimes,2)) = sacdurs;
                
                saccade_rate{nvr,set}(which_images,1:size(fixationtimes,2)-1) = diff(fixationtimes(1,:));
            end
        end
    end
end
%% Fixuation Durations by Ordinal Fixation
all_durs = cell(1,2);
for set = 1:size(fixation_durations,2);
    all_durs{1} = [all_durs{1}; fixation_durations{1,set}(:,1:40)];
    all_durs{2} = [all_durs{2}; fixation_durations{2,set}(:,1:40)];
end
% all_durs{1}(all_durs{1} > 500) = NaN;
% all_durs{2}(all_durs{2} > 500) = NaN;
%%
% 
% figure
subplot(1,2,1)
hold on
plot(nanmean(all_durs{1}(:,1:20)))
errorb(1:20,nanmean(all_durs{1}(:,1:20)),nanstd(all_durs{1}(:,1:20))./sqrt(sum(~isnan(all_durs{1}(:,1:20)))),'color','b')
plot(nanmean(all_durs{2}(:,1:20)),'r')
errorb(1:20,nanmean(all_durs{2}(:,1:20)),nanstd(all_durs{2}(:,1:20))./sqrt(sum(~isnan(all_durs{2}(:,1:20)))),'color','r')
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

