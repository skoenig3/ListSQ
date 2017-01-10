%written by Seth Konig September 14, 2016
%code extacts fixation durations, saccade rate, saccade_durations, and saccade amplitudes
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
saccade_durations = cell(2,85);

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
        saccade_durations{1,set} = NaN(96,40);
        saccade_durations{2,set} = NaN(96,40);
        
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
                
                sacdurs = diff(saccadetimes)+1; %calculate sacccade duration
                saccade_durations{nvr,set}(which_image,1:size(saccadetimes,2)) = sacdurs;
                
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
all_fix_durs = []; %all fixaiton durations
nov_fix_durs = NaN(size(fixation_durations,2),20); %median by set novel fixation durations
rep_fix_durs = NaN(size(fixation_durations,2),20); %median by set repeat fixation durations
for set = 1:size(fixation_durations,2);
    nov_durs = fixation_durations{1,set}; %novel images
    rep_durs = fixation_durations{2,set}; %repeat images
    
    %remove fixations shorter than 100 msin duration since removed from analysis
    nov_durs(nov_durs < 100) = NaN; 
    rep_durs(rep_durs < 100) = NaN;
        
    nov_fix_durs(set,:) = nanmedian(nov_durs(:,1:20)); %novel
    rep_fix_durs(set,:) = nanmedian(rep_durs(:,1:20)); %repeat
    
    all_fix_durs = [all_fix_durs; nov_durs; rep_durs]; %all fixation durations
end

%---Stats Test---%
p_wilx = signrank(median(nov_fix_durs),median(rep_fix_durs)); %Wilcoxon signed rank test for zero median between novel and repeated images
%nonparametric, repeated measures (across multiple fixations), analysiss
p_vals = [];
for f = 1:size(nov_fix_durs,2)
    [~,p_vals(f)] = ttest(nov_fix_durs(:,f),rep_fix_durs(:,f)); %paired ttest, not sure if really valid but easy to interpret
end

%---Plot Results---%
figure
hold on
plot(nanmean(nov_fix_durs))
errorb(1:20,nanmean(nov_fix_durs),nanstd(nov_fix_durs)./sqrt(sum(~isnan(nov_fix_durs))),'color','b')
plot(nanmean(rep_fix_durs),'r')
errorb(1:20,nanmean(rep_fix_durs),nanstd(rep_fix_durs)./sqrt(sum(~isnan(rep_fix_durs))),'color','r')
for f = 1:size(nov_fix_durs,2)
   if p_vals(f) < 0.05/size(nov_fix_durs,2) %Bonferroni correction
      plot(f,215,'k*')
   end
end
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
xlim([0 21])
ylim([160 220])
axis square
legend('Novel','Repeat')
title(['PW and TO n_{sessions} = ' num2str(size(fixation_durations,2)) ', p_{Wilcoxon} = ' num2str(p_wilx,3)])

%% Distribution of All Fixation Durations
durs = all_fix_durs;
durs(durs > 400) = [];
figure
hist(durs,301)
hold on
plot([round(nanmean(durs)) round(nanmean(durs))],[0 2000],'k--')
hold off
xlabel('Fixation Duration (ms)')
ylabel('Count')
title(sprintf(['Distribution of Fixation Durations \n'...
    'mean = ' num2str(round(nanmean(durs))) ' ms median = ' num2str(round(nanmedian(durs))) ' ms']))
axis square

%% Distribution of All Saccade Durations
%novel/differences really don't exist (effect ~1-2 ms)
all_sac_durs = [];
for set = 1:size(saccade_durations,2);
    all_sac_durs = [all_sac_durs; saccade_durations{1,set}; saccade_durations{2,set}];
end

durs =  all_sac_durs;
durs(durs > 80) = [];
figure
hist(durs,70)
hold on
plot([round(nanmean(durs)) round(nanmean(durs))],[0 14000],'k--')
hold off
xlim([10 80])
xlabel('saccade Duration (ms)')
ylabel('Count')
title(sprintf(['Distribution of Sacccade Durations \n'...
    'mean = ' num2str(round(nanmean(durs))) ' ms median = ' num2str(round(nanmedian(durs))) ' ms']))
axis square
box off
%% Distribution of Saccade Ampltiudes

nov_sac_amps = NaN(size(saccade_amplitudes,2),20); %median by set novel saccade_amplitude
rep_sac_amps = NaN(size(saccade_amplitudes,2),20); %median by set repeat saccade_amplitude
for set = 1:size(saccade_amplitudes,2);
    nov_sac = saccade_amplitudes{1,set}/24; %novel  images and convert from pixel to dva
    rep_sac = saccade_amplitudes{2,set}/24; %repeat images and convert from pixel to dva
    
    %remove saccades smaller than 2 dva 
    nov_sac(nov_durs < 2) = NaN; 
    rep_sac(rep_durs < 2) = NaN;
        
    nov_sac_amps(set,:) = nanmedian(nov_sac(:,1:20)); %novel
    rep_sac_amps(set,:) = nanmedian(rep_sac(:,1:20)); %repeat
end


%---Stats Test---%
p_wilx = signrank(median(nov_sac_amps),median(rep_sac_amps)); %Wilcoxon signed rank test for zero median between novel and repeated images
%nonparametric, repeated measures (across multiple fixations), analysiss
p_vals = [];
for f = 2:size(nov_sac_amps,2)
    [~,p_vals(f)] = ttest(nov_sac_amps(:,f),rep_sac_amps(:,f)); %paired ttest, not sure if really valid but easy to interpret
end

%---Plot Results---%
figure
hold on
plot(nanmean(nov_sac_amps))
errorb(1:20,nanmean(nov_sac_amps),nanstd(nov_sac_amps)./sqrt(sum(~isnan(nov_sac_amps))),'color','b')
plot(nanmean(rep_sac_amps),'r')
errorb(1:20,nanmean(rep_sac_amps),nanstd(rep_sac_amps)./sqrt(sum(~isnan(rep_sac_amps))),'color','r')
for f = 1:size(nov_sac_amps,2)
   if p_vals(f) < 0.05/size(nov_sac_amps,2) %Bonferroni correction
      plot(f,8.5,'k*')
   end
end
hold off
xlabel('Ordinal Saccade #')
ylabel('Saccade Amplitude (dva)')
xlim([0 21])
ylim([5 9])
axis square
legend('Novel','Repeat')
title(['PW and TO n_{sessions} = ' num2str(size(saccade_amplitudes,2)) ', p_{Wilcoxon} = ' num2str(p_wilx,3)])
