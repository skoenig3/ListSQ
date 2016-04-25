function make_rasters_and_plot_waveformsV2(data_dir,figure_dir,session_data,task)
% written by Seth Konig August 2014. Updated January 2016 by SDK
% 1) Plots spike trains aligned to trial start. Below is the PSTH.
% 2) Plots waveform divided into 1/4 of the task at time
% 3) Plots firing rate over time/trials
% 4) defines which trials should be used for data analysis if firing rate
% is not stable

%aligns data to item 1 on for sequence task
%aligns data to image onset for image task
%aligns data to dot onset for cvtnew task
%also aligns data to ITI since there seems to be a lot of reward/ITI neurons

%don't want to waste processing power on units that are a) poorly isolated,
% b) don't have enough spikes, and c) aren't stable for long enough to do
% any analysis with sufficient confidence. But should plot the rasters
% nevertheless.
%
% !!!WARNING CODE WILL SET VIABLE TRILAS TO NULL FOR THESE UNITS TO SAVE PROCESSING TIME!!!!
%
too_few_folder = [figure_dir '\Rasters Too Few Trials or Spikes\']; %where to put rasters see line 47
PoorIsolationFolder = [figure_dir '\PoorIsolation\']; %where to put poorly isolated units
MultiUnit_folder = [figure_dir '\MultiUnit\'];

reward_code = 3;
ITI_code = 15;
dot_on_code = 25;%for cvtnew
trial_start_code = 23; %image on and item 1 on
binsize = 35;%ms probably want 25-100 ms
trialbinsize = 6;%averging data by trial number

%for neurons with too few spikes or too few stable trials. Hard to properly
%isolate but also over isolation period has a very low average firing rate.
% Prior experinece suggest these neurons typically don't fire except
% sporadically in burst or on specifc fixations. May be useful to analyze
% in the future.
minimum_spikes = 100; %currently setup for all trials excludes before and after task was running
poor_isolation_threshold = 2;%0,1,2 on cluster cutting quality
confidence_threshold = 0.79;% percent condifence that I'm willing to take anything over 80%

%get important task related data
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
    get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'cfg','data','hdr','fixationstats','waveforms');


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

%get unit data
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);

num_trials = length(cfg.trl);
switch task
    case 'ListSQ'
        [itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        
        valid_trials = NaN(2,num_units);%trials we want to keep based on firing rate
        stability_attribute = ones(1,num_units); %1 stable and using, 2 insufficient number of spikes,
        %3 not stable for sufficient trials, 4 poorly isolated/not confident unit is real
        for unit = 1:num_units
            
            type = NaN(1,length(cfg.trl));
            allspikes = NaN(length(cfg.trl),7000);
            allspikesITI = NaN(length(cfg.trl),1750);
            allsaccadespikes = NaN(length(cfg.trl),1000); %looking at eye related activity mostly for low firing rate neurons
            
            avg_waveform_amplitude = max(waveforms{1,unit})-min(waveforms{1,unit}); 
            spike_count = length(avg_waveform_amplitude);
            trial_averaged_amplitude = NaN(1,num_trials);
            trial_num = NaN(1,num_trials);
            for t = 1:num_trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITI_code); %ms sample index
                event = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code)-trial_start; %event start within trial time
                
                if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(1) % for the sequence trials
                    if  any(cfg.trl(t).allval == reward_code)
                        type(t) = 1;
                    else
                        continue;
                    end
                elseif  itmlist(cfg.trl(t).cnd-1000) <= sequence_items(2)
                    if  any(cfg.trl(t).allval == reward_code)
                        type(t) = 2;
                    else
                        continue;
                    end
                else
                    if  any(cfg.trl(t).allval == trial_start_code) %if image was shown
                        type(t) = 3;
                    else
                        continue;
                    end
                end
                
                trial_num(t) = t;
                
                %locked to the ITI
                spikes =data(unit).values{t};
                allspikesITI(t,:) = spikes(1:1750);
                
                %locked to main event
                spikes = find(spikes);
                
                waveform_ind = spikes+trial_start-1;
                
                ind = NaN(1,length(waveform_ind));
                for ii = 1:length(waveform_ind)
                    %give 1 ms quantization error buffer 1 ms for spike
                    %times + 1 ms for trial start = 2 ms
                    temp = find(waveforms{2,unit} <= waveform_ind(ii)+10 & waveforms{2,unit} >= waveform_ind(ii)-10);
                    temp = temp(1); %only take the first one in case several spikes occur within a short period of time
                    ind(ii) = temp;
                end
                if ~isempty(ind)
                    trial_averaged_amplitude(t) =  mean(avg_waveform_amplitude(ind)); %grab waveforms from trial and average       
                else
                    trial_averaged_amplitude(t) = 0;  
                end
                
                spikes = spikes-event;
                spikes(spikes < 1) = [];
                spikes(spikes > 7000) = [];
                spks = zeros(1,7000);
                spks(spikes) = 1;
                allspikes(t,:) = spks;
                
                
                if type(t) == 3 %if image type look at saccade aligned activity
                    saccadetimes = fixationstats{t}.saccadetimes;
                    invalid= find(saccadetimes(1,:) < event+500); %ignore first 500 ms
                    saccadetimes(:,invalid) = [];
                    
                    saccade_algined = zeros(1,1000);
                    
                    %going to average across all saccades mostly for low firing rate neurons
                    for s = 1:size(saccadetimes,2)
                        spks = spikes-saccadetimes(1,s);
                        spks(spks < -499) =[];
                        spks(spks > 500) = [];
                        saccade_algined(spks+500) = 1; %not an error only care when spikes not amount
                    end
                    allsaccadespikes(t,:) = saccade_algined;
                end
            end
            
            %cleanup by remove unsuccessful trial #'s
            type = laundry(type);
            allspikes = laundry(allspikes,1);
            allspikesITI = laundry(allspikesITI,1);
            allsaccadespikes = laundry(allsaccadespikes,1);
            trial_averaged_amplitude = laundry(trial_averaged_amplitude);
            trial_averaged_amplitude(trial_averaged_amplitude == 0) = NaN; 
            %replace 0's with NaNs didn't want to erase the data before
            trial_num = laundry(trial_num);
            
            %determine spike per groups of trials to determine if firing rate is approximately stable over time
            spikespertrial = bin1(allspikes',trialbinsize)'./trialbinsize;
            ITIspikespertrial = bin1(allspikesITI',trialbinsize)'./trialbinsize;
            allspikespertrial = spikespertrial+ITIspikespertrial; %not mutually exclusive but doens't matter for this
            trial_averaged_amplitude = bin1(trial_averaged_amplitude,trialbinsize,'mean');
            %have same number as bins as trial data though time of spikes
            %may be slightly off
            
            title_str = unit_stats{1,unit};
            if  multiunit(unit)
                title_str = ['Multiunit '  title_str];
            end
            title_str = [title_str ' Unit Confidence ' num2str(100*unit_stats{2,unit}) '% ' ...
                'Cutting Quality ' num2str(unit_stats{3,unit})];
            
            if spike_count < minimum_spikes %too few spikes
                listsqplot(allspikes,type,allspikesITI,allspikespertrial,spikespertrial,...
                    ITIspikespertrial,allsaccadespikes,trialbinsize,binsize,trial_averaged_amplitude,trial_num)
                subtitle(title_str);
                save_and_close_fig(too_few_folder,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                
                stability_attribute(unit) = 2;
                valid_trials(:,unit) = 0; %don't want to process further so don't analyze any trials for this unit
                continue %make raster then go on to next neuron
            elseif unit_stats{3,unit} <= poor_isolation_threshold || unit_stats{2,unit} <= confidence_threshold
                %not confident it's a unit or good isolation
                listsqplot(allspikes,type,allspikesITI,allspikespertrial,spikespertrial,...
                    ITIspikespertrial,allsaccadespikes,trialbinsize,binsize,trial_averaged_amplitude,trial_num)
                subtitle(title_str);
                save_and_close_fig(PoorIsolationFolder,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                
                stability_attribute(unit) = 4;
                valid_trials(:,unit) = 0; %don't want to process further so don't analyze any trials for this unit
                continue %make raster then go on to next neuron
            end
            
            listsqplot(allspikes,type,allspikesITI,allspikespertrial,spikespertrial,...
                ITIspikespertrial,allsaccadespikes,trialbinsize,binsize,trial_averaged_amplitude,trial_num)
            subtitle(title_str);
            
            %do manual check of which trials to take. Humans are just
            %really good at seeing pattens....
            try
                reply = input(['You can keep all trials (max trial #' num2str(num_trials) '). Is this Ok? [Y/#s]: \n' ...
                    'If NO please state [trial start and trial end].']);
            catch
                reply = input(['You can keep all trials (max trial #' num2str(num_trials) '). Is this Ok? [Y/#s]: \n' ...
                    'If NO please state [trial start and trial end].']);
            end
            if isnumeric(reply)
                [start_end] = listsqTrialCutoffplots(reply,allspikes,type,allspikesITI,allspikespertrial,...
                    spikespertrial,ITIspikespertrial,allsaccadespikes,trialbinsize,binsize,trial_averaged_amplitude,trial_num);
            else
                start_end = [NaN NaN]; %then take all trials
            end
            
            %determine if the number of trials is actually sufficient to
            %use for data analysis
            if any(~isnan(start_end)) && all(start_end ~= 0); %if all NaNs take all trials
                if isnan(start_end(2))
                    max_trial = cfg.trl(t).cnd-1000; %get the last condition number
                    min_trial = cfg.trl(start_end(1)).cnd-1000; %get the first condition number
                    if max_trial >  minimum_trials_1 %if last trial was greater than the end of the first novel/repeat block
                        if max_trial-min_trial < minimum_trials_2 %not enough trials didn't get through 1 novel/repeat block
                            start_end = [0 0];
                        end
                    else %didn't even get through first novel/repeat block
                        start_end = [0 0];
                    end
                elseif isnan(start_end(1)) %then take first trial to max_trial
                    max_trial = cfg.trl(start_end(2)).cnd-1000; %get the last condition number
                    if max_trial < minimum_trials_1 %didn't finish fam block + at least 1 novel/repeat block
                        start_end = [0 0];
                    end
                else %both trials numeric
                    min_trial = cfg.trl(start_end(1)).cnd-1000; %get the last condition number
                    max_trial = cfg.trl(start_end(2)).cnd-1000; %get the last condition number
                    if max_trial < minimum_trials_1 %didn't finish fam block + at least 1 novel/repeat block
                        start_end = [0 0];
                    elseif max_trial-min_trial < minimum_trials_2 %not enough trials didn't get through 1 novel/repeat block
                        start_end = [0 0];
                    end
                end
            end
            
            valid_trials(:,unit) = start_end';
            
            subtitle(title_str);
            
            if all(valid_trials(:,unit) == 0)%not enough trials to do data analysis on
                save_and_close_fig(too_few_folder,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                stability_attribute(unit) = 3;
            else
                if multiunit(unit) == 1 %then treat as multiunit
                    save_and_close_fig(MultiUnit_folder,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                else
                    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                end
            end
        end
    case {'cvtnew','CVTNEW'}
        
        valid_trials = NaN(2,num_units);%trials we want to keep based on firing rate
        for unit = 1:num_units
            allspikes = NaN(length(cfg.trl),7000);
            allspikesITI = NaN(length(cfg.trl),2000);
            for t = 1:num_trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITI_code);
                event = cfg.trl(t).alltim(cfg.trl(t).allval == dot_on_code)-trial_start;
                
                spikes =data(unit).values{t};
                %locked to the ITI
                allspikesITI(t,:) = spikes(1:2000);
                
                %locked to main event
                spikes = find(spikes);
                spikes = spikes-event;
                spikes(spikes < 1) = [];
                spikes(spikes > 3500) = [];
                allspikes(t,spikes) = 1;
            end
            
            
            %determine spike per groups of trials to determine if firing rate is approximately stable over time
            spikespertrial = bin1(allspikes',trialbinsize)'./trialbinsize;
            ITIspikespertrial = bin1(allspikesITI',trialbinsize)'./trialbinsize;
            allspikespertrial = spikespertrial+ITIspikespertrial; %not mutually exclusive but doens't matter for this
            
            title_str = unit_stats{1,unit};
            if  multiunit(unit)
                title_str = ['Multiunit '  title_str];
            end
            title_str = [title_str ' Unit Confidence ' num2str(100*unit_stats{2,unit}) '% ' ...
                'Cutting Quality ' num2str(unit_stats{3,unit})];
            
            cvtnewplot(allspikes,allspikesITI,allspikespertrial,spikespertrial,...
                ITIspikespertrial,trialbinsize,binsize,trial_averaged_amplitude,trial_num)
            subtitle(title_str);
            
            %do manual check of which trials to take. Humans are just
            %really good at seeing pattens....
            reply = input(['You can keep all trials (max trial #' num2str(num_trials) '). Is this Ok? [Y/#s]: \n' ...
                'If NO please state [trial start and trial end].']);
            if isnumeric(reply)
                [start_end] = listsqTrialCutoffplots(reply,allspikes,type,allspikesITI,allspikespertrial,...
                    spikespertrial,ITIspikespertrial,trialbinsize,binsize,trial_averaged_amplitude,trial_num);
            end
            valid_trials(:,unit) = start_end';
            
            subtitle(title_str);
            save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
        end
end

%add the valid trials variable to preprocess file
save([data_dir task_file(1:end-11) '-preprocessed.mat'],'-append','valid_trials','stability_attribute')
end

function listsqplot(allspikes,type,allspikesITI,allspikespertrial,...
    spikespertrial,ITIspikespertrial,allsaccadespikes,trialbinsize,binsize,trial_averaged_amplitude,trial_num)
%plot the listsq

if ~isempty(findall(0,'Type','Figure'))
    g = gcf;
    if g.Number == 101;
        close
    end
end

% Rasters from Main Event
maxtime = 5000*ones(1,3);

figure(101);

screen_size = get(0, 'ScreenSize');
set(gcf, 'Position', [0 0 screen_size(3) screen_size(4)]);
pause(0.5) %give time for plot to reach final size

for i = 1:3;
    sub_trial_num =	trial_num(find(type == i));
    subplot(3,3,i)
    [trial,time] = find(allspikes(type == i,:)== 1);
    plot(time,sub_trial_num(trial),'.k')
    ylim([0 max(sub_trial_num(trial))+5])
    ylabel('Trial #')
    if ~isempty(time)
        maxtime(i) = ceil(max(time)/1000)*1000;
    end
end

%just raster from saccade aligned activtiy. Any spike aligned to a saccade
%mostly for low firing rate neurons to determine stability for these when
%not obvious
subplot(3,3,8)
[trial,time] = find(allsaccadespikes == 1);
if ~isempty(trial)
    plot(time-500,sub_trial_num(trial),'.k')
    ylim([0 max(sub_trial_num(trial))+5])
end
ylabel('Trial #')
xlabel('Time from Saccade Start (ms)')

%PSTHs from main event
for i = 1:3;
    subplot(3,3,3+i)
    asp = bin1(allspikes(type == i,:),binsize,'lower','sum');
    bar(binsize:binsize:binsize*length(asp),asp,'k');
    box off
    xlim([0 maxtime(i)])
    ylabel('Count')
    if i < 3
        xlabel('Time from Trial Start (ms)')
    else
        xlabel('Time from Image On (ms)')
    end
    title('PSTH')
end


%Just raster from ITI start
subplot(3,3,7)
[trial,time] = find(allspikesITI == 1);
if ~isempty(trial)
    plot(time,trial_num(trial),'.k')
    ylim([0 max(trial_num(trial))+5])
end
ylabel('Trial #')
xlabel('Time from ITI start (ms)')


allspikespertrial(allspikespertrial == 0) = NaN;
subplot(3,3,9)
hold on
b = bar([spikespertrial ITIspikespertrial],'stacked');
set(b,'edgecolor','none','FaceAlpha',0.5)
if ~isempty(trial)
    xlim([floor(min(trial)/6)-1 ceil(max(trial)/6)+1]);
end
xl = xlim;
plot([xl(1) xl(2)],[nanmedian(allspikespertrial) nanmedian(allspikespertrial)],'k-','linewidth',5)
plot([xl(1) xl(2)],[nanmedian(allspikespertrial)-nanstd(allspikespertrial) nanmedian(allspikespertrial)-nanstd(allspikespertrial)],'k--')
plot([xl(1) xl(2)],[nanmedian(allspikespertrial)+nanstd(allspikespertrial) nanmedian(allspikespertrial)+nanstd(allspikespertrial)],'k--')

yl = ylim;
ylim([0 yl(2)]);

%scale waveform amplitude so sits nicely on plot
trial_averaged_amplitude = trial_averaged_amplitude-mean(trial_averaged_amplitude); %zero
trial_averaged_amplitude = trial_averaged_amplitude/(max(abs(trial_averaged_amplitude))); %scale to 1
trial_averaged_amplitude = yl(2)/3*trial_averaged_amplitude; %rescale scale
avg_waveform_ampltiude = trial_averaged_amplitude+yl(2)/2;%set level to median
plot(avg_waveform_ampltiude,'r')

hold off
xlabel(['Groups of Trials (' num2str(trialbinsize) 'trials/group)'])
ylabel('Average Spikes per trial')

end

function cvtnewplot(allspikes,allspikesITI,allspikespertrial,...
    spikespertrial,ITIspikespertrial,trialbinsize,binsize,trial_averaged_amplitude,trial_num)
%plot the cvtnew

if ~isempty(findall(0,'Type','Figure'))
    g = gcf;
    if g.Number == 101;
        close
    end
end


figure(101);

screen_size = get(0, 'ScreenSize');
set(gcf, 'Position', [0 0 screen_size(3) screen_size(4)]);
pause(0.5) %give time for plot to reach final size

% Rasters from Main Event
subplot(3,2,1)
[trial,time] = find(allspikesITI == 1);
plot(time,trial,'.k')
ylim([0 max(trial)])
ylabel('Trial #')
xlim([0 2000])

subplot(3,2,2)
[trial,time] = find(allspikes == 1);
plot(time,trial,'.k')
ylim([0 max(trial)])
ylabel('Trial #')
xlim([0 3500])


%PSTHs from main event
subplot(3,2,3)
asp = bin1(allspikesITI,binsize,'lower','sum');
bar(binsize:binsize:binsize*length(asp),asp,'k');
box off
xlim([0 2000])
ylabel('Count')
xlabel('Time from ITI Start (ms)')
title('PSTH')

subplot(3,2,4)
asp = bin1(allspikes,binsize,'lower','sum');
bar(binsize:binsize:binsize*length(asp),asp,'k');
box off
xlim([0 3500])
ylabel('Count')
xlabel('Time from Dot On (ms)')
title('PSTH')

subplot(3,2,5:6)
hold on
b = bar([spikespertrial ITIspikespertrial],'stacked');
set(b,'edgecolor','none','FaceAlpha',0.5)
xl = xlim;
plot([xl(1) xl(2)],[median(allspikespertrial) median(allspikespertrial)],'k-','linewidth',5)
plot([xl(1) xl(2)],[median(allspikespertrial)-std(allspikespertrial) median(allspikespertrial)-std(allspikespertrial)],'k--')
plot([xl(1) xl(2)],[median(allspikespertrial)+std(allspikespertrial) median(allspikespertrial)+std(allspikespertrial)],'k--')
hold off
xlabel(['Groups of Trials (' num2str(trialbinsize) 'trials/group)'])
ylabel('Average Spikes per trial')

end

function [start_end] = listsqTrialCutoffplots(reply,allspikes,type,allspikesITI,allspikespertrial,...
    spikespertrial,ITIspikespertrial,allsaccadespikes,trialbinsize,binsize,trial_averaged_amplitude,trial_num)
%plot the cutoffs from replys
while isnumeric(reply)
    listsqplot(allspikes,type,allspikesITI,allspikespertrial,spikespertrial,...
        ITIspikespertrial,allsaccadespikes,trialbinsize,binsize,trial_averaged_amplitude,trial_num)
    
    start_end = reply;
    
    %for plotting and visualization should put a line down
    nano = 0;
    if isnan(reply(1));
        reply(1) = 0;
    elseif isnan(reply(2))
        nano = 1;
    end
    for sb = [1:3 7:9];
        if sb == 9
            subplot(3,3,9)
            yl = ylim;
            hold on
            plot([reply(1)/trialbinsize reply(1)/trialbinsize],[0 yl(2)],'r');
            plot([reply(2)/trialbinsize reply(2)/trialbinsize],[0 yl(2)],'r');
            hold off
        else
            subplot(3,3,sb)
            xl = xlim;
            if nano
                yl = ylim;
                hold on
                plot([0 xl(2)],[reply(1) reply(1)],'r');
                plot([0 xl(2)],[yl(2) yl(2)],'r');
                hold off
            else
                hold on
                plot([0 xl(2)],[reply(1) reply(1)],'r');
                plot([0 xl(2)],[reply(2) reply(2)],'r');
                hold off
            end
        end
    end
    
    try
        reply = input(['You can keep these trials. Is this Ok? [Y/#s]: \n' ...
            'If NO please state [trial start and trial end].']);
    catch
        reply = input(['You can keep these trials. Is this Ok? [Y/#s]: \n' ...
            'If NO please state [trial start and trial end].']);
    end
    
end
end

function [start_end] = cvtnewTrialCutoffplots(reply,allspikes,allspikesITI,allspikespertrial,...
    spikespertrial,ITIspikespertrial,trialbinsize,binsize,trial_averaged_amplitude,trial_num)
%plot the cutoffs from replys
while isnumeric(reply)
    cvtnewplot(allspikes,allspikesITI,allspikespertrial,spikespertrial,...
        ITIspikespertrial,trialbinsize,binsize,trial_averaged_amplitude,trial_num)
    
    start_end = reply;
    
    %for plotting and visualization should put a line down
    nano = 0;
    if isnan(reply(1));
        reply(1) = 0;
    elseif isnan(reply(2))
        nano = 1;
    end
    
    subplot(3,2,5:6)
    yl = ylim;
    hold on
    plot([reply(1)/trialbinsize reply(1)/trialbinsize],[0 yl(2)],'r');
    plot([reply(2)/trialbinsize reply(2)/trialbinsize],[0 yl(2)],'r');
    hold off
    
    subplot(3,2,1)
    xl = xlim;
    if nano
        yl = ylim;
        hold on
        plot([0 xl(2)],[reply(1) reply(1)],'r');
        plot([0 xl(2)],[yl(2) yl(2)],'r');
        hold off
    else
        hold on
        plot([0 xl(2)],[reply(1) reply(1)],'r');
        plot([0 xl(2)],[reply(2) reply(2)],'r');
        hold off
    end
    
    subplot(3,2,2)
    xl = xlim;
    if nano
        yl = ylim;
        hold on
        plot([0 xl(2)],[reply(1) reply(1)],'r');
        plot([0 xl(2)],[yl(2) yl(2)],'r');
        hold off
    else
        hold on
        plot([0 xl(2)],[reply(1) reply(1)],'r');
        plot([0 xl(2)],[reply(2) reply(2)],'r');
        hold off
    end
    
    
    reply = input(['You can keep all trials. Is this Ok? [Y/#s]: \n' ...
        'If NO please state [trial start and trial end].']);
    
end
end

