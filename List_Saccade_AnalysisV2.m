function List_Saccade_AnalysisV2(data_dir,figure_dir,session_data)
% written by Seth Konig December, 2014
% updated by Seth Konig April 7, 2015 to include estimating temporal
% information locked to saccades and fixations on each item.
% updated SDK 1/11/17 to handlde new format and partial session data for
% vaild trials only.
%
% Function analyizes spike times correlated with eye movements in the
% List image portion of the ListSQ task.
%
% Inputs:
%   1) data_dir: directory where preprocessed data is located
%   2) figure_dir: location of where to put save figures
%   3) session_data
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-Eyemovement_Locked_List_results'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set important Analysis parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
twin = 750;% how much time to take before and after saccade.
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
trial_start_code = 15;
trial_end_code = 20;
imgon_code = 23;
imgoff_code = 24;
smval = 60 ;%gaussian 1/2 width for smoothing
numshuffs = 500; %recommend this is between 100 & 1000, for bootstrapping to
% dtermine significance
info_type = 'temporal';
Fs = 1000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Load important Session Data and Information---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task = 'ListSQ';
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
if isempty(task_file)
    warning('No ListSQ file could be found. Exiting function...')
    return
end
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials','hdr','fixationstats');
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);

num_trials = length(cfg.trl);
%NaNs are for start and end trials otherwise cut
valid_trials(1,isnan(valid_trials(1,:))) = 1;
valid_trials(2,isnan(valid_trials(2,:))) = num_trials;

%set the image duration
if str2double(task_file(3:8)) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Get successful trials---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get important task specific information
[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);

%preallocate space and parallel structure of cfg
num_trials = length(cfg.trl);
image_trials = zeros(1,num_trials);
for t = 1:num_trials %only take trials in which image was actually shown
    if any(cfg.trl(t).allval == 23) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end);
        image_trials(t) = 1;
    end
end

%remove excess data associated with non image trials
image_trials = find(image_trials);
num_trials = length(image_trials);

%only take the eye data from successful trials
fixationstats = fixationstats(image_trials);

%vivian makes about 3.5 sacs/sec and tobii makes 4.25 sacs/sec so 25 is a reasonable number
% and should have enough samples for each fixation/saccade #
fixation_locked_data = cell(1,25);
saccade_locked_data = cell(1,25);
for fixation = 1:25
    fixation_locked_data{fixation} = NaN(num_trials,1000);
    saccade_locked_data{fixation} = NaN(num_trials,1000);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Process eye data locked to Eye Movements---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saccade_start_time = NaN(length(fixationstats),25);%when did saccade to item start
fixation_start_time = NaN(length(fixationstats),25);%when did fixation on item start
for t = 1:num_trials
    fixationtimes = fixationstats{t}.fixationtimes;
    saccadetimes = fixationstats{t}.saccadetimes;
    
    trial_start = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == trial_start_code);
    trial_end = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == trial_end_code);
    imgon = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == imgon_code)-trial_start;
    imgoff = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == imgoff_code)-trial_start;
    
    
    % if monkey isn't paying attention data isn't probably
    % worth much plus have to cut off somewhere
    if imgoff-imgon > 1.5*imgdur-1
        imgoff = imgon+1.5*imgdur-1;
    end
    
    %find fiations and saccades that did not occur during the image period;
    %should also take care of the 1st fixation on the crosshair
    
    %fixation started before image turned on 
    invalid= find(fixationtimes(1,:) < imgon); 
    fixationtimes(:,invalid) = [];
    
    %fixation started after the image turned off and/or firing rate could corrupted by image turning off
    invalid= find(fixationtimes(1,:) > imgoff+twin); 
    fixationtimes(:,invalid) = [];
    
    %saccade started before image turned on 
    invalid= find(saccadetimes(1,:) < imgon); 
    saccadetimes(:,invalid) = [];
    
    %saccade started after the image turned off and/or firing rate could corrupted by image turning off
    invalid= find(saccadetimes(1,:) > imgoff+twin); 
    saccadetimes(:,invalid) = [];
    
    for s = 1:size(saccadetimes,2);
        saccade_start_time(t,s) = saccadetimes(1,s);
    end
    for f = 1:size(fixationtimes,2);
        fixation_start_time(t,f) = fixationtimes(1,f);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Firing Rate Locked to Eye Movements---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saccade_locked_firing = cell(25,num_units);
fixation_locked_firing = cell(25,num_units);
for c = 1:25
    for unit = 1:num_units
        saccade_locked_firing{c,unit} = NaN(num_trials,twin*2);
        fixation_locked_firing{c,unit} = NaN(num_trials,twin*2);
    end
end
for trial = 1:num_trials
    for unit = 1:num_units;
        if image_trials(trial) >= valid_trials(1,unit) && image_trials(trial) <= valid_trials(2,unit)
            spikes = find(data(unit).values{image_trials(trial)});
            for c = 1:25;
                fixt = fixation_start_time(trial,c);
                sact = saccade_start_time(trial,c);
                if ~isnan(fixt)
                    fix_spikes = spikes(spikes > fixt-twin & spikes <= fixt+twin)-fixt+twin;
                    temp = zeros(1,twin*2);
                    temp(fix_spikes) = 1;
                    fixation_locked_firing{c,unit}(trial,:) = temp;
                end
                if ~isnan(sact)
                    sac_spikes = spikes(spikes > sact-twin & spikes <= sact+twin)-sact+twin;
                    temp = zeros(1,twin*2);
                    temp(sac_spikes) = 1;
                    saccade_locked_firing{c,unit}(trial,:) = temp;
                end
            end
        end
    end
end

%combine data across all fixation/saccade numbers
all_saccade_locked_firing = cell(1,num_units);
all_fixation_locked_firing = cell(1,num_units);
for unit = 1:num_units
    for c = 5:25; %ignore the first 4 ordinal saccades as these might be 
        all_saccade_locked_firing{unit} = [all_saccade_locked_firing{unit}; saccade_locked_firing{c,unit}];
        all_fixation_locked_firing{unit} = [all_fixation_locked_firing{unit}; fixation_locked_firing{c,unit}];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Mutual Info for Eye Movements and Firing Rate---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fixation_info = NaN(1,num_units);
fixation_info_90 = NaN(1,num_units);
fixation_info_95 = NaN(1,num_units);
fixation_shuffled_info = cell(1,num_units);
saccade_info = NaN(1,num_units);
saccade_info_90 = NaN(1,num_units);
saccade_info_95 = NaN(1,num_units);
saccade_shuffled_info = cell(1,num_units);
for unit = 1:num_units
    [fixation_info(unit),fixation_shuffled_info{unit}]...
        = estimated_mutual_information(all_fixation_locked_firing{unit},numshuffs,info_type,smval,Fs);
    [saccade_info(unit),saccade_shuffled_info{unit}]...
        = estimated_mutual_information(all_saccade_locked_firing{unit},numshuffs,info_type,smval,Fs);
    fixation_info_90(unit) = prctile(fixation_shuffled_info{unit},90);
    fixation_info_95(unit) = prctile(fixation_shuffled_info{unit},95);
    saccade_info_90(unit) = prctile(saccade_shuffled_info{unit},90);
    saccade_info_95(unit) = prctile(saccade_shuffled_info{unit},95);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot Firing Rates Locked to Eye Movements---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%only plot for fixations 1-16 though 
t = -twin:twin-1;
unit_names = unit_stats(1,:);
for unit = 1:num_units
    figure
    ylims = NaN(2,16);
    for c = 1:16
        if ~all(all(isnan(fixation_locked_firing{c,unit}))) && ~all(all(isnan(saccade_locked_firing{c,unit})))
            subplot(4,4,c)
            hold on
            [~,~,~,y1] = dofill(t,fixation_locked_firing{c,unit},'blue',1,smval);
            [~,~,~,y2] = dofill(t,saccade_locked_firing{c,unit},'red',1,smval);
            ylims(1,c) = 0.8*min(min(y1),min(y2));
            ylims(2,c) = 1.2*max(max(y1),max(y2));
            hold off
            xlabel(['Time from Eye Movement ' num2str(c) ' (ms)'])
            ylabel('Firing Rate (Hz)')
        end
    end
    
    ymin = min(ylims(1,:));
    ymax = max(ylims(2,:));
    for c = 1:16;
        subplot(4,4,c)
        ylim([ymin ymax]);
        hold on
        plot([0 0],[ymin ymax],'k')
        hold off
    end
    
    if multiunit(unit)
        subtitle(['Saccade-Locked Multiunit ' unit_names{unit} ' n =' num2str(size(fixation_locked_firing{1,unit},1))]);
    else
        subtitle(['Saccade-Locked ' unit_names{unit} ' n ='  num2str(size(fixation_locked_firing{1,unit},1))]);
    end
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_List_Eye_Locked_analysis']);
    
    figure
    hold on
    dofill(t,all_fixation_locked_firing{unit},'blue',1,smval);
    dofill(t,all_saccade_locked_firing{unit},'red',1,smval);
    plot([0 0],[ylim],'k')
    hold off
    xlabel('Time from Eye Movement (ms)')
    ylabel('Firing Rate (Hz)')
    legend('All fixations','All Saccades')
    
    bitstr = [];
    if fixation_info(unit) > fixation_info_95(unit)
        bitstr = ['fix_{95} = ' num2str(fixation_info(unit))];
    elseif fixation_info(unit) > fixation_info_90(unit)
        bitstr = ['fix_{90} = ' num2str(fixation_info(unit))];
    end
    if saccade_info(unit) > saccade_info_95(unit)
        bitstr = [bitstr ' sac_{95} = ' num2str(saccade_info(unit))];
    elseif saccade_info(unit) > saccade_info_90(unit)
        bitstr = [bitstr ' sac_{90} = ' num2str(saccade_info(unit))];
    end
    
    num_trial_str = [' n_{fix} =' num2str(sum(~isnan(all_fixation_locked_firing{unit}(:,1)))) ...
        ' n_{sac} =' num2str(sum(~isnan(all_saccade_locked_firing{unit}(:,1)))) ' '];
    if multiunit(unit)
        title(['Saccade-Locked Multiunit ' unit_names{unit} num_trial_str bitstr]);
    else
        title(['Saccade-Locked ' unit_names{unit} num_trial_str ' ' bitstr]);
    end
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_List_All_Eye_Locked_analysis']);
end

save([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat'],...
    'twin','smval','all_saccade_locked_firing','all_fixation_locked_firing',...
    'fixation_info','fixation_info_90','fixation_info_95','fixation_shuffled_info',...
    'saccade_info','saccade_info_90','saccade_info_95','saccade_shuffled_info');
end