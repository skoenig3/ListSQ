function List_Saccade_Analysis(data_dir,preprocessed_data_file,figure_dir)
% written by Seth Konig December, 2014
% updated by Seth Konig April 7, 2015 to include estimating temporal
% information locked to saccades and fixations on each item. 
%
% Function analyizes spike times correlated with eye movements in the
% List image portion of the task.
%
% Inputs:
%   1) data_dir: direction where preprocessed data is located
%   2) preprocessed_data_file: preprocessed ListSQ data file with LFPs
%   divided by channel and trial.
%   3) figure_dir: location of where to put save figures
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-Eyemovement_Locked_List_results'


twin = 500;% how much time to take before and after saccade.
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
trial_start_code = 15;
imgon_code = 23;
imgoff_code = 24;
smval = 60 ;%gaussian 1/2 width for smoothing
numshuffs = 100; %recommend this is between 100 & 1000, for bootstrapping to
% dtermine significance
info_type = 'temporal';
Fs = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import data & get successful trials---%%%

load([data_dir preprocessed_data_file]);

%get important task specific information
[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
[which_img,novel_vs_repeat] = get_image_numbers(cfg,itmlist,sequence_items,23);

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

saccade_start_time = NaN(length(fixationstats),25);%when did saccade to item start
fixation_start_time = NaN(length(fixationstats),25);%when did fixation on item start
for t = 1:num_trials
    fixationtimes = fixationstats{t}.fixationtimes;
    saccadetimes = fixationstats{t}.saccadetimes;
    
    trial_start = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == trial_start_code);
    imgon = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == imgon_code)-trial_start;
    imgoff = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == imgoff_code)-trial_start;
    
    %find fiations and saccades that did not occur during the image period;
    %should also take care of the 1st fixation on the crosshair
    invalid= find(fixationtimes(1,:) < imgon);
    fixationtimes(:,invalid) = [];
    invalid= find(fixationtimes(1,:) > imgoff);
    fixationtimes(:,invalid) = [];
    invalid= find(saccadetimes(1,:) < imgon);
    saccadetimes(:,invalid) = [];
    invalid= find(saccadetimes(1,:) > imgoff);
    saccadetimes(:,invalid) = [];
    
    for s = 1:size(saccadetimes,2);
        saccade_start_time(t,s) = saccadetimes(1,s);
    end
    for f = 1:size(fixationtimes,2);
        fixation_start_time(t,f) = fixationtimes(1,f);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Firing Rate Locked to Eye Data---%%%
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

all_saccade_locked_firing = cell(1,num_units);
all_fixation_locked_firing = cell(1,num_units);
for unit = 1:num_units
    for c= 1:25;
        all_saccade_locked_firing{unit} = [all_saccade_locked_firing{unit}; saccade_locked_firing{c,unit}];
        all_fixation_locked_firing{unit} = [all_fixation_locked_firing{unit}; fixation_locked_firing{c,unit}];
    end
end

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


t = -twin:twin-1;
unit_names = cfg.channel(1:num_units);
for unit = 1:num_units
    figure
    ylims = NaN(1,25);
    for c = 1:25
        if ~all(all(isnan(fixation_locked_firing{c,unit}))) && ~all(all(isnan(saccade_locked_firing{c,unit})))
            subplot(5,5,c)
            hold on
            dofill(t,fixation_locked_firing{c,unit},'blue',1,smval);
            dofill(t,saccade_locked_firing{c,unit},'red',1,smval);
            hold off
            xlabel('Time from Eye Movement (ms)')
            ylabel('Firing Rate (Hz)')
            yl = ylim;
            ylims(c) = yl(2);
        end
    end
    for c = 1:25;
        subplot(5,5,c)
        ylim([0 max(ylims)]);
        hold on
        plot([0 0],[0 max(ylims)],'k')
        hold off
    end
    if multiunit(unit)
        subtitle(['Saccade-Locked Multiunit ' unit_names{unit} ' n =' num2str(num_trials)]);
    else
        subtitle(['Saccade-Locked ' unit_names{unit} ' n =' num2str(num_trials)]);
    end
    save_and_close_fig(figure_dir,[preprocessed_data_file(1:end-11) '_' unit_names{unit} '_List_Eye_Locked_analysis']);
    
    figure
    hold on
    dofill(t,all_fixation_locked_firing{unit},'blue',1,smval);
    dofill(t,all_saccade_locked_firing{unit},'red',1,smval);
    yl = ylim;
    plot([0 0],[0 yl(2)],'k')
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
        bitstr = [bitstr 'sac_{95} = ' num2str(saccade_info(unit))];
    elseif saccade_info(unit) > saccade_info_90(unit)
        bitstr = [bitstr 'sac_{90} = ' num2str(saccade_info(unit))];
    end
    if multiunit(unit)
        title(['Saccade-Locked Multiunit ' unit_names{unit} ' n =' num2str(num_trials) ' ' bitstr]);
    else
        title(['Saccade-Locked ' unit_names{unit} ' n =' num2str(num_trials) ' ' bitstr]);
    end
    save_and_close_fig(figure_dir,[preprocessed_data_file(1:end-11) '_' unit_names{unit} '_List_All_Eye_Locked_analysis']);
end

save([data_dir preprocessed_data_file(1:8) '-Eyemovement_Locked_List_results.mat'],...
    'twin','smval','all_saccade_locked_firing','all_fixation_locked_firing',...
    'fixation_info','fixation_info_90','fixation_info_95','fixation_shuffled_info',...
    'saccade_info','saccade_info_90','saccade_info_95','saccade_shuffled_info');
end