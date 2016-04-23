function List_Saccade_Direction_Analysis(data_dir,preprocessed_data_file,figure_dir)
% written by Seth Konig January 18, 2016 based on List_Saccade_Analysis
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


twin = 100;% how much time to take before and after saccade.
bin_deg =4; %number of degrees per bin
trial_start_code = 15;
imgon_code = 23;
imgoff_code = 24;
smval = 4;% moving average filter width
% dtermine significance
info_type = 'directional';
Fs = 1000;
numshuffs = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import data & get successful trials---%%%

load([data_dir preprocessed_data_file]);

%get important task specific information

if strcmpi(preprocessed_data_file(1:2),'TO')
    if strcmpi('ListSQ04.itm',item_set)
        cnd_file = 'ListSQ53.cnd';
    elseif strcmpi('ListSQ05.itm',item_set)
        cnd_file = 'ListSQ10.cnd';
    elseif strcmpi('ListSQ06.itm',item_set)
        cnd_file = 'ListSQ41.cnd';
    else
        cnd_file = [item_set(1:end-3) 'cnd'];
    end
    [itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_set,cnd_file);
else
    [itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_set);
    
end
[which_img,~] = get_image_numbers(cfg,itmlist,sequence_items,23);

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

sac_direction_pre = NaN(length(fixationstats),50);%direction from previous fixation to current fixation
sac_direction_post = NaN(length(fixationstats),50);%direction from current fixation to future fixation
sac_distance_pre = NaN(length(fixationstats),50);%amplitude from previous fixation to current fixation
sac_distance_post = NaN(length(fixationstats),50);%amplitude from current fixation to future fixation

saccade_start_time = NaN(length(fixationstats),50);%when did saccade to item start
fixation_start_time = NaN(length(fixationstats),50);%when did fixation on item start
for t = 1:num_trials
    fixationtimes = fixationstats{t}.fixationtimes;
    saccadetimes = fixationstats{t}.saccadetimes;
    fixations = fixationstats{t}.fixations;
    
    trial_start = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == trial_start_code);
    imgon = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == imgon_code)-trial_start;
    imgoff = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == imgoff_code)-trial_start;
    
    if fixationtimes(1,1) == 1 %session starts with fixation
        sacind = 0;
    else
        sacind = -1;
    end
    
    %find fiations and saccades that did not occur during the image period;
    %should also take care of the 1st fixation on the crosshair
    fixnum = 1:size(fixations,2);
    invalid= find(fixationtimes(1,:) < imgon);
    fixationtimes(:,invalid) = NaN;
    invalid= find(fixationtimes(1,:) > imgoff);
    fixationtimes(:,invalid) = NaN;
    invalid= find(saccadetimes(1,:) < imgon);
    saccadetimes(:,invalid) = NaN;
    invalid= find(saccadetimes(1,:) > imgoff);
    saccadetimes(:,invalid) = NaN;
    
    for s = 1:size(saccadetimes,2);
        if ~isnan(saccadetimes(1,s))
            saccade_start_time(t,s) = saccadetimes(1,s);
            
            %doing pre direction/amplitude locked to saccades upcomming movement direction
            f = s+sacind;
            if f+1 <= length(fixnum)
                sac_direction_pre(t,s) = atan2d(fixations(2,fixnum(f+1))-fixations(2,fixnum(f)),fixations(1,fixnum(f+1))-fixations(1,fixnum(f)));
                sac_distance_pre(t,s) = sqrt((fixations(2,fixnum(f+1))-fixations(2,fixnum(f))).^2+(fixations(1,fixnum(f+1))-fixations(1,fixnum(f))).^2);
            end
        end
    end
    
    for f = 1:size(fixationtimes,2);
        if ~isnan(fixationtimes(1,f))
            fixation_start_time(t,f) = fixationtimes(1,f);
            
            %doing post direction/amplitude locked to fixations previous movment direction
            sac_direction_post(t,f) = atan2d(fixations(2,fixnum(f))-fixations(2,fixnum(f-1)),fixations(1,fixnum(f))-fixations(1,fixnum(f-1)));
            sac_distance_post(t,f) = sqrt((fixations(2,fixnum(f))-fixations(2,fixnum(f-1))).^2+(fixations(1,fixnum(f))-fixations(1,fixnum(f-1))).^2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Firing Rate Locked to Eye Data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saccade_direction_pre = cell(1,num_units);
saccade_direction_post = cell(1,num_units);
saccade_amplitude_pre = cell(1,num_units);
saccade_amplitude_post = cell(1,num_units);
saccade_direction_firing_pre = cell(1,num_units);
saccade_direction_firing_post = cell(1,num_units);
saccade_amplitude_firing_pre = cell(1,num_units);
saccade_amplitude_firing_post = cell(1,num_units);

for unit = 1:num_units
    saccade_direction_pre{unit} = NaN(num_trials,size(saccade_start_time,2));
    saccade_direction_post{unit} = NaN(num_trials,size(fixation_start_time,2));
    saccade_amplitude_pre{unit} = NaN(num_trials,size(saccade_start_time,2));
    saccade_amplitude_post{unit} = NaN(num_trials,size(fixation_start_time,2));
    saccade_direction_firing_pre{unit} = NaN(num_trials,size(saccade_start_time,2));
    saccade_direction_firing_post{unit} = NaN(num_trials,size(fixation_start_time,2));
    saccade_amplitude_firing_pre{unit} = NaN(num_trials,size(saccade_start_time,2));
    saccade_amplitude_firing_post{unit} = NaN(num_trials,size(fixation_start_time,2));
end

for trial = 1:num_trials
    for unit = 1:num_units;
        spikes = find(data(unit).values{image_trials(trial)});
        for c = 1:size(saccade_start_time,2);
            
            sact = saccade_start_time(trial,c);
            
            if ~isnan(sact)
                sac_spikes = find(spikes >= sact-twin & spikes < sact);
                
                if sac_distance_pre(trial,c) < 2*24;%hard to conistently detect < 1.5, and within fovea
                    continue
                end
                
                %pre firing
                saccade_amplitude_pre{unit}(trial,c) = sac_distance_pre(trial,c)/24; %amplitude
                saccade_direction_pre{unit}(trial,c) = sac_direction_pre(trial,c);%direction
                saccade_amplitude_firing_pre{unit}(trial,c) = length(sac_spikes)*(Fs/twin); %firing rate
                saccade_direction_firing_pre{unit}(trial,c) = length(sac_spikes)*(Fs/twin);%firing rate
            end
        end
        
        for c = 1:size(fixation_start_time,2);
            fixt = fixation_start_time(trial,c);
            if ~isnan(fixt)
                fix_spikes = find(spikes > fixt & spikes <= fixt+twin);
                
                if sac_distance_post(trial,c) < 2*24 %hard to conistently detect < 1.5, and within fovea
                    continue
                end
                
                %post firing
                saccade_amplitude_post{unit}(trial,c) = sac_distance_post(trial,c)/24; %amplitude
                saccade_direction_post{unit}(trial,c) = sac_direction_post(trial,c);%direction
                saccade_amplitude_firing_post{unit}(trial,c) = length(fix_spikes)*(Fs/twin);%firing rate
                saccade_direction_firing_post{unit}(trial,c) = length(fix_spikes)*(Fs/twin);%firing rate
                
            end
        end
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Plots for Saccade Direction/Ampltiude Activity---%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
degrees = [0:bin_deg:360]-180;
amplitude = 2:1:25;
directional_firing_pre = cell(length(degrees),num_units);
directional_firing_post = cell(length(degrees),num_units);
amplitude_firing_pre = cell(length(amplitude)-1,num_units);
amplitude_firing_post = cell(length(amplitude)-1,num_units);
uniformity_pval = NaN(2,num_units);

for unit = 1:num_units
    
    for bin = 1:length(amplitude)-1
        pre_ind = find(saccade_amplitude_pre{unit} >= amplitude(bin) & saccade_amplitude_pre{unit} < amplitude(bin+1));
        amplitude_firing_pre{bin,unit} = [amplitude_firing_pre{bin,unit}; saccade_amplitude_firing_pre{unit}(pre_ind)];
        
        post_ind = find(saccade_amplitude_post{unit} >= amplitude(bin) & saccade_amplitude_post{unit} < amplitude(bin+1));
        amplitude_firing_post{bin,unit} = [amplitude_firing_post{bin,unit}; saccade_amplitude_firing_post{unit}(post_ind)];
        
    end
    
    trial_data.polar_firing = saccade_direction_pre{unit};
    trial_data.degrees = saccade_direction_firing_pre{unit};
    %uniformity_pval(1,unit) = estimated_mutual_information(trial_data,numshuffs,info_type,NaN,NaN);
    
    trial_data.polar_firing = saccade_direction_post{unit};
    trial_data.degrees = saccade_direction_firing_post{unit};
    %uniformity_pval(2,unit) = estimated_mutual_information(trial_data,numshuffs,info_type,NaN,NaN);
    
    
    
    for bin = 2:length(degrees)
        pre_ind = find(saccade_direction_pre{unit} < degrees(bin) & saccade_direction_pre{unit} >= degrees(bin-1));
        directional_firing_pre{bin,unit} = [directional_firing_pre{bin,unit}; saccade_direction_firing_pre{unit}(pre_ind)];
        
        post_ind = find(saccade_direction_post{unit} < degrees(bin) & saccade_direction_post{unit} >= degrees(bin-1));
        directional_firing_post{bin,unit} = [directional_firing_post{bin,unit}; saccade_direction_firing_post{unit}(post_ind)];
    end
end
%%
degrees = degrees(2:end);
degrees = degrees*pi/180;
degrees = [degrees degrees(1)];
polar_firing = cell(2,num_units);
%uniformity_pval = NaN(2,num_units);
for unit = 1:num_units;
    means = cellfun(@mean,directional_firing_pre(2:end,unit))';
    means =  [means(end-6:end) means means(1:7)];
    means = filtfilt(1/smval*ones(1,smval),1,means);
    means = means(8:end-7);
    
    trial_data.polar_firing = means;
    trial_data.degrees = degrees(1:end-1);
    uniformity_pval(1,unit) = estimated_mutual_information(trial_data,numshuffs,info_type,NaN,NaN);
    polar_firing{1,unit} = [means means(1)];
    
    
    means = cellfun(@mean,directional_firing_post(2:end,unit))';
    means =  [means(end-6:end) means means(1:7)];
    means = filtfilt(1/smval*ones(1,smval),1,means);
    means = means(8:end-7);
    
    trial_data.polar_firing = means;
    trial_data.degrees = degrees(1:end-1);
    uniformity_pval(2,unit) = estimated_mutual_information(trial_data,numshuffs,info_type,NaN,NaN);
    polar_firing{2,unit} = [means means(1)];
end




unit_names = cfg.channel(1:num_units);
for unit = 1:num_units
    
    plot_lim = max(max(polar_firing{1,unit}),max(polar_firing{2,unit}));
    plot_lim(plot_lim < .025) = 0.025;
    
    figure
    subplot(2,2,1)
    bar(amplitude(1:end-1),cellfun(@nanmean,amplitude_firing_pre(:,unit)));
    xlabel('Saccade Amplitude (dva)')
    ylabel('Firing Rate (Hz)')
    title('Pre-Saccade: Saccade Amplitude')
    
    subplot(2,2,2)
    bar(amplitude(1:end-1),cellfun(@nanmean,amplitude_firing_pre(:,unit)));
    xlabel('Saccade Amplitude (dva)')
    ylabel('Firing Rate (Hz)')
    title('Post-Saccade: Saccade Amplitude')
    
    subplot(2,2,3)
    polar(pi,plot_lim,'w') %invisible point so scales are the same
    hold on
    polar(degrees,polar_firing{1,unit},'b')
    hold off
    if uniformity_pval(1,unit) < 0.05
        title(['Pre-Saccade: Saccade Direction, p = ' num2str(uniformity_pval(1,unit))])
    else
        title('Pre-Saccade: Saccade Direction')
    end
    
    
    subplot(2,2,4)
    polar(pi,plot_lim,'w') %invisible point so scales are the same
    hold on
    polar(degrees,polar_firing{2,unit},'b')
    hold off
    if uniformity_pval(2,unit) < 0.05
        title(['Post-Saccade: Saccade Direction, p = ' num2str(uniformity_pval(2,unit))])
    else
        title('Post-Saccade: Saccade Direction')
    end
    
    bitstr = [];
    
    if multiunit(unit)
        subtitle(['Multiunit ' unit_names{unit} ' n =' num2str(sum(cellfun(@numel,amplitude_firing_post(:,unit)))) ' ' bitstr]);
    else
        subtitle([unit_names{unit} ' n =' num2str(sum(cellfun(@numel,amplitude_firing_post(:,unit)))) ' ' bitstr]);
    end
    save_and_close_fig(figure_dir,[preprocessed_data_file(1:end-11) '_' unit_names{unit} '_List_Direction_analysis']);
    
end

save([data_dir preprocessed_data_file(1:8) '-List_Direction_results.mat'],...
    'twin','smval','bin_deg','amplitude_firing_pre','amplitude_firing_post',...
    'directional_firing_pre','directional_firing_post','uniformity_pval');
