%plots show fixation (motor) vs event (visual) aligned rasters

%% For Raster Aligned to CrossHair on Image Trials---%

% data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
% %data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
%
% task_file = 'TO151208';
% unit_name = 'sig001a';
%
% load([data_dir  task_file(1:8) '-ListSQ-Visual_Response_results.mat']);
% load([data_dir  task_file(1:8) '_3-preprocessed.mat']);
% task = 'ListSQ';
%
% %grab unit data
% [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
%     multiunit,unit_confidence,sorting_quality);
% clear unit_names
%
% if num_units == 0
%     return %if no units exit function
% end
% num_trials = length(cfg.trl); %number of trials
%
% %get trials in which spikes are valid i.e. trials for which neuron is
% %stable for
% min_blks = 2;
% valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
%
% %get important task specific information
% [itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
% [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,23);
%
%
% this_unit = [];
% for unit = 1:size(unit_stats,2)
%     if strcmpi(unit_stats{1,unit},unit_name)
%         this_unit = unit;
%         break
%     end
% end
% unit = this_unit;
%
%
% fix_rt = cell(1,num_units); %how long does it take on average to fixate cross once displayed
% for unit = 1:num_units
%     fix_rt{unit}=  NaN(1,192); %how long does it take on average to fixate cross once displayed
% end
%
% event_codes = [35 8 23 24 23];
% for t = 1:num_trials
%     if any(cfg.trl(t).allval == 23); %in which image was displayed or 1st item in sequence was displayed
%         if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(end) %then sequence trial
%             continue % go to the next trial
%         end
%
%         trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == 15); %start of ITI
%         trial_end = cfg.trl(t).alltim(end)-trial_start; %end of trial after image off
%         crosson = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(1))-trial_start; %cross hair on
%         fix_cross = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(2))-trial_start; %fixation on cross hair according to cortex
%         fix_cross = fix_cross(1);%since fixation on image also counts as event 8
%         imgon = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(3))-trial_start; %when image turns on
%
%         for unit = 1:num_units
%             if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only put data in for valid trials
%
%                 %---image info---%
%                 img_index = find(cfg.trl(t).cnd == img_cnd); %image index
%                 if any(isnan(which_img(img_index)))
%                     continue
%                 end
%                 fix_rt{unit}(img_index) = fix_cross-crosson; %reaction time to fixate on cross
%             end
%         end
%     end
% end
% fix_rt = laundry(fix_rt);
%
%
% unit = this_unit;
% t12 = -twin1:twin2-1;
% t13 = -twin1:twin3-1;
% t133 = -twin1:twin3-1+250;
% t3 = -twin3:twin3-1;
% t14 = -twin1:twin4-1;
% unit_names = unit_stats(1,:);
%
% figure
%
%
% %---Plot Firing rate curve locked to Cross and Fixation on Cross---%
% subplot(2,2,1)
% hold on
% [~,~,~,y] =dofill(t13-round(nanmean(fix_rt{unit})),time_lock_firing{unit,1},'black',1,smval);%smoothed cross onset curve
% enough_samples = find(sum(~isnan(time_lock_firing{unit,2})) > 40); %variable event length so cut time points without enough samples
% [~,~,~, y1,~] =dofill2(enough_samples-twin1,time_lock_firing{unit,2}(:,enough_samples),'green',1,smval); %smoothed fixation on cross curve
% plot([-twin3 twin3+250],[baseline_firing_rate(1,unit) baseline_firing_rate(1,unit)],'k--') %pre-image "baseline" line
% hold off
% xlim([-twin3 twin3+250])
% set(gca,'Xtick',[-500 -250 0 250 500 750])
% ylabel('Firing Rate (Hz)')
% xlabel('Time from Fixation on Cross (ms)')
% legend('Cross On','Fix on Cross','"Baseline"')
% title('Firing Rate Curves Aligned by average Reaction Time')
% box off
%
% %---Plot Raster Aligned to Fixation Sorted by Reaction Time---%
% subplot(2,2,2)
% [~,si] = sort(fix_rt{unit});
% [trial,time] = find(time_lock_firing{unit,2}(si,:) == 1);
%  plot(time-twin1,trial,'.k')
% xlim([-250 750])
% set(gca,'Xtick',[-250 0 250 500 750])
% ylabel('Ranked Trial #')
% xlabel('Time from Fixation Start (ms)')
% box off
% title('Reaction Time Ranked Fixation Raster')
%
%
% %---Plot Raster aligned to fixation on cross sorted by  Reaction Time---%
% subplot(2,2,3)
% [trial,time] = find(time_lock_firing{unit,1}(si,:) == 1);
% plot(time-250,trial,'.k')
% ylim([0 max(trial)+1]);
% xlim([-250  500])
% set(gca,'Xtick',[-250 0 250 500])
% ylabel('Ranked Trial #')
% xlabel('Time from Fixation on Cross (ms)')
% box off
% title('Reaction Time Ranked Cross On Raster')
%
% %---Plot Raster aligned to fixation on cross sorted by  Reaction Time---%
% subplot(2,2,4)
% [trial,time] = find(time_lock_firing{unit,1} == 1);
% plot(time-250,trial,'.k')
% ylim([0 max(trial)+1]);
% xlim([-250  500])
% set(gca,'Xtick',[-250 0 250 500])
% ylabel('Trial #')
% xlabel('Time from Fixation on Cross (ms)')
% box off
% title('Trial Ranked Cross On Raster')
%
% subtitle(['Fixation on Cross: ' task_file(1:8) ' ' unit_names{unit} n_str]);


%% For Sequence Task aligned to fixation on shapes/shape appearance

clar
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';

task_file = 'PW141013';
unit_name = 'sig001c';


load([data_dir  task_file(1:8) '-Fixation_locked_Sequence_results.mat']);
load([data_dir  task_file(1:8) '_3-preprocessed.mat']);
task = 'ListSQ';

%grab unit data
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

if num_units == 0
    return %if no units exit function
end
num_trials = length(cfg.trl); %number of trials

%get trials in which spikes are valid i.e. trials for which neuron is
%stable for
min_blks = 2;
valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);

%get important task specific information
[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,23);


this_unit = [];
for unit = 1:size(unit_stats,2)
    if strcmpi(unit_stats{1,unit},unit_name)
        this_unit = unit;
        break
    end
end
unit = this_unit;
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---All Fixation Aligned Activity---%
fixation_firing = fixation_locked_firing{unit};
event_firing = event_locked_firing{unit};

ff = nandens(fixation_firing,50,'gauss',1000);

yls = NaN(2,2);

t = -twin+1:twin;
figure

subplot(2,2,1)
dofill(t,fixation_firing,'black',1,smval);
yls(1,:) = ylim;
xlabel('Time from Fixation Start (ms)')
ylabel('Firing Rate (Hz)')
xlim([-twin twin])
title('Fixation Aligned')


[~,rti] = sort(reaction_time{unit});
subplot(2,2,3)
[trial,time] = find(fixation_firing(rti,:) == 1);
plot(time-twin,(trial),'.k')
ylim([0 max(trial)])
ylabel('Occurence #')
xlabel('Time from fixation Start (ms)')
title('Fixation Aligned')
xlim([-twin twin])
ylim([0 max(trial)+1])
box off

subplot(2,2,2)
hold on
dofill(t,event_firing,'blue',1,smval);
yls(2,:) = ylim;
xlabel('Time from Item On (ms)')
ylabel('Firing Rate (Hz)')
xlim([-twin twin])
title('Event Aligned')

subplot(2,2,4)
[trial,time] = find(event_firing(rti,:) == 1);
plot(time-twin,(trial),'.k')
if ~isempty(trial)
    ylim([0 max(trial)+1])
end
ylabel('Occurence #')
xlabel('Time from Item On (ms)')
title('Event Aligned')
xlim([-twin twin])
ylim([0 max(trial)+1])
box off

ymin = min(yls(:,1));
ymin(ymin < 0) = 0;
ymax = max(yls(:,2));

subplot(2,2,1)
ylim([ymin ymax])
hold on
plot([0 0],[ymin ymax],'k--')
hold off

subplot(2,2,2)
ylim([ymin ymax])
hold on
plot([0 0],[ymin ymax],'k--')
hold off


subtitle([task_file(1:8) ' ' unit_stats{1,unit}]);

