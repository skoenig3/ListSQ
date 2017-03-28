% %% recreate plots so can pull them into illustrator
% 
%% Just List


% %plots fixation aligned rasters so that can export to eps
%
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';

task_file = 'TO160318_3';
unit_name = 'sig001a';

load([data_dir  task_file(1:8) '-Place_Cell_Analysis.mat']);

this_unit = [];
for unit = 1:size(unit_stats,2)
    if strcmpi(unit_stats{1,unit},unit_name)
        this_unit = unit;
        break
    end
end

fix_locked_firing = list_fixation_locked_firing{this_unit};
fix_in_out = in_out{this_unit};
t = -twin1:twin2-1;

num_in = sum(fix_in_out == 1);
num_out = sum(fix_in_out == 4);
downsample = floor(num_out/num_in)/2; %many more fixations outside so downsample to show ~equal number
figure
%---Fixations in->out vs out->out---%

out_matrix = fix_locked_firing(fix_in_out == 4,:);
out_matrix = out_matrix(1:downsample:end,:);
[trial,time] = find(out_matrix == 1);
plot(time-twin1,(trial),'.b')
hold on
if ~isempty(trial)
    b4 = max(trial);
else
    b4 = 0;
end
[trial,time] = find(fix_locked_firing(fix_in_out == 1,:) == 1);
trial = trial+b4;
plot(time-twin1,(trial),'.r')
if ~isempty(trial)
    ylim([0 max(trial)])
else
    ylim([0 b4])
end
box off
plot([0 0],[0 max(trial)+1],'k--')
ylabel('Occurence')
set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
xlabel('Time from Fixation Start (ms)')
title('Fixation Aligned Rasters')
axis square

figure
hold on
[~,~,~,y_list,~] = dofill(t,fix_locked_firing(fix_in_out == 1,:),'red',1,smval);%out-> in
dofill(t,fix_locked_firing(fix_in_out == 4,:),'blue',1,smval);%out->out
%     plot(t,list_95_curve{1,unit},'k','linewidth',2);%95% confidence interval
[pks,locs] = findpeaks(y_list,'MinPeakWidth',40);
locs(pks < 0.66*max(y_list)) = [];
pks(pks < 0.66*max(y_list)) = [];
plot(locs-twin1,pks,'*k')
yl = ylim;
if yl(1) < 0
    yl(1) = 0;
    ylim(yl);
end
plot([0 0],[yl(1) yl(2)],'k--')
% gaps = findgaps(sig_ind);
% if ~isempty(gaps)
%     for g = 1:size(gaps,1)
%         gp = gaps(g,:);
%         gp(gp == 0) = [];
%         if length(gp) > 40
%             h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
%                 [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
%             uistack(h,'down')
%             set(h,'facealpha',.25,'EdgeColor','None')
%         end
%     end
% end
xlim([-twin1 twin2]);
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Firing Rate (Hz)')
legend('out->in','out->out','Location','NorthWest')
set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
% title(sprintf(['n_{out->in} = ' num2str(sum(fix_in_out == 1))...
%     ', n_{out->out} = ' num2str(sum(fix_in_out == 4))]));
title(['Peak Firing Rate of ' num2str(pks,3) ' Hz at ' num2str(locs-twin1) ' ms'])
axis square
%%
%% List vs Sequence Analysis
% 
% clar
% data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
% %data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
% 
% task_file = 'TO160318_3';
% unit_name = 'sig003a';
% 
% % load([data_dir  task_file '-preprocessed.mat'],'cfg','item_file','cnd_file');
% load([data_dir  task_file(1:8) '-Place_Cell_Analysis.mat']);
% load([data_dir  task_file(1:8) '_3-spatial_analysis_results.mat']);
% load([data_dir  task_file(1:8) '_3-preprocessed.mat'],'item_file','cnd_file');
% [~,~,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
% 
% H = define_spatial_filter(filter_width);
% 
% 
% this_unit = [];
% for unit = 1:size(unit_stats,2)
%     if strcmpi(unit_stats{1,unit},unit_name)
%         this_unit = unit;
%         break
%     end
% end
% imageX = 800;
% imageY = 600;
% 
% firing_rate_map = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all'); %get firing rate map
% maxfr = prctile(firing_rate_map(:),97.5);
% 
% figure
% %---plot firing rate map for all images---%
% h = imagesc(firing_rate_map);
% set(h,'alphadata',~isnan(firing_rate_map));
% title('All images')
% axis off
% axis equal
% colorbar
% colormap('jet')
% clim = caxis;
% caxis([clim(1) maxfr])
% 
% this_unit = [];
% for unit = 1:size(unit_stats,2)
%     if strcmpi(unit_stats{1,unit},unit_name)
%         this_unit = unit;
%         break
%     end
% end
% 
% fix_locked_firing = list_fixation_locked_firing{this_unit};
% fix_in_out = in_out{this_unit};
% t = -twin1:twin2-1;
% 
% num_in = sum(fix_in_out == 1);
% num_out = sum(fix_in_out == 4);
% downsample = floor(num_out/num_in)/2; %many more fixations outside so downsample to show ~equal number
% 
% place_field_matrix = all_place_field_matrix{unit};
% %Determine if any of the items are inside the matrix or not?
% sequence_inside = NaN(2,4);
% for c = 1:4
%     for seq = 1:2
%         yloc = imageY-sequence_locations{seq}(2,c);
%         xloc = sequence_locations{seq}(1,c);
%         
%         %---Simple version...is item in or out of field---%
%         %             if place_field_matrix(yloc,xloc) == 1 %item is in field
%         %                 sequence_inside(seq,c) = 1;
%         %             elseif place_field_matrix(yloc,xloc) == 0 %item is out of the field
%         %                 sequence_inside(seq,c) = 0;
%         %             else %coverage unknown
%         %                 sequence_inside(seq,c) = NaN;
%         %             end
%         
%         %---Less Simple Version...is item on border of field---%
%         if place_field_matrix(yloc,xloc) == 1 %item is in field
%             %then check if item is on border of field, if yes don't
%             %count. Check if field extends at least 0.5 dva, should be
%             %effected by coverage on edges items are at least 3.5 dva
%             %away from border of image. Median eye tracking error (fixation accuracy)
%             %on this task ~0.56 dva so should be ok with 0.5 and most
%             %fixations within 1 dva
%             if place_field_matrix(yloc-12,xloc) == 1 && place_field_matrix(yloc-12,xloc-12) == 1 &&...
%                     place_field_matrix(yloc-12,xloc+12) == 1 && place_field_matrix(yloc+12,xloc) == 1 && ...
%                     place_field_matrix(yloc+12,xloc-12) == 1 && place_field_matrix(yloc+12,xloc+12) == 1 && ...
%                     place_field_matrix(yloc,xloc+12) == 1 && place_field_matrix(yloc,xloc-12) == 1
%                 sequence_inside(seq,c) =1;
%             else
%                 sequence_inside(seq,c) = NaN; %don't want to use border for any category
%             end
%         else %check if item outside is also close to the border
%             if place_field_matrix(yloc-12,xloc) == 1 || place_field_matrix(yloc-12,xloc-12) == 1 ||...
%                     place_field_matrix(yloc-12,xloc+12) == 1 || place_field_matrix(yloc+12,xloc) == 1 || ...
%                     place_field_matrix(yloc+12,xloc-12) == 1 || place_field_matrix(yloc+12,xloc+12) == 1 || ...
%                     place_field_matrix(yloc,xloc+12) == 1 || place_field_matrix(yloc,xloc-12) == 1
%                 sequence_inside(seq,c) =NaN; %don't want to use border for any category
%             else
%                 sequence_inside(seq,c) = 0;
%             end
%         end
%     end
% end
% 
% figure
% hold on
% [~,~,~,y_list,~] = dofill(t,fix_locked_firing(fix_in_out == 1,:),'red',1,smval);%out-> in
% dofill(t,fix_locked_firing(fix_in_out == 4,:),'blue',1,smval);%out->out
% %     plot(t,list_95_curve{1,unit},'k','linewidth',2);%95% confidence interval
% [pks,locs] = findpeaks(y_list,'MinPeakWidth',40);
% locs(pks < 0.66*max(y_list)) = [];
% pks(pks < 0.66*max(y_list)) = [];
% plot(locs-twin1,pks,'*k')
% yl = ylim;
% if yl(1) < 0
%     yl(1) = 0;
%     ylim(yl);
% end
% plot([0 0],[yl(1) yl(2)],'k--')
% %
% xlim([-twin1 twin2]);
% hold off
% xlabel('Time from Fixation Start (ms)')
% ylabel('Firing Rate (Hz)')
% legend('out->in','out->out','Location','NorthWest')
% set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
% % title(sprintf(['n_{out->in} = ' num2str(sum(fix_in_out == 1))...
% %     ', n_{out->out} = ' num2str(sum(fix_in_out == 4))]));
% title(['Peak Firing Rate of ' num2str(pks,3) ' Hz at ' num2str(locs-twin1) ' ms'])
% axis square
% 
% %---Plot Place Field---%
% colors = 'rg';
% shapes = 'xo';
% 
% figure
% imagesc(all_place_field_matrix{unit});
% colormap('gray')
% hold on
% for c = 1:4
%     for seq = 1:2
%         if ~isnan(sequence_inside(seq,c))
%             plot(sequence_locations{seq}(1,c),imageY-sequence_locations{seq}(2,c),[colors(seq) shapes(seq)],'markersize',16)
%         end
%     end
% end
% hold off
% xlim([0 800])
% ylim([0 240])
% axis equal
% axis off
% title(sprintf(['Place Field Location \n Area = ' num2str(area(unit),2) '%%']));
% 
% 
% %---Plot Firing Rate Curves for Suquence Trials---%
% which_sequence = all_which_sequence{unit};
% seq_in_out = [];
% fixation_firing = [];
% for c = 1:4
%     for seq = 1:2
%         fixation_firing = [fixation_firing; ...
%             sequence_fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == seq,:)];
%         if  sequence_inside(seq,c) == 1;
%             seq_in_out = [ seq_in_out ones(1,sum(which_sequence(trial_nums{c,unit}) == seq))];
%         elseif isnan(sequence_inside(seq,c))
%             seq_in_out = [ seq_in_out NaN(1,sum(which_sequence(trial_nums{c,unit}) == seq))];
%         else
%             seq_in_out = [ seq_in_out zeros(1,sum(which_sequence(trial_nums{c,unit}) == seq))];
%         end
%     end
% end
% 
% figure
% hold on
% dofill(t,fixation_firing(in_out_sequence{unit} == 1,:),'red',1,smval);
% dofill(t,fixation_firing(in_out_sequence{unit} == 0,:),'blue',1,smval);
% yl = ylim;
% if yl(1) < 0
%     yl(1) = 0;
%     ylim(yl);
% end
% plot([0 0],[yl(1) yl(2)],'k--')
% xlabel('Time from Fixation Start (ms)')
% ylabel('Firing Rate (Hz)')
% legend('Items Inside','Items Outside','Location','NorthWest')
% if ~isnan(stats_across_tasks(3,unit))
%     plot(stats_across_tasks(3,unit)-twin1,stats_across_tasks(4,unit),'*k')
%     title(['Sequence Trials: peak of ' num2str(stats_across_tasks(4,unit),3) 'Hz @ ' ...
%         num2str(stats_across_tasks(3,unit)-twin1) ' ms']);
% else
%     title(['Sequence Trials: No peak, max firing of ' num2str(max(seq_in_curve),3) 'Hz']);
% end
% hold off
% 
% 
% 
% %---Plot Firing Rate Curves for Suquence vs Image Trials---%
% figure
% hold on
% dofill(t(1:twin1+twin2),list_fixation_locked_firing{unit}(in_out{unit} == 1,:),'black',1,smval);%list out-> in
% dofill(t,fixation_firing(in_out_sequence{unit} == 1,:),'green',1,smval);%sequence
% yl = ylim;
% if yl(1) < 0
%     yl(1) = 0;
%     ylim(yl);
% end
% plot([0 0],[yl(1) yl(2)],'k--')
% hold off
% xlabel('Time from Fixation Start (ms)')
% ylabel('Firing Rate (Hz)')
% legend('Images','Sequences','Location','NorthWest')
% 
% if stats_across_tasks(2,unit) > stats_across_tasks(4,unit)
%     title(['Contextual Gain: ' num2str(100*(stats_across_tasks(2,unit)/stats_across_tasks(4,unit)),3) '%'])
% else
%     title(['Contextual Gain: ' num2str(-100*(stats_across_tasks(4,unit)/stats_across_tasks(2,unit)),3) '%'])
% end

%% Saccade Direction Stuff

%plots fixation aligned rasters so that can export to eps

clar
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
%data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';

% direction tuned place cell 1
% task_file = 'TO151208_3';
% unit_name = 'sig001a';

%diretion tuned place cell 2
% task_file = 'TO160108_3';
% unit_name = 'sig003b';

%non direction tuned place cell
% task_file = 'TO160209_3';
% unit_name = 'sig001a';

%non-place cell direction tuning 1
task_file = 'TO160105';
unit_name = 'sig004b';

%non-place cell direction tuning 2
% task_file = 'TO151204_3';
% unit_name = 'sig002b';




load([data_dir  task_file(1:8) '_3-spatial_analysis_results.mat']);
load([data_dir  task_file(1:8) '-Saccade_Direction_Analysis.mat']);
load([data_dir  task_file(1:8) '-Place_Cell_Analysis.mat']);

H = define_spatial_filter(filter_width);

this_unit = [];
for unit = 1:size(unit_stats,2)
    if strcmpi(unit_stats{1,unit},unit_name)
        this_unit = unit;
        break
    end
end

imageX = 800;
imageY = 600;
bin_deg = 1; %number of degrees per bin
bin_deg2 = 45; %initial degree bins to determine time of peak direction modualtion
smval_deg = 18; %9 degrees std
smval = 30; %smoothing in time

degrees = [0:bin_deg:360]-180;
egrees = degrees(2:end);
degrees = degrees*pi/180;
% degrees = [degrees degrees(1)];

fix_aligned = list_fixation_locked_firing{unit}; %fixation aligned firing rate
sac_dirs = saccade_direction{unit}; %saccade directions organized the same way as fix_algined
fr = nandens(fix_aligned,smval,'gauss',1000,'nanflt'); %firing rate curve aligned to fixations


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Method 2: Determine time window of interest based on Time of Peak Firing Rate Modulation---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prefered since we currently assume most units show a preference for
%certain stimulus features at a particular time relative to eye
%movements. Method simply finds the maximum firing rate deviation from
%the mean when we bin the firing rate curves into arbitrary directions.
%Note: the exact bins are not important for this method just to remove
%the anti-prefered direction since this seems to often cause inhbition.
degrees2 = [0:bin_deg2:360]-180;
curves = NaN(length(degrees2),twin1+twin2);
for bin = 2:length(degrees2)
    these_directions = (sac_dirs < degrees2(bin) & sac_dirs >= degrees2(bin-1));
    curves(bin,:) = nandens(fix_aligned(these_directions,:),30,'gauss',Fs);
end
mean_of_curves = nanmean(fr);
crude_direction_curves = curves; %save for later
curves = abs(curves-mean_of_curves); %negatvie and positive work
%don't include more than 100 ms before fixation and 200 ms after, have
%yet to see a good example with prominent tuning outside this anyway
time_of_max_modulation = find(curves(:,100:400) == max(max(curves(:,100:400))));
[~,time_of_max_modulation] = ind2sub(size(curves(:,100:400)),time_of_max_modulation); %got to convert back into time index
time_of_max_modulation = time_of_max_modulation+100; %since ingored the first 100 ms before
window = time_of_max_modulation-((window_width/2)-1):time_of_max_modulation+window_width/2; %take window around peak

if (spatial_info.shuffled_rate_prctile(unit) > 95) || (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
    fix_aligned = fix_aligned(in_out{unit} == 2 | in_out{unit} == 4,:); %fixations in2in or out2out
    sac_dirs = sac_dirs(in_out{unit} == 2 | in_out{unit} == 4); %directions for fixations in2in or out2out
end


this_unit = [];
for unit = 1:size(unit_stats,2)
    if strcmpi(unit_stats{1,unit},unit_name)
        this_unit = unit;
        break
    end
end


firing_rate_map = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all'); %get firing rate map
maxfr = prctile(firing_rate_map(:),97.5);
figure
%---plot firing rate map for all images---%
h = imagesc(firing_rate_map);
set(h,'alphadata',~isnan(firing_rate_map));
title_str = ['\\rho_{1/2} = ' num2str(spatial_info.spatialstability_halves(1,unit),2)];
title(sprintf(title_str))
axis off
axis equal
colorbar
colormap('jet')
clim = caxis;
caxis([clim(1) maxfr])

bin_deg = 12;
  [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,fix_aligned,window,sac_dirs);
    binned_firing_rate_curves{1,unit} = mean_binned_firing_rate; %binned firing rates
    
    [prefered_dirs,anti_prefered_dirs,smoothed_direction_curve] = ...
    select_prefred_indeces(binned_firing_rate_curves{1,unit},degrees,sac_dirs,3);

figure
t = -twin1:twin2-1; %time for plotting
%---Fixation Aligned Acitivity out2out & in2in Plotted by Prefered and AntiPrefered Direction-All Fixations---%
[prefered_dirs,anti_prefered_dirs,smoothed_direction_curve] = ...
    select_prefred_indeces(binned_firing_rate_curves{1,unit},degrees,sac_dirs,2);
hold on
dofill(t,fix_aligned(prefered_dirs,:),'green',1,smval); %plot fixation aligned activity for prefered direction
dofill(t,fix_aligned(anti_prefered_dirs,:),'black',1,smval); %plot fixation aligned activity for anti-prefered direction
yl = ylim;
if yl(1) < 0;
    yl(1) = 0;
    ylim([0 yl(2)]);
end
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-twin1 twin2]);
set(gca,'Xtick',[-twin1 0 twin1 twin2])
xlabel('Time from Fixation Start (ms)')
ylabel('Firing Rate (Hz)')
title('All Fixation Aligned Activity')
legend('Prefered','Anti-Preferred')
axis square

figure
polarplot(degrees,[smoothed_direction_curve smoothed_direction_curve(1)],'b')
title(sprintf(['All mrl: ' num2str(mrls.all_fixations(unit),2) ' (' num2str(mrls.all_fixations_shuffled_prctile(unit),3) '%%)']))
rl = max(smoothed_direction_curve);
rlim([0 rl])




%---Fixation Aligned Raster Sorted by Saccade Direction for out2out & in2in fixations---%
figure
[~,si] = sort(sac_dirs);
fix_aligned_sorted = fix_aligned(si,:);
[trial,time] = find(fix_aligned_sorted == 1);
plot(time-twin1,trial,'.k')
xlim([-twin1 twin2])
if ~isempty(trial)
    ylim([0 max(trial)+1]);
    hold on
    h = fill([window(1) window(end) window(end) window(1) window(1)]-twin1,...
        [0 0 max(trial)+1 max(trial)+1 0],'r');
    %     set(h,'facealpha',.25,'EdgeColor','None')
    hold off
end
xlim([-twin1 twin2]);
set(gca,'Xtick',[-twin1 0 twin1 twin2])
ylabel('Ranked Sac. Dir.')
xlabel('Time from Fixation Start (ms)')
title('Fixation Aligned-Saccade Direction')
box off
axis square

%---Plots for Fixations Out2Out Only---%
%not plotting In2In since this is less interesting and coverage isn't
%great for all neurons anyway. I'm more intesteted in showing direction
%tuning out of the field
if (spatial_info.shuffled_rate_prctile(unit) > 95) || (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
    %only run if spatially modulated or possibly modulated (95% skaggs OR 95% corr 1/2)
    
       %---Out2Out Fixations---%
        fix_out_out = list_fixation_locked_firing{unit}(in_out{unit} == 4,:); %fixation aligned firing rate
        dir_out_out = saccade_direction{unit}(in_out{unit} == 4); %saccade directions
        
    
    %---Smoothed Plot of Firing Rate by Saccade Direction Out2Out ixations---%
%     binned_firing_rate_curves{3,unit}(isnan(binned_firing_rate_curves{3,unit})) = nanmean(binned_firing_rate_curves{3,unit});
% %     [prefered_dirs,anti_prefered_dirs,smoothed_direction_curve] = ...
% %         select_prefred_indeces(binned_firing_rate_curves{3,unit},degrees,dir_out_out,smval_deg);

    bin_deg = 12;
  [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,fix_out_out,window,dir_out_out);
    binned_firing_rate_curves{3,unit} = mean_binned_firing_rate; %binned firing rates
    
    [prefered_dirs,anti_prefered_dirs,smoothed_direction_curve] = ...
    select_prefred_indeces(binned_firing_rate_curves{3,unit},degrees,sac_dirs,3);


    figure
    polarplot(degrees,[smoothed_direction_curve smoothed_direction_curve(1)],'b')
    title(['Out2Out mrl: ' num2str(mrls.out2out(unit),2) ' (' num2str(mrls.out2out_shuffled_prctile(unit),3) '%)'])
    rlim([0 rl])
  
%     %---Fixation Aligned Acitivity Plotted by Prefered and AntiPrefered Direction-All Fixations---%
%     figure
%     hold on
%     dofill(t,fix_out_out(prefered_dirs,:),'green',1,smval);
%     dofill(t,fix_out_out(anti_prefered_dirs,:),'black',1,smval);
%     yl = ylim;
%     if yl(1) < 0;
%         yl(1) = 0;
%         ylim([0 yl(2)]);
%     end
%     plot([0 0],[yl(1) yl(2)],'k--')
%     hold off
%     xlabel('Time from Fixation Start (ms)')
%     ylabel('Firing Rate (Hz)')
%     title('Out2Out Fixation Aligned Activity')
%     legend('Prefered','Anti-Preferred')%,'Neutral1','Neutral2')
%     ylim([0 16])
%     axis square


end
% %% Sacade Amplitude Stuff
% 
% %plots fixation aligned rasters so that can export to eps
% 
% clar
% data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
% %data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
% 
% %non-place cell amplitude tuning 1
% task_file = 'TO151208_3';
% unit_name = 'sig003b';
% 
% bin_amplitude = 2; %dva for calculating saccade amplitude tuning
% bin_amplitude2 = 2;%dva for estimating window of interest
% max_amplitude = 16;%dva, approximately 95th percentile, don't have many large ones so remove
% 
% 
% load([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat']);
% load([data_dir task_file(1:end) '-spatial_analysis_results.mat'],...
%     'spatial_info','eyepos','spike_times','binsize','Fs','filter_width')
% load([data_dir task_file(1:8) '-Saccade_amplitude_Analysis.mat']);
% 
% H = define_spatial_filter(filter_width);
% 
% 
% unit_stats = unit_names;
% this_unit = [];
% for unit = 1:size(unit_stats,2)
%     if strcmpi(unit_stats{1,unit},unit_name)
%         this_unit = unit;
%         break
%     end
% end
% 
% imageX = 800;
% imageY = 600;
% 
% 
% sac_amps = fixation_information{unit}(:,6); %fixation aligned firing rate
% sac_amps = sac_amps/24; %convert to dva
% fix_aligned = fixation_locked_firing{unit}; %saccade amplitudes
% 
% %remove fixation within the first 500 ms of images on
% %amplitude is affected by time within image period
% time_from_image_on = fixation_information{unit}(:,4); %fixation aligned firing rate
% too_early = find(time_from_image_on < image_on_twin);
% fix_aligned(too_early,:) = [];
% sac_amps(too_early,:) = [];
% 
% window_width = 100;
% 
% %remove saccades that are too big since wont have many anyway
% too_large = find(sac_amps > max_amplitude);
% fix_aligned(too_large,:) = [];
% sac_amps(too_large)=[];
% 
% fr = nandens(fix_aligned,smval,'gauss',1000,'nanflt'); %firing rate curve aligned to fixations
% window = all_windows{unit};
% 
% 
%  firing_rates = sum(fix_aligned(:,window),2)*1000/window_width;   
% 
% 
% t = -twin1:twin2-1; %time for plotting
% 
% 
% %---Plot Firing Rate Curves for Large vs Small Amplitude Saccades---%
% small = prctile(sac_amps,25);
% large = prctile(sac_amps,75);
% 
% figure
% hold on
% dofill(t,fix_aligned(sac_amps <= small,:),'black',1,smval); %smoothed fiirng rate curve
% dofill(t,fix_aligned(sac_amps >= large,:),'green',1,smval); %smoothed fiirng rate curve
% yl = ylim;
% if yl(1) < 0
%     ylim([0 yl(2)]);
%     yl(1) = 0;
% end
% plot([0 0],[yl(1) yl(2)],'k--')
% hold off
% xlim([-twin1 twin2]);
% set(gca,'Xtick',[-twin1 0 twin1 twin2])
% xlabel('Time From Fixation Start')
% ylabel('Firing Rate')
% legend('Small','Large')
% title(['Firing Rate Curves for Small and Large Saccades']);
% yl = ylim;
% 
% 
% figure
% 
% %---Fixation Aligned Firing Rate Curve--%
% dofill(t,fix_aligned,'black',1,smval); %smoothed fiirng rate curve
% hold on
% h = fill([window(1) window(end) window(end) window(1) window(1)]-twin1,...
%     [yl(1) yl(1) yl(2) yl(2) yl(1)],'r'); %window of interest
% % set(h,'facealpha',.25,'EdgeColor','None')
% plot([0 0],[yl(1) yl(2)],'k--')
% hold off
% xlim([-twin1 twin2]);
% set(gca,'Xtick',[-twin1 0 twin1 twin2])
% xlabel('Time from Fixation Start (ms)')
% ylabel('Firing Rate (Hz)')
% title('All Fixation Aligned Activity')
% ylim(yl)
% 
% 
% %---Fixation Aligned Raster Sorted by Saccade amplitude for out2out & in2in fixations---%
% figure
% [~,si] = sort(sac_amps);
% fix_aligned_sorted = fix_aligned(si,:);
% [trial,time] = find(fix_aligned_sorted == 1);
% plot(time-twin1,trial,'.k')
% xlim([-twin1 twin2])
% if ~isempty(trial)
%     ylim([0 max(trial)+1]);
%     hold on
%     h = fill([window(1) window(end) window(end) window(1) window(1)]-twin1,...
%         [0 0 max(trial)+1 max(trial)+1 0],'r');
%     set(h,'facealpha',.25,'EdgeColor','None')
%     hold off
% end
% xlim([-twin1 twin2]);
% set(gca,'Xtick',[-twin1 0 twin1 twin2])
% ylabel('Ranked Sac. Amp.')
% xlabel('Time from Fixation Start (ms)')
% title('Fixation Aligned-Saccade amplitude')
% box off
% 
% %---Firing Rate Curve for Saccade Amplitude---%
% amps = [2:bin_amplitude:max_amplitude];
% figure
% plot(amps,binned_firing_rate_curves{unit})
% xlabel('Saccade Ampltiude')
% ylabel('Firing Rate')
% title(sprintf(['\\rho_{amp} = '  num2str(amplitude_correlations(unit),3) ' (' num2str(amplitude_correlations_percentile(unit),3) '%%)']))
% box off
% 
% %---Firing Rate Map---%
% figure
% firing_rate_map = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all'); %get firing rate map
% maxfr = prctile(firing_rate_map(:),97.5); %set 97.5 percentile of firing rate map to help visualize
% h = imagesc(firing_rate_map);
% set(h,'alphadata',~isnan(firing_rate_map));
% axis off
% axis equal
% colormap('jet')
% colorbar
% clim = caxis;
% caxis([clim(1) maxfr])
% title_str = ['All images, peak rate = ' num2str(max(firing_rate_map(:)),3) ' Hz\n'];
% if spatial_info.shuffled_rate_prctile(unit) > 95;
%     title_str = [title_str 'Bits = ' num2str(spatial_info.rate(unit),3) ...
%         '(' num2str(spatial_info.shuffled_rate_prctile(unit),3) '%%) '];
% end
% if (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
%     title_str = [title_str '\n \\rho_{1/2} = ' num2str(spatial_info.spatialstability_halves(1,unit),2) ...
%         '(' num2str(spatial_info.spatialstability_halves_prctile(1,unit),3) '%%)'];
% end
% if (spatial_info.spatialstability_even_odd_prctile(1,unit) > 95)
%     title_str = [title_str ' \\rho_{e/o} = ' num2str(spatial_info.spatialstability_even_odd(1,unit),2) ...
%         '(' num2str(spatial_info.spatialstability_even_odd_prctile(1,unit),3) '%%)'];
% end
% title(sprintf(title_str));
