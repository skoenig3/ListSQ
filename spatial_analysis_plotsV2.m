function spatial_analysis_plotsV2(figure_dir,task_file,position,spike_times,spatial_info,...
    task_type,unit_names,smval,imageX,imageY,trial_type,Fs,peak_firing_rate,fr_threshold,...
    min_bin_dur)
% written by Seth Konig August, 2014
% generates and save plots for spatial spike analysis
% % updated SDK 1/11/17 to handlde new format and partial session data for
% vaild trials only. Updated CVTNEW section on 1/19/16
%
% Inputs:
%   1) figure_dir: directory where to save figures
%   2) task_file: name of the preprocessed_file
%   3) position: of attention by trial x in odd rows y in even rows
%   4) spike_times: aligned to eye by col and trial in row
%   6) spatial_info: spatial information
%   5) task_type: type of plot to create by task and/or subtask...
%       a) 'List_spatial': plot firing rate over space
%       b) 'List_spatial_time_shifted': time shift variations of 'List_spatial'
%       c) 'cvtnew_spatial': plot firing rate over space
%   6) unit_names...
%       a) unit_names.name: name of unit e.g. sig001a
%       b) unit_names.multiunit: 1 for multiunit 0 for single unit
%   7) smval: smoothing parameters of varying lengths depending on task
%   8) imageX: horizontal size of the image
%   9) imageY: vertical size of the image
%   10) trial_type: for List task only e.g. novel vs repeat
%   11) Sampling rate
%
% Outputs:
%   1) saved figures into figure_dir

figure_dir = [figure_dir 'Spatial Analysis\'];
tooslow_firing_folder = [figure_dir '\FiringRateTooSlow_Spatial\'];
num_units = length(position);


switch task_type
    case 'List_spatial'
        binsize = smval(1); %number of pixels per spatial bing
        filter_width = smval(2); %std of 2D gaussian smoothing filter
        filter_size = filter_width*10;
        H = fspecial('gaussian',filter_size,filter_width);
        
        for unit = 1:num_units
            if ~isnan(peak_firing_rate(1,unit))
                clims = NaN(2,3);
                figure
                
                
                %Calculated smoothed firing rate for all images
                filtered_time = filter_time(position{unit},imageX,imageY,Fs,binsize,H);
                filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
                filtered_space = filter_space(position{unit},spike_times{unit},imageX,imageY,binsize,H);
                firing_rate = filtered_space./filtered_time;
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %---Define Place Field---%
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                fr = sort(firing_rate(1:end));
                fr(isnan(fr)) = [];
                sil = zeros(1,4); %determines the number of clusters by comparing the ratio
                %of intercluster and intracluster distances, faster mod of silhouette
                for numclusts = 2:4
                    T = kmeans(fr',numclusts,'replicate',5);
                    [silh] = InterVSIntraDist(fr',T);
                    sil(numclusts) = mean(silh);
                end
                numclusters = find(sil == max(sil));
                T = kmeans(fr',numclusters(end),'replicate',25);
                
                region_mean_fr = zeros(1,numclusters);
                for t = 1:numclusters
                    region_mean_fr(t) = mean(fr(T == t));
                end
                [~,place_field_cluster] = max(region_mean_fr);
                min_firing_rate = min(fr(T == place_field_cluster));

                [r,c] = find(firing_rate > min_firing_rate);
                firing_ind = sub2ind(size(firing_rate),r,c);
                threshold_matrix = zeros(size(firing_rate));
                threshold_matrix(firing_ind) = 1;%want opposite to remove outside field
                threshold_matrix(isnan(firing_rate)) = NaN; %remove locations with too little data from further analysis
                if sum(sum(~isnan(threshold_matrix))) < 20
                    continue %way too little area to process
                end
                threshold_matrix = imresize(threshold_matrix,[imageY,imageX],'method','nearest');%upsample to images size in pixel coordinates
                
                make_spike_jittered_plot_threshold(position{unit},spike_times{unit},threshold_matrix,[2 3],1)
                title('Thresholded: All images')
                set(gca,'Xcolor','w')
                set(gca,'Ycolor','w')
                
                %draw raw spikes on all images all locations
                make_spike_jittered_plot(position{unit},spike_times{unit},[2 3],2)
                title('Raw: All images')
                set(gca,'Xcolor','w')
                set(gca,'Ycolor','w')
                
                %plot color coded by portion of session/by time
                make_spike_jittered_colored_plot(position{unit},spike_times{unit},[2 3],3)
                set(gca,'Xcolor','w')
                set(gca,'Ycolor','w')
                if any(spatial_info.shuffled_spatialstability_prctile(unit) >= 90)
                    title(sprintf(['All images color-coded by time \n' ....
                        'r = ' num2str(spatial_info.spatialstability(unit),2) ' ' ...
                        num2str(spatial_info.shuffled_spatialstability_prctile(unit),2) '%%']))
                else
                    title('All images color-coded by time')
                end
                
                subplot(2,3,4)
                h = imagesc(filtered_space./filtered_time);
                set(h,'alphadata',~isnan(filtered_time));
                axis off
                axis equal
                
                fr = sort(firing_rate(1:end));
                fr(isnan(fr)) = [];
                clims(:,1) = caxis;
                if length(fr) > 20
                    clims(2,1) = fr(round(0.99*length(fr)));% the ~99%-tile
                end
                
                if spatial_info.shuffled_rate_prctile(unit) >= 90;
                    title(sprintf(['All images, peak rate = ' num2str(max(fr),3) ...
                        ' Hz \n Bits = ' num2str(spatial_info.rate(unit),3) ' ' ...
                        num2str(spatial_info.shuffled_rate_prctile(unit),3) '%%']))
                else
                    title(['All images, peak rate = ' num2str(peak_firing_rate(3,unit),3) ' Hz'])
                end
                
                subplot(2,3,5)
                filtered_time = filter_time(select_eyepos(position{unit},trial_type{unit} ==1),imageX,imageY,Fs,binsize,H);
                filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
                filtered_space = filter_space(select_eyepos(position{unit},trial_type{unit} ==1),spike_times{unit}(trial_type{unit} ==1,:),imageX,imageY,binsize,H);
                h = imagesc(filtered_space./filtered_time);
                set(h,'alphadata',~isnan(filtered_time));
                axis off
                axis equal
                
                
                fr = sort(filtered_space(1:end)./filtered_time(1:end));
                fr(isnan(fr)) = [];
                clims(:,2) = caxis;
                if length(fr) > 20
                    clims(2,2) = fr(round(0.99*length(fr)));% the ~99%-tile
                end

                title(sprintf(['Novel images, peak rate = '  num2str(peak_firing_rate(1,unit),3) ' Hz']))
                
                
                subplot(2,3,6)
                filtered_time = filter_time(select_eyepos(position{unit},trial_type{unit} ==2),imageX,imageY,Fs,binsize,H);
                filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
                filtered_space = filter_space(select_eyepos(position{unit},trial_type{unit} ==2),spike_times{unit}(trial_type{unit} ==2,:),imageX,imageY,binsize,H);
                h = imagesc(filtered_space./filtered_time);
                set(h,'alphadata',~isnan(filtered_time));
                axis off
                axis equal
                
                
                fr = sort(filtered_space(1:end)./filtered_time(1:end));
                fr(isnan(fr)) = [];
                clims(:,3) = caxis;
                if length(fr) > 20
                    clims(2,3) = fr(round(0.99*length(fr)));% the ~99%-tile
                end
                
                title(sprintf(['Repeat images, peak rate = '  num2str(peak_firing_rate(2,unit),3) ' Hz']))
                
                minc = nanmin(clims(1,:));
                maxc = nanmax(clims(2,:));
                diffc = maxc-minc;
                minc = minc - 0.15*diffc;
                for sp = 4:6
                    subplot(2,3,sp)
                    colormap('jet')
                    shading interp
                    caxis([minc maxc])
                    
                    if sp == 4;
                        colorbar
                    end
                    
                    c = colormap;
                    c(1,:) = [1 1 1];%turn nans in to white pixels
                    colormap(c);
                end
                
                num_trials_str = [' n_ = ' num2str(size(spike_times{unit},1))];
                if unit_names.multiunit(unit)
                    multi_str = 'Multiunit ';
                else
                    multi_str = ' ';
                end
                
                subtitle(['Spatial Plots' num_trials_str multi_str unit_names.name{unit}]);

                if all(peak_firing_rate(:,unit) < fr_threshold) %peak rate must be greater than 1 Hz to process
                    save_and_close_fig(tooslow_firing_folder,[task_file(1:end-11) '_' unit_names.name{unit} '_List_spatial_analysis']);
                else
                    if unit_names.multiunit(unit)
                        save_and_close_fig([figure_dir '\MultiUnit\'],[task_file(1:end-11) '_' unit_names.name{unit} '_List_spatial_analysis']);
                    else
                        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names.name{unit} '_List_spatial_analysis']);
                    end
                end
            end
        end

    case 'cvtnew_spatial'
        binsize = smval(1); %number of pixels per spatial bing
        filter_width = smval(2); %std of 2D gaussian smoothing filter
        filter_size = filter_width*10;
        H = fspecial('gaussian',filter_size,filter_width);
        
        for unit = 1:num_units
            if ~isempty(position{unit})
                figure
                
                make_spike_jittered_plot(position{unit},spike_times{unit},[1 2],1)
                set(gca,'Xcolor','w')
                set(gca,'Ycolor','w')
                xlim([100 700])
                ylim([0 600])
                
                subplot(1,2,2)
                filtered_time = filter_time(position{unit},imageX,imageY,Fs,binsize,H);
                filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
                filtered_space = filter_space(position{unit},spike_times{unit},imageX,imageY,binsize,H);
                imagesc(filtered_space./filtered_time);
                
                fr = sort(filtered_space(1:end)./filtered_time(1:end));
                fr(isnan(fr)) = [];
                [cmin,cmax] = caxis;
                if length(fr) > 20
                    cmax = fr(round(0.99*length(fr)));% the ~99%-tile
                end
                
                colormap('jet')
                c = colormap;
                diffc = cmax-cmin;
                cmin = cmin-0.15*diffc;
                caxis([cmin cmax]);
                c(1,:) = [1 1 1];%turn nans in to white pixels
                colormap(c);
                
                axis off
                axis square
                colorbar
                
                title_str = ['peak rate = ' num2str(max(fr),3) ' Hz'];
                if spatial_info.shuffled_rate_prctile(unit) >= 90;
                    title_str = [title_str ' \n Bits = ' num2str(spatial_info.rate(unit),3) ' ' ...
                        num2str(spatial_info.shuffled_rate_prctile(unit),3) '%%'];
                end
                if spatial_info.shuffled_spatialstability_prctile(unit) >= 90;
                    title_str = [title_str  '\n r = ' num2str(spatial_info.spatialstability(unit),2) ' ' ...
                        num2str(spatial_info.shuffled_spatialstability_prctile(unit),2) '%%'];
                end
                title(sprintf(title_str))
                
                
                if unit_names.multiunit(unit)
                    subtitle(['Multiunit ' unit_names.name{unit} ]);
                else
                    subtitle(unit_names.name{unit});
                end
                save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names.name{unit} '_cvtnew_spatial_analysis']);
            end
        end
end
end