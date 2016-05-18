function spatial_analysis_plotsV2(figure_dir,task_file,position,spike_times,spatial_info,...
    task_type,unit_names,smval,imageX,imageY,trial_type,Fs,peak_firing_rate)
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

tooslow_firing_folder = [figure_dir '\FiringRateTooSlow_Spatial\'];
fr_threshold = 1; %peak rate must be greater than 1 Hz to process
min_bin_dur = 0.100; %minimum of 500 ms in each bin to use so no outlier
num_units = length(position);

switch task_type
    case 'List_spatial'
        binsize = smval(1); %number of pixels per spatial bing
        filter_width = smval(2); %std of 2D gaussian smoothing filter
        filter_size = filter_width*10;
        H = fspecial('gaussian',filter_size,filter_width);
        
        for unit = 1:num_units
            if ~all(isnan(peak_firing_rate(:,unit)))
                clims = NaN(2,3);
                figure
                
                
                %Calculated smoothed firing rate for all images
                filtered_time = filter_time(position{unit},imageX,imageY,Fs,binsize,H);
                filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
                filtered_space = filter_space(position{unit},spike_times{unit},imageX,imageY,binsize,H);
                firing_rate = filtered_space./filtered_time;
                
                %draw spikes in locations above threshold (currently median firing rate)
                median_firing_rate = nanmedian(firing_rate(:))+0.6*nanstd(firing_rate(:));
                [r,c] = find(firing_rate < median_firing_rate);
                firing_ind = sub2ind(size(firing_rate),r,c);
                threshold_matrix = zeros(size(firing_rate));
                threshold_matrix(firing_ind) = 1;
                threshold_matrix = imresize(threshold_matrix,[600,800],'method','nearest');%upsample to images size in pixel coordinates
                threshold_matrix = threshold_matrix(end:-1:1,:);%flip to work with eye data
                
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
                if any(spatial_info.shuffled_spatialstability_prctile(3,unit) >= 90)
                    title(sprintf(['All images color-coded by time \n' ....
                        'r = ' num2str(spatial_info.spatialstability(3,unit),2) ' ' ...
                        num2str(spatial_info.shuffled_spatialstability_prctile(3,unit),2) '%%']))
                else
                    title(sprintf(['All images color-coded by time']))
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
                    clims(2,1) = fr(round(0.95*length(fr)));% the ~95%-tile
                end
                
                if spatial_info.shuffled_rate_prctile(3,unit) >= 90;
                    title(sprintf(['All images, peak rate = ' num2str(max(fr),3) ...
                        ' Hz \n Bits = ' num2str(spatial_info.rate(3,unit),3) ' ' ...
                        num2str(spatial_info.shuffled_rate_prctile(3,unit),3) '%%']))
                else
                    title(['All images, peak rate = ' num2str(max(fr),3) ' Hz'])
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
                    clims(2,2) = fr(round(0.95*length(fr)));% the ~95%-tile
                end
                
                
                if spatial_info.shuffled_rate_prctile(1,unit) >= 90;
                    title(sprintf(['Novel images, peak rate = ' num2str(max(fr),3) ...
                        ' Hz \n Bits = ' num2str(spatial_info.rate(1,unit),3) ' ' ...
                        num2str(spatial_info.shuffled_rate_prctile(1,unit),3) '%%']))
                else
                    title(['Novel Images, peak rate = ' num2str(max(fr),3) ' Hz'])
                end
                
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
                    clims(2,3) = fr(round(0.95*length(fr)));% the ~95%-tile
                end
                
                if spatial_info.shuffled_rate_prctile(2,unit) >= 90;
                    title(sprintf(['Repeat images, peak rate = ' num2str(max(fr),3) ...
                        ' Hz \n Bits = ' num2str(spatial_info.rate(2,unit),3) ' ' ...
                        num2str(spatial_info.shuffled_rate_prctile(2,unit),3) '%%']))
                else
                    title(['Repeat images, peak rate = ' num2str(max(fr),3) ' Hz'])
                end
                
                minc = nanmin(clims(1,:));
                maxc = nanmax(clims(2,:));
                diffc = maxc-minc;
                minc = minc - 0.15*diffc;
                for sp = 4:6
                    subplot(2,3,sp)
                    colormap('jet')
                    shading interp
                    caxis([minc maxc])
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
        
    case 'List_spatial_time_shifted'
        binsize = smval(1); %number of pixels per spatial bing
        filter_width = smval(2); %std of 2D gaussian smoothing filter
        filter_size = filter_width*10;
        H = fspecial('gaussian',filter_size,filter_width);
        
        temporal_shifts = spatial_info.temporal_shifts;
        
        %---plot all novel and repeat images---%
        for unit = 1:num_units
            clims = NaN(2,length(temporal_shifts));
            
            figure
            for ts = 1:1:length(temporal_shifts)
                
                make_spike_jittered_plot(position{unit},spike_times{ts,unit},[4,4],2*ts-1)
                title(['Time shift: ' num2str(temporal_shifts(ts)) ' ms'])
                set(gca,'Xcolor','w')
                set(gca,'Ycolor','w')
                
                subplot(4,4,2*ts)
                filtered_time = filter_time(position{unit},imageX,imageY,Fs,binsize,H);
                filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
                filtered_space = filter_space(position{unit},spike_times{ts,unit},imageX,imageY,binsize,H);
                imagesc(filtered_space./filtered_time);
                axis off
                axis equal
                clims(:,ts) = caxis;
                
                shift_str = ['TS ' num2str(temporal_shifts(ts)) ' ms '];
                if spatial_info.rate(3,unit,ts) > spatial_info.shuffled_95_percentile(3,unit);
                    bit_str = [' Bits = ' num2str(spatial_info.rate(3,unit))];
                else
                    bit_str = [];
                end
                title([shift_str 'Novel, peak = ' num2str(max(max(filtered_space./filtered_time))) 'Hz ' bit_str]);
                
            end
            
            cmin = min(clims(1,:));
            cmax = max(clims(2,:));
            diffc = cmax-cmin;
            cmin = cmin-0.15*diffc;
            for sp = 1:length(temporal_shifts)
                subplot(4,4,2*ts)
                caxis([cmin cmax])
                c = colormap;
                c(1,:) = [1 1 1];%turn nans in to white pixels
                colormap(c);
            end
            
            num_trials_str = [' n_{nov} = ' num2str(sum(trial_type{unit} == 1)) ' n_{rep} = ' num2str(sum(trial_type{unit} == 2))];
            if unit_names.multiunit(unit)
                multi_str = 'Multiunit ';
            else
                multi_str = ' ';
            end
            subtitle(['Time shifted ' num_trials_str multi_str unit_names.name{unit}]);
            
            
            save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names.name{unit} '_List_time_shift_spatial_analysis']);
        end
        
        %---Plot novel and repeats seperately---%
        climsn = NaN(2,length(temporal_shifts)); %for novel only
        climsr = NaN(2,length(temporal_shifts)); %for repeats only
        figure
        for ts = 1:1:length(temporal_shifts)
            
            subplot(4,4,2*ts-1)
            filtered_time = filter_time(select_eyepos(position{unit},trial_type{unit} ==1),imageX,imageY,Fs,binsize,H);
            filtered_space = filter_space(select_eyepos(position{unit},trial_type{unit} ==1),spike_times{ts,unit}(trial_type{unit} ==1,:),imageX,imageY,binsize,H);
            h = imagesc(filtered_space./filtered_time);
            set(h,'alphadata',~isnan(filtered_time));
            axis off
            climsn(:,ts) = caxis;
            
            shift_str = ['TS ' num2str(temporal_shifts(ts)) ' ms '];
            if spatial_info.rate(1,unit,ts) > spatial_info.shuffled_95_percentile(1,unit);
                bit_str = [' Bits = ' num2str(spatial_info.rate(1,unit))];
            else
                bit_str = [];
            end
            title([shift_str 'Novel, peak = ' num2str(max(max(filtered_space./filtered_time))) 'Hz ' bit_str]);
            
            subplot(4,4,2*ts)
            filtered_time = filter_time(select_eyepos(position{unit},trial_type{unit} ==2),imageX,imageY,Fs,binsize,H);
            filtered_space = filter_space(select_eyepos(position{unit},trial_type{unit} ==2),spike_times{ts,unit}(trial_type{unit} ==2,:),imageX,imageY,binsize,H);
            h = imagesc(filtered_space./filtered_time);
            set(h,'alphadata',~isnan(filtered_time));
            axis off
            climsr(:,ts) = caxis;
            
            
            shift_str = ['TS ' num2str(temporal_shifts(ts)) ' ms '];
            if spatial_info.rate(2,unit,ts) > spatial_info.shuffled_95_percentile(2,unit);
                bit_str = [' Bits = ' num2str(spatial_info.rate(2,unit))];
            else
                bit_str = [];
            end
            title([shift_str 'Repeat, peak = ' num2str(max(max(filtered_space./filtered_time))) 'Hz ' bit_str]);
        end
        
        cmin = min(clims(1,:));
        cmax = max(clims(2,:));
        diffc = cmax-cmin;
        cmin = cmin-0.15*diffc;
        for sp = 1:length(temporal_shifts)
            subplot(4,4,2*ts)
            caxis([cmin cmax])
            c = colormap;
            c(1,:) = [1 1 1];%turn nans in to white pixels
            colormap(c);
        end
        
        num_trials_str = [' n_ = ' num2str(size(spike_times{ts,unit},1))];
        if unit_names.multiunit(unit)
            multi_str = 'Multiunit ';
        else
            multi_str = ' ';
        end
        subtitle(['Time shifted ' num_trials_str multi_str unit_names.name{unit}]);
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names.name{unit} '_List_novelty_time_shift_spatial_analysis']);
        
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
                    cmax = fr(round(0.95*length(fr)));% the ~95%-tile
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