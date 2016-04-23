function spatial_analysis_plots(figure_dir,preprocessed_data_file,position,spike_times,spatial_info,...
    task_type,unit_names,smval,imageX,imageY,trial_type,Fs)
% written by Seth Konig August, 2014
% generates and save plots for spatial spike analysis
%
% Inputs:
%   1) figure_dir: directory where to save figures
%   2) preprocessed_data_file: name of the preprocessed_file
%   2) time_lock_firing: matrix containing spike times locked to events by trial
%   3) position: of attention by trial x in odd rows y in even rows
%   4) spike_times: aligned to eye by col and trial in row
%   6) spatial_info: spatial information
%   5) task_type: type of plot to create by task and/or subtask...
%       a) 'List_spatial': plot firing rate over space
%       b) 'cvtnew_spatial': plot firing rate over space
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

num_units = length(spike_times);


switch task_type
    case 'List_spatial'
        binsize = smval(1); %number of pixels per spatial bing
        filter_width = smval(2); %std of 2D gaussian smoothing filter
        filter_size = filter_width*10;
        H = fspecial('gaussian',filter_size,filter_width);
        
        for unit = 1:num_units
            clims = NaN(2,3);
            figure
            
            make_spike_jittered_plot(position,spike_times{unit},[2 2],1)
            title('All images')
            set(gca,'Xcolor','w')
            set(gca,'Ycolor','w')
            
            
            subplot(2,2,2)
            filtered_time = filter_time(position,imageX,imageY,Fs,binsize,H);
            filtered_time(filtered_time < 0.01) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
            filtered_space = filter_space(position,spike_times{unit},imageX,imageY,binsize,H);
            h = imagesc(filtered_space./filtered_time);
            set(h,'alphadata',~isnan(filtered_time));
            axis off
            if spatial_info.rate(1,unit) > spatial_info.shuffled_95_percentile(1,unit);
                title(['All images, max firing rate = ' num2str(max(max(filtered_space./filtered_time))) ...
                    ' Bits = ' num2str(spatial_info.rate(1,unit))])
            else
                title(['All images, max firing rate = ' num2str(max(max(filtered_space./filtered_time)))])
            end
            clims(:,1) = caxis;
            
            subplot(2,2,3)
            filtered_time = filter_time(select_eyepos(position,trial_type == 1),imageX,imageY,Fs,binsize,H);
            filtered_time(filtered_time < 0.01) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
            filtered_space = filter_space(select_eyepos(position,trial_type == 1),spike_times{unit}(trial_type == 1,:),imageX,imageY,binsize,H);
            h = imagesc(filtered_space./filtered_time);
            set(h,'alphadata',~isnan(filtered_time));
            axis off
            if spatial_info.rate(2,unit) > spatial_info.shuffled_95_percentile(2,unit);
                title(['Novel images, max firing rate = ' num2str(max(max(filtered_space./filtered_time))) ...
                    ' Bits = ' num2str(spatial_info.rate(2,unit))])
            else
                title(['Novel images, max firing rate = ' num2str(max(max(filtered_space./filtered_time)))])
            end
            clims(:,2) = caxis;
            
            subplot(2,2,4)
            filtered_time = filter_time(select_eyepos(position,trial_type == 2),imageX,imageY,Fs,binsize,H);
            filtered_time(filtered_time < 0.01) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
            filtered_space = filter_space(select_eyepos(position,trial_type == 2),spike_times{unit}(trial_type == 2,:),imageX,imageY,binsize,H);
            h = imagesc(filtered_space./filtered_time);
            set(h,'alphadata',~isnan(filtered_time));
            axis off
            if spatial_info.rate(3,unit) > spatial_info.shuffled_95_percentile(3,unit);
                title(['Repeat images, max firing rate = ' num2str(max(max(filtered_space./filtered_time))) ...
                    ' Bits = ' num2str(spatial_info.rate(3,unit))])
            else
                title(['Repeat images, max firing rate = ' num2str(max(max(filtered_space./filtered_time)))])
            end
            clims(:,3) = caxis;
            
            minc = min(clims(1,:));
            maxc = max(clims(2,:));
            diffc = maxc-minc;
            minc = minc - 0.1*diffc; 
            for sp = 2:4
                subplot(2,2,sp)
                colormap('jet')
                caxis([minc maxc])
                c = colormap;
                c(1,:) = [1 1 1];%turn nans in to white pixels
                colormap(c);
            end
            
            if unit_names.multiunit(unit)
                subtitle(['Multiunit ' unit_names.name{unit} ]);
            else
                subtitle(unit_names.name{unit});
            end
            save_and_close_fig(figure_dir,[preprocessed_data_file(1:end-11) '_' unit_names.name{unit} '_List_spatial_analysis']);
        end
    case 'cvtnew_spatial'
        binsize = smval(1); %number of pixels per spatial bing
        filter_width = smval(2); %std of 2D gaussian smoothing filter
        filter_size = filter_width*10;
        H = fspecial('gaussian',filter_size,filter_width);
        
        for unit = 1:num_units
            figure

            make_spike_jittered_plot(position,spike_times{unit},[2 1],1)
            set(gca,'Xcolor','w')
            set(gca,'Ycolor','w')

            subplot(2,1,2)
            filtered_time = filter_time(position,imageX,imageY,Fs,binsize,H);
            filtered_time(filtered_time < 0.025) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
            filtered_space = filter_space(position,spike_times{unit},imageX,imageY,binsize,H);
            h = imagesc(filtered_space./filtered_time);
            set(h,'alphadata',~isnan(filtered_time));
            colormap('jet');
            [minc,maxc] = caxis;
            diffc = maxc-minc;
            minc = minc - 0.1*diffc;
            caxis([minc maxc]);
            c = colormap;
            c(1,:) = [1 1 1];%turn nans in to white pixels
            colormap(c);
            
            axis off
            if spatial_info.rate(1,unit) > spatial_info.shuffled_95_percentile(1,unit);
                title(['Max firing rate = ' num2str(max(max(filtered_space./filtered_time))) ...
                    ' Bits = ' num2str(spatial_info.rate(1,unit))])
            else
                title(['Max firing rate = ' num2str(max(max(filtered_space./filtered_time)))])
            end
         
            if unit_names.multiunit(unit)
                subtitle(['Multiunit ' unit_names.name{unit} ]);
            else
                subtitle(unit_names.name{unit});
            end
            save_and_close_fig(figure_dir,[preprocessed_data_file(1:end-11) '_' unit_names.name{unit} '_cvtnew_spatial_analysis']);
        end
end
end

function make_spike_jittered_plot(position,spike_times,src,subnum)
%src: subplot row and subplot col
%subnum: subplot number
jitter = 6; %1/4 dva

x = position(1:2:end,:);
y = position(2:2:end,:);

[trial,time] = find(spike_times);
spikeind = sub2ind(size(x),trial,time);

xs = x(spikeind);
ys = y(spikeind);
[xs,ys] = remove_nans(xs,ys);

if size(xs,1) == 1; %should only happen if theres 1 trial 
    xs = xs';
    ys = ys';
end

xs = xs+randi(jitter,length(xs),1);
ys = ys+randi(jitter,length(ys),1);

subplot(src(1),src(2),subnum)
hold on
plot(x',y','color',[0.8 0.8 0.8])
plot(xs,ys,'.r')
hold off
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
set(gca,'Xcolor',[1 1 1]);
set(gca,'Ycolor',[1 1 1]);

end

function [filtered_time] = filter_time(position,imageX,imageY,Fs,binsize,H)
%calculate the total time spent at any locaitons in binned pixels
spatial_time = time_per_pixel(position,imageX,imageY,Fs);
filtered_time = bin2(spatial_time,binsize,binsize);
filtered_time = imfilter(filtered_time,H,'replicate');
filtered_time = filtered_time(end:-1:1,:); %flip so rightside up
end

function [filtered_space] = filter_space(position,spike_times,imageX,imageY,binsize,H)
%caluclate total spikes over space
[firing_location] = pixel_spike_location(position,spike_times,imageX,imageY);
filtered_space = bin2(firing_location,binsize,binsize);
filtered_space = imfilter(filtered_space,H,'replicate');
filtered_space = filtered_space(end:-1:1,:); %flip so rightside up
end

function [selected] = select_eyepos(eyepos,select_rows)
%x eye position is in odd rows and y ye position is in even rows
%select rows should be a logical index
x = eyepos(1:2:end,:);
y = eyepos(2:2:end,:);
x = x(select_rows,:);
y = y(select_rows,:);
selected = NaN(sum(select_rows)*2,size(x,2));
selected(1:2:end,:) = x;
selected(2:2:end,:) = y;
end
