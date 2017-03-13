function ListSQ_Saccade_Amplitude_Analysis(data_dir,figure_dir,session_data)
%modified from ListSQ_Saccade_Direction_Analysis on 1/31/2017 SDK
%Code determines whether units are modulation by saccade ampltiude. Currently
%rasters are aligned to fixation but could be aligned to saccade, but since
%saccade duration is failry stereotype it's probably not a big deal.
%I hypothosize that there is no amplitude encoding in the hippocampus.

figure_dir = [figure_dir 'Saccade Amplitude\'];
task = 'ListSQ';
window_width = 100; %time around peak to take for firing rate, based on FWHM of place cell fixation analysis
imageX = 800; %horizontal image size
imageY = 600; %vertical image size
numshuffs = 10000; %number of times to shuffled/permutate saccade amplitude analysis

bin_amplitude = 2; %dva for calculating saccade amplitude tuning
bin_amplitude2 = 2;%dva for estimating window of interest
max_amplitude = 16;%dva, approximately 95th percentile, don't have many large ones so remove

%load important task parameters
task_file = get_task_data(session_data,task);
if isempty(task_file)
    return
end
if exist([data_dir task_file(1:8) '-Place_Cell_Analysis.mat'],'file') %want to remove later
    load([data_dir task_file(1:8) '-Place_Cell_Analysis.mat']);
    load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
        'spatial_info','eyepos','spike_times','binsize','Fs','filter_width')
    H = define_spatial_filter(filter_width);
else%no file found
    return
end

num_units = size(unit_stats,2); %number of units
t = -twin1:twin2-1; %time for plotting

%---Where to Store Summary Data---%
amplitude_correlations = NaN(1,num_units);%observed correlation between firing rate and saccade amplitude
shuffled_amplitude_correlations = cell(1,num_units); %shuffled correlations
amplitude_correlations_percentile = NaN(1,num_units); %observed correlation percentile
binned_firing_rate_curves = cell(1,num_units); %"firing rate curves"
all_windows = cell(1,num_units);%window used to calculate amplitude tuning

for unit = 1:num_units
   if isempty(list_fixation_locked_firing{unit}) %no valid trials so go to next unit
        continue
    end
    
    fix_aligned = list_fixation_locked_firing{unit}; %fixation aligned firing rate
    sac_amps = saccade_amplitude{unit}; %fixation aligned firing rate
    sac_amps = sac_amps/24; %convert to dva    
    
    %remove saccades that are too big since wont have many anyway
    too_large = find(sac_amps > max_amplitude);
    fix_aligned(too_large,:) = [];
    sac_amps(too_large)=[];
    
    fr = nandens(fix_aligned,smval,'gauss',1000,'nanflt'); %firing rate curve aligned to fixations
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine Time window of Interest---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Method 1: Determine time window of interest based on mrl---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %NOT rechecked for bugs on 1/13/2017
    %Finds the peaks and troughs in the firing rate curves and asks which
    %peak/trough has the maximum mrl. Not always the best since many
    %neurons have different shaped peak and troughs. Doesn't seem to work
    %particularly well in gener and not for neurons with wider curves.
    %
    %             [~,LOCS]= findpeaks(abs(fr),'MinPeakWidth',25);
    %             LOCS(LOCS < 49) = [];
    %             LOCS(LOCS > twin1+twin2-50) = [];
    %
    %             mrls = NaN(1,length(LOCS));
    %             degrees = [0:bin_deg:360]-180; %binned degrees
    %             for w = 1:length(LOCS)
    %                 window = LOCS(w)-49:LOCS(w)+50;
    %                 degrees = [0:bin_deg:360]-180; %binned degrees
    %                 amplitude_firing_rate = cell(1,length(degrees));
    %                 for bin = 2:length(degrees)
    %                     these_amplitudes = find(saccade_amplitude{unit} < degrees(bin) & saccade_amplitude{unit} >= degrees(bin-1));
    %                     amplitude_firing_rate{bin} = 1000*nansum(fix_aligned(these_amplitudes,window),2)/length(window);
    %                 end
    %                 degrees = degrees(2:end);
    %                 degrees = degrees*pi/180;
    %                 degrees = [degrees degrees(1)];
    %                 means = cellfun(@nanmean,amplitude_firing_rate(2:end));
    %                 mrls(w) = circ_r(degrees(1:end-1),means,bin_deg/180*pi);
    %             end
    %             w = find(mrls == max(mrls)); %best amplitude tuning window
    %             window = LOCS(w)-49:LOCS(w)+50;
    %     %             subplot(3,3,2)
    %     %             plot(mrls)
    %     %             xlabel('Bin #')
    %     %             ylabel('amplitude tunining (mrl)')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Method 2: Determine time window of interest based on mrl---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %NOT rechecked for bugs on 1/13/2017
    %Similar to Method 1, but looks at all time points with a sliding
    %window to determine time point of maximum mrl. Seems too sensitive
    %and often takes the time in which neuron is firing  the least rather
    %than the peak; could be important but not what I wanted for first pass.
    %
    %             mrls = NaN(1,50);
    %             mean_frs = NaN(1,50);
    %             for w = 10:40;
    %                 window = (step*(w-1)+1):(step*w+window_width);
    %                 degrees = [0:bin_deg:360]-180; %binned degrees
    %                 mean_frs(w) = nanmean(fr(window));
    %
    %                 amplitude_firing_rate = cell(1,length(degrees));
    %                 for bin = 2:length(degrees)
    %                     these_amplitudes = find(saccade_amplitude{unit} < degrees(bin) & saccade_amplitude{unit} >= degrees(bin-1));
    %                     amplitude_firing_rate{bin} = 1000*nansum(fix_aligned(these_amplitudes,window),2)/length(window);
    %                 end
    %                 degrees = degrees(2:end);
    %                 degrees = degrees*pi/180;
    %                 degrees = [degrees degrees(1)];
    %                 means = cellfun(@nanmean,amplitude_firing_rate(2:end));
    %
    %                 mrls(w) = circ_r(degrees(1:end-1),means,bin_deg/180*pi);
    %             end
    %             w = find(mrls == max(mrls)); %best amplitude tuning window
    %             window = ((step-1)*w+1):(step*w+window_width);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Method 2: Determine time window of interest based on Time of Peak Firing Rate Modulation---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %prefered since we currently assume most units show a preference for
    %certain stimulus features at a particular time relative to eye
    %movements. Method simply finds the maximum firing rate deviation from
    %the mean when we bin the firing rate curves into arbitrary amplitudes.
    %Note: the exact bins are not important for this method just to remove
    %the anti-prefered amplitude since this seems to often cause inhbition.
    amps = [0:bin_amplitude2:18]+2;
    curves = NaN(length(amps),twin1+twin2);
    for bin = 1:length(amps)-1
        these_amplitudes = (sac_amps < amps(bin+1) & sac_amps >= amps(bin));
        curves(bin,:) = nandens(fix_aligned(these_amplitudes,:),smval,'gauss',Fs);
    end
    mean_of_curves = nanmean(fr);
    crude_amplitude_curves = curves; %save for later
    curves = abs(curves-mean_of_curves); %negatvie and positive work
    %don't include more than 100 ms before fixation and 200 ms after, have
    %yet to see a good example with prominent tuning outside this anyway
    time_of_max_modulation = find(curves(:,100:400) == max(max(curves(:,100:400))));
    [~,time_of_max_modulation] = ind2sub(size(curves(:,100:400)),time_of_max_modulation); %got to convert back into time index
    time_of_max_modulation = time_of_max_modulation+100; %since ingored the first 100 ms before
    window = time_of_max_modulation-((window_width/2)-1):time_of_max_modulation+window_width/2; %take window around peak
    all_windows{unit} = window;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Determine if Unit is Significantly Modulated by Saccade Amplitude--%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %take eye data from fixations in2in or out2out since in2out or out2in
    %could be biased by field location creating artificial amplitude tuning
    % but only do this if spatial (both criterion) or possibly spatial (either criterion)
    if (spatial_info.shuffled_rate_prctile(unit) > 95) || (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
        fix_aligned = list_fixation_locked_firing{unit}; %fixation aligned firing rate
        sac_amps = saccade_amplitude{unit}; %fixation aligned firing rate
        
        fix_aligned = fix_aligned(in_out{unit} == 2 | in_out{unit} == 4,:); %fixations in2in or out2out
        sac_amps = sac_amps(in_out{unit} == 2 | in_out{unit} == 4); %directions for fixations in2in or out2out
        
        
        %remove saccades that are too big since wont have many anyway
        sac_amps = sac_amps/24; %convert to dva
        too_large = find(sac_amps > max_amplitude);
        fix_aligned(too_large,:) = [];
        sac_amps(too_large)=[];
    end
    
    %---Shuffle Data---%
    firing_rates = sum(fix_aligned(:,window),2)*1000/window_width;    
    [amps,binned_firing_rate_curves{unit}] = bin_firing(sac_amps,firing_rates,bin_amplitude,max_amplitude);
    shuffled_corrs = NaN(1,numshuffs); %shuffled correlations
    parfor shuff = 1:numshuffs
        ind = randperm(length(sac_amps)); %shuffle indeces
        shuffled_sac_amps = sac_amps(ind); %shuffled saccade amplitudes
        [amps,bfrc] = bin_firing(shuffled_sac_amps,firing_rates,bin_amplitude,max_amplitude);
        shuffled_corrs(shuff) = corr(bfrc(2:end)',amps(2:end)','row','pairwise','type','Spearman');
    end
    
    amplitude_correlations(unit) = corr(binned_firing_rate_curves{unit}',amps','row','pairwise','type','Spearman');%observed value
    shuffled_amplitude_correlations{unit} = shuffled_corrs; %shuffled correlations
    if  amplitude_correlations(unit) < 0
        amplitude_correlations_percentile(unit)= 100*sum(amplitude_correlations(unit) <  shuffled_amplitude_correlations{unit})/numshuffs;
    else
        amplitude_correlations_percentile(unit)= 100*sum(amplitude_correlations(unit) >  shuffled_amplitude_correlations{unit})/numshuffs;
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    %%%---Make Plot---%%%
    %%%%%%%%%%%%%%%%%%%%%
    yls = NaN(3,2);
    
    figure
    
    %---Fixation Aligned Firing Rate Curve--%
    subplot(2,3,1)
    dofill(t,fix_aligned,'black',1,smval); %smoothed fiirng rate curve
    xlim([-twin1 twin2]);
    set(gca,'Xtick',[-twin1 0 twin1 twin2])
    xlabel('Time from Fixation Start (ms)')
    ylabel('Firing Rate (Hz)')
    title('All Fixation Aligned Activity')
    yls(1,:) = ylim;
    
    %---Plot Firing Rate Curves for Large vs Small Amplitude Saccades---%
    small = prctile(sac_amps,25);
    large = prctile(sac_amps,75);
    subplot(2,3,2)
    hold on
    dofill(t,fix_aligned(sac_amps <= small,:),'black',1,smval); %smoothed fiirng rate curve
    dofill(t,fix_aligned(sac_amps >= large,:),'green',1,smval); %smoothed fiirng rate curve
    hold off
    xlim([-twin1 twin2]);
    set(gca,'Xtick',[-twin1 0 twin1 twin2])
    xlabel('Time From Fixation Start')
    ylabel('Firing Rate')
    legend('Small','Large')
    title(['Firing Rate Curves for Small and Large Saccades']);
    yls(2,:) = ylim;

    
    %---Plot Firing Rate Curves for  saccades of different amplitudes---%
    subplot(2,3,3)
    plot(t,crude_amplitude_curves')
    xlim([-twin1 twin2]);
    set(gca,'Xtick',[-twin1 0 twin1 twin2])
    xlabel('Time From Fixation Start')
    ylabel('Firing Rate')
    title(['Firing Rate Curves by Sac. Ampl. in 2 dva bins']);
    box off
    yls(3,:) = ylim;

    ymin = min(yls(:,1));
    ymin(ymin < 0) = 0;
    ymax = max(yls(:,2));
    
    for sb = 1:3
        subplot(2,3,sb)
        ylim([ymin ymax])
        
        hold on
        plot([0 0],[ymin ymax],'k--')
        
        if sb == 1
            h = fill([window(1) window(end) window(end) window(1) window(1)]-twin1,...
                [ymin ymin ymax ymax ymin],'r'); %window of interest
            set(h,'facealpha',.25,'EdgeColor','None')
        end
        hold off
    end
    
    
    %---Fixation Aligned Raster Sorted by Saccade amplitude for out2out & in2in fixations---%
    subplot(2,3,4)
    [~,si] = sort(sac_amps);
    fix_aligned_sorted = fix_aligned(si,:);
    [trial,time] = find(fix_aligned_sorted == 1);
    plot(time-twin1,trial,'.k')
    xlim([-twin1 twin2])
    if ~isempty(trial)
        ylim([0 max(trial)+1]);
        hold on
        h = fill([window(1) window(end) window(end) window(1) window(1)]-twin1,...
            [0 0 max(trial)+1 max(trial)+1 0],'r');
        set(h,'facealpha',.25,'EdgeColor','None')
        hold off
    end
    xlim([-twin1 twin2]);
    set(gca,'Xtick',[-twin1 0 twin1 twin2])
    ylabel('Ranked Sac. Amp.')
    xlabel('Time from Fixation Start (ms)')
    title('Fixation Aligned-Saccade amplitude')
    box off
    
    %---Firing Rate Curve for Saccade Amplitude---%
    [amps,binned_firing_rate_curves{unit}] = bin_firing(sac_amps,firing_rates,bin_amplitude,max_amplitude);
    subplot(2,3,5)
    plot(amps,binned_firing_rate_curves{unit})
    xlabel('Saccade Ampltiude')
    ylabel('Firing Rate')
    title(sprintf(['\\rho_{amp} = '  num2str(amplitude_correlations(unit),2) ' (' num2str(amplitude_correlations_percentile(unit),3) '%%)']))
    box off 
    
    %---Firing Rate Map---%
    subplot(2,3,6)
    firing_rate_map = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all'); %get firing rate map
    maxfr = prctile(firing_rate_map(:),97.5); %set 97.5 percentile of firing rate map to help visualize
    h = imagesc(firing_rate_map);
    set(h,'alphadata',~isnan(firing_rate_map));
    axis off
    axis equal
    colormap('jet')
    colorbar
    clim = caxis;
    caxis([clim(1) maxfr])
    title_str = ['All images, peak rate = ' num2str(max(firing_rate_map(:)),3) ' Hz\n'];
    if spatial_info.shuffled_rate_prctile(unit) > 95;
        title_str = [title_str 'Bits = ' num2str(spatial_info.rate(unit),3) ...
            '(' num2str(spatial_info.shuffled_rate_prctile(unit),3) '%%) '];
    end
    if (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
        title_str = [title_str '\n \\rho_{1/2} = ' num2str(spatial_info.spatialstability_halves(1,unit),2) ...
            '(' num2str(spatial_info.spatialstability_halves_prctile(1,unit),3) '%%)'];
    end
    if (spatial_info.spatialstability_even_odd_prctile(1,unit) > 95)
        title_str = [title_str ' \\rho_{e/o} = ' num2str(spatial_info.spatialstability_even_odd(1,unit),2) ...
            '(' num2str(spatial_info.spatialstability_even_odd_prctile(1,unit),3) '%%)'];
    end
    title(sprintf(title_str));
    
  
    subtitle([task_file(1:8) ' ' unit_stats{1,unit}]);
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Saccade_amplitude_Analysis']);
end


save([data_dir task_file(1:8) '-Saccade_amplitude_Analysis.mat'],...
    'window_width','numshuffs','bin_amplitude','bin_amplitude2','max_amplitude',...
    'amplitude_correlations','shuffled_amplitude_correlations','amplitude_correlations_percentile',...
    'binned_firing_rate_curves','smval','all_windows','unit_stats','all_windows','unit_stats')
end

function [amps,binned_firing_rates] = bin_firing(sacamps,firing_rates,bin_amplitude,max_amplitude)
%%
amps = [2:bin_amplitude:max_amplitude];
binned_firing_rates = NaN(1,length(amps));
for bin = 1:length(amps)-1
    these_amplitudes = (sacamps < amps(bin+1) & sacamps >= amps(bin));
    binned_firing_rates(bin)= mean(firing_rates(these_amplitudes));
end
end