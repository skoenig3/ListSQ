function ListSQ_Saccade_Direction_Analysis(data_dir,figure_dir,session_data)
%written by Seth Konig 11/28/16
%code determines whether units are modulation by saccade direction. Currently
%rasters are aligned to fixation but could be aligned to saccade, but since
%saccade duration is failry stereotype it's probably not a big deal. For units
%that are spatially modulated (95% corr 1/2 &  95% skaggs) or possibly modulated
%(95% skaggs OR 95% corr 1/2) the data is split into high firing and low firing
% regions. Code calculates significance with circ_rtest and by
% bootstrapping (permutation) of mrl. mrl seems more accurate in the sense
% that it taskes into account sample size better and should be more
% sensitive to low firing rate neurons.
%
% code rechecked for bugs on 1/13/2017 by SDK

figure_dir = [figure_dir 'Saccade Direction\'];
task = 'ListSQ';
window_width = 100; %time around peak to take for firing rate, based on FWHM of place cell fixation analysis
imageX = 800; %horizontal image size
imageY = 600; %vertical image size
numshuffs = 1000; %number of times to shuffled/permutate saccade direction analysis

bin_deg = 1; %number of degrees per bin
bin_deg2 = 45; %initial degree bins to determine time of peak direction modualtion
smval_deg = 18; %9 degrees std

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
mrls.all_fixations = NaN(1,num_units); %MRLs all fixations ignoring out2in and in2out
mrls.all_fixations_shuffled = cell(1,num_units); %shuffled MRLs all fixations ignoring out2in and in2out
mrls.all_fixations_shuffled_prctile = NaN(1,num_units); %observed MRL percentile
mrls.in2in = NaN(1,num_units);  %MRLs for in2in fixations ignoring all others
mrls.in2in_shuffled = cell(1,num_units);  %shuffled MRLs for in2in fixations ignoring all others
mrls.in2in_shuffled_prctile = NaN(1,num_units);  %observed MRLs percentile for in2in fixations ignoring all others
mrls.out2out = NaN(1,num_units);  %MRLs for out2out fixations ignoring all others
mrls.out2out_shuffled = cell(1,num_units); %shuffled MRLs for out2out fixations ignoring all others
mrls.out2out_shuffled_prctile = NaN(1,num_units);%MRLs percentile for out2out fixations ignoring all others

uniformity_pvalue = NaN(3,num_units); %circ_rtest...row 1: all fixations, row 2: in2in, row3: out2out
binned_firing_rate_curves = cell(3,num_units); %row 1: all fixations, row 2: in2in, row3: out2out

for unit = 1:num_units
    if isempty(list_fixation_locked_firing{unit}) %no valid trials so go to next unit
        continue
    end
    
    fix_aligned = list_fixation_locked_firing{unit}; %fixation aligned firing rate
    sac_dirs = saccade_direction{unit}; %saccade directions organized the same way as fix_algined
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
    %                 direction_firing_rate = cell(1,length(degrees));
    %                 for bin = 2:length(degrees)
    %                     these_directions = find(saccade_direction{unit} < degrees(bin) & saccade_direction{unit} >= degrees(bin-1));
    %                     direction_firing_rate{bin} = 1000*nansum(fix_aligned(these_directions,window),2)/length(window);
    %                 end
    %                 degrees = degrees(2:end);
    %                 degrees = degrees*pi/180;
    %                 degrees = [degrees degrees(1)];
    %                 means = cellfun(@nanmean,direction_firing_rate(2:end));
    %                 mrls(w) = circ_r(degrees(1:end-1),means,bin_deg/180*pi);
    %             end
    %             w = find(mrls == max(mrls)); %best direction tuning window
    %             window = LOCS(w)-49:LOCS(w)+50;
    %     %             subplot(3,3,2)
    %     %             plot(mrls)
    %     %             xlabel('Bin #')
    %     %             ylabel('Direction tunining (mrl)')
    
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
    %                 direction_firing_rate = cell(1,length(degrees));
    %                 for bin = 2:length(degrees)
    %                     these_directions = find(saccade_direction{unit} < degrees(bin) & saccade_direction{unit} >= degrees(bin-1));
    %                     direction_firing_rate{bin} = 1000*nansum(fix_aligned(these_directions,window),2)/length(window);
    %                 end
    %                 degrees = degrees(2:end);
    %                 degrees = degrees*pi/180;
    %                 degrees = [degrees degrees(1)];
    %                 means = cellfun(@nanmean,direction_firing_rate(2:end));
    %
    %                 mrls(w) = circ_r(degrees(1:end-1),means,bin_deg/180*pi);
    %             end
    %             w = find(mrls == max(mrls)); %best direction tuning window
    %             window = ((step-1)*w+1):(step*w+window_width);
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Determine if Unit is Significantly Modulated by Saccade Direction--%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %take eye data from fixations in2in or out2out since in2out or out2in
    %could be biased by field location creating artificial direction tuning
    % but only do this if spatial (both criterion) or possibly spatial (either criterion)
    if (spatial_info.shuffled_rate_prctile(unit) > 95) || (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
        fix_aligned = fix_aligned(in_out{unit} == 2 | in_out{unit} == 4,:); %fixations in2in or out2out
        sac_dirs = sac_dirs(in_out{unit} == 2 | in_out{unit} == 4); %directions for fixations in2in or out2out
    end
    
    %---Saccade Direction across all Fixations Regardless of Fixation Location---%
    [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,fix_aligned,window,sac_dirs);
    mrls.all_fixations(unit) = mrl; %observed MRL
    binned_firing_rate_curves{1,unit} = mean_binned_firing_rate; %binned firing rates
    uniformity_pvalue(1,unit) = circ_rtest(sac_dirs*pi/180,sum(fix_aligned(:,window),2)*1000/length(window)); 
    %circ p_values may not be that useful since depends on firing rate and
    %sample size
    
    shuffled_mrls = NaN(1,numshuffs); %shuffled mrls
    parfor shuff = 1:numshuffs
        ind = randperm(length(sac_dirs)); %shuffle indeces
        dirs = sac_dirs(ind); %shuffled saccade directions
        mrl = circ_r(dirs*pi/180,sum(fix_aligned(:,window),2)/window_width); %MRL for unbinned data
        shuffled_mrls(shuff) = mrl;
    end
    mrls.all_fixations_shuffled{unit} = shuffled_mrls;
    mrls.all_fixations_shuffled_prctile(unit) = 100*sum(mrls.all_fixations(unit) > shuffled_mrls)/numshuffs; %shuffled percentile
    
    
    %---Saccade Direction across Fixations In2In or Out2Out---%
    %only run if spatial (both criterion) or possibly spatial (either criterion)
    if (spatial_info.shuffled_rate_prctile(unit) > 95) || (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
        
        %---In2In Fixations---%
        fix_in_in = list_fixation_locked_firing{unit}(in_out{unit} == 2,:); %fixation aligned firing rate
        dir_in_in = saccade_direction{unit}(in_out{unit} == 2); %saccade directions
        
        [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,fix_in_in,window,dir_in_in);
        mrls.in2in(unit) = mrl; %observed MRL
        if isnan(mrl);
            disp('now')
        end
        binned_firing_rate_curves{2,unit} = mean_binned_firing_rate; %binned firing rate
        uniformity_pvalue(2,unit) = circ_rtest(dir_in_in*pi/180,sum(fix_in_in(:,window),2)*1000/length(window));
        %circ p_values may not be that useful since depends on firing rate and
        %sample size
        
        shuffled_mrls = NaN(1,numshuffs);
        for shuff = 1:numshuffs
            ind = randperm(length(dir_in_in)); %shuffle indeces
            dirs = dir_in_in(ind); %shuffled directions
            mrl = circ_r(dirs*pi/180,sum(fix_in_in(:,window),2)/window_width); %MRL for unbinned data
            shuffled_mrls(shuff) = mrl;
        end
        mrls.in2in_shuffled{unit} = shuffled_mrls;
        mrls.in2in_shuffled_prctile(unit) = 100*sum(mrls.all_fixations(unit) > shuffled_mrls)/numshuffs; %shuffled percentile
        
        %---Out2Out Fixations---%
        fix_out_out = list_fixation_locked_firing{unit}(in_out{unit} == 4,:); %fixation aligned firing rate
        dir_out_out = saccade_direction{unit}(in_out{unit} == 4); %saccade directions
        
        [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,fix_out_out,window,dir_out_out);
        mrls.out2out(unit) = mrl; %observed MRL
        binned_firing_rate_curves{3,unit} = mean_binned_firing_rate; %binned firing rate
        uniformity_pvalue(3,unit) = circ_rtest(dir_out_out*pi/180,sum(fix_out_out(:,window),2)*1000/length(window));
        %circ p_values may not be that useful since depends on firing rate and
        %sample size
        
        shuffled_mrls = NaN(1,numshuffs);
        parfor shuff = 1:numshuffs
            ind = randperm(length(dir_out_out)); %shuffle indeces
            dirs = dir_out_out(ind);%shuffled saccade directions
            mrl = circ_r(dirs*pi/180,sum(fix_out_out(:,window),2)/window_width); %MRL for unbinned data
            shuffled_mrls(shuff) = mrl;
        end
        mrls.out2out_shuffled{unit} = shuffled_mrls;
        mrls.out2out_shuffled_prctile(unit) = 100*sum(mrls.all_fixations(unit) > shuffled_mrls)/numshuffs;%shuffled perecentile
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    %%%---Make Plot---%%%
    %%%%%%%%%%%%%%%%%%%%%
    
    figure
    
    %---Fixation Aligned Activity out2out & in2in---%
    subplot(3,3,1)
    dofill(t,fix_aligned,'black',1,smval); %smoothed fiirng rate curve
    hold on
    yl = ylim;
    if yl(1) < 0
        ylim([0 yl(2)]);
        yl(1) = 0;
    end
    h = fill([window(1) window(end) window(end) window(1) window(1)]-twin1,...
        [yl(1) yl(1) yl(2) yl(2) yl(1)],'r'); %window of interest
    set(h,'facealpha',.25,'EdgeColor','None')
    plot([0 0],[yl(1) yl(2)],'k--')
    hold off
    xlim([-twin1 twin2]);
    set(gca,'Xtick',[-twin1 0 twin1 twin2])
    xlabel('Time from Fixation Start (ms)')
    ylabel('Firing Rate (Hz)')
    title('All Fixation Aligned Activity')
    
    %---Fixation Aligned Acitivity out2out & in2in Plotted by Prefered and AntiPrefered Direction-All Fixations---%
    [prefered_dirs,anti_prefered_dirs,smoothed_direction_curve] = ...
        select_prefred_indeces(binned_firing_rate_curves{1,unit},degrees,sac_dirs,smval_deg);
    subplot(3,3,2)
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
    
    %---Smoothed Polar Plot of Prefered Firing Direction-Fixations out2out & in2in---%
    subplot(3,3,3)
    polar(degrees,[smoothed_direction_curve(end) smoothed_direction_curve],'b')
    if (spatial_info.shuffled_rate_prctile(unit) > 95) || (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
        %only ran if spatially modulated or possibly modulated (95% skaggs OR 95% corr 1/2)
        title(sprintf(['All mrl: ' num2str(mrls.all_fixations(unit),2) ' (' num2str(mrls.all_fixations_shuffled_prctile(unit),3) '%%)'...
            ' \n In2In mrl: ' num2str(mrls.in2in(unit),2) ' (' num2str(mrls.in2in_shuffled_prctile(unit),3) '%%)']))
    else
        title(sprintf(['All mrl: ' num2str(mrls.all_fixations(unit),2) ' (' num2str(mrls.all_fixations_shuffled_prctile(unit),3) '%%)']))
    end
    
    %---Plot Firing Rate Curves for each of the cardinal and intermediate directions---%
    subplot(3,3,4)
    plot(t,crude_direction_curves')
    hold on
    yl = ylim;
    plot([0 0],[yl(1) yl(2)],'k--')
    hold off
    xlim([-twin1 twin2]);
    set(gca,'Xtick',[-twin1 0 twin1 twin2])
    xlabel('Time From Fixation Start')
    ylabel('Firing Rate')
    title(['Firing Rate Curves by Crude 45' char(176) ' Direction']);
    box off

    %---Fixation Aligned Raster Sorted by Saccade Direction for out2out & in2in fixations---%
    subplot(3,3,5)
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
        set(h,'facealpha',.25,'EdgeColor','None')
        hold off
    end
    xlim([-twin1 twin2]);
    set(gca,'Xtick',[-twin1 0 twin1 twin2])
    ylabel('Ranked Sac. Dir.')
    xlabel('Time from Fixation Start (ms)')
    title('Fixation Aligned-Saccade Direction')
    box off
    
    %---Firing Rate Map---%
    subplot(3,3,6)
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
    
    %---Plots for Fixations Out2Out Only---%
    %not plotting In2In since this is less interesting and coverage isn't
    %great for all neurons anyway. I'm more intesteted in showing direction
    %tuning out of the field
    if (spatial_info.shuffled_rate_prctile(unit) > 95) || (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
        %only run if spatially modulated or possibly modulated (95% skaggs OR 95% corr 1/2)
        
        %---Smoothed Plot of Firing Rate by Saccade Direction Out2Out ixations---%
        [prefered_dirs,anti_prefered_dirs,smoothed_direction_curve] = ...
            select_prefred_indeces(binned_firing_rate_curves{3,unit},degrees,dir_out_out,smval_deg);
        subplot(3,3,7)
        polar(degrees,[smoothed_direction_curve(end) smoothed_direction_curve],'b')
        title(['Out2Out mrl: ' num2str(mrls.out2out(unit),2) ' (' num2str(mrls.out2out_shuffled_prctile(unit),3) '%)'])
        
        %---Out2Out Fixation Aligned Raster Sorted by Saccade Direction---%
        subplot(3,3,8)
        [~,si] = sort(dir_out_out);
        fix_aligned_sorted = fix_out_out(si,:);
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
        ylabel('Ranked Sac. Dir.')
        xlabel('Time from Fixation Start (ms)')
        title('Out2Out Fixation Aligned-Saccade Direction')
        box off
        
        %---Fixation Aligned Acitivity Plotted by Prefered and AntiPrefered Direction-All Fixations---%
        subplot(3,3,9)
        hold on
        dofill(t,fix_out_out(prefered_dirs,:),'green',1,smval);
        dofill(t,fix_out_out(anti_prefered_dirs,:),'black',1,smval);
        yl = ylim;
        if yl(1) < 0;
            yl(1) = 0;
            ylim([0 yl(2)]);
        end
        plot([0 0],[yl(1) yl(2)],'k--')
        hold off

        xlabel('Time from Fixation Start (ms)')
        ylabel('Firing Rate (Hz)')
        title('Out2Out Fixation Aligned Activity')
        legend('Prefered','Anti-Preferred')%,'Neutral1','Neutral2')
    end
    
    subtitle([task_file(1:8) ' ' unit_stats{1,unit}]);
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Saccade_Direction_Analysis']);
end


save([data_dir task_file(1:8) '-Saccade_Direction_Analysis.mat'],...
    'window_width','numshuffs','bin_deg','bin_deg2','mrls','uniformity_pvalue',...
    'binned_firing_rate_curves','unit_stats','smval')
end

function [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,fixation_aligned_firing,time_window,saccade_directions)
%smooth data, binned into 1 degree bins
degrees = [0:bin_deg:360]-180;
binned_firing_rate = cell(1,length(degrees));
for bins = 2:length(degrees)
    these_dirs = (saccade_directions < degrees(bins) & saccade_directions >= degrees(bins-1));
    binned_firing_rate{bins} = 1000*nansum(fixation_aligned_firing(these_dirs,time_window),2)/length(time_window);
end
degrees = degrees(2:end);
degrees = degrees*pi/180;
degrees = [degrees degrees(1)];
mean_binned_firing_rate = cellfun(@nanmean,binned_firing_rate(2:end));
mrl = circ_r(saccade_directions*pi/180,sum(fixation_aligned_firing(:,time_window),2)/length(time_window)); %for unbinned data

% NOT rechecked for bugs 1/17/2017 SDK
% binned data
% degrees = [0:bin_deg:360]-180; %binned degrees
% binned_firing_rate = cell(1,length(degrees));
% for bins = 2:length(degrees)
%     these_dirs = (saccade_directions < degrees(bins) & saccade_directions >= degrees(bins-1));
%     binned_firing_rate{bins} = 1000*nansum(fixation_aligned_firing(these_dirs,time_window),2)/length(time_window);
% end
% degrees = degrees(2:end);
% degrees = degrees*pi/180;
% degrees = [degrees degrees(1)];
% mean_binned_firing_rate = cellfun(@nanmean,binned_firing_rate(2:end));
% mbfr = mean_binned_firing_rate;
% mbfr(isnan(mbfr)) = 0;
% mrl = circ_r(degrees(1:end-1),mbfr,bin_deg/180*pi); %calculate for binned data
end

function [prefered_dirs,anti_prefered_dirs,binned_firing] = select_prefred_indeces(binned_firing,binned_directions,dirs,smval)
bf = [binned_firing(end-(3*smval):end) binned_firing binned_firing(1:3*smval)];%so don't get edege artifacts
binned_firing = nandens(bf,smval,'gauss',1); %smooth to display firing rate by saccade direction
binned_firing = binned_firing(3*smval+2:end-3*smval);%remove buffers

%smoothed more so small peaks don't interfer detecting prefered direction
binned_firing2 = nandens(bf,2*smval,'gauss',1);
binned_firing2 = binned_firing2(3*smval+2:end-3*smval);%remove buffers
prefered_direction = 180/pi*binned_directions(binned_firing2(1:end-1) == max(binned_firing2(1:end-1)));
prefered_direction = prefered_direction(1);%if multiple just pick 1
anti_prefered_direction = prefered_direction-180;
dirs2 = dirs; %store for later since may be modifying dirs variable

if prefered_direction > 135 %will have to look at negative angles too
    dirs(dirs < 0) = dirs(dirs < 0)+360;
elseif prefered_direction < -135
    prefered_direction = prefered_direction+360;
    dirs(dirs < 0) = dirs(dirs < 0)+360;
else
    dirs = dirs;
end
prefered_dirs = find((dirs <= prefered_direction+30)& (dirs >= prefered_direction-30));

dirs = dirs2;%restore dirs variable since may have just modifid above
if anti_prefered_direction > 135 %will have to look at negative angles too
    dirs(dirs < 0) = dirs(dirs < 0)+360;
elseif anti_prefered_direction < -135
    anti_prefered_direction = anti_prefered_direction+360;
    dirs(dirs < 0) = dirs(dirs < 0)+360;
else
    dirs = dirs;
end
anti_prefered_dirs = find((dirs <= anti_prefered_direction+30)& (dirs >= anti_prefered_direction-30));

end