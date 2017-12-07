% Population Saccadic Modulation Phase
%draft written 10/11/16
%
% code imports data from List_Saccade_AnalysisV2 results and looks at how
% modulated the whole population is by saccades. Analysis includes...

%set(0,'DefaultFigureVisible','off');

clear,clc
task = 'ListSQ';
Fs = 1000;%Hz
smval = 60;
imageX = 800;
imageY = 600;
max_time = 250; %250 ms after saccade start

figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Population Figures\Saccade Spike Phase\';

smooth_deg = 3;% moving average filter width for saccade angle
bin_deg = 12;
degrees = 0:bin_deg:360; %binned degrees

%spike times vs spike phase
num_frequencies = 32;
all_slopes = NaN(350,num_frequencies);
all_correlations = NaN(350,num_frequencies);
all_phase_offsets = NaN(350,num_frequencies);
all_spike_counts = NaN(350,num_frequencies);
all_saccade_counts = NaN(350,num_frequencies);

%phase clustering within fixation period
all_MLR = NaN(350,num_frequencies);%mean resultant vector
all_mean_phase = NaN(350,num_frequencies);%mean prefered phase
NaN(350,num_frequencies);%p-value of MLR
all_median_phase = NaN(350,num_frequencies);%median prefered phase
all_std_phase = NaN(350,num_frequencies);%std around prefered phase
all_MLR_p = NaN(350,num_frequencies);%p-value of MLR
all_peak_phase=NaN(350,num_frequencies); %maximum in distribution
peak_firing_time = NaN(1,350);%time of peak firing aligned to saccade
peak_firing_time_limited = NaN(1,350);%time of peak firing aligned to saccade within 250 ms of saccade start
all_which_monkey = NaN(1,350);


cells = 1;

for monkey = 2:-1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
%         listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 155;%155.85 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 135;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
%         listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for session = 1:length(session_data)
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~,lfp_quality]=...
            get_task_data(session_data{session},task);
        if isempty(task_file)
            continue
        elseif all(lfp_quality == 0) %all LFPs are bad
            continue
        end
        disp(['Session #' num2str(session)])
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials','cfg','hdr','data');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        
        if exist([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat']);
        else
            continue
        end
        
        %code to find spike/LFP channels
        LFPchannels = find_desired_channels(cfg,'LFP');
        if length(LFPchannels) < 4
           for channel = 1:4
               if isempty(cell2mat(strfind(hdr.label,['AD0' num2str(channel)]))) %channel wasn't recorded
                   LFPchannels = [LFPchannels(1:channel-1) NaN LFPchannels(channel:end)]; %add NaN for formatting
                   break 
               end
           end
        end
        bad_channels = [];
        for channel = 1:length(LFPchannels)
            if cell2mat(strfind(hdr.label,['AD0' num2str(channel)])) %make sure have recorded channel
                if  lfp_quality(channel) == 0; %if it is bad
                    bad_channels = [bad_channels channel];
                end
            end
        end
        LFPchannels(bad_channels) = NaN;
        
        for unit = 1:num_units
            if nansum(nansum(fixation_locked_firing{unit})) > 100
                if temporal_info.fixation.shuffled_temporalstability_prctile(1,unit) > 95  %1st 2nd half corr
                             
                    lfp_channel_num = str2double(unit_stats{1,unit}(6));
                    if isnan(LFPchannels(lfp_channel_num)) %no LFP or bad lfp for channel spike was recorded on
                        lfp_channel_num = find(~isnan(LFPchannels));
                        lfp_channel_num =  lfp_channel_num(1);%just take the first one
                    end
                    desired_channel = LFPchannels(lfp_channel_num);
                    
                    %fixation_information{unit}(1,sac_ind) %trial #
                    %fixation_information{unit}(2,sac_ind) %saccade start time from trial start
                    
                    trial_nums = unique(fixation_information{unit}(:,1));
                    
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
                    image_trials = find(image_trials);
                    
                    LFPs = cell(1,length(image_trials));
                    for t = 1:length(image_trials)
%                         for  chan = 1:length(desired_channels);
                            LFPs{t} = [LFPs{t}; data(desired_channel).values{image_trials(t)}];
%                         end
                    end
                    
                    LFP_phase = cell(1,length(image_trials));
                    LFP_power = cell(1,length(image_trials));
                    [~,~,wfq]=waveletanalysis(LFPs{1},'n_fq',num_frequencies);
                    parfor t = 1:length(image_trials)
                        [trialpower,trialphase]=waveletanalysis(LFPs{t},'n_fq',num_frequencies);
                        LFP_phase{t} = trialphase;
                        LFP_power{t} = trialpower;
                    end
                    
                    saccade_firing = fixation_locked_firing{unit}(fixation_information{unit}(:,4) > image_on_twin,:);
                    saccade_firing(nansum(saccade_firing,2) == 0,:) = [];%remove saccades without spikes
                    [firing_rate,~]= nandens(saccade_firing,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    
                    spike_phase = cell(1,twin1+twin2);
                    saccade_num = cell(1,twin1+twin2);
                    saccade_aligned_phase = cell(1,32);
                    for freq = 1:32
                        saccade_aligned_phase{freq} = NaN(2000,twin1+twin2);
                    end
                    saccade_aligned_LFP = NaN(2000,twin1+twin2);
                    count = 0;
                    fixcount = 1;
                    for f = 1:length(fixation_information{unit})
                        if fixation_information{unit}(f,4) > image_on_twin;
                            if sum(fixation_locked_firing{unit}(f,:)) > 0 %there were spikes
                                trial = fixation_information{unit}(f,1);
                                spks = find(fixation_locked_firing{unit}(f,:));
                                sac_num = fixation_information{unit}(f,3);
                                t0 = fixation_information{unit}(f,2)-twin1-1;
                                fixt = fixation_information{unit}(f,2);
                                if fixt+twin2 < length(LFPs{trial})
                                    saccade_aligned_LFP(fixcount,:) = LFPs{trial}(fixt-twin1:fixt+twin2-1);
                                    for freq = 1:32
                                        saccade_aligned_phase{freq}(fixcount,:) =LFP_phase{trial}(freq,fixt-twin1:fixt+twin2-1);
                                    end
                                    fixcount = fixcount+1;
                                end
                                for s = 1:length(spks);
                                    tind = spks(s)+t0;
                                    spike_phase{spks(s)} = [spike_phase{spks(s)} LFP_phase{trial}(:,tind,:)];
                                    saccade_num{spks(s)} = [saccade_num{spks(s)} sac_num];
                                    count = count+1;
                                end
                            else
                                continue
                            end
                        else
                            continue
                        end
                    end
                    
                    time_phase = cell(1,length(wfq));
                    phase_sac_num = cell(1,length(wfq));
                    for freq = 1:length(wfq)
                        for t = 1:length(spike_phase)
                            dat = spike_phase{t};
                            sacnums = saccade_num{t};
                            if isempty(dat)
                                continue
                            end
                            if size(dat,3) == 1
                                sacnums = sacnums;
                            elseif size(dat,3) == 2
                                sacnums = [sacnums sacnums];
                            elseif size(dat,3) == 3
                                sacnums = [sacnums sacnums sacnums];
                            elseif size(dat,3) == 4
                                sacnums = [sacnums sacnums sacnums sacnums];
                            end
                            samples = size(dat,2)*size(dat,3);
                            phase = dat(freq,:,:);
                            phase = reshape(phase,1,numel(phase));
                            time_phase{freq} = [time_phase{freq} [t*ones(1,samples);phase(1:end)+pi]];
                            phase_sac_num{freq} = [phase_sac_num{freq} sacnums];
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%---Spike Times vs Spike Phase---%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    tm1 = -twin1+1:twin2;
                    figure
                    subplot(4,4,1)
                    plot(tm1,firing_rate)
                    hold on
                    yl = ylim;
                    if yl(1) < 0
                        yl(1) = 0;
                    end
                    yl(2) = 1.1*max(firing_rate);
                    ylim(yl);
                    plot([0 0],[yl(1) yl(2)],'k--')
                    hold off
                    xlabel('Time from saccade Start (ms)')
                    ylabel('Firing Rate (Hz)')
                    title('Saccade Aligned Firing Rate')
                    
                    %%
                    time_zero = twin1; %alinged to start of saccade
                    for freq = 1:length(wfq);
                       time_period = 1000/wfq(freq); %time period for 1 cycle
                        time_period = time_period*2; %2 cycles so doesn't go all the way around

                        use_these = find(time_phase{freq}(1,:) > time_zero & time_phase{freq}(1,:) <= time_zero+time_period);
                        
                        x = time_phase{freq}(1,use_these); %time
                        alpha = time_phase{freq}(2,use_these);%phase
                        
                        a=-2/time_period:0.001:2/time_period; %should adjust based on frequency
                        cosinepart=zeros(length(a),1);
                        sinepart=zeros(length(a),1);
                        R=zeros(length(a),1);
                        
                        phi = [alpha alpha+2*pi];
                        theta = [x x];
                        parfor i=1:length(a)
                            cosinepart(i)=sum(cos(phi-(2*pi*a(i)*theta)));
                            sinepart(i)=sum(sin(phi-(2*pi*a(i)*theta)));
                            firstterm=(cosinepart(i)/length(phi))^2;
                            secondterm=(sinepart(i)/length(phi))^2;
                            R(i)=sqrt(firstterm+secondterm);
                        end
                        
                        tm = twin2:twin2+max_time;
                        slope=a(R==max(R));
                        if length(slope) > 1 %usually occurs when no good fit
                            slope = median(slope);
                        elseif isempty(slope)
                            continue
                        end
                        
                        offset = atan2(sum(sin(phi-2*pi*slope*theta)),sum(cos(phi-2*pi*slope*theta)));
                        %from "Quantifying phiular–thetaear associations: Hippocampal phase
                        %precession" by Kempter et al
                        
                        
                        mean_phase = atan2(sum(sin(alpha)),sum(cos(alpha)))+pi;
                        phi_time = mod(2*pi*slope*x+offset,2*pi);
                        mean_phi_time = atan2(sum(sin(phi_time)),sum(cos(phi_time)))+pi;
                        
                        correlation = sum(sin(alpha-mean_phase).*sin(phi_time-mean_phi_time))./...
                            sqrt(sum(sin(alpha-mean_phase).^2)*sum(sin(phi_time-mean_phi_time).^2));
                        correlation = abs(correlation); %added not originaly part of equation since doing 2 period direction can be wrong
                        
                        if mod(freq,2) == 0 && freq > 1%plot every other
                            subplot(4,4,(freq)/2)
                            plot(time_phase{freq}(1,:),time_phase{freq}(2,:),'.k')
                            set(gca,'Xtick',0:100:600);
                            xlim([100 600])
                            set(gca,'XtickLabel',{'-200','-100','0','100','200','300','400'});
                            ylim([0 2*pi])
                            
                            tm = min(x):max(x);
                            y = 2*pi*slope*tm+offset;
                            yb = mod(y,2*pi) ;
                            hold on
                            plot(tm,yb,'r.')
                            
                             title([num2str(wfq(freq),3) ' Hz, r = ' num2str(correlation,3)]);% ', p = ' num2str(pval,3)])
                        end
                        
                        all_slopes(cells,freq) = slope;
                        all_correlations(cells,freq) = correlation;
                        all_phase_offsets(cells,freq) = offset;
                        all_spike_counts(cells,freq) = length(alpha);
                    end
                    
                    subtitle([task_file(1:8) ' ' unit_names{unit} ', ' num2str(count) ' saccades']);
                    %save_and_close_fig(figure_dir,[task_file(1:8) '_' unit_names{unit} '_Fit_Spike_Phase'])
                    %%
                        
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%---Spike Prefered Phase---%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %sort of PPC
                    
                    tm1 = -twin1+1:twin2;
                    figure
                    subplot(4,4,1)
                    plot(tm1,firing_rate)
                    hold on
                    yl = ylim;
                    if yl(1) < 0
                        yl(1) = 0;
                    end
                    yl(2) = 1.1*max(firing_rate);
                    ylim(yl);
                    plot([0 0],[yl(1) yl(2)],'k--')
                    hold off
                    xlabel('Time from saccade Start (ms)')
                    ylabel('Firing Rate (Hz)')
                    title('Saccade Aligned Firing Rate')
                    time_zero = twin1; %alinged to start of saccade
                    
                    for freq = 1:length(wfq);
                        time_period = max_time;%sacadde + avg fixation duration
                        
                        use_these = find(time_phase{freq}(1,:) > time_zero & time_phase{freq}(1,:) <= time_zero+time_period);
                                                if isempty(use_these)
                            continue
                        end
                        
                        
                        alpha = time_phase{freq}(2,use_these);%phases from -pi to pi
                        

                        MLR = circ_r(alpha');
                        mean_phase = circ_mean(alpha');
                        %stats = circ_stats(alpha'); %mean,median,etc...
                        [pval, z] = circ_rtest(alpha');
                        
                     
                            
                            alpha = alpha*180/pi;
                            phase_distribution = zeros(1,length(degrees));
                            for bin = 2:length(degrees)
                                phase_distribution(bin) = sum(alpha > degrees(bin-1) & alpha < degrees(bin));
                            end
                            phase_distribution(1) = [];
                            means =  [phase_distribution(end-6:end) phase_distribution phase_distribution(1:7)];
                            means = filtfilt(1/smooth_deg*ones(1,smooth_deg),1,means);
                            means = means(8:end-7);
                            means = means/sum(means);
                            
                      if mod(freq,2) == 1 && freq > 1%plot every other    
                            subplot(4,4,(freq+1)/2)
                            polarplot(degrees*pi/180,[means means(1)]);
                            hold on
                           %polarplot([stats.median stats.median],[0 max(means)])
                            hold off
                            set(gca,'RTick',[]);
                            set(gca,'ThetaTick',[0 45 90 135 180 225 270 315])

                            %title([num2str(wfq(freq),3) ' Hz, phase: ' num2str(stats.median*180/pi,3)]);
                      end
                        
                        maxpp = degrees(find(means == max(means)))*pi/180;
                        all_peak_phase(cells,freq) = maxpp(1); %maximum in distribution
                        all_MLR(cells,freq) = MLR;%mean resultant vector
                        all_mean_phase(cells,freq) = mean_phase;%mean prefered phase
                        all_MLR_p(cells,freq) = pval;%p-value of uniformity
                    end
                    
                    subtitle([task_file(1:8) ' ' unit_names{unit} ', ' num2str(count) ' saccades']);
                    %save_and_close_fig(figure_dir,[task_file(1:8) '_' unit_names{unit} '_Phase_Distribution'])
                    
                    %---Define Peak Times to Associated with phase---%
                    firing_rate = firing_rate-mean(firing_rate);%normalize firing rate
%                     if max(firing_rate) > abs(min(firing_rate))
                         max_peak_time = find(firing_rate == max(firing_rate))-twin1;
%                     else %trough is stronger than peak so take min
%                          max_peak_time = find(firing_rate == min(firing_rate))-twin1;
%                     end
                    peak_firing_time(cells) = max_peak_time;%time of peak firing aligned to saccade
                    
                    %take ~saccade + ~fixation duration and define peak
                    firing_rate = firing_rate(twin1:twin1+max_time);
%                     if max(firing_rate) > abs(min(firing_rate))
                        max_peak_time = find(firing_rate == max(firing_rate));
%                     else %trough is stronger than peak so take min
%                         max_peak_time = find(firing_rate == min(firing_rate));
%                     end                    
                    peak_firing_time_limited(cells) =  max_peak_time;%time of peak firing aligned to saccade within 250 ms of saccade start
                    
                    all_which_monkey(cells) = monkey;
                    cells = cells+1;
                end
            end
        end
    end
end
save
%%
%
% figure
% polarplot(alpha,x-min(x),'.')
% hold on
% polarplot(y,tm-min(x),'k')
% hold off
% %%
% figure
% polarplot(alpha(1:2000),phi_time(1:2000),'k.')

figure
plot(nanmean(all_MLR))
subplot(2,2,1)
plot(wfq,nanmean(all_MLR))
xlabel('Frequency (Hz)')
ylabel('MLR')
title('Stength of Phase Clustering: MLR')

subplot(2,2,3)
plot(wfq,nanmedian(all_MLR_p))
hold on
plot([wfq(1) wfq(end)],[0.05 0.05],'k--')
hold off
xlabel('Frequency (Hz)')
ylabel('Median MLR p-values')
title('Phase Clustering: Median MLR P-values')

subplot(2,2,2)
plot(all_mean_phase(:,3),peak_firing_time_limited,'k.')
xlabel('Theta @ 4.8 Hz Phase')
ylabel('Time @ Saccade Aligned Peak Modulation')
title('Theta Phase and Peak Firing Time 0:250 ms, r = 0.10')

subplot(2,2,4)
plot(all_mean_phase(:,3),peak_firing_time,'k.')
xlabel('Theta @ 4.8 Hz Phase')
ylabel('Time @ Saccade Aligned Peak Modulation')
title('Theta Phase and Peak Firing Time -200:400 ms, r = 0.19')

subtitle('Spike Phase Properties: Saccade Modulated Neurons')
 %%
 nans = find(isnan(peak_firing_time));
 peak_firing_time(nans) = [];
 peak_firing_time_limited(nans) = [];
 all_peak_phase(nans,:) = [];
 all_mean_phase(nans,:) = [];
 

all_mean_phase(all_mean_phase < 0) =  all_mean_phase(all_mean_phase < 0) + 2*pi;
all_peak_phase(all_peak_phase < 0) =  all_peak_phase(all_peak_phase < 0) + 2*pi;

corrs1 = [];
corrs2 = [];
corrs3 = [];
corrs4 = [];
 for f = 1:11;
    r = corrcoef(all_mean_phase(:,f),peak_firing_time);
    corrs1 = [corrs1 r(2)];
    r2 =  corrcoef(all_mean_phase(:,f),peak_firing_time_limited); 
    corrs2 = [corrs2 r2(2)];
    
    
    r = corrcoef(all_peak_phase(:,f),peak_firing_time);
    corrs3 = [corrs3 r(2)];
    r2 =  corrcoef(all_peak_phase(:,f),peak_firing_time_limited);
    corrs4 = [corrs4 r2(2)];
 end
 
 figure
 hold on
 plot(wfq(1:11),corrs1)
 plot(wfq(1:11),corrs2)
 plot(wfq(1:11),corrs3)
 plot(wfq(1:11),corrs4)
 hold off
 xlabel('Frequency')
 ylabel('Correlation')
 legend('Mean Phase/Peak Time','Mean Phase/Peak Time Limited','Peak Phase/Peak Time','Peak Phase/Peak time limited')
 
 %%
 figure
 hold on
 plot(wfq,nanmedian(all_correlations(all_which_monkey == 1,:)))
 plot(wfq,nanmedian(all_correlations(all_which_monkey == 2,:)))
 hold off
 xlabel('Frequency')
 ylabel('Spike Field Correlation')
 legend('Vivian','Tobii')
 xlim([4 80])
 title('Median Correlation-all units')