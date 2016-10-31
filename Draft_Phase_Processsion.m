% NO CLEAR Phase Precession :(
%% LFP = data(9).values{12};
% plot(LFP);
% [trialpower,trialphase,wfq]=waveletanalysis(LFP,'n_fq',100);
% %%
% t = 1:length(trialphase);
%
% imagesc(t,wfq,trialphase)
% axis xy
%%
% Population Saccadic Modulation
%draft written 10/10/16-10/11/16
%
% code imports data from List_Saccade_AnalysisV2 results and looks at how
% modulated the whole population is by saccades. Analysis includes...
% 1) Average saccade/saccade-locked firing rate
% 2) AP Differences
%   a) % of neurons significantly modulated by eye movements
%   b) degree of modulation e.g. % change in firing rate?
% 3) Relative timing of modulation i.e. lag
%   a) on average
%   b) distribution
% clar

clear,clc
task = 'ListSQ';
Fs = 1000;%Hz
smval = 60;
imageX = 800;
imageY = 600;
max_time = 250; %250 ms after saccade start

figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Population Figures\Phase Procession\';

all_unit_positions = zeros(1,length(locations));
saccade_cell_poisitions = zeros(1,length(locations));
avg_saccade_firing_rates = [];
all_unit_names = {};

num_frequencies = 32;
all_slopes = NaN(350,num_frequencies);
all_correlations = NaN(350,num_frequencies);
all_phase_offsets = NaN(350,num_frequencies);
all_spike_counts = NaN(350,num_frequencies);
all_saccade_counts = NaN(350,num_frequencies);
cells = 1;

for monkey = 2:-1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for session = 1:length(session_data)
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{session},task);
        if isempty(task_file)
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
        
        for unit = 1:num_units
            if nansum(nansum(saccade_locked_firing{unit})) > 100
                if 1%temporal_info.saccade.shuffled_temporalstability_prctile(1,unit) > 95  %1st 2nd half corr
                    
                    [desired_channels] = find_desired_channels(cfg,'LFP');
                    
                    %saccade_information{unit}(1,sac_ind) %trial #
                    %saccade_information{unit}(2,sac_ind) %saccade start time from trial start
                    
                    trial_nums = unique(saccade_information{unit}(:,1));
                    
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
                        for  chan = 1:length(desired_channels);
                            LFPs{t} = [LFPs{t}; data(desired_channels(chan)).values{image_trials(t)}];
                        end
                    end
                    
                    LFP_phase = cell(1,length(image_trials));
                    LFP_power = cell(1,length(image_trials));
                    [~,~,wfq]=waveletanalysis(LFPs{1},'n_fq',num_frequencies);
                    parfor t = 1:length(image_trials)
                        [trialpower,trialphase]=waveletanalysis(LFPs{t},'n_fq',num_frequencies);
                        LFP_phase{t} = trialphase;
                        LFP_power{t} = trialpower;
%                         disp(num2str(t));
                    end
                    
                    saccade_aligned_LFP = zeros(length(desired_channels),2*twin2);
                    saccade_aligned_power = zeros(size(LFP_power{1},1),2*twin2);
                    saccade_aligned_phase = zeros(size(LFP_power{1},1),2*twin2);
                    count = 0;
                    for f = 1:length(saccade_information{unit})
                        if saccade_information{unit}(f,4) > image_on_twin;
                            tind = saccade_information{unit}(f,2)-(twin2-1):saccade_information{unit}(f,2)+twin2;
                            trial = saccade_information{unit}(f,1);
                            saccade_aligned_LFP = [saccade_aligned_LFP+LFPs{trial}(:,tind)];
                            count = count+1;
                            saccade_aligned_power = saccade_aligned_power+mean(LFP_power{trial}(:,tind,:),3);
                            saccade_aligned_phase = saccade_aligned_phase+mean(LFP_phase{trial}(:,tind,:),3);
                        else
                            continue
                        end
                    end
                    saccade_aligned_LFP =  saccade_aligned_LFP/count;
                    saccade_aligned_power = saccade_aligned_power/count;
                    saccade_aligned_phase = saccade_aligned_phase/count;
                    
                    saccade_firing = saccade_locked_firing{unit}(saccade_information{unit}(:,4) > image_on_twin,:);
                    saccade_firing(nansum(saccade_firing,2) == 0,:) = [];%remove saccades without spikes
                    [firing_rate,~]= nandens(saccade_firing,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    
                    tm1 = -twin1+1:twin2;
                    tm2 = -twin2+1:twin2;
                    figure
                    subplot(2,2,1)
                    plot(tm1,firing_rate)
                    hold on
                    yl = ylim;
                    plot([0 0],[yl(1) yl(2)],'k--')
                    hold off
                    xlabel('Time from saccade Start (ms)')
                    ylabel('Firing Rate (Hz)')
                    title('saccade Aligned Firing Rate')
                    
                    subplot(2,2,2)
                    plot(tm2,saccade_aligned_LFP')
                    hold on
                    yl = ylim;
                    plot([0 0],[yl(1) yl(2)],'k--')
                    hold off
                    xlabel('Time from saccade Start (ms)')
                    ylabel('LFP (uV)')
                    title('saccade Aligned LFP')
                    
                    low_ind = find(wfq <= 30);
                    low_freq = wfq(wfq <= 30);
                    high_ind = find(wfq > 30);
                    high_freq = wfq(wfq > 30);
                    
                    subplot(4,4,[9 10])
                    imagesc(tm2,high_freq,saccade_aligned_power(high_ind,:))
                    hold on
                    plot([0 0],[high_freq(1) high_freq(end)],'w--')
                    hold off
                    axis xy
                    colormap('jet')
                    ylabel('LFP power')
                    title('saccade Aligned LFP power')
                    
                    subplot(4,4,[13 14])
                    imagesc(tm2,low_freq,saccade_aligned_power(low_ind,:))
                    hold on
                    plot([0 0],[low_freq(1) low_freq(end)],'w--')
                    hold off
                    axis xy
                    colormap('jet')
                    xlabel('Time from saccade Start (ms)')
                    ylabel('LFP power')
                    
                    
                    %                     subplot(2,2,4)
                    %                     imagesc(tm2,wfq,saccade_aligned_phase)
                    %                     hold on
                    %                     plot([0 0],[wfq(1) wfq(end)],'w--')
                    %                     hold off
                    %                     axis xy
                    %                     colormap('jet')
                    %                     xlabel('Time from saccade Start (ms)')
                    %                     ylabel('Phase')
                    %                     title('saccade Aligned LFP Phase')
                    
                    subplot(4,4,[11 12])
                    imagesc(tm2,high_freq,saccade_aligned_phase(high_ind,:))
                    hold on
                    plot([0 0],[high_freq(1) high_freq(end)],'w--')
                    hold off
                    axis xy
                    colormap('jet')
                    ylabel('LFP power')
                    title('saccade Aligned LFP Phase')
                    
                    subplot(4,4,[15 16])
                    imagesc(tm2,low_freq,saccade_aligned_phase(low_ind,:))
                    hold on
                    plot([0 0],[low_freq(1) low_freq(end)],'w--')
                    hold off
                    axis xy
                    colormap('jet')
                    xlabel('Time from saccade Start (ms)')
                    ylabel('LFP power')
                    
                    subtitle([task_file(1:8) ' ' unit_names{unit} ':' num2str(count) ' Saccades'])
                    
                    save_and_close_fig(figure_dir,[task_file(1:8) '_' unit_names{unit} '_Saccade_Triggered_Average_LFP'])
                    
                    
                    STA_LFP = zeros(length(desired_channels),2*twin2);
                    spike_phase = cell(1,twin1+twin2);
                    saccade_duration = cell(1,twin1+twin2);
                    count = 0;
                    for f = 1:length(saccade_information{unit})
                        if saccade_information{unit}(f,4) > image_on_twin;
                            if sum(saccade_locked_firing{unit}(f,:)) > 0 %there were spikes
                                trial = saccade_information{unit}(f,1);
                                spks = find(saccade_locked_firing{unit}(f,:));
                                sacdur = saccade_information{unit}(f,5);
                                
                                t0 = saccade_information{unit}(f,2)-twin1-1;
                                
                                for s = 1:length(spks);
                                    tind = (spks(s)-(twin2-1):spks(s)+twin2)+t0;
                                    if tind(end) > length(LFPs{trial})
                                        continue
                                    end
                                    STA_LFP = [STA_LFP+LFPs{trial}(:,tind)];
                                    
                                    tind = spks(s)+t0;
                                    spike_phase{spks(s)} = [spike_phase{spks(s)} LFP_phase{trial}(:,tind,:)];
                                    saccade_duration{spks(s)} = [saccade_duration{spks(s)} sacdur];
                                    count = count+1;
                                end
                            else
                                continue
                            end
                        else
                            continue
                        end
                    end
                    STA_LFP = STA_LFP/count;
                    
                    
                    %                     figure
                    %                     subplot(2,2,1)
                    %                     plot(tm2,STA_LFP')
                    %                     xlabel('Time from saccade Start (ms)')
                    %                     ylabel('STA LFP')
                    %                     title('Spike Triggered Average LFP')
                    
                    
                    time_phase = cell(1,length(wfq));
                    phase_sac_durs = cell(1,length(wfq));
                    for freq = 1:length(wfq)
                        for t = 1:length(spike_phase)
                            dat = spike_phase{t};
                            sacdurs = saccade_duration{t};
                            if isempty(dat)
                                continue
                            end
                            if size(dat,3) == 1
                                sacdurs = sacdurs;
                            elseif size(dat,3) == 2
                                sacdurs = [sacdurs sacdurs];
                            elseif size(dat,3) == 3
                                sacdurs = [sacdurs sacdurs sacdurs];
                            elseif size(dat,3) == 4
                                sacdurs = [sacdurs sacdurs sacdurs sacdurs];
                            end
                            samples = size(dat,2)*size(dat,3);
                            phase = dat(freq,:,:);
                            phase = reshape(phase,1,numel(phase));
                            time_phase{freq} = [time_phase{freq} [t*ones(1,samples);phase(1:end)+pi]];
                            phase_sac_durs{freq} = [phase_sac_durs{freq} sacdurs];
                        end
                    end
                    
                    figure
                    for i = 1:2:num_frequencies
                        subplot(4,4,(i+1)/2)
                        plot(time_phase{i}(1,:),time_phase{i}(2,:),'.k')
                        set(gca,'Xtick',0:100:600);
                        xlim([100 600])
                        set(gca,'XtickLabel',{'-200','-100','0','100','200','300','400'});
                        title(num2str(wfq(i)))
                        ylim([-pi pi])
                    end
                    
                    subtitle([task_file(1:8) ' ' unit_names{unit} ':' num2str(count) ' saccades'])
                    save_and_close_fig(figure_dir,[task_file(1:8) '_' unit_names{unit} '_Spike_Phase'])
                    
                    mean_sacdur = mean(phase_sac_durs{1});
                    
                    sb = [1 2 3 4 9 10 11 12];
                    figure
                    for i = 1:8
                        subplot(4,4,sb(i));
                        plot(time_phase{i}(1,phase_sac_durs{freq} <mean_sacdur),time_phase{i}(2,phase_sac_durs{freq} <mean_sacdur),'.k')
                        set(gca,'Xtick',0:100:600);
                        xlim([100 600])
                        set(gca,'XtickLabel',{'-200','-100','0','100','200','300','400'});
                        title([num2str(wfq(i)) '-short'])
                        ylim([-pi pi])
                        
                        subplot(4,4,4+sb(i));
                        plot(time_phase{i}(1,phase_sac_durs{freq} >=mean_sacdur),time_phase{i}(2,phase_sac_durs{freq} >=mean_sacdur),'.k')
                        set(gca,'Xtick',0:100:600);
                        xlim([100 600])
                        set(gca,'XtickLabel',{'-200','-100','0','100','200','300','400'});
                        title([num2str(wfq(i)) '-long'])
                        ylim([-pi pi])
                        
                    end
                    subtitle([task_file(1:8) ' ' unit_names{unit} ':' num2str(count) ' saccades, median dur: ' num2str(mean_sacdur) ' ms']);
                    save_and_close_fig(figure_dir,[task_file(1:8) '_' unit_names{unit} '_Spike_Phase_sacDur'])
                    
                    time_zero = twin1; %alinged to start of saccade
                    rhos = [];
                    pvals = [];
                    Rs = [];
                    for freq = 1:length(wfq);
                        time_period = 1000/wfq(freq); %time period for 1 cycle
                        time_period = time_period*3/4; %0.75 cycles so doesn't go all the way around
                        
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
                        
                        tm = min(x):max(x);
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
                        [rho pval] = circ_corrcc(phi_time,alpha);
                        
                        if rho-correlation > 0.01 %rounding error
                            disp('error wft')
                        end
                        
                        if mod(freq,2) == 1 %plot every other
                            subplot(4,4,(freq+1)/2)
                            plot(x,alpha,'k.')
                            xlim([time_zero  round(time_period)+time_zero])
                            ylim([0 2*pi])
                            
                            
                            tm = min(x):max(x);
                            y = 2*pi*slope*tm+offset;
                            yb = mod(y,2*pi) ;
                            hold on
                            plot(tm,yb,'r.')
                            
                            title([num2str(wfq(freq)) ' Hz, r = ' num2str(correlation,3)]);% ', p = ' num2str(pval,3)])
                        end
                        
                        
                        if correlation < -0.2 && slope > 0
                            disp('what')
                        end

                        
                        all_slopes(cells,freq) = slope;
                        all_correlations(cells,freq) = rho;
                        all_phase_offsets(cells,freq) = offset;
                        all_spike_counts(cells,freq) = length(alpha);
                        all_saccade_counts(cells,freq) = sum(cellfun(@numel,saccade_duration));
                    end
                    subtitle([task_file(1:8) ' ' unit_names{unit} ':' num2str(count) ' saccades']);
                    save_and_close_fig(figure_dir,[task_file(1:8) '_' unit_names{unit} '_Fit_Spike_Phase'])
                    
                    cells = cells+1;
                end
            end
        end
    end
end

%%
% tp = time_phase{10};
% tp(:,tp(1,:) < 200) = [];
% tp(:,tp(1,:) > 400) = [];
% 
% figure
% rose(tp(1,:),30)
% %%
% tp = tp(2,:)*180/pi;
% %%
% 
% bin_deg = 6 ; %number of degrees per bin
% moving_avg = 4;% moving average filter width for saccade angle
% degrees = [0:bin_deg:360]-180; %binned degrees
% degrees = degrees*pi/180;
% 
% figure
% for i = 10:18
%     subplot(3,3,i-9)
%     
%     tp = time_phase{i};
%     tp(:,tp(1,:) < 200) = [];
%     tp(:,tp(1,:) > 400) = [];
%     tp = tp(2,:);
%     
%     
%     
%     avg_phase =  zeros(1,length(degrees));
%     for bin = 2:length(degrees)
%         avg_phase(bin) = sum(tp < degrees(bin) & tp >= degrees(bin-1));
%     end
%     avg_phase(1) = [];
%     means =  [avg_phase(end-6:end) avg_phase avg_phase(1:7)];
%     means = filtfilt(1/moving_avg*ones(1,moving_avg),1,means);
%     avg_phase = means(8:end-7);
%     
%     polarplot(degrees,[avg_phase(end) avg_phase],'b')
%     title(num2str(wfq(i)))
% end
% %%
% time_zero = 150; %50 ms before saccade starts
% for freq = 1:length(wfq);
%     time_period = 1000/wfq(freq); %time period for 1 cycle
%     time_period = time_period*2; %2 cycles
%     
%     use_these = find(time_phase{freq}(1,:) > time_zero & time_phase{freq}(1,:) <= time_zero+time_period);
%     
%     x = time_phase{freq}(1,use_these); %time
%     alpha = time_phase{freq}(2,use_these);%phase
%     [rho pval ts] = phi_corrcl(alpha,x);
%     
%     
% end
%%

figure
polarplot(alpha,x-min(x),'.')
hold on
polarplot(y,tm-min(x),'k')
hold off
%%
figure
polarplot(alpha(1:2000),phi_time(1:2000),'k.')

