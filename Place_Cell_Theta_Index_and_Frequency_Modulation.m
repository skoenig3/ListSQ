% Code below creates population summary for Significnat Place Cells
% Written by Seth Konig
% Code does the following
% 1) Summarizes spatial correlations for place cell and non-place cells
% 2) Calculates fixation aligned firing rate curves for place cells
% 3) Calculates eye coverage and place field coverage for place cells
% 4) Tracks AP location, unit counts, and which monkey (not currently used)
% 5) Contextual differences between list and sequence task
% 6) Copies relevant figures for place cells to summary directory


clar %clear,clc
set(0,'DefaultFigureVisible','OFF');

%where to store spatial analysis figure copies
summary_directory = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Population Figures\Place Cell Frequency Modulation\';
if ~isdir(summary_directory)
    mkdir(summary_directory)
end

task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size

%---Misc. Parameters (Mostly Place Cells)---%
monkey_all_unit_count = zeros(2,2);%row 1 place row 2 non-place, column by monkey
all_multi_unit_count = zeros(1,2); %all_multi_unit_countunits
all_place_cell_unit_names = {}; %place cell unit names
all_place_cell_monkeys = []; %1s and 2s
all_place_cell_spatial_corrs = []; %place cell spatial correlations

theta_index_fixation_aligned = [];
theta_index_image_viewing_period = [];
thetaindex2 = [];
SS_PSD_theta_index = [];
all_SS_PSDs = [];
all_PSD_est = [];
all_PSD_est_rel = [];

monkeys = {'Vivian','Tobii'};
figure_dir = {};

all_PSD_est_rel = [];
for monk =2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir{1} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir{2} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    end
    
    for sess = 1:length(session_data)
        %read in task file data
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        %load task file data
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials',...
            'stability_attribute','cfg','hdr','data','fixationstats');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        if num_units == 0
            continue
        end
        
        num_trials = length(cfg.trl); %number of trials
        valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
        if all(isnan(valid_trials(1,:)))
            continue
        end
        
        %load spatial analysis data
        disp(task_file(1:8))
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
            'spatial_info','spike_times')
        
        %load Place Cell Fixation Analysis data
        load([data_dir task_file(1:8) '-Place_Cell_Analysis.mat'])
        if numshuffs < 1000
            error('Should have 1000 shuffles')%make sure file has right number of shuffles
        end
        if smval ~= 30
            error('Smoothing Value (2xStd) does not match expectations!')
        end
        
        for unit = 1:num_units
            if ~isempty(spatial_info.shuffled_info_rate{unit})&& length(spatial_info.shuffled_info_rate{unit}) < 1000
                error('Should have 1000 shuffles') %make sure file has right number of shuffles
            elseif isempty(spatial_info.shuffled_info_rate{unit})
                continue
            end
            
            if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                    && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                
                %---Misc Parameters---%
                monkey_all_unit_count(1,monk) = monkey_all_unit_count(1,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                all_place_cell_unit_names = [all_place_cell_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}]; %unit name
                all_place_cell_monkeys = [all_place_cell_monkeys monk]; %1s and 2s for monkey
                all_place_cell_spatial_corrs = [all_place_cell_spatial_corrs spatial_info.spatialstability_halves(unit)];
                
                %firing rate out-> in
                in_firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 1,:); %get spike trains
                in_xc = NaN(size(in_firing_rate,1),1001);
                for f = 1:size(in_firing_rate,1);
                    xc = xcorr(in_firing_rate(f,:),in_firing_rate(f,:));
                    in_xc(f,:) = xc(600-500:600+500);
                end
                in_xc(:,501) = 0;
                
                %firing rate out-> in & in ->in
                in_in_firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 1 |in_out{unit} == 2,:); %get spike trains
                in_in_xc = NaN(size(in_in_firing_rate,1),1001);
                for f = 1:size(in_in_firing_rate,1);
                    xc = xcorr(in_in_firing_rate(f,:),in_in_firing_rate(f,:));
                    in_in_xc(f,:) = xc(600-500:600+500);
                end
                in_in_xc(:,501) = 0;
                
                %firing rate for all fixations
                all_firing_rate = list_fixation_locked_firing{unit};%get spike trains
                all_xc = NaN(size(all_firing_rate,1),1001);
                for f = 1:size(all_firing_rate,1);
                    xc = xcorr(all_firing_rate(f,:),all_firing_rate(f,:));
                    all_xc(f,:) = xc(600-500:600+500);
                end
                all_xc(:,501) = 0;
                
                %firing rate across trial
                spikes_tms = spike_times{unit};
                trial_xc =  NaN(size(spikes_tms,1),1001);
                for tr = 1:size(spikes_tms,1);
                    spk_trl = spikes_tms(tr,1:4500);
                    xc = xcorr(spk_trl,spk_trl);
                    trial_xc(tr,:) = xc(4500-500:4500+500);
                end
                trial_xc(:,501) = 0;
                
                %---Plot Fourier Spectrum---%
                Fs = 1000;            % Sampling frequency
                L = 1001;             % Length of signal
                n = 2*2^nextpow2(L);
                f = Fs*(0:(n/2))/n;
                
                tm = -500:500;
                figure
                subplot(2,3,1)
                hold on
                dofill(tm,in_xc/sum(in_xc(:)),'red',1,20);
                dofill(tm,in_in_xc/sum(in_xc(:)),'blue',1,20);
                dofill(tm,all_xc/sum(in_xc(:)),'black',1,20);
                hold off
                yl = ylim;
                yl(1) = 0;
                ylim(yl);
                xlim([-500 500])
                ylabel('Correlation (a.u.)')
                xlabel('lag (ms)')
                legend('Out2In','?2In','All')
                title('AutoCorrelation for Fixations')
                
                subplot(2,3,2)
                [trial,time] = find(in_in_firing_rate == 1);
                plot(time-200,trial,'.r')
                xlabel('Time Fixation Onset (ms)')
                ylabel('Trial')
                xlim([-200 400])
                ylim([0 max(trial)+1])
                title('?2in fixations')
                box off
                
                b = mean(in_in_xc);
                Fs = 1000;
                NFFT = 2^16;
                f2 = Fs/2*linspace(0,1,NFFT/2+1);
                Y2 = fft(b-mean(b),NFFT);
                Y2 = abs(Y2(1:NFFT/2+1)).^2;
                Y2 = conv(Y2,ones(1,round(2/mode(diff(f)))),'same')';
                [~,pk] = nanmax((f2(:)>3&f2(:)<12).*Y2(:));
                pk0 = pk;
                thetaindex = mean(Y2(abs(f2-f2(pk))<=1))/mean(Y2(f2<125));
                thetaindex2 = [thetaindex2 thetaindex];
                
                subplot(2,3,3)
                hold on
                Y = fft(mean(in_xc)-mean(in_xc(:)),n);
                P = abs(Y/n);
                plot(f,P(1:n/2+1)/max(P),'red')
                Y = fft(mean(in_in_xc)-mean(in_in_xc(:)),n);
                P = abs(Y/n);
                plot(f,P(1:n/2+1)/max(P),'blue')
                Y = fft(mean(all_xc)-mean(all_xc(:)),n);
                P = abs(Y/n);
                plot(f,P(1:n/2+1)/max(P),'black')
                plot(f2,Y2/max(Y2),'green')
                xlabel('f (Hz)')
                xlim([0 31])
                ylabel('Power')
                
                Y = fft(mean(in_in_xc),n);
                P = abs(Y/n);
                [PKS,LOCS] = findpeaks(P(1:n/2+1));
                peak_freq = f(LOCS);
                peak_theta_freq = find(peak_freq >= 3 & peak_freq <= 12);
                if ~isempty(peak_theta_freq)
                    [~,max_theta_peak] = max(PKS(peak_theta_freq));
                    peak_loc = LOCS(peak_theta_freq(max_theta_peak));
%                     plot(f(peak_loc),P(peak_loc),'k*')
                    max_theta_peak_power = mean(P(peak_loc-1:peak_loc+1));
                    freq_range = find(f < 125);
                    theta_index = max_theta_peak_power/mean(P(freq_range));
                    title(['?2In Theta Index: ' num2str(theta_index,2)])
                else
                    theta_index = NaN;
                end
                  title(['?2In Theta Index: ' num2str(thetaindex,2)])
                hold off
                theta_index_fixation_aligned = [theta_index_fixation_aligned theta_index];
                
                
                subplot(2,3,4)
                dofill(tm,trial_xc/sum(trial_xc(:)),'black',1,20);
                yl = ylim;
                yl(1) = 0;
                ylim(yl);
                xlim([-500 500])
                ylabel('Correlation (a.u.)')
                xlabel('lag (ms)')
                title('Correlation for Whole Image Trial')
                
                %---Power Spectrum Modulation Emery Brown---%
                % SS-PSD Method
                
                N = 500; % number of spectrum samples in [0,fs/2) determining the
                % frequency spacing fs/(2N) between samples
                N_max = 140; % number of desired spectrum samples which tends to be much
                % lower then N in neural signals because of oversampling
                iter_Newton = 4; % number of newton iterations per EM iteration
                iter_EM = 130; % number of EM iterations
                gamma = 1e-4; % cross-validated value for sparsity parameter gamma
                
                spike_raster = in_in_firing_rate';
                K = 600;
                fs = 1000;
                n = 0:K-1;
                % calculating the PSTH
                PSTH = mean(spike_raster,2);
                PSTH_rate = length(find(PSTH>0))/K;
                                
                % initializing parameters for the gaussian state space method
                EM_iterations = 8;
                sig_EM_init = 0.5;
                x_init = -5.7;
                sig_x_init = 1;
                
                x_k_k = zeros(K+1,1);
                x_k_k1 = zeros(K+1,1);
                x_k_K = zeros(K+1,1);
                sig_k_k = zeros(K+1,1);
                sig_k_k1 = zeros(K+1,1);
                sig_k_K = zeros(K+1,1);
                sig_k_k1_K = zeros(K+1,1);
                
                x_k_k(1) = x_init;
                sig_k_k(1) = sig_x_init;
                
                sig = zeros(EM_iterations,1);
                sig(1) = sig_EM_init;
                
                for i = 1:EM_iterations-1
                    for j = 2:K+1
                        x_k_k1(j) = x_k_k(j-1);
                        sig_k_k1(j) = sig_k_k(j-1)+sig(i);
                        x_k_k(j) = x_k_k1(j) + sig_k_k1(j) * ( PSTH(j-1)-exp( x_k_k1(j) ) );
                        sig_k_k(j) = 1/( exp(x_k_k1(j)) + 1/sig_k_k1(j) );
                    end
                    
                    x_k_K(K+1) = x_k_k(K+1);
                    sig_k_K(K+1) = sig_k_k(K+1);
                    
                    for z = K:-1:1
                        A = sig_k_k(z)/sig_k_k1(z+1);
                        x_k_K(z) = x_k_k(z)+A*(x_k_K(z+1)-x_k_k1(z+1));
                        sig_k_K(z) = sig_k_k(z)+(A^2)*(sig_k_K(z+1)-sig_k_k1(z+1));
                        sig_k_k1_K(z) = A*sig_k_K(z+1);
                    end
                    
                    x_k_k(1) = x_k_K(1);
                    sig_k_k(1) = sig_k_K(1);
                    
                    B = sum(sig_k_K(2:end))+sum(sig_k_K(1:end-1))+sum((x_k_K(2:end)).^2)+sum((x_k_K(1:end-1)).^2)-2*sum(x_k_K(2:end).*x_k_K(1:end-1)+sig_k_k1_K(1:end-1));
                    sig(i+1) = B/K;
                end
                
                x_padded = zeros(2*N,1);
                x_padded(1:K) = x_k_K(2:end) - mean(x_k_K(2:end));
                x_k_K_PSD = abs(fft(x_padded)).^2;
                x_k_K_PSD = x_k_K_PSD(1:N_max)/max(x_k_K_PSD(1:N_max));
                
                all_SS_PSDs =[all_SS_PSDs; x_k_K_PSD'];

                ind_f = 0:0.5*fs/N:0.5*fs*(N_max-1)/N;
                subplot(2,3,5)
                plot(ind_f,x_k_K_PSD,'LineWidth',1.13);
                hold on
                [PKS,LOCS] = findpeaks(x_k_K_PSD);
                peak_freq = ind_f(LOCS);
                peak_theta_freq = find(peak_freq >= 3 & peak_freq <= 12);
                if ~isempty(peak_theta_freq)
                    [~,max_theta_peak] = max(PKS(peak_theta_freq));
                    peak_loc = LOCS(peak_theta_freq(max_theta_peak));
                    %plot(ind_f(peak_loc),x_k_K_PSD(peak_loc),'k*')
                    max_theta_peak_power = mean(x_k_K_PSD(peak_loc-1:peak_loc+1));
                    freq_range = find(ind_f < 125);
                    theta_index = max_theta_peak_power/mean(x_k_K_PSD(freq_range));
                    title(['SS-PSD, theta-index ' num2str(theta_index)]);
                else
                    theta_index = NaN;
                    title('SS-PSD');
                end
                hold off
                xlabel('frequency(Hz)');
                box off;
                xlim([0 30])
                SS_PSD_theta_index = [SS_PSD_theta_index theta_index];
                
                
                %% initializing the algorithm
                % constructing matrix A
                A=zeros(K,2*N_max);
                for i=1:K
                    for j=1:N_max
                        A(i,2*j-1)=cos(i*pi*(j-1)/N);
                        A(i,2*j)=-sin(i*pi*(j-1)/N);
                    end
                end
                A=2*pi*A/N;
                A(:,2)=[];
                
                mu_v = zeros(2*N_max-1,1);
                theta = zeros(2*N_max-1,iter_EM);
                PSD_est = zeros(N_max,iter_EM);
                
                % initialize theta for EM
                theta(:,1) = ones(2*N_max-1,1);
                PSD_est(2:end,1) = theta(2:2:end,1) + theta(3:2:end,1);
                
                % initialize Newton
                x = zeros(K,1);
                alpha = 0.5*ones(K,1);
                g = L*A'*(PSTH-alpha);
                H = -L*A'*diag( alpha.*(1-alpha) )*A - diag( 1./theta(:,1) );
                
                % start of Alg. 2
                
                % EM iterations
                for i=2:iter_EM
                    
                    % Newton iterations
                    
                    for j=1:iter_Newton
                        
                        mu_v = mu_v - H\g;
                        x = A*mu_v;
                        alpha = 1./(1+exp(-x));
                        g = L*A'*(PSTH-alpha) - mu_v./theta(:,i-1);
                        H = -L*A'*diag( alpha.*(1-alpha) )*A - diag( 1./theta(:,i-1) );
                        
                    end
                    
                    E = diag(-H\eye(2*N_max-1)) + mu_v.^2;
                    theta(:,i) = (-1+sqrt(1+8*E*gamma))/(4*gamma);
                    PSD_est(2:end,i) = theta(2:2:end,i)+theta(3:2:end,i);
                    
                end
                
                m = max(max(PSD_est));
                ind_EM = 1:iter_EM;
                ind_f = 0:0.5*fs/N:0.5*fs*(N_max-1)/N;
                
                subplot(2,3,6);
                plot(ind_f,PSD_est(:,iter_EM)/m,'LineWidth',1.13);
                axis tight;
                title('Output of EM at Iteration 130')
                xlabel('frequency(Hz)');
                ylabel('PSD MLE estimate');
                grid on;
                title('Emery Brown Method')
                xlim([0 30])
                box off
                
                all_PSD_est = [all_PSD_est; PSD_est(:,iter_EM)'];
                all_PSD_est_rel = [all_PSD_est_rel; PSD_est(:,iter_EM)'/m];

                %%
                
                %close
                save_and_close_fig(summary_directory,[task_file(1:end-11) '_' unit_stats{1,unit} '_SpikeAutoCorrelation'])
            end
        end
    end
end
set(0,'DefaultFigureVisible','ON');
%%
figure
plot(ind_f,mean(all_SS_PSDs))
xlim([0 30])
xlabel('Frequncy (Hz)')
ylabel('Relative Power')

figure
plot(ind_f,mean(all_PSD_est))
xlim([0 30])
xlabel('Frequncy (Hz)')
ylabel('Relative Power')

figure
plot(ind_f,mean(all_PSD_est_rel))
xlim([0 30])
xlabel('Frequncy (Hz)')
ylabel('Relative Power')
