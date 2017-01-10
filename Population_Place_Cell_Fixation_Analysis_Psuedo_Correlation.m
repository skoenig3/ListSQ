% %get recording locations for view cells
% clar
% monkeys = {'Vivian','Tobii'};
%
% task = 'ListSQ';
% min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
%
% all_spatial_corrs = []; %place cell spatial correlations
% all_not_spatial_corrs = [];%non-place cell spatial correlations
%
% %---Firing Rate Curves for List---%
% all_firing_rates = []; %all fixations
% all_in_rates = [];%fixations out-> in
% all_out_rates = []; %fixaiton out-> out
% all_in_out_rates = [];%fixaions out-> in
% all_peaks = [];%peak firing time of place cells for out-> in fixations
%
% %---Other Parameters---%
% monkey_count = zeros(2,2);%row 1 place row 2 non-place
% all_unit_names = {}; %place cell unit names
% all_monkeys = []; %1s and 2s
% AP_location = []; %AP location of recorded place cell
%
%
% Fs = 1000;
% imageX = 800;
% imageY = 600;
%
% in_field_spike_times = cell(1,107);
% place_cell_index = 1;
%
% for monk = 2:-1:1
%     monkey = monkeys{monk};
%
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %---Read in Excel Sheet for Session data---%%%
%     %only need to run when somethings changed or sessions have been added
%     if strcmpi(monkey,'Vivian')
%         excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
%         excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
%         data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
%         figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
%
%
%         listsq_read_excel(data_dir,excel_file);
%         load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
%
%         predict_rt = 156;%156 ms prediction 5-percentile
%         chamber_zero = [13.5 -11]; %AP ML
%
%     elseif strcmpi(monkey,'Tobii')
%         excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
%         excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
%         data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
%         figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
%
%         predict_rt = 138;%ms prediction 5-percentile
%         chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
%
%         listsq_read_excel(data_dir,excel_file);
%         load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
%         session_data(end) = [];%last file doesn't have strobe signal working on importing the data
%     end
%
%     for sess = 1:length(session_data)
%         [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
%             get_task_data(session_data{sess},task);
%         if isempty(task_file)
%             continue
%         end
%
%         load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials','stability_attribute','cfg','hdr','data');
%
%         %get unit data
%         [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
%             multiunit,unit_confidence,sorting_quality);
%
%         if num_units == 0
%             continue
%         end
%
%         num_trials = length(cfg.trl);
%         valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
%         if all(isnan(valid_trials(1,:)))
%             continue
%         end
%
%         disp(task_file(1:8))
%         load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
%             'spatial_info','spike_times','eyepos','binsize','filter_width')
%
%         load([data_dir task_file(1:8) '-Place_Cell_Analysis.mat'],...
%             'twin1','twin2','smval','list_fixation_locked_firing','in_out',...
%             'task_file','area','unit_stats','all_place_field_matrix')
%
%         for unit = 1:num_units
%             if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
%                     && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
%
%                 firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 1,:);
%                 [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
%
%                 %find peaks
%                 [PKS,LOCS]= findpeaks(firing_rate,'MinPeakWidth',50);
%
%                 if ~isempty(LOCS)
%                     %remove peaks less than 1/2 the max
%                     LOCS(PKS < 0.5) = [];
%                     PKS(PKS < 0.5) = [];
%                 end
%
%                 if isempty(LOCS)
%                     continue %ignore neuron
%                 end
%
%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %---Firing Rate Curves for List---%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 %1) first fixation in: out-> in
%                 firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 1,:);
%                 in_field_spike_times{place_cell_index} = firing_rate;
%                 place_cell_index = place_cell_index+1;
%                 [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
%                 firing_rate = firing_rate-nanmean(firing_rate(1:twin1));%normalize to firing rate before fixation inside field
%                 firing_rate = firing_rate/max(firing_rate);
%
%                 firing_rate = firing_rate/firing_rate(LOCS(1));
%                 all_in_rates = [all_in_rates; firing_rate];
%
%
%
%                 all_spatial_corrs = [all_spatial_corrs spatial_info.spatialstability_halves(unit)];
%
%                 %---Other Parameters---%
%                 monkey_count(1,monk) = monkey_count(1,monk)+1;
%                 all_unit_names = [all_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}];
%                 all_monkeys = [all_monkeys monk]; %1s and 2s
%                 AP_location = [AP_location chamber_zero(1)+ session_data{sess}.location(1)]; %AP location of recorded place cell
%             end
%         end
%     end
% end
% %% Plot Population Firing Rate Curves
% figure
% subplot(2,2,1)
% [m,i] = max(all_in_rates(:,twin1:end),[],2);
% [mm,ii] = sort(all_peaks);
% [mm,ii] = sort(i);
% imagesc([-twin1:twin2-1],[1:size(all_in_rates,1)],all_in_rates(ii,:))
% colormap('jet')
% hold on
% plot([0 0],[1 size(all_in_rates,1)],'w--');
% hold off
% xlabel('Time from Fixation Start')
% ylabel('Neuron #')
% title('out -> in')
% caxis([0 1])

%%
% base_unit = 75;%putative inhbitory unit peaks just after fixaiton starts, TO160325_3_sig001a

lag = -599:599;
relative_lags = NaN(length(in_field_spike_times),length(in_field_spike_times));
for i = 1:length(in_field_spike_times)
    for ii = 1:length(in_field_spike_times);
        if i <= ii
           continue 
        end
        disp([num2str(i) ' ' num2str(ii)])
        
        numfix1 = size(in_field_spike_times{i},1);
        numfix2 = size(in_field_spike_times{ii},1);
        
        all_xc = NaN(numfix1*numfix2,2*size(in_field_spike_times{i},2)-1);
        indeces = reshape(1:numfix1*numfix2,numfix1,numfix2);
        
        for n1 = 1:numfix1;
            for n2 = 1:numfix2;
                
                all_xc(indeces(n1,n2),:) = xcorr(in_field_spike_times{i}(n1,:),in_field_spike_times{ii}(n2,:));
            end
        end
        
        smoothed_xc = nandens(all_xc,10,'gauss',1,'nanflt');
        
        relative_lags(i,ii) = lag(find(smoothed_xc == max(smoothed_xc)));
    end
end

%%
% lag = -599:599;
% smoothed_xc = nandens(all_xc,10,'gauss',1,'nanflt');
% figure
% plot(lag,smoothed_xc)
% hold on
% yl = ylim;
% plot([0 0],[0 yl(2)],'k--')
% ylim([0 yl(2)])
% hold off
% xlabel('Lag (ms)')
% ylabel('Coincendence/spike')
%%
for i = 1:length(in_field_spike_times)
    for ii = 1:length(in_field_spike_times);
        if i >  ii
            relative_lags(ii,i) = relative_lags(i,ii);
        end
    end
end