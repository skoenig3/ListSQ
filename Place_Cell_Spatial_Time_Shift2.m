% Population Saccadic Modulation
% written Seth Konig 8/12/16
%
% code imports data from List_Saccade_AnalysisV2 results and looks at how
% modulated the whole population is by saccades. Analysis includes...
% 1) Average saccade/fixation-locked firing rate
% 2) AP Differences
%   a) % of neurons significantly modulated by eye movements
%   b) degree of modulation e.g. % change in firing rate?
% 3) Relative timing of modulation i.e. lag
%   a) on average
%   b) distribution
clar
task = 'ListSQ';
Fs = 1000;%Hz
smval = 60;
imageX = 800;
imageY = 600;
min_shift_dur = 40; % ms intervals relative to fixation peak 250/8

t = -500:499;

best_delay = [];
fix_delays = [];
worst_delay = [];
best_bit_increase = [];
worst_bit_decrease = [];

figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Population Figures\Place Cell Eye Time Shift2\';
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
        
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials','cfg','hdr','data','fixationstats');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        
        if exist([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat']);
            load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat']);
            filter_size = filter_width*10;
            H = fspecial('gaussian',filter_size,filter_width);
        else
            continue
        end
        
        for unit = 1:num_units
            if ~isempty(spike_times{unit})
                if ((spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.shuffled_spatialstability_prctile(unit) > 95)) 
                        if ((fixation_info.rate(unit) >  fixation_info_95.rate(unit)) &&  (saccade_info.rate(unit) >  saccade_info_95.rate(unit)))...
                            || ((fixation_info2.rate(unit) >  fixation_info2_95.rate(unit)) &&  (saccade_info2.rate(unit) >  saccade_info2_95.rate(unit)))
                        col = 1;
                        else
                            col = 0;
                        end
                        saccade_firing = saccade_locked_firing{unit}(saccade_information{unit}(:,4) > 2*twin,:);
                        fixation_firing = fixation_locked_firing{unit}(fixation_information{unit}(:,4) > 2*twin,:);
                        saccade_firing(nansum(saccade_firing,2) == 0,:) = [];%remove saccades without spikes
                        fixation_firing(nansum(fixation_firing,2) == 0,:) = [];%remove fixations without spikes
                        
                        
                        figure
                        %draw raw spikes on all images all locations
                        make_spike_jittered_plot(eyepos{unit},spike_times{unit},[4 5],1)
                        title(['Raw: All images, n = ' num2str(size(spike_times{unit},1))])
                        set(gca,'Xcolor','w')
                        set(gca,'Ycolor','w')
                        
                        subplot(4,5,6)
                        %Calculated smoothed firing rate for all images
                        filtered_time = filter_time(eyepos{unit},imageX,imageY,Fs,binsize,H);
                        filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
                        filtered_space = filter_space(eyepos{unit},spike_times{unit},imageX,imageY,binsize,H);
                        firing_rate = filtered_space./filtered_time;
                        
                        h = imagesc(firing_rate);
                        set(h,'alphadata',~isnan(filtered_time));
                        axis off
                        axis equal
                        fr = sort(firing_rate(1:end));
                        fr(isnan(fr)) = [];
                        clims(:,1) = caxis;
                        if length(fr) > 20
                            clims(2,1) = fr(round(0.99*length(fr)));% the ~99%-tile
                        end
                        colorbar
                        colormap('jet')
                        
                        title(sprintf(['All images, peak rate = ' num2str(max(fr),3) ...
                            ' Hz \n Bit: ' num2str(spatial_info.shuffled_rate_prctile(unit),3) '%%,'...
                            ' r = ' num2str(spatial_info.spatialstability(unit),2) ' ' ...
                            num2str(spatial_info.shuffled_spatialstability_prctile(unit),2) '%%']))
                        
                        
                        subplot(4,5,11)
                        hold on
                        % [yfix,~]= nandens(fixation_firing,smval,'gauss',Fs,'nanflt');
                        % [ysac,~]= nandens(saccade_firing,smval,'gauss',Fs,'nanflt');
                        [~,~,~,yfix,~] = dofill(t,fixation_firing,'blue',1,smval);
                        [~,~,~,ysac,~] = dofill(t,saccade_firing,'red',1,smval);
                        ylim([0.8*min([yfix ysac]),1.2*max([yfix ysac])]);
                        yl = ylim;
                        plot([0 0],[yl(1) yl(2)],'k--')
                        hold off
                        xlabel('Time from Eye Movement')
                        ylabel('FR (Hz)')
                        %                         legend('+Fixations','+Saccades')
                        if col == 1
                            title('Significant Eye Movement Modulation!')
                        end
                        
                        %find max modulation time so take abs
                        [~,fix_delay] = max(yfix(400:800));
                        [~,sac_delay] = max(yfix(400:800));
                        shift_dur = fix_delay-sac_delay; %pic shift relative to difference in fix and sac peaks s
                        if shift_dur < min_shift_dur
                            shift_dur = min_shift_dur;
                        end
                        fix_delay = fix_delay-100;
                        
                        shifts = fix_delay-3*shift_dur:shift_dur:fix_delay+4*shift_dur;
                        sbnums = [2 4 7 9 12 14 17 19];
                        shift_bits = NaN(1,length(shifts));
                        
                        for shift = 1:8
                            %shift eye movements relative to spikes
                            shift_eye = eye_time_shift(eyepos{unit},shifts(shift));
                            
                            %draw raw spikes on all images all locations
                            make_spike_jittered_plot(shift_eye,spike_times{unit},[4 5],sbnums(shift))
                            title(['Shift: ' num2str(shifts(shift)) 'ms'])
                            set(gca,'Xcolor','w')
                            set(gca,'Ycolor','w')
                            
                            filtered_time = filter_time(shift_eye,imageX,imageY,Fs,binsize,H);
                            filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
                            filtered_space = filter_space(shift_eye,spike_times{unit},imageX,imageY,binsize,H);
                            firing_rate = filtered_space./filtered_time;
                            
                            [shift_bits(shift),~] = estimated_mutual_information({filtered_time,filtered_space},1,'spatial_noshuff',smval,Fs);
                            
                            subplot(4,5,sbnums(shift)+1)
                            h = imagesc(firing_rate);
                            set(h,'alphadata',~isnan(filtered_time));
                            axis off
                            axis equal
                            caxis(clims)
                            colormap('jet')
                            title(['Shift: ' num2str(shifts(shift)) 'ms, max ' num2str(max(firing_rate(1:end)),3) ' Hz'])
                        end
                        
                        subplot(4,5,16)
                        hold on
                        bar(1,spatial_info.rate(unit),'r')
                        bar(2:9,shift_bits)
                        hold off
                        xlabel('Shift (ms)')
                        ylabel('Skaggs (bits)')
                        set(gca,'Xtick',1:9)
                        set(gca,'XtickLabel',[{'None'} num2cell(shifts) ])
                        set(gca,'XTickLabelRotation',45)
                        
                        [~,best] = max(shift_bits);
                        [~,worst] = min(shift_bits);
                        
                        best_delay = [best_delay shifts(best)];
                        fix_delays = [fix_delays fix_delay];
                        worst_delay = [worst_delay shifts(worst)];
                        best_bit_increase = [best_bit_increase (shift_bits(best)-spatial_info.rate(unit))/spatial_info.rate(unit)];
                        worst_bit_decrease = [worst_bit_decrease (shift_bits(worst)-spatial_info.rate(unit))/spatial_info.rate(unit)];
                        
                        if multiunit(unit)
                            multi_str = 'Multiunit ';
                        else
                            multi_str = ' ';
                        end
                        subtitle(['Spatial Plots' multi_str unit_names{unit}]);
                        
                        save_and_close_fig(figure_dir,[task_file(1:end-11) '-' unit_names{unit} '_Place_Cell_TimeShift'])
                end
            end
        end
    end
end