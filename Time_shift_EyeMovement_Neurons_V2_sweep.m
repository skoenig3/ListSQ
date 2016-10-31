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
min_shift_dur = 30; % ms intervals relative to fixation peak 250/8

binsize = 12; %pixels per bin spatial bin in either dimension 1/2 dva
filter_width = 6; %std of 2D guassian filter ~ 3 dva
filter_size = filter_width*10;
H = fspecial('gaussian',filter_size,filter_width);
%10 ms seems too low and 100 ms may be a little high

min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)

best_delay = [];
fix_delays = [];
worst_delay = [];
best_bit_increase = [];
worst_bit_decrease = [];
new_spatial_cells = 0;
total_cells = 0;

trough_to_best = [];

base = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Population Figures\Spatial Time Shift\';
figure_dir = [base 'Eye Movement\'];
figure_dir1 = [base 'Eye Movement\Place Cells\'];
figure_dir2 = base;
figure_dir3 = [base 'Place Cells\'];

for monkey = 1:2
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
        else
            continue
        end
        
        for unit = 1:num_units
            if ~isempty(spike_times{unit})
                if temporal_info.fixation.shuffled_temporalstability_prctile(1,unit) > 95 ... %1st 2nd half corr
                        col = 1;
                else
                    col = 2;
                end
                saccade_firing = saccade_locked_firing{unit}(saccade_information{unit}(:,4) > 2*twin,:);
                fixation_firing = fixation_locked_firing{unit}(fixation_information{unit}(:,4) > 2*twin,:);
                saccade_firing(nansum(saccade_firing,2) == 0,:) = [];%remove saccades without spikes
                fixation_firing(nansum(fixation_firing,2) == 0,:) = [];%remove fixations without spikes
                
                
                figure
                
                %plot color coded by portion of session/by time
                make_spike_jittered_colored_plot(eyepos{unit},spike_times{unit},[3 3],1)
                set(gca,'Xcolor','w')
                set(gca,'Ycolor','w')
                title(sprintf([    '\\rho_{1/2} = ' num2str(spatial_info.spatialstability_halves(1,unit),2) ...
                    ' ' num2str(spatial_info.spatialstability_halves_prctile(1,unit),3) '%%']))
                
                %draw raw spikes on all images all locations
                make_spike_jittered_plot(eyepos{unit},spike_times{unit},[3 3],4)
                title(['Raw: All images, n = ' num2str(size(spike_times{unit},1))])
                set(gca,'Xcolor','w')
                set(gca,'Ycolor','w')
                
                subplot(3,3,7)
                %Calculated smoothed firing rate for all images
                filtered_time = filter_time(eyepos{unit},imageX,imageY,Fs,binsize,H);
                filtered_time(filtered_time == 0) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
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
                
                title(sprintf(['max FR = ' num2str(max(fr),3) ...
                    ' Hz \n Bit: ' num2str(spatial_info.shuffled_rate_prctile(unit),3) '%%']))
                
                
                t = -twin1:twin2-1;
                subplot(3,3,2)
                hold on
                [~,~,~,yfix,~] = dofill(t,fixation_firing,'blue',1,smval);
                yl = ylim;
                plot([0 0],[yl(1) yl(2)],'k--')
                hold off
                xlabel('Time from Fixation Start (ms)')
                ylabel('FR (Hz)')
                xlim([-twin1 twin2])
                
                
                %find max modulation time so take abs
                [~,fix_delay] = min(yfix);
                shift_dur = min_shift_dur;
                fix_delay = fix_delay-twin1;
 
                
                shifts = fix_delay-20*shift_dur:shift_dur:fix_delay+20*shift_dur;
                shifts(shifts < -twin1) = [];
                shifts(shifts > twin2) = [];
                
                
                shift_bits = NaN(1,length(shifts));
                shift_corrs = NaN(2,length(shifts));
                
                
                for shift = 1:length(shifts)
                    %shift eye movements relative to spikes
                    shift_eye = eye_time_shift(eyepos{unit},shifts(shift));
                    
                    trial_data{1} = shift_eye;
                    trial_data{2} = spike_times{unit};
                    trial_data{3} = [imageX,imageY];
                    spatial_smval = [binsize filter_width];
                    
                    [spatial_shift_info,~] = estimated_mutual_information(trial_data,1,'spatial_noshuff',spatial_smval,Fs);
                    
                    shift_bits(shift) = spatial_shift_info.skaggs;
                    shift_corrs(:,shift) = spatial_shift_info.spatialstability';
                end
                
                bit_95 = prctile(spatial_info.shuffled_info_rate{unit},95);
                subplot(3,3,9)
                hold on
                bar(1,spatial_info.rate(unit),'r')
                bar(2:length(shift_bits)+1,shift_bits)
                plot([0 length(shift_bits)+1],[bit_95 bit_95],'k--')
                hold off
                xlabel('Shift (ms)')
                ylabel('Skaggs (bits\\spike)')
                set(gca,'Xtick',1:3:length(shift_bits))
                set(gca,'XtickLabel',[{'None'} num2cell(shifts(1:3:end))])
                set(gca,'XTickLabelRotation',45)
                xlim([0 length(shift_bits)+2])
                
                if all(spatial_info.rate(unit) > shift_bits)
                    title('Skaggs Info, Best: No Shft')
                else
                    title(['Skaggs Info, Best: ' num2str(shifts(shift_bits == max(shift_bits)))])
                end
                
                corr_95_half = prctile(spatial_info.shuffled_spatialstability_halves{unit}(1,:),95);
                subplot(3,3,3)
                hold on
                bar(1,spatial_info.spatialstability_halves(1,unit),'r')
                bar(2:length(shift_corrs)+1,shift_corrs(1,:))
                plot([0 length(shift_corrs)+1],[corr_95_half corr_95_half],'k--')
                hold off
                xlabel('Shift (ms)')
                ylabel('Spatial corr 1/2')
                set(gca,'Xtick',1:3:length(shift_corrs(1,:)))
                set(gca,'XtickLabel',[{'None'} num2cell(shifts(1:3:end))])
                set(gca,'XTickLabelRotation',45)
                xlim([0 length(shift_bits)+2])
                
                if all(spatial_info.spatialstability_halves(1,unit) > shift_corrs(1,:))
                    title('Spatial Corr 1/2, Best: No Shft')
                else
                    title(['Spatial Corr 1/2: ' num2str(shifts(shift_corrs(1,:) == max(shift_corrs(1,:))))])
                end
                
                
                corr_95_even_odds = prctile(spatial_info.shuffled_spatialstability_even_odd{unit}(1,:),95);
                subplot(3,3,6)
                hold on
                bar(1,spatial_info.spatialstability_even_odd(1,unit),'r')
                bar(2:length(shift_corrs)+1,shift_corrs(2,:))
                plot([0 length(shift_corrs)+1],[corr_95_even_odds corr_95_even_odds],'k--')
                hold off
                xlabel('Shift (ms)')
                ylabel('Spatial corr even/odd')
                set(gca,'Xtick',1:3:length(shift_corrs(2,:)))
                set(gca,'XtickLabel',[{'None'} num2cell(shifts(1:3:end))])
                set(gca,'XTickLabelRotation',45)
                xlim([0 length(shift_bits)+2])
                
                if all(spatial_info.spatialstability_even_odd(1,unit) > shift_corrs(2,:))
                    title('Spatial Corr E/O, Best: No Shft')
                else
                    title(['Spatial Corr E/O, Best: ' num2str(shifts(shift_corrs(2,:) == max(shift_corrs(2,:))))])
                end
                
                
                %choose best average correlation....cuz I don't know what's
                %best for now
                
                best_shift = find(shift_bits == max(shift_bits));
                
                %shift eye movements relative to spikes
                shift_eye = eye_time_shift(eyepos{unit},shifts(best_shift));
                
                %draw raw spikes on all images all locations
                make_spike_jittered_plot(shift_eye,spike_times{unit},[3 3],5)
                title(['Best Avg Shift: ' num2str(shifts(best_shift)) 'ms'])
                set(gca,'Xcolor','w')
                set(gca,'Ycolor','w')
                
                
                [ratemaps,timemaps] = get_firing_rate_map({shift_eye,spike_times{unit}},imageX,imageY,binsize,H,Fs,'all');
                
                subplot(3,3,8)
                h = imagesc(ratemaps);
                set(h,'alphadata',~isnan(timemaps));
                axis off
                axis equal
                caxis(clims);
                colormap('jet')
                title(['Best Shift: ' num2str(shifts(best_shift)) ', max ' num2str(max(ratemaps(:)),3) ' Hz'])
                
                
                
                trough_to_best = [trough_to_best best_shift-fix_delay];
                
                if multiunit(unit)
                    multi_str = 'Multiunit ';
                else
                    multi_str = ' ';
                end
                subtitle(['Spatial Plots' multi_str unit_names{unit}]);
                
                if col == 1 %eye movement modulated eurons
                    if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                            && (spatial_info.spatialstability_even_odd_prctile(2,unit) > 95) ... %spatial consistency
                            && (spatial_info.spatialstability_halves_prctile(2,unit) > 95) %spatial stability
                        save_and_close_fig(figure_dir3,[task_file(1:end-11) '-' unit_names{unit} '_Place_Cell_TimeShift'])
                    else
                        save_and_close_fig(figure_dir,[task_file(1:end-11) '-' unit_names{unit} '_Place_Cell_TimeShift'])
                    end
                else%not eye movement modulated
                    if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                            && (spatial_info.spatialstability_even_odd_prctile(2,unit) > 95) ... %spatial consistency
                            && (spatial_info.spatialstability_halves_prctile(2,unit) > 95) %spatial stability
                        save_and_close_fig(figure_dir3,[task_file(1:end-11) '-' unit_names{unit} '_Place_Cell_TimeShift'])
                    else
                        save_and_close_fig(figure_dir2,[task_file(1:end-11) '-' unit_names{unit} '_Place_Cell_TimeShift'])
                    end
                end
            end
        end
    end
end
