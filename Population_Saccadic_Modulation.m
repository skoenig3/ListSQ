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
locations = [2:17];
Fs = 1000;%Hz
smval = 60;
imageX = 800;
imageY = 600;
all_unit_positions = zeros(1,length(locations));
saccade_cell_poisitions = zeros(1,length(locations));

avg_saccade_firing_rates = cell(1,3);%colums: sig not sig
avg_fixation_firing_rates = cell(1,3);%colums: sig not sig
avg_fixation_firing_rates_limited =[];%for significant ones limted to duration of prececeding saccade+fixation

t = -500:499;
avg_saccade_firing_rates_place_cells = cell(1,3);
avg_fixation_firing_rates_place_cells = cell(1,3);

sig_count = zeros(2,3);
place_sig_count = zeros(2,3);

max_firing = NaN(347,2);
un = 1;
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Population Figures\Place Cell Eye Movements\';
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
            if ~isempty(fixation_locked_firing{unit})
                if (fixation_info.rate(unit) >  fixation_info_95.rate(unit)) &&  (saccade_info.rate(unit) >  saccade_info_95.rate(unit))
                    col = 1;%significant modulation
                    sig_count(1,1) = sig_count(1,1)+1;
                elseif (fixation_info2.rate(unit) >  fixation_info2_95.rate(unit)) &&  (saccade_info2.rate(unit) >  saccade_info2_95.rate(unit))
                    col = 1;%significant modulation
                    sig_count(2,1) = sig_count(2,1)+1;
                elseif (fixation_info.rate(unit) >  fixation_info_95.rate(unit)) || (saccade_info.rate(unit) >  saccade_info_95.rate(unit))
                    sig_count(1,3) = sig_count(1,3)+1;
                    col = 3;
                elseif (fixation_info2.rate(unit) >  fixation_info2_95.rate(unit)) || (saccade_info2.rate(unit) >  saccade_info2_95.rate(unit))
                    col = 3;
                    sig_count(2,3) = sig_count(2,3)+1;
                else
                    sig_count(1,2) = sig_count(1,2)+1;
                    col = 2;%not saccade modulated
                end
                
                if ((spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.shuffled_spatialstability_prctile(unit) > 95)) ...
                        
                if (fixation_info.rate(unit) >  fixation_info_95.rate(unit)) &&  (saccade_info.rate(unit) >  saccade_info_95.rate(unit))
                    place_sig_count(1,1) = place_sig_count(1,1)+1;
                elseif (fixation_info2.rate(unit) >  fixation_info2_95.rate(unit)) &&  (saccade_info2.rate(unit) >  saccade_info2_95.rate(unit))
                    place_sig_count(2,1) = place_sig_count(2,1)+1;
                elseif (fixation_info.rate(unit) >  fixation_info_95.rate(unit)) || (saccade_info.rate(unit) >  saccade_info_95.rate(unit))
                    place_sig_count(1,3) = place_sig_count(1,3)+1;
                elseif (fixation_info2.rate(unit) >  fixation_info2_95.rate(unit)) || (saccade_info2.rate(unit) >  saccade_info2_95.rate(unit))
                    place_sig_count(2,3) = place_sig_count(2,3)+1;
                else
                    place_sig_count(1,2) = place_sig_count(1,2)+1;
                end
                end
                
                %                     figure
                %                     %draw raw spikes on all images all locations
                %                     make_spike_jittered_plot(eyepos{unit},spike_times{unit},[2 2],1)
                %                     title(['Raw: All images, n = ' num2str(size(spike_times{unit},1))])
                %                     set(gca,'Xcolor','w')
                %                     set(gca,'Ycolor','w')
                %
                %                     subplot(2,2,2)
                %                     %Calculated smoothed firing rate for all images
                %                     filtered_time = filter_time(eyepos{unit},imageX,imageY,Fs,binsize,H);
                %                     filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
                %                     filtered_space = filter_space(eyepos{unit},spike_times{unit},imageX,imageY,binsize,H);
                %                     firing_rate = filtered_space./filtered_time;
                %
                %                     h = imagesc(firing_rate);
                %                     set(h,'alphadata',~isnan(filtered_time));
                %                     axis off
                %                     axis equal
                %                     fr = sort(firing_rate(1:end));
                %                     fr(isnan(fr)) = [];
                %                     clims(:,1) = caxis;
                %                     if length(fr) > 20
                %                         clims(2,1) = fr(round(0.99*length(fr)));% the ~99%-tile
                %                     end
                %                     colorbar
                %                     colormap('jet')
                %
                %                     title(sprintf(['All images, peak rate = ' num2str(max(fr),3) ...
                %                         ' Hz \n Bit: ' num2str(spatial_info.shuffled_rate_prctile(unit),3) '%%,'...
                %                         ' r = ' num2str(spatial_info.spatialstability(unit),2) ' ' ...
                %                         num2str(spatial_info.shuffled_spatialstability_prctile(unit),2) '%%']))
                %
                %
                %                     subplot(2,2,4)
                %                     %plot raster over time by spikes/saccade epoch aka by average firing rate
                saccade_firing = saccade_locked_firing{unit}(saccade_information{unit}(:,4) > 2*twin,:);
                fixation_firing = fixation_locked_firing{unit}(fixation_information{unit}(:,4) > 2*twin,:);
                %                     f1 = sum( saccade_firing,2);
                %                     [~,fi] = sort(f1);
                %                     [trial,time] = find(saccade_firing(fi,:) == 1);
                %                     plot(time-500,(trial),'.k')
                %                     ylim([0 max(trial)])
                %                     set(gca,'Xtick',[-500 -250 0 250 500])
                %                     xlabel('Time from Saccade Start')
                %                     ylabel('Ranked Spike Count')
                %
                saccade_firing(nansum(saccade_firing,2) == 0,:) = [];%remove saccades without spikes
                fixation_firing(nansum(fixation_firing,2) == 0,:) = [];%remove fixations without spikes
                
                if col == 1
                    info = fixation_information{unit};
                    %want to look only looking at window around fixation/fixation
                    limited_firing = NaN(size(fixation_firing));
                    if sum(isnan(info(:,6))) ~= 0
                        info(isnan(info(:,6)),6) = round(nanmean(info(:,6)));%go back on average if no saccade before
                    end
                    info(info(:,5) > twin,5) = twin; %fixations longer than twin should be set to twin
                    for f = 1:size(limited_firing,1)
                        ind = twin-info(f,6):twin+info(f,5);
                        limited_firing(f,ind) = fixation_firing(f,ind);
                    end
                    
                    %         %make raster plot for spikes limited time period of 1 fixation around fixation
                    %         subplot(3,3,9)
                    %         [trial,time] = find(limited_firing == 1);
                    %         plot(time-500,trial,'.k')
                    %         hold off
                    %         ylim([0 max(trial)])
                    %         set(gca,'Xtick',[-500 -250 0 250 500])
                    %         ylabel('Occurence #')
                    %         xlabel('Time from fixation Start')
                    %         title('Time Period 1 saccade before fixation')
                    
                    %plot firing rate curve over time for spikes limited time period of 1 fixation around fixation
                    num_not_nans = sum(~isnan(limited_firing));%will be fewer sikes
                    not_enough_samples = find(num_not_nans < .5*size(limited_firing,1));
                    
                    [firing_rate,~]= nandens(limited_firing,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate(:,not_enough_samples) = NaN;  %remove any time points with less than 1/4 of the samples on avg ~< 500 trials
                    firing_rate = firing_rate-nanmean(firing_rate);
                    firing_rate = firing_rate/max(abs(firing_rate));
                    
                    avg_fixation_firing_rates_limited = [avg_fixation_firing_rates_limited; firing_rate];
                end
                
                
                %
                %                     subplot(2,2,3)
                %                     hold on
                %                     dofill(t,fixation_locked_firing{unit}(fixation_information{unit}(:,4) > 2*twin,:),'blue',1,smval);
                %                     dofill(t,saccade_locked_firing{unit}(saccade_information{unit}(:,4) > 2*twin,:),'red',1,smval);
                %                     dofill(t,fixation_firing,'green',1,smval);
                %                     dofill(t,saccade_firing,'magenta',1,smval);
                %                     yl = ylim;
                %                     plot([0 0],[yl(1) yl(2)],'k--')
                %                     hold off
                %                     xlabel('Time from Eye Movement Start')
                %                     ylabel('Firing Rate (Hz)')
                %                     legend('All Fixations','All Saccades','+Fixations','+Saccades')
                %
                %
                %                     fix_bit_prctile = 100*sum(fixation_info.rate(unit) >  fixation_shuffled_info.rate{unit})/length(fixation_shuffled_info.rate{unit});
                %                     sac_bit_prctile = 100*sum(saccade_info.rate(unit) >  saccade_shuffled_info.rate{unit})/length(saccade_shuffled_info.rate{unit});
                %                     %                     fix_r_prctile = 100*sum(fixation_info.temporalstability(unit) > fixation_shuffled_info.temporalstability{unit})/length(fixation_shuffled_info.temporalstability{unit});
                %                     %                     sac_r_prctile = 100*sum(saccade_info.temporalstability(unit) > saccade_shuffled_info.temporalstability{unit})/length(saccade_shuffled_info.temporalstability{unit});
                %                     fix_bit_prctile2 = 100*sum(fixation_info2.rate(unit) >  fixation_shuffled_info2.rate{unit})/length(fixation_shuffled_info2.rate{unit});
                %                     sac_bit_prctile2 = 100*sum(saccade_info2.rate(unit) >  saccade_shuffled_info2.rate{unit})/length(saccade_shuffled_info2.rate{unit});
                %
                % %                     if fix_bit_prctile > 95 && sac_bit_prctile > 95
                % %                        if fix_r_prctile < 95 || sac_r_prctile < 95
                % %                            disp('now')
                % %                        end
                % %                     end
                %
                %                     title(sprintf(['fix_{bit} = ' num2str(fix_bit_prctile,3)...
                %                         '%% sac_{bit} = ' num2str(sac_bit_prctile,3) '%% ' ...
                %                         'fix+_{bit} = ' num2str(fix_bit_prctile2,3)...
                %                         '%% sac+_{bit} = ' num2str(sac_bit_prctile2,3) '%%']));
                %                         %%fix_r = ' ...
                % %                         num2str(fixation_info.temporalstability(unit),2)  ...
                % %                         ' ' num2str(fix_r_prctile,3) '%% sac_r = ' ...
                % %                         num2str(saccade_info.temporalstability(unit),2) ...
                % %                         ' ' num2str(sac_r_prctile,3) '%%']))
                %
                %                     if multiunit(unit)
                %                         multi_str = 'Multiunit ';
                %                     else
                %                         multi_str = ' ';
                %                     end
                %                     subtitle(['Spatial Plots' multi_str unit_names{unit}]);
                %
                %                     if col == 1;
                %                         disp('now')
                %                     end
                %
                %                     save_and_close_fig(figure_dir,[task_file(1:end-11) '-' unit_names{unit} '_Place_Cell_EyeMovements'])
                [saccade_firing,~]= nandens(saccade_firing,smval,'gauss',Fs,'nanflt');%nandens made by nathan killian
                
                saccade_firing = saccade_firing-mean(saccade_firing);
                max_firing(un,1) = max(saccade_firing);
                
                saccade_firing = saccade_firing/max(abs(saccade_firing));
                avg_saccade_firing_rates_place_cells{col} = [avg_saccade_firing_rates_place_cells{col};  saccade_firing];
                [fixation_firing,~]= nandens(fixation_firing,smval,'gauss',Fs,'nanflt');%nandens made by nathan killian
                
                fixation_firing = fixation_firing-mean(fixation_firing);
                max_firing(un,2) = max(fixation_firing);
                un = un+1;
                fixation_firing = fixation_firing/max(abs(fixation_firing));
                avg_fixation_firing_rates_place_cells{col} = [avg_fixation_firing_rates_place_cells{col};  fixation_firing];
                
                %for fixations
                [firing_rate,~]= nandens(fixation_locked_firing{unit}(fixation_information{unit}(:,4) > 2*twin,:)...
                    ,smval,'gauss',Fs,'nanflt');%nandens made by nathan killian
                firing_rate = firing_rate-mean(firing_rate); %zero/remove average firing rate
                firing_rate = firing_rate/max(abs(firing_rate)); %normalize to 1
                avg_fixation_firing_rates{col} = [ avg_fixation_firing_rates{col}; firing_rate];
                
                %for saccades
                [firing_rate,~]= nandens(saccade_locked_firing{unit}(saccade_information{unit}(:,4) > 2*twin,:)...
                    ,smval,'gauss',Fs,'nanflt');%nandens made by nathan killian
                
                firing_rate = firing_rate-mean(firing_rate);  %zero/remove average firing rate
                firing_rate = firing_rate/max(abs(firing_rate)); %normalize to 1
                avg_saccade_firing_rates{col} = [ avg_saccade_firing_rates{col}; firing_rate];
                
                
                AP_location = round(chamber_zero(1)+ session_data{session}.location(1));
                location_index = find(locations == AP_location);
                
                all_unit_positions(location_index) = all_unit_positions(location_index)+1;
                if col == 1
                    saccade_cell_poisitions(location_index) = saccade_cell_poisitions(location_index)+1;
                end
            end
        end
    end
end
%%
figure
plot(mean( avg_fixation_firing_rates{1}))
hold on
plot(mean(avg_saccade_firing_rates{1}),'r')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Saccade  Locked Activity---%
[U,S,V] = pca(avg_saccade_firing_rates{1});
T = kmeans(U,4);
t = -500:499;
figure
hold all
for i = 1:max(T)
    if sum(T == i) > 1
        plot(t,mean(avg_saccade_firing_rates{1}(T == i,:)))
    else
        plot(t,avg_saccade_firing_rates{1}(T == i,:));
    end
end
plot([0 0],[-1 1],'k--')
plot([-500 500],[0 0],'k')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Fixation  Locked Activity---%
[U,S,V] = pca(avg_fixation_firing_rates{1});
T = kmeans(U,4);
t = -500:499;
figure
hold all
for i = 1:max(T)
    if sum(T == i) > 1
        plot(t,mean(avg_fixation_firing_rates{1}(T == i,:)))
    else
        plot(t,avg_fixation_firing_rates{1}(T == i,:));
    end
end
plot([0 0],[-1 1],'k--')
plot([-500 500],[0 0],'k')

%%
plot(mean(avg_saccade_firing_rates{1,1}))
%%
figure
[m,i] = max(avg_saccade_firing_rates{1}(:,350:850),[],2);
[mm,ii] = sort(i);
imagesc([-500:499],[1:size(avg_saccade_firing_rates{1},1)],avg_saccade_firing_rates{1}(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(avg_saccade_firing_rates{1},1)],'w--');
hold off
xlabel('Time from Saccade Start')
ylabel('Neuron #')

%% All Cells
figure
[m,i] = max(avg_fixation_firing_rates{1}(:,350:850),[],2);
[mm,ii] = sort(i);
imagesc([-500:499],[1:size(avg_fixation_firing_rates{1},1)],avg_fixation_firing_rates{1}(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(avg_fixation_firing_rates{1},1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')

%% Place Cells
figure
[m,i] = max(avg_fixation_firing_rates_place_cells{1}(:,350:850),[],2);
[mm,ii] = sort(i);
imagesc([-500:499],[1:size(avg_fixation_firing_rates_place_cells{1},1)],avg_fixation_firing_rates_place_cells{1}(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(avg_fixation_firing_rates{1},1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')

figure
[m,i] = max(avg_saccade_firing_rates_place_cells{1}(:,350:850),[],2);
[mm,ii] = sort(i);
imagesc([-500:499],[1:size(avg_saccade_firing_rates_place_cells{1},1)],avg_saccade_firing_rates_place_cells{1}(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(avg_saccade_firing_rates{1},1)],'w--');
hold off
xlabel('Time from saccade Start')
ylabel('Neuron #')
%%
figure
[m,i] = max(avg_fixation_firing_rates_place_cells{2}(:,250:850),[],2);
[mm,ii] = sort(i);
imagesc([-500:499],[1:size(avg_fixation_firing_rates_place_cells{2},1)],avg_fixation_firing_rates_place_cells{2}(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(avg_fixation_firing_rates{1},1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')

figure
[m,i] = max(avg_saccade_firing_rates_place_cells{2},[],2);
[mm,ii] = sort(i);
imagesc([-500:499],[1:size(avg_saccade_firing_rates_place_cells{2},1)],avg_saccade_firing_rates_place_cells{2}(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(avg_saccade_firing_rates{1},1)],'w--');
hold off
xlabel('Time from saccade Start')
ylabel('Neuron #')
%%
figure
[m,i] = nanmax(avg_fixation_firing_rates_limited,[],2);
[mm,ii] = sort(i);
imagesc([-500:499],[1:size(avg_fixation_firing_rates_limited,1)],avg_fixation_firing_rates_limited(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(avg_fixation_firing_rates{1},1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
