% % View Cell Psuedo-Population Decoding Analysis
% % written by Seth Konig 6/14/16
% %code only analyzes units with significant (95+ percentile) skaggs information
% %score and Miriam's spatial stability score.
% 
% all_firing_rate_maps = cell(1,50);
% rate_map_ind = 1;
% 
% all_time = [];
% all_n = [];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%---Set important task parameters---%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% twin = 500;% here how much after image onset to ignore data
% imageX = 800;
% imageY = 600;
% img_on_code = 23;
% img_off_code = 24;
% ITIstart_code = 15;
% smval = 60;
% Fs = 1000;
% min_bin_dur = 0.100; %minimum of 500 ms in each bin to use so no outlier
% smval = 60;%gaussian 1/2 width for smoothing
% task = 'ListSQ';
% 
% for monk = 1:2
%     if monk == 1
%         %---load Vivian's data---%
%         excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
%         excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
%         data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
%         load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
%     else
%         %---Load Tobiis data---%
%         excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
%         excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
%         data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
%         load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
%         session_data(end-2:end) = [];%last file doesn't have strobe signal working on importing the data
%     end
%     
%     for sess = 1:length(session_data)
%         
%         task_file = get_task_data(session_data{sess},task);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%---import task and unit data---%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if isempty(task_file)
%             warning('No file could be found for specificed task. Exiting function...')
%             continue
%         end
%         load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat']);
%         load([data_dir task_file(1:end-11) '-preprocessed.mat']);
%         absolute_fixationstats = fixationstats;
%         absolute_cfg = cfg;
%         
%         [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
%             = get_task_data(session_data{sess},task);
%         [multiunit,~,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
%             multiunit,unit_confidence,sorting_quality);
%       
%         %filter parameters
%         filter_size = filter_width*10;
%         H = fspecial('gaussian',filter_size,filter_width);
%         
%         for unit = 1:num_units
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %%%---Calculate Firing Rate Maps---%%%
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%             %draw spikes in locations above threshold (currently > 20% of max)
%             which_condition_skaggs = find(spatial_info.shuffled_spatialstability_prctile(:,unit) >= 95);
%             which_condition_stable = find(spatial_info.shuffled_spatialstability_prctile(:,unit) >= 95);
%             which_condition_both = which_condition_skaggs & which_condition_stable;
%             if any(which_condition_skaggs == 3) && any(which_condition_stable == 3);
%                 condition = 3;                
%             else
%                continue 
%             end
%             
%             if condition == 3 %if condition is 3 then all images
%                 filtered_time = filter_time(eyepos{unit},imageX,imageY,Fs,binsize,H);
%                 filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
%                 filtered_space = filter_space(eyepos{unit},spike_times{unit},imageX,imageY,binsize,H);
%                 
%                 all_firing_rate_maps{rate_map_ind} = filtered_space./filtered_time;
%                 rate_map_ind = rate_map_ind+1;
%                 if isempty(all_time)
%                     all_n = zeros(size(filtered_time));
%                     all_n(~isnan(filtered_time)) = 1;
%                     filtered_time(isnan(filtered_time)) = 0;
%                     all_time = filtered_time;
%                     
%                 else
%                     all_n(~isnan(filtered_time)) =  all_n(~isnan(filtered_time))+1;
%                     filtered_time(isnan(filtered_time)) = 0;
%                     all_time = all_time+filtered_time;
%                 end
%             end
%         end
%     end
% end
% %%
% save('ViewCellPopulation','all_time','all_n','all_firing_rate_maps')
%%

average_prob = zeros(size(all_firing_rate_maps{1}));
firing_rate_prob = [];
for map = 1:length(all_firing_rate_maps)
    this_map = all_firing_rate_maps{map}./nansum(nansum(all_firing_rate_maps{map}));
    this_map(isnan(this_map)) = 0;
    average_prob = average_prob+this_map;
    firing_rate_prob = cat(3,firing_rate_prob,this_map);
    this_map = all_firing_rate_maps{map}./nansum(nansum(all_firing_rate_maps{map}));
end

figure
imagesc(100*average_prob./all_n)
axis equal
box off
axis off
title('Average Probability')
colormap('jet')
colorbar
%%
distance_percentile_5 = NaN(size(average_prob));
corr_percentile_95 = NaN(size(average_prob));
for row = 1:size(average_prob,1)
    for col = 1:size(average_prob,2)
        distance_matrix = NaN(size(average_prob)); 
        corr_matrix = NaN(size(average_prob)); 
        for rrow = 1:size(average_prob,1)
            for ccol = 1:size(average_prob,2)
                a = squeeze(firing_rate_prob(row,col,:));
                b = squeeze(firing_rate_prob(rrow,ccol,:));
                b(isnan(a)) = [];
                a(isnan(a)) = [];
                a(isnan(b)) = [];
                b(isnan(b)) = [];
                distance_matrix(rrow,ccol) = sqrt(sum((a-b).^2));
                corr_matrix(rrow,ccol) = corr(a,b);
            end
        end
        distance_prct5 = prctile(distance_matrix(1:end),5);
        [i,j] = find(distance_matrix <= distance_prct5);
        dist = sqrt((i-row).^2+(j-col).^2);
        distance_percentile_5(row,col) = max(dist);
        corr_prct95 = prctile(corr_matrix(1:end),95);
        [i,j] = find(corr_matrix >= corr_prct95);
        dist = sqrt((i-row).^2+(j-col).^2);
        corr_percentile_95(row,col) = max(dist);
    end
end
%%
figure
subplot(1,2,1)
imagesc(distance_percentile_5)
axis equal
box off
axis off
title('5% of Distance from Origin')
colormap('jet')
colorbar

subplot(1,2,2)
imagesc(corr_percentile_95)
axis equal
box off
axis off
title('95% of Correlation from Origin')
colormap('jet')
colorbar
%%

input = 