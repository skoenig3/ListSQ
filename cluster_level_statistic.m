function [cluster_level_percentile,significant_indeces] = cluster_level_statistic(observed_diff,shuffled_diff_curves,tails,min_cluster_width)
%written by Seth Konig January 18, 2017
%code implements method described by Maris and Oostenveld, 2007 to correct
%for multiple comparisons during time series analysis between 2 difference
%conditions.
% min_cluster_width should be scaled with smoothing std, if unknown try 50

%rechecked for bugs Feburary 3, 2017

alpha = 5; %significance level, p < 0.05
%min_cluster_width = 50; %minimum size of contiguous significant time points
numshuffs = size(shuffled_diff_curves,1);

if tails == 1 %1 tail test under assumption that 1 curve will always be higher than the other
    
    %---Determine Upper Bounds of Shuffled Data---%
    upper_shuffled = prctile(shuffled_diff_curves,100-alpha,1);%95% of shuffled differences
    
    %---Determine time points that are significantly different---%
    sig_times = zeros(1,length(observed_diff));
    sig_times (observed_diff > upper_shuffled) = 2; %positive differences
    
    %---Find Contigous Time Points that Show Significat Differences---%
    gaps = findgaps(find(sig_times)); %find contigous significant regions
    if ~isempty(gaps)
        cluster_values = NaN(1,size(gaps,1)); %observed sum of diff in cluster
        shuffled_cluster_values = cell(1,size(gaps,1)); %shuffled sum of diff in cluster
        cluster_widths = NaN(1,size(gaps,1));%size of clusters
        for g = 1:size(gaps,1)
            gp = gaps(g,:);
            gp(gp == 0) = [];
            if length(gp) >= min_cluster_width
                cluster_values(g) = sum(observed_diff(gp)); %sum of the diff in the curve
                shuffled_cluster_values{g} = sum(shuffled_diff_curves(:,gp),2); %sum of the shuffled diff in the curve
                cluster_widths(g) = length(gp); %cluster width
            end
        end
    else %no differences observed so exit function
        cluster_level_percentile = zeros(1,length(observed_diff));
        significant_indeces = zeros(1,length(observed_diff));%no significant differences
        return
    end
    
    if all(isnan(cluster_values)) %significant clusters were too short
        cluster_level_percentile = zeros(1,length(observed_diff));
        significant_indeces = zeros(1,length(observed_diff));%no significant differences
        return
    end
    
    %---Compare Cluster-level Statistic to Cluster-level Distribution of Biggest Difference---%
    biggest_difference = max(abs(cluster_values))== abs(cluster_values); %find cluster with biggest avg/sum difference
    new_shuffled_distribution = shuffled_cluster_values{biggest_difference}; %distribution of biggest cluster
    cluster_level_percentile = NaN(1,size(gaps,1)); %percentile at cluster level
    significant_indeces = zeros(1,length(observed_diff));%FWER corrected significant differences
    for g = 1:size(gaps,1)
        %cluster value should be positive, average diff should be 0
        cluster_level_percentile(g) = 100*sum(cluster_values(g) > new_shuffled_distribution)/numshuffs;
        if cluster_level_percentile(g) > (100-alpha) %if pass criterion set as significant indeces
            gp = gaps(g,:);
            gp(gp == 0) = [];
            significant_indeces(gp) = 2;
        end
    end
    
elseif tails == 2 %two tailed test assuming difference is 0
    
    %---Determine Upper and Lower Bounds of Shuffled Data---%
    upper_shuffled = prctile(shuffled_diff_curves,100-alpha/2,1);%97.5% of shuffled differences
    lower_shuffled = prctile(shuffled_diff_curves,alpha/2,1);%2.5% of shuffled differences
    
    %---Determine time points that are significantly different---%
    sig_times = zeros(1,length(observed_diff));
    sig_times ((observed_diff > upper_shuffled) | ... %positive differences
        (observed_diff <lower_shuffled)) = 2; %negative differences
    
    %---Find Contigous Time Points that Show Significat Differences---%
    gaps = findgaps(find(sig_times)); %find contigous significant regions
    if ~isempty(gaps)
        cluster_values = NaN(1,size(gaps,1)); %observed sum of diff in cluster
        shuffled_cluster_values = cell(1,size(gaps,1)); %shuffled sum of diff in cluster
        cluster_widths = NaN(1,size(gaps,1));%size of clusters
        for g = 1:size(gaps,1)
            gp = gaps(g,:);
            gp(gp == 0) = [];
            if length(gp) >= min_cluster_width
                cluster_values(g) = sum(observed_diff(gp)); %sum of the diff in the curve
                shuffled_cluster_values{g} = sum(shuffled_diff_curves(:,gp),2); %sum of the shuffled diff in the curve
                cluster_widths(g) = length(gp);%cluster width
            end
        end
    else %no differences observed so exit function
        cluster_level_percentile = zeros(1,length(observed_diff));
        significant_indeces = zeros(1,length(observed_diff));%no significant differences
        return
    end
    
    if all(isnan(cluster_values)) %significant clusters were too short
        cluster_level_percentile = zeros(1,length(observed_diff));
        significant_indeces = zeros(1,length(observed_diff));%no significant differences
        return
    end
    
    %---Compare Cluster-level Statistic to Cluster-level Distribution of Biggest Difference---%
    biggest_difference = max(abs(cluster_values))== abs(cluster_values); %find cluster with biggest avg/sum difference
    new_shuffled_distribution = shuffled_cluster_values{biggest_difference}; %distribution of biggest cluster
    cluster_level_percentile = NaN(1,size(gaps,1)); %percentile at cluster level
    significant_indeces = zeros(1,length(observed_diff));%FWER corrected significant differences
    for g = 1:size(gaps,1)
        if cluster_values(g) < 0 %negative value, average diff should be 0
            cluster_level_percentile(g) = 100*sum(cluster_values(g) < new_shuffled_distribution)/numshuffs;
        elseif cluster_values(g) > 0 %positive value, average diff should be 0
            cluster_level_percentile(g) = 100*sum(cluster_values(g) > new_shuffled_distribution)/numshuffs;
        end
        if cluster_level_percentile(g) > (100-alpha/2) %if pass criterion set as significant indeces
            gp = gaps(g,:);
            gp(gp == 0) = [];
            significant_indeces(gp) = 2;
        end
    end
else
    error('need to specify number of tails')
end
