%get Cell counts SDK written 4/25/16

%monkey = 'Vivian';
monkey = 'Tobii';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Read in Excel Sheet for Session data---%%%
%only need to run when somethings changed or sessions have been added
task = 'ListSQ';
if strcmpi(monkey,'Vivian')
    data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
    load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
    figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
elseif strcmpi(monkey,'Tobii')
    data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
    load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
    session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
end

%---Variables for General Spike info---%
%column 1 single units, 2 column 2 multiunits
total_recorded = [0 0];
total_used = [0 0];
num_blocks = zeros(3,2); %row 1, 1+ blocks, row 2 2+ blocks, row 3 3+ blocks
reason_unused = zeros(4,2); %row 1 insufficeint spike count, row 2,
% not stable for sufficeint trials, row 3 poorly isolated/not confident
% unit is real, row 4 sum of rows 1-3
barely_active_count = [0 0];%column 1 total average firing rate < 1 but not all, total with sig skagg and stability socres
low_sig = 0;

%---Variables for Spatial analysis---%
spatial_total_used = [0 0]; %all usable units should be same as total_used
spatial_total_low_firing = [0 0]; %units with < 1 Hz peak firing rate under all conditions
spatial_skaggs = [0 0]; %units > 95% shuffled for skaggs infomration score
spatial_stability = [0 0];%units > 95% shuffled for Miriam's spatial stability
spatial_skaggs_stability = [0 0]; %units with > 95% shuffled for skaggs AND spatial stability

count2 = 0;

for sess = 1:length(session_data)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import General Data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get important task related data
    [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
        get_task_data(session_data{sess},task);
    if isempty(task_file)
        disp('No file could be found for specificed task. Exiting function...')
        continue
    end
    try
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials','stability_attribute','cfg','hdr','data');
    catch
        disp(['Cant find ' task_file(1:end-11)])
        continue
    end
    
    %get unit data
    [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
        multiunit,unit_confidence,sorting_quality);
    
    
    %these are the absolute minimum data required to do data analysis may want
    %to be more strigent later but not worth doing analysis (especially
    %shuffling) on such few trials for these neurons
    if str2double(task_file(3:8)) < 140805 %7 second viewing with fewere sequence trials
        minimum_trials_1 = 101; %for session start minimum number of trials: includes fam block and 16 novel/repeat images + sequence trials
        minimum_trials_2 =  80;%for rest of session: includes 16 novel/repeat images + sequence trials
    else
        minimum_trials_1 = 117; %for session start minimum number of trials: includes fam block and 16 novel/repeat images + sequence trials
        minimum_trials_2 =  96;%for rest of session: includes 16 novel/repeat images + sequence trials
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Get General Unit Attributes---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for unit = 1:num_units
        if multiunit(unit) %then multiunit
            su_mu_ind = 2;
        else %single unit
            su_mu_ind = 1;
        end
        
        total_recorded(su_mu_ind) =   total_recorded(su_mu_ind)+1;
        
        if stability_attribute(unit) == 1 %stable usable cell
            %double check unit is good
            if unit_stats{2,unit} < 0.8 || unit_stats{3,unit} < 2.5
                error('Poor quality unit, it should not have been used')
            end
            total_used(su_mu_ind) = total_used(su_mu_ind)+1;
            
            %---determine number of blocks unit was stable for
            start_end = valid_trials(:,unit);
            if isnan(start_end(1))
                start_end(1) = 1;
            end
            if isnan(start_end(2))
                start_end(2) = length(cfg.trl);
            end
            start_end(start_end == 0) = 1;
            min_trial = cfg.trl(start_end(1)).cnd-1000; %get the first condition number
            max_trial = cfg.trl(start_end(2)).cnd-1000; %get the last condition number
            
            if max_trial < minimum_trials_1
                error('Unit should have been stable for more trials to be processed')
            end
            
            if min_trial < 22 %includes fam block
                min_trial = 22; %remove then count from there
            end
            num_blks = floor((max_trial-min_trial+1)/minimum_trials_2);
            
            if num_blks >= 1
                num_blocks(1,su_mu_ind) = num_blocks(1,su_mu_ind)+1;
            end
            if num_blks >= 2
                num_blocks(2,su_mu_ind) = num_blocks(2,su_mu_ind)+1;
            end
            if num_blks >= 3
                num_blocks(3,su_mu_ind) = num_blocks(3,su_mu_ind)+1;
            end
            
        elseif stability_attribute(unit) == 2 %insufficeint spike count
            reason_unused(1,su_mu_ind) = reason_unused(1,su_mu_ind)+1;
        elseif stability_attribute(unit) == 3 %not stable for sufficient trials
            reason_unused(2,su_mu_ind) = reason_unused(2,su_mu_ind)+1;
        elseif stability_attribute(unit) == 4 %poorly isolated/not confident unit is real
            if unit_stats{2,unit} >= 0.8 && unit_stats{3,unit} >= 3%should have been used but wasn't
                error('Unit should have been good enought to use')
            end
            reason_unused(3,su_mu_ind) = reason_unused(3,su_mu_ind)+1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import Spatial Analysis Data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get important task related data
    try
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'peak_firing_rate','spatial_info');
    catch
        if num_units == 0 %then no units so shouldn't have anything
            continue
        else
            error('should have run spatial analysis had units')
        end
    end
    
    if num_units ~= size(peak_firing_rate,2)
        error('# of units mismatch: Wrong file imported or data not processed properly')
    end
    
    for unit = 1:num_units
        %unstable/usable units
        if stability_attribute(unit) ~=  1 %not used unit
            if any(~isnan(peak_firing_rate(:,unit)))
                error('Unit was stable and should have been processed')
            end
            continue %move to the next unit
        end
        
        %usable units
        if multiunit(unit) %then multiunit
            su_mu_ind = 2;
        else %single unit
            su_mu_ind = 1;
        end
        
%         if all(isnan(peak_firing_rate(:,unit)))
%             error('Stable unit but was not procssed for spatial analyis')
%         end
        if any(isnan(spatial_info.rate(:,unit)) & peak_firing_rate(:,unit) > 1)
            error('Firing rate for unit is > 1 Hz but not processed')
        end
        
        %---determine number of blocks unit was stable for
        start_end = valid_trials(:,unit);
        if isnan(start_end(1))
            start_end(1) = 1;
        end
        if isnan(start_end(2))
            start_end(2) = length(cfg.trl);
        end
        start_end(start_end == 0) = 1;
        min_trial = cfg.trl(start_end(1)).cnd-1000; %get the first condition number
        max_trial = cfg.trl(start_end(2)).cnd-1000; %get the last condition number
        
        if max_trial < minimum_trials_1
            error('Unit should have been stable for more trials to be processed')
        end
        
        if min_trial < 22 %includes fam block
            min_trial = 22; %remove then count from there
        end
        num_blks = floor((max_trial-min_trial+1)/minimum_trials_2);
        
        if num_blks < 2
            continue
        end
        
        spatial_total_used(su_mu_ind) = spatial_total_used(su_mu_ind)+1;
        if all(peak_firing_rate(:,unit) < 1)
            spatial_total_low_firing(su_mu_ind) =  spatial_total_low_firing(su_mu_ind)+1;
        end
        
        go_no_go = 0;
        if any(spatial_info.shuffled_rate_prctile(:,unit) > 95) %significant skaggs score
            spatial_skaggs(su_mu_ind) = spatial_skaggs(su_mu_ind)+1;
        end
        if any(spatial_info.shuffled_spatialstability_prctile(:,unit) > 95) %significant spatial stability score
            spatial_stability(su_mu_ind) = spatial_stability(su_mu_ind)+1;
        end
        if any((spatial_info.shuffled_rate_prctile(:,unit) > 95) & (spatial_info.shuffled_spatialstability_prctile(:,unit) > 95)) %both skagg and spatially stable
            spatial_skaggs_stability(su_mu_ind) = spatial_skaggs_stability(su_mu_ind)+1;
            go_no_go = 1;
        end
        
        if all(peak_firing_rate(:,unit) < 1) && (spatial_info.shuffled_rate_prctile(:,unit) > 95) && (spatial_info.shuffled_spatialstability_prctile(:,unit) > 95)
           low_sig = low_sig+1; 
        end
        
        if peak_firing_rate(3,unit) < 1 && any(peak_firing_rate(:,unit) > 1)
            barely_active_count(1) =   barely_active_count(1)+1;
            if any((spatial_info.shuffled_rate_prctile(:,unit) > 95) & (spatial_info.shuffled_spatialstability_prctile(:,unit) > 95))
                 barely_active_count(2) =  barely_active_count(2)+1;
            end
        end
        
        
        if all(isnan(peak_firing_rate(:,unit))) || all(peak_firing_rate(:,unit) < 1)
            continue %unit doesn't fire enough go to next unit
        elseif all((spatial_info.shuffled_rate_prctile(:,unit) < 95) & (spatial_info.shuffled_spatialstability_prctile(:,unit) < 95))
            continue %unit likely not spatial
        elseif  ~any((spatial_info.shuffled_rate_prctile(:,unit) > 95) & (spatial_info.shuffled_spatialstability_prctile(:,unit) > 95))
            continue %only process cells that pass both criterion under at least 1 condition
        end
        if go_no_go == 0
            disp('Why')
        end
        count2 = count2+1;
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Format Generic Data into a "Table"---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---Calculate Percent of Total for Variables of Interest---%
reason_unused(4,:) = sum(reason_unused(1:3,:));  %calculate total unused number of cells
percent_total_usable = round(100*total_used./total_recorded);
percent_numblks = round(100*num_blocks./((ones(3,1)*total_used)));
percent_unused = round(100*reason_unused(4,:)./total_recorded);
percent_reason_unused = round(100*reason_unused(1:3,:)./(ones(3,1)*reason_unused(4,:)));

%---Make "Table"--%
fprintf('\n\n\n')
fprintf('%%----------Cell Counts----------%% \n')
fprintf('\t\t\t\t\t SingleUnit \t\t MultiUnit \n')
fprintf(['Total Recorded Units \t ' num2str(total_recorded(1)) '\t\t\t\t '  num2str(total_recorded(2)) '\n'])
fprintf('\n')

%usable units
fprintf(['Total Usable(%% total) \t  '  num2str(total_used(1)) '(' num2str(percent_total_usable(1)) '%%) \t\t\t '...
    num2str(total_used(2)) '(' num2str(percent_total_usable(2)) '%%)\n'])
fprintf(['\t 1+ Blks(%% usable) \t  '  num2str(num_blocks(1,1)) '(' num2str(percent_numblks(1,1)) '%%) \t\t ' ...
    num2str(num_blocks(1,2)) '(' num2str(percent_numblks(1,2)) '%%) \n'])
fprintf(['\t 2+ Blks(%% usable) \t  '  num2str(num_blocks(2,1)) '(' num2str(percent_numblks(2,1)) '%%) \t\t\t ' ...
    num2str(num_blocks(2,2)) '(' num2str(percent_numblks(2,2)) '%%) \n'])
fprintf(['\t 3+ Blks(%% usable) \t  '  num2str(num_blocks(3,1)) '(' num2str(percent_numblks(3,1)) '%%) \t\t\t ' ...
    num2str(num_blocks(3,2)) '(' num2str(percent_numblks(3,2)) '%%) \n'])

%unusable units
fprintf('\n')
fprintf(['Total Unsable(%% total) \t  '  num2str(reason_unused(4,1)) '(' num2str(percent_unused(1)) '%%) \t\t\t '...
    num2str(reason_unused(4,2)) '(' num2str(percent_unused(2)) '%%)\n'])

fprintf(['\t Spk Count(%% unusable) \t'  num2str(reason_unused(1,1)) '(' num2str(percent_reason_unused(1,1)) '%%) \t\t ' ...
    num2str(reason_unused(1,2)) '(' num2str(percent_reason_unused(1,2)) '%%) \n'])
fprintf(['\t Unstable(%% unusable) \t'  num2str(reason_unused(2,1)) '(' num2str(percent_reason_unused(2,1)) '%%) \t\t ' ...
    num2str(reason_unused(2,2)) '(' num2str(percent_reason_unused(2,2)) '%%) \n'])
fprintf(['\t Quality(%% unusable) \t '  num2str(reason_unused(3,1)) '(' num2str(percent_reason_unused(3,1)) '%%) \t\t ' ...
    num2str(reason_unused(3,2)) '(' num2str(percent_reason_unused(3,2)) '%%) \n'])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Format Spatial Data into a "Table"---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if all(spatial_total_used ~= total_used)
    warning('Mismatch in total usable units and units used for spatial analysis')
end

%---Calculate Percent of Total for Variables of Interest---%
percent_low_firing = round(100*spatial_total_low_firing./spatial_total_used);
percent_skaggs = round(100*spatial_skaggs./spatial_total_used);
percent_skaggs_active = round(100*spatial_skaggs./(spatial_total_used-spatial_total_low_firing));
percent_stability = round(100*spatial_stability./spatial_total_used);
percent_stability_active = round(100*spatial_stability./(spatial_total_used-spatial_total_low_firing));
percent_skaggs_stability = round(100*spatial_skaggs_stability./spatial_total_used);
percent_skaggs_stability_active = round(100*spatial_skaggs_stability./(spatial_total_used-spatial_total_low_firing));
%%
fprintf('\n\n\n')
fprintf('%%------------------Spatial Analysis------------------%% \n')
fprintf('\t\t\t\t\t\t SingleUnit \t\t MultiUnit \n')
fprintf(['Total Units \t\t\t\t ' num2str(spatial_total_used(1)) '\t\t\t\t '  num2str(spatial_total_used(2)) '\n'])
fprintf('\n')
fprintf(['Inactive(%% total) \t\t\t ' num2str(spatial_total_low_firing(1)) '(' num2str(percent_low_firing(1)) '%%) \t\t\t '...
    num2str(spatial_total_low_firing(2)) '(' num2str(percent_low_firing(2)) '%%) \n'])
fprintf('\n')

fprintf(['Skaggs(%%total/%%active) \t\t ' num2str(spatial_skaggs(1))...
    '(' num2str(percent_skaggs(1)) '%%/' num2str(percent_skaggs_active(1)) '%%) \t\t'...
    num2str(spatial_skaggs(2)) '(' num2str(percent_skaggs(2)) '%%/' num2str(percent_skaggs_active(2)) '%%)\n'])
fprintf(['Stability(%%total/%%active) \t ' num2str(spatial_stability(1))...
    '(' num2str(percent_stability(1)) '%%/' num2str(percent_stability_active(1)) '%%) \t\t'...
    num2str(spatial_stability(2)) '(' num2str(percent_stability(2)) '%%/' num2str(percent_stability_active(2)) '%%)\n'])
fprintf(['Skaggs & Stability \t\t\t ' num2str(spatial_skaggs_stability(1))...
    '(' num2str(percent_skaggs_stability(1)) '%%/' num2str(percent_skaggs_stability_active(1)) '%%) \t\t'...
    num2str(spatial_skaggs_stability(2)) '(' num2str(percent_skaggs_stability(2)) '%%/' num2str(percent_skaggs_stability_active(2)) '%%)\n'])
fprintf(['(%%total/%%active) \n']')
