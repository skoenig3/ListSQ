clar
task = 'ListSQ';
window_width = 100;
step = 10;
imageX = 800;
imageY = 600;
numshuffs = 100;

all_mlrs = [];
all_mlr_pctiles = [];
all_pvals = [];
% all_mlrs_in2in = [];
% all_mlrs_in2in_prctile = [];
% all_mlrs_out2out = [];
% all_mlrs_out2out_prctile = [];
% all_pvals_in2in = [];
% all_pvals_out2out = [];
max_fr = [];
spatialness = [];
all_unit_names = {};
for monkey = 2:-1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML[
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for session =1:length(session_data)
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{session},task);
        if isempty(task_file)
            continue
        end
        disp(num2str(session))
        if exist([data_dir task_file(1:end-13) '-Saccade_Direction_Analysis.mat'],'file') %want to remove later
            load([data_dir task_file(1:end-13) '-Saccade_Direction_Analysis.mat'])
            load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'spatial_info')
        else
            continue
        end
        
        num_units = size(unit_stats,2);
        for unit = 1:num_units
            if ~isnan(mlrs.all_fixations(unit))
                all_mlrs = [all_mlrs mlrs.all_fixations(unit)];
                all_mlr_pctiles = [all_mlr_pctiles mlrs.all_fixations_shuffled_prctile(unit)];
                all_pvals = [all_pvals uniformity_pvalue(1,unit)];
                max_fr = [max_fr max([max(binned_firing_rate_curves{1,unit}),...
                    max(binned_firing_rate_curves{2,unit}),max(binned_firing_rate_curves{3,unit})])];
                
%                 all_mlrs_in2in = [all_mlrs_in2in mlrs.in2in(unit)];
%                 all_mlrs_in2in_prctile = [all_mlrs_in2in_prctile mlrs.in2in_shuffled_prctile(unit)];
%                 all_pvals_in2in = [all_pvals_in2in uniformity_pvalue(2,unit)];
%                 
%                 all_mlrs_out2out = [all_mlrs_out2out mlrs.out2out(unit)];
%                 all_mlrs_out2out_prctile = [all_mlrs_out2out_prctile mlrs.out2out_shuffled_prctile(unit)];
%                 all_pvals_out2out = [all_pvals_out2out uniformity_pvalue(3,unit)];
                
                if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                    spatialness = [spatialness 1];
                elseif (spatial_info.shuffled_rate_prctile(unit) > 95)
                    spatialness = [spatialness 2];
                elseif  (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                    spatialness = [spatialness 3];
                else
                    spatialness = [spatialness 0];
                end
                
                all_unit_names = [all_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}];
            end
        end
    end
end