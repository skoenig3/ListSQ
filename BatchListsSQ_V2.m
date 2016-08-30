% written by Seth Konig August 2014. Updated to V2 by SDK 1/7/16. Updated
% to import files, mutliunit stats, etc. from excel file. V2 also analyzes
% firing rates over time and uses only stable portion of the task.
% code runs all the other code preprocess then process recording data
clar

% monkey = 'Vivian';
% monkey = 'Tobii';

% parpool
task = 'ListSQ';
for monkey = 1:2
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import and Pre-Process Recording Data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    for session = 1:length(session_data)
    %         ImportListSQRecordingDataV2(data_dir,figure_dir,session_data{session})
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Plot Waveforms and Rasters to Determine Firing Stability---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %addpath('C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PlexonSDK_v1_6\')
    %     for session = 1:length(session_data)
    %         disp(['#' num2str(session) ' ' session_data{session}.task1_file])
    %         make_rasters_and_plot_waveformsV2(data_dir,figure_dir,session_data{session},task)
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Spatially Modulated---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for session = 1:length(session_data)
        disp(['#' num2str(session) ' ' session_data{session}.task1_file])
        spatial_analysisV2(data_dir,figure_dir,session_data{session},task)
    end
    emailme('Done spatial')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Temporally Modulated---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session = 1:length(session_data)
    %         time_locked_analysisV2(data_dir,figure_dir,session_data{session},task)
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Saccade Modulated---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session = 2%1:length(session_data)
    %         List_Saccade_AnalysisV2(data_dir,figure_dir,session_data{session})
    %     end
end
emailme('Batch Done Running')
