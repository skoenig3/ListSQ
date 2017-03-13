% written by Seth Konig August 2014. Updated to V2 by SDK 1/7/16. Updated
% to import files, mutliunit stats, etc. from excel file. V2 also analyzes
% firing rates over time and uses only stable portion of the task.
% code runs all the other code preprocess then process recording data
clar

task = 'ListSQ';
set(0,'DefaultFigureVisible','OFF');
for monkey = 2:-1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
%         listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 155;%155.85 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 135;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
%         listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import and Pre-Process Recording Data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        for session = 1:length(session_data)
    %             ImportListSQRecordingDataV2(data_dir,figure_dir,session_data{session})
    %         end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Plot Waveforms and Rasters to Determine Firing Stability---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %addpath('C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PlexonSDK_v1_6\')
    %         for session = 1:length(session_data)
    %             disp(['#' num2str(session) ' ' session_data{session}.task1_file])
    %             make_rasters_and_plot_waveformsV2(data_dir,figure_dir,session_data{session},task)
    %         end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Spatially Modulated---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session = 1:length(session_data)
    %         disp(['#' num2str(session) ' ' session_data{session}.task1_file])
    %         spatial_analysisV2(data_dir,figure_dir,session_data{session},task)
    %     end
    %     emailme('Done spatial')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine Place Cell Reliablity Across Fixations in Field vs out Field---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session =1:length(session_data)
    %         disp(['#' num2str(session) ' ' session_data{session}.task1_file])
    %         Place_Cell_Fixation_Analysis(data_dir,figure_dir,session_data{session},predict_rt)
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Modulated by Saccade Direction---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session = 1:length(session_data)
    %         disp(['#' num2str(session) ' ' session_data{session}.task1_file])
    %         ListSQ_Saccade_Direction_Analysis(data_dir,figure_dir,session_data{session})
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Modulated by Saccade Amplitude---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session = 23%1:length(session_data)
    %         disp(['#' num2str(session) ' ' session_data{session}.task1_file])
    %         ListSQ_Saccade_Amplitude_Analysis(data_dir,figure_dir,session_data{session})
    %     end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Determine if Neurons are Fixation/Saccade Modulated in List Task---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         for session = 1:length(session_data)
    %             List_Fixation_Analysis(data_dir,figure_dir,session_data{session})
    %             %List_Saccade_AnalysisV2(data_dir,figure_dir,session_data{sesseion}) %not rechecked for bugs
    %         end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Visually reponsive/like Cross hair---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session = 1:length(session_data)
    %         Visual_Response_AnalysisV2(data_dir,figure_dir,session_data{session});
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Visually Responsive Neurons Show Novel/Repeat Differences---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session = 1 :length(session_data)
    %         Visual_Response_Memory(data_dir,figure_dir,session_data{session});
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Temporally Modulated---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session = 1:length(session_data)
    %         time_locked_analysisV2(data_dir,figure_dir,session_data{session},task)
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Fixaiton Modulated in Sequence Task---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for session = 1:length(session_data)
%         Sequence_Fixation_AnalysisV2(data_dir,figure_dir,session_data{session})
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Contextually Modulated in Sequence Task---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session = 1:length(session_data)
    %         Sequence_Contextual_Fixation_Analysis(data_dir,figure_dir,session_data{session})
    %     end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Modulated by Reaction Time in Sequence Task---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session = 1:length(session_data)
    %         Sequence_ReactionTime_Fixation_Analysis(data_dir,figure_dir,session_data{session},predict_rt)
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Firing During ITI Period---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for session =1:length(session_data)
            ITI_analysis(data_dir,figure_dir,session_data{session})
        end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Firing During ITI Period---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session = 1:length(session_data)
    %         Reward_analysis(data_dir,figure_dir,session_data{session})
    %     end
end
emailme('Batch Done Running')
set(0,'DefaultFigureVisible','ON');