% written by Seth Konig August 2014. Updated to V2 by SDK 1/7/16. Updated
% to import files, mutliunit stats, etc. from excel file. V2 also analyzes
% firing rates over time and uses only stable portion of the task. 
% code runs all the other code preprocess then process recording data
cla
excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';

predict_rt = 156;%156 ms for Vivian, 138 ms for Tobii


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Read in Excel Sheet for Session data---%%%
%only need to run when somethings changed or sessions have been added
excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx'];
listsq_read_excel(data_dir,excel_file);

load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Preprocess all the ListSQ data---%%%
% for sess = 1:length(session_data) %last session imported #27
%     ImportListSQRecordingDataV2(data_dir,figure_dir,session_data{sess})
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%---Automatically plot rasters---%%%
% for sess = 1:length(session_data)
%     task = 'ListSQ';
%     make_rasters_and_plot_waveformsV2(data_dir,figure_dir,session_data{sess},task)
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%---Autmoatically analyze for event locked cells for ListSQ---%%%
% for sess = length(session_data)
%     task = 'ListSQ';
%     time_locked_analysisV2(data_dir,figure_dir,session_data{sess},task)
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%---Autmoatically analyze ITI cells for ListSQ---%%%
% for sess = 10:length(session_data)
%     task = 'ListSQ';
%     ITI_analysis(data_dir,figure_dir,session_data{sess},task)
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%---Autmoatically analyze Sequence Eye Movements---%%%
% for sess = length(session_data)
%     task = 'ListSQ';
%     Sequence_Saccade_AnalysisV2(data_dir,figure_dir,session_data{sess},task,predict_rt)
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%---Autmoatically analyze List Eye Movements---%%%
% for sess = length(session_data)
%     List_Saccade_AnalysisV2(data_dir,figure_dir,session_data{sess})
% end
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Autmoatically analyze ListSQ for Spatially Modulated Cells---%%%
for sess = 1:length(session_data)
    spatial_analysisV2(data_dir,figure_dir,session_data{sess},'ListSQ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Autmoatically analyze ListSQ for Saccade Direction Cells---%%%
% for sess = 1:length(session_data)
%    List_Saccade_Direction_AnalysisV2(data_dir,figure_dir,session_data{sess})
% end
