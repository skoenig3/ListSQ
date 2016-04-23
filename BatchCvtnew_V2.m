% written by Seth Konig August 2014. Updated to V2 by SDK 1/19/16. Updated
% to import files, mutliunit stats, etc. from excel file. V2 also analyzes
% firing rates over time and uses only stable portion of the task. 

clc

excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Read in Excel Sheet for Session data---%%%
%only need to run when somethings changed or sessions have been added
excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx'];
listsq_read_excel(data_dir,excel_file);

load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---Preprocess all the cvtnew data---%%%
for sess = 12:length(session_data)
    ImportCVTNEWRecordingDataV2(data_dir,session_data{sess})
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Automatically plot rasters---%%%
% for sess = length(session_data)
%     task = 'cvtnew';
%     make_rasters_and_plot_waveformsV2(data_dir,figure_dir,session_data{sess},task)
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Autmoatically analyze for event locked cells for CVTNEW--%%%
% for sess = length(session_data)
%     task = 'cvtnew';
%     time_locked_analysisV2(data_dir,figure_dir,session_data{sess},task)
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Autmoatically analyze for Spatial cells for CVTNEW---%%%
% for sess = 1:length(session_data)
%     spatial_analysisV2(data_dir,figure_dir,session_data{sess},'cvtnew')
% end