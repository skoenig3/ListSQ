% written by Seth Konig August 2014. Updated to V2 by SDK 1/19/16. Updated
% to import files, mutliunit stats, etc. from excel file. V2 also analyzes
% firing rates over time and uses only stable portion of the task. 

monkey = 'Vivian'; 
% monkey = 'Tobii'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Read in Excel Sheet for Session data---%%%
%only need to run when somethings changed or sessions have been added
if strcmpi(monkey,'Vivian')
    excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
    excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
    data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
    figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW CVTNEW Figures\';
    
    
    listsq_read_excel(data_dir,excel_file);
    load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
    
    chamber_zero = [13.5 -11]; %AP ML

elseif strcmpi(monkey,'Tobii')
    excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
    excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
    data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
    figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO CVTNEW Figures\';
    
    chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
    
    listsq_read_excel(data_dir,excel_file);
    load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
    session_data(end-1:end) = [];%last file doesn't have strobe signal working on importing the data
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Preprocess all the cvtnew data---%%%
% for sess = 1:length(session_data)
%     ImportCVTNEWRecordingDataV2(data_dir,session_data{sess})
% end

% emailme('Done')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Automatically plot rasters---%%%
for sess = 1:length(session_data)
    task = 'cvtnew';
    make_rasters_and_plot_waveformsV2(data_dir,figure_dir,session_data{sess},task)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Autmoatically analyze for Spatial cells for CVTNEW---%%%
% for sess = 1:length(session_data)
%     spatial_analysisV2(data_dir,figure_dir,session_data{sess},'cvtnew')
% end
% emailme('Done processing spatial')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Autmoatically analyze for event locked cells for CVTNEW--%%%
% for sess = 1:length(session_data)
%     task = 'cvtnew';
%     time_locked_analysisV2(data_dir,figure_dir,session_data{sess},task)
% end
% emailme('Done processing temporal')