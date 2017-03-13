% written by Seth Konig August 2014. Updated to V2 by SDK 1/19/16. Updated
% to import files, mutliunit stats, etc. from excel file. V2 also analyzes
% firing rates over time and uses only stable portion of the task.

clar
task = 'cvtnew';
set(0,'DefaultFigureVisible','OFF');
for monkey = 2%:-1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        listsq_read_excel(data_dir,excel_file);
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
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Preprocess all the cvtnew data---%%%
%         for sess = 11%1:length(session_data)
%             ImportCVTNEWRecordingDataV2(data_dir,session_data{sess})
%         end
    
    % emailme('Done')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Automatically plot rasters---%%%
    %     for sess = 1:length(session_data)
    %         task = 'cvtnew';
    %         make_rasters_ and_plot_waveformsV2(data_dir,figure_dir,session_data{sess},task)
    %         disp('------------------------------')
    %         disp(['Sess #: ' num2str(sess)]);
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Autmoatically analyze for Spatial cells for CVTNEW---%%%
    for sess = 46:length(session_data)
        spatial_analysisV2(data_dir,figure_dir,session_data{sess},'cvtnew')
    end
    emailme('Done processing spatial')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%---Autmoatically analyze for event locked cells for CVTNEW--%%%
    %     for sess = 1:length(session_data)
    %         cvtnew_time_locked_analysis(data_dir,figure_dir,session_data{sess})
    %     end
    %     emailme('Done processing temporal')
    
end
set(0,'DefaultFigureVisible','ON');