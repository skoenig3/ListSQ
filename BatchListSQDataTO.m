% written by Seth Konig August 2014
% code runs all the other code preprocess then process recording data

clc
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Figures2\';

predicted_rt = 135;%maximum "reaction time" for what constitutes as predictive, everything else is reactive

multiunits = {zeros(1,16)};
%           
cch25_files = {'TO160401_2.nex'};
listsq_files = {'TO160401_3-sorted.nex'};
item_sets =  {'ListSQ52.itm'};
 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %---Preprocess all the ListSQ data---%%%
% for ls =1:length(listsq_files)
%     ImportListSqRecordingData(data_dir,cch25_files{ls},listsq_files{ls},item_sets{ls},multiunits{ls},figure_dir)
%     %ImportListSqRecordingDataNoStrobe(data_dir,cch25_files{ls},listsq_files{ls},item_sets{ls},multiunits{ls},figure_dir)
% end
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%---Automatically plot rasters---%%%
% for ls =1:length(listsq_files)
%     task = 'ListSQ';
%     load([data_dir,listsq_files{ls}(1:end-11)  '-preprocessed'],'cfg','data','item_set');
%     make_rasters_and_plot_waveforms(listsq_files{ls},cfg,data,task,figure_dir,item_set,multiunits{ls})
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %---Autmoatically analyze for event locked cells for ListSQ---%%%
% for ls =1:length(listsq_files)
%     time_locked_analysis(data_dir,[listsq_files{ls}(1:end-11) '-preprocessed'],figure_dir,'ListSQ')
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %---Autmoatically analyze for Spatial cells for ListSQ---%%%
% for ls =1:length(listsq_files)
%     spatial_analysis(data_dir,[listsq_files{ls}(1:end-11) '-preprocessed'],figure_dir,'ListSQ')
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %---Autmoatically analyze Sequence Eye Movements---%%%
% for ls =1:length(listsq_files)
%     Sequence_Saccade_Analysis(data_dir,[listsq_files{ls}(1:end-11) '-preprocessed'],figure_dir,'ListSQ',predicted_rt)
% end
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %---Autmoatically analyze List Eye Movements---%%%
% for ls =1:length(listsq_files)
%     List_Saccade_Analysis(data_dir,[listsq_files{ls}(1:end-11) '-preprocessed'],figure_dir)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Autmoatically analyze List Direction Specific---%%%
for ls =1:length(listsq_files)
    List_Saccade_Direction_Analysis(data_dir,[listsq_files{ls}(1:end-11) '-preprocessed'],figure_dir)
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Autmoatically analyze Eye movement Locked Sequence Data with T-test---%%%
% for ls = [1 3:13 15:length(listsq_files)]
%     Sequence_Saccade_TtestAnalysis(data_dir,[listsq_files{ls}(1:end-11) '-preprocessed'],figure_dir,'ListSQ')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Autmoatically analyze Visual Response to Sequence Items using T-test---%%%
% for ls = [1 3:13 15:length(listsq_files)]
%     Sequence_VisualResponse_TtestAnalysis(data_dir,[listsq_files{ls}(1:end-11) '-preprocessed'],figure_dir,'ListSQ')
% end