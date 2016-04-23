% written by Seth Konig August 2014
% code runs all the other code preprocess then process recording data

clc
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Figures\'; 
task = 'Sequence';

multiunits = {0,[zeros(1,9) 1 zeros(1,7) 1],[0 0]}; %for plotting and other purpose 1 mulitunit 0 single unit
            
cch25_files = {'PW141030_1.nex','PW140801_2.nex','PW140827_1.nex'};
  
sequence_files = {'PW141030_4-sorted.nex','PW140801_4-sorted.nex','PW140827_4-sorted.nex'};
item_sets =  {'ckwenzrx.itm','ckwenzrr.itm','ckwenzrw.itm'};

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%---Preprocess all the Sequence data---%%%
% for ls =2%:length(sequence_files)
%     ImportListSqRecordingData(data_dir,cch25_files{ls},sequence_files{ls},item_sets{ls},multiunits{ls})
%     save_and_close_fig(figure_dir,[sequence_files{ls}(1:10) '-calibration data from all trials'])
%     save_and_close_fig(figure_dir,[sequence_files{ls}(1:10) '-accuracy of  calibration'])
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%---Automatically plot rasters---%%%
% for ls = 2%:length(sequence_files)
%     load([data_dir,sequence_files{ls}(1:10) '-preprocessed'],'cfg','data','item_set','num_units');
%     make_rasters_sequence(sequence_files{ls},cfg,data,figure_dir,item_set,multiunits{ls})
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%---Autmoatically analyze for event locked cells for Sequence---%%%
% for ls = 2%:length(sequence_files)
%     time_locked_analysis(data_dir,[sequence_files{ls}(1:10) '-preprocessed'],figure_dir,task)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Autmoatically analyze Sequence Eye Movements---%%%
for ls = 3%1:length(sequence_files)
    Sequence_Saccade_Analysis(data_dir,[sequence_files{ls}(1:end-11) '-preprocessed'],figure_dir,task)
end