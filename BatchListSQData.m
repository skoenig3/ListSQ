% written by Seth Konig August 2014
% code runs all the other code preprocess then process recording data

clc
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Figures\'; 

multiunits = {[0 0 1],... %for plotting and other purpose 1 mulitunit 0 single unit
            [0 0 0 1 1 1 0 0 0],[zeros(1,9) 1],[zeros(1,5) 1],[0 1 zeros(1,5) 1 zeros(1,5) 1],...
            [0 0 0 1 0 0],[0 0 0],[0 0 0 0],[ 0 1 0 0 0 0],[zeros(1,7)],[0 0 0],[0 0 0 0],...
              [zeros(1,4) 1],0,[1 0 0], [0 1 0],zeros(1,7),zeros(1,6),[zeros(1,8) 1 0],...
              zeros(1,11),zeros(1,10),zeros(1,3),[zeros(1,4) ones(1,2) zeros(1,7)],...
              zeros(1,5),zeros(1,7),zeros(1,8),[1 0 0 1],zeros(1,9),[0 0 0 1 0 0 0],...
            zeros(1,8),zeros(1,4),[zeros(1,5) 1 0],zeros(1,3),[1 0],[0 0],...
            zeros(1,9),zeros(1,8),zeros(1,3),[0 1 1 1],zeros(1,20)};
          
cch25_files = {'PW140725_2.nex','PW140728_2_2-sorted.nex','PW140729_2-sorted.nex',...
  'PW140730_2-sorted.nex','PW140801_2.nex','PW140805_2.nex','PW140806_2.nex',...
  'PW140825_1.nex','PW140826_1.nex','PW140829_1.nex','PW140908_1.nex','PW140910_1.nex',...
  'PW140915_1.nex','PW140917_1.nex','PW140101_1.nex','PW141006_1.nex',...
  'PW141007_1.nex','PW141013_1.nex','PW141014_1.nex','PW141009_1.nex','PW141015_1.nex',...
  'PW141008_1.nex','PW141010_1.nex','PW141024_1.nex','PW141027_1.nex','PW141028_1.nex',...
  'PW141029_1.nex','PW141023_1.nex','PW141031_1.nex','PW141103_1.nex','PW141105_1.nex',...
  'PW141106_1.nex','PW141110_1.nex','PW150121_1.nex','PW150122_1.nex','PW150123_1.nex',...
  'PW150126_1.nex','PW150127_1.nex','PW150204_1.nex','PW150205_1.nex'};
listsq_files = {'PW140725_3-sorted.nex','PW140728_2_3-sorted.nex','PW140729_3-sorted.nex',...
  'PW140730_3-sorted.nex','PW140801_3-sorted.nex','PW140805_3-sorted.nex','PW140806_3-sorted.nex',...
  'PW140825_3-sorted.nex','PW140826_3-sorted.nex','PW140829_3-sorted.nex','PW140908_3-sorted.nex',...
  'PW140910_3-sorted.nex','PW140915_3-sorted.nex','PW140917_3-sorted.nex','PW140101_3-sorted.nex',...
  'PW141006_3-sorted.nex','PW141007_3-sorted.nex','PW141013_3-sorted.nex','PW141014_3-sorted.nex',...
  'PW141009_3-sorted.nex','PW141015_3-sorted.nex','PW141008_3-sorted.nex','PW141010_3-sorted.nex',...
  'PW141024_3-sorted.nex','PW141027_3-sorted.nex','PW141028_3-sorted.nex','PW141029_3-sorted.nex',...
  'PW141023_3-sorted.nex','PW141031_3-sorted.nex','PW141103_3-sorted.nex','PW141105_3-sorted.nex',...
  'PW141106_3-sorted.nex','PW141110_3-sorted.nex','PW150121_3-sorted.nex','PW150122_3-sorted.nex',...
  'PW150123_3-sorted.nex','PW150126_3-sorted.nex','PW150127_3-sorted.nex','PW150204_3-sorted.nex',...
  'PW150205_3-sorted.nex'};
item_sets =  {'ListSQ02.itm','ListSQ03.itm','ListSQ04.itm','ListSQ05.itm','ListSQ06.itm',...
    'ListSQ07.itm','ListSQ08.itm','ListSQ09.itm','ListSQ10.itm','ListSQ12.itm','ListSQ13.itm',...
    'ListSQ14.itm','ListSQ15.itm','ListSQ16.itm','ListSQ17.itm','ListSQ19.itm',...
    'ListSQ20.itm','ListSQ23.itm','ListSQ24.itm','ListSQ21.itm','ListSQ25.itm',...
    'ListSQ21.itm','ListSQ22.itm','ListSQ29.itm','ListSQ30.itm','ListSQ31.itm',...
    'ListSQ32.itm','ListSQ28.itm','listSQ33.itm','listsq34.itm','Listsq35.itm',...
    'ListSQ37.itm','ListSQ39.itm','ListSQVA.itm','ListSQVB.itm','ListSQVC.itm',...
    'ListSQVD.itm','ListSQVE.itm','ListSQVH.itm','ListSQVI.itm'};
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Preprocess all the ListSQ data---%%%
for ls = 2%length(listsq_files)
    ImportListSqRecordingData(data_dir,cch25_files{ls},listsq_files{ls},item_sets{ls},multiunits{ls})
    save_and_close_fig(figure_dir,[listsq_files{ls}(1:10) '-calibration data from all trials'])
    save_and_close_fig(figure_dir,[listsq_files{ls}(1:10) '-accuracy of  calibration'])
end
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%---Automatically process ListSQ LFP data---%%%
% for ls = length(listsq_files)
%     ListSQ_wfANOVA_LFP_analysis(data_dir,[listsq_files{ls}(1:10) '-preprocessed'],figure_dir)
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%---Automatically plot rasters---%%%
% for ls = 12:13%:length(listsq_files)
%     task = 'List';
%     load([data_dir,listsq_files{ls}(1:end-11)  '-preprocessed'],'cfg','hdr','data','item_set','num_units');
%     make_rasters_and_plot_waveforms(listsq_files{ls},cfg,hdr,data,task,figure_dir,item_set,multiunits{ls})
% end
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%---Autmoatically analyze for event locked cells for ListSQ---%%%
% for ls = 11:13%31:length(listsq_files)
%     time_locked_analysis(data_dir,[listsq_files{ls}(1:end-11) '-preprocessed'],figure_dir,'ListSQ')
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%---Autmoatically replot all ListSQ time analysis--%%%
% for ls =  5%[1 3:13 15:length(listsq_files)]
%     preprocessed_data_file = [listsq_files{ls}(1:end-11) '-preprocessed.mat'];
%     load([data_dir,preprocessed_data_file]);
%     load([data_dir listsq_files{ls}(1:end-13) '-ListSQ-time_locked_results.mat']);
%     
%     [itmlist,sequence_items] = read_ListSQ_itm_and_cnd_files(item_set);
%     [~,novel_vs_repeat] = get_image_numbers(cfg,itmlist,sequence_items,23);
%      
%     unit_names.name = cfg.channel(1:num_units);
%     unit_names.multiunit = multiunit;
%     smval = 120; 
%     task_data.task_type = 'ListSQ_Sequence';
%     task_data.which_sequence = which_sequence;
%     task_data.novel_vs_repeat = novel_vs_repeat;
%     time_locked_plots(figure_dir,preprocessed_data_file,time_lock_firing,epoch_data,temporal_info,task_data,unit_names,smval)
%     task_data.task_type = 'ListSQ_List';
%     time_locked_plots(figure_dir,preprocessed_data_file,time_lock_firing,epoch_data,temporal_info,task_data,unit_names,[smval,smval_novrep]) 
% end
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %---Autmoatically analyze for Spatial cells for ListSQ---%%%
% for ls = 1:length(listsq_files)
%     spatial_analysis(data_dir,[listsq_files{ls}(1:end-11) '-preprocessed'],figure_dir,'ListSQ')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Autmoatically analyze Sequence Eye Movements---%%%
% for ls = 1:length(listsq_files)
%     Sequence_Saccade_Analysis(data_dir,[listsq_files{ls}(1:end-11) '-preprocessed'],figure_dir,'ListSQ')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Autmoatically analyze List Eye Movements---%%%
% for ls = 15:length(listsq_files)
%     List_Saccade_Analysis(data_dir,[listsq_files{ls}(1:end-11) '-preprocessed'],figure_dir)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Autmoatically analyze for Spatial cells for ListSQ using ANOVA analysis---%%%
% for ls = [1 3:13 15:length(listsq_files)]
%     spatial_ANOVA_analysis(data_dir,[listsq_files{ls}(1:end-11) '-preprocessed'],figure_dir)
% end

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