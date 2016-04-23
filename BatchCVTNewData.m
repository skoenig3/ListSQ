% written by Seth Konig August 2014
% code runs all the other code. 1st preprocess the data i.e. imports data 
% into matlab. 2nd runs temporal anylsis (time/event-locked) and spatial
% analysis. 

clc
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\'; %where to get data from
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Figures\'; %where to save figures

cvtnew_files = {'PW140728_2_4-sorted.nex','PW140729_4-sorted.nex','PW140730_4-sorted.nex',...
    'PW140806_4-sorted.nex','PW140807_3-sorted.nex','PW140829_4-sorted.nex',...
    'PW140903_3-sorted.nex','PW140909_3-sorted.nex','PW140916_3-sorted.nex',...
    'PW141007_4-sorted.nex','PW141015_4-sorted.nex','PW141009_4-sorted.nex',...
    'PW141008_4-sorted.nex','PW141024_4-sorted.nex','PW141028_4-sorted.nex',...
    'PW141105_4-sorted.nex','PW150127_4-sorted.nex','PW150130_3-sorted.nex',...
    'PW150205_4-sorted.nex'};
multiunits = {[1 0 0 0 1 1 1 0 0],... %for plotting and other purpose 1 mulitunit 0 single unit
            [1 0 0 1 0 1],[0 0 0 0 0 0],[zeros(1,5)],[0],[zeros(1,6)],[0 0],zeros(1,4),...
            zeros(1,7),[1 0 0 0 0 0 1 0],zeros(1,9),zeros(1,14),zeros(1,2),...
            [1 zeros(1,4)],zeros(1,6),zeros(1,3),zeros(1,3),0,zeros(1,17)};


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%---Preprocess all the cvtnew data---%%%
% for cvt = length(cvtnew_files)
%     ImportCVTNEWRecordingData(data_dir,cvtnew_files{cvt},multiunits{cvt})
% end
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---Autmoatically analyze for event locked cells for CVTNEW--%%%
for ls = 1:length(cvtnew_files)
    time_locked_analysis(data_dir,[cvtnew_files{ls}(1:end-11) '-preprocessed'],figure_dir,'cvtnew')
end
%%
% %%%---Automatically plot rasters---%%%
% for ls = 2:length(cvtnew_files)
%     task = 'CVTNEW';
%     if ls == 1
%         load([data_dir,'PW140728_2_4-sorted' '-preprocessed'],'cfg','hdr','data','num_units');
%     else
%         load([data_dir,cvtnew_files{ls}(1:10) '-preprocessed'],'cfg','hdr','data','num_units');
%     end
%       item_set = [];
%     make_rasters_and_plot_waveforms(cvtnew_files{ls},cfg,hdr,data,task,figure_dir,item_set,multiunits{ls})
% end
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Autmoatically analyze for Spatial cells for CVTNEW---%%%
for ls = 1:length(cvtnew_files)
    spatial_analysis(data_dir,[cvtnew_files{ls}(1:end-11) '-preprocessed'],figure_dir,'cvtnew')
end