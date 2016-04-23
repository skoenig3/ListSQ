% written by Seth Konig August 2014
% code runs all the other code. 1st preprocess the data i.e. imports data 
% into matlab. 2nd runs temporal anylsis (time/event-locked) and spatial
% analysis. 

clc
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\'; %where to get data from
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Figures2\'; %where to save figures

cvtnew_files = {'TO160331_5-sorted.nex'};
multiunits = {zeros(1,25)};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---Preprocess all the cvtnew data---%%%
for cvt = length(cvtnew_files)
    ImportCVTNEWRecordingData(data_dir,cvtnew_files{cvt},multiunits{cvt})
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Autmoatically analyze for event locked cells for CVTNEW--%%%
for ls = 1:length(cvtnew_files)
    time_locked_analysis(data_dir,[cvtnew_files{ls}(1:end-11) '-preprocessed'],figure_dir,'cvtnew')
end
%%
% % %%---Automatically plot rasters---%%%
for ls = 1:length(cvtnew_files)
    task = 'CVTNEW';
    load([data_dir,cvtnew_files{ls}(1:end-11) '-preprocessed'],'cfg','data','num_units');
    item_set = [];
    make_rasters_and_plot_waveforms(cvtnew_files{ls},cfg,data,task,figure_dir,item_set,multiunits{ls})
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Autmoatically analyze for Spatial cells for CVTNEW---%%%
for ls = 1:length(cvtnew_files)
    spatial_analysis(data_dir,[cvtnew_files{ls}(1:end-11) '-preprocessed'],figure_dir,'cvtnew')
end
%%
% for ls = 1%:length(cvtnew_files);
%     cvtnew_LFPanalysis(data_dir,[cvtnew_files{ls}(1:end-11) '-preprocessed'],figure_dir)
% end