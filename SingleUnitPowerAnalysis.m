% attempt at a "power analysis" for number of images to show reliable
% difference between conditions. 
% Written Seth Konig January 18, 2015

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\'
%%
% %---novel vs repeat analysis---%
% files = {'PW141013','PW140908'};
% cellnum = [3,1];%refers to cell index not actual # such as sig001c
% time_range = {[500 1500],[200 500]}; %time where it looks like there's a
% % difference in firing rate between novel and repeat presentations
% 
% N = NaN(2,length(files));%row 1 power 0.8, row 2 power 0.9
% novrep = [ones(1,16),zeros(1,16),ones(1,16),zeros(1,16),ones(1,16),zeros(1,16),...
%     ones(1,16),zeros(1,16),ones(1,16),zeros(1,16),ones(1,16),zeros(1,16)]; 
% for f = 1:length(files);
%     load([data_dir files{f} '-ListSQ-time_locked_results.mat']);
%     data = time_lock_firing{16,cellnum(f)};
%     
%     num_spikes = sum(data(:,500+time_range{f}(1):500+time_range{f}(2))'); 
%     nov = novrep(1:length(num_spikes)); 
%     novspikes = num_spikes(nov==1);
%     repspikes = num_spikes(nov == 0);
%     N(1,f) = sampsizepwr('t',[mean(novspikes),std(novspikes)],mean(repspikes),0.8);
%     N(2,f) = sampsizepwr('t',[mean(novspikes),std(novspikes)],mean(repspikes),0.9);
% end
%%
% %---Saccade Locked Sequence analysis---% 
% % sequence 1 vs sequence 2 results locked to saccaes
% % most are perisaccadic activity and these are good representatives
% files = {'PW140908','PW140908','PW141010'}; %cell 1 2x on purporse
% cellnum = [3,3,11];%refers to cell index not actual # such as sig001c
% time_range = {[-200 200],[-200 200],[-150 50]}; %time where it looks like there's a
% item_num = [1 4 2];
% N = NaN(2,length(files));%row 1 power 0.8, row 2 power 0.9
% % difference in firing rate between novel and repeat presentations
% 
% for f = 3%:length(files);
%     load([data_dir files{f} '-ListSQ-time_locked_results.mat'],'which_sequence');
%     load([data_dir files{f} '-Eyemovement_Locked_Sequence_results.mat']);
%     
%     data = saccade_locked_firing{item_num(f),cellnum(f)};
%     
%     num_spikes = sum(data(:,500+time_range{f}(1):500+time_range{f}(2))'); 
%     % not sure but some of the trials have NaNs probably becuase couldn't
%     % properly detect a saccade or something 
%     
%     seq1 = num_spikes(which_sequence == 1);
%     seq2 = num_spikes(which_sequence == 2);
%     N(1,f) = sampsizepwr('t',[nanmean(seq1),nanstd(seq1)],nanmean(seq2),0.8);
%     N(2,f) = sampsizepwr('t',[nanmean(seq1),nanstd(seq1)],nanmean(seq2),0.9);
% end
%%
%---Sequence Task analysis---% 
% sequence 1 vs sequence 2 results locked to particular events
files = {'PW140801-Sequence','PW140825_3','PW141023-ListSQ'}; %cell 1 2x on purporse
cellnum = [1 1 2];%refers to cell index not actual # such as sig001c
time_range = {[50 350],[0 500],[-100 200]}; %time where it looks like there's a
epic_num = [3 13 10];
N = NaN(2,length(files));%row 1 power 0.8, row 2 power 0.9
% difference in firing rate between novel and repeat presentations

for f = 1%:length(files);
    load([data_dir files{f} '-time_locked_results.mat']);
    
    data = time_lock_firing{epic_num(f),cellnum(f)};
    
    num_spikes = sum(data(:,500+time_range{f}(1):500+time_range{f}(2))'); 
    % not sure but some of the trials have NaNs probably becuase couldn't
    % properly detect a saccade or something 
    
    seq1 = num_spikes(which_sequence == 1);
    seq2 = num_spikes(which_sequence == 2);
    N(1,f) = sampsizepwr('t',[nanmean(seq1),nanstd(seq1)],nanmean(seq2),0.8);
    N(2,f) = sampsizepwr('t',[nanmean(seq1),nanstd(seq1)],nanmean(seq2),0.9);
end