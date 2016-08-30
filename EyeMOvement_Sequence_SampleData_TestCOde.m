% %look at relationship between 6 simulatenously recorded eye movement
% %modulated neurons
% 
% Fs = 1000;
% smval = 60;
% 
% load('C:\Users\seth.koenig\Documents\MATLAB\ListSQ\EyeMovement_Sequence_SampleData.mat');
% 
% num_units = length(unit_names);
% saccade_firing_rate = cell(1,length(num_units));
% fixation_firing_rate =  cell(1,length(num_units));
% 
% for unit = 1:num_units
%     
%     saccade_firing = saccade_locked_firing{unit}(saccade_info(4,:) > 2*twin,:);
%     saccade_firing(sum(isnan(saccade_firing),2) == 2*twin,:) = [];
%     saccade_firing(nansum(saccade_firing,2) == 0,:) = [];%remove saccades without spikes
%     [saccade_firing,~]= nandens(saccade_firing,smval,'gauss',Fs,'nanflt');%nandens made by nathan killian
% %     saccade_firing = saccade_firing-mean(saccade_firing);
% %     saccade_firing = saccade_firing/max(abs(saccade_firing));
%     saccade_firing_rate{unit} = saccade_firing;
%     
%     fixation_firing = fixation_locked_firing{unit}(fixation_info(4,:) > 2*twin,:);
%     fixation_firing(nansum(fixation_firing,2) == 0,:) = [];%remove fixations without spikes
%     fixation_firing(sum(isnan(fixation_firing),2) == 2*twin,:) = [];
%     [fixation_firing,~]= nandens(fixation_firing,smval,'gauss',Fs,'nanflt');%nandens made by nathan killian
% %     fixation_firing = fixation_firing-mean(fixation_firing);
% %     fixation_firing = fixation_firing/max(abs(fixation_firing));
%     fixation_firing_rate{unit} = fixation_firing;
%     
% end
%%
figure
subplot(1,2,1)
hold on
for unit = 1:2%num_units
   plot(saccade_firing_rate{unit}) 
end
xlabel('Time from Saccade Start (ms)')
ylabel('Firing Rate (Hz)')


subplot(1,2,2)
hold on
for unit = 1:num_units
   plot(fixation_firing_rate{unit}) 
end
xlabel('Time from Fixation Start (ms)')
ylabel('Firing Rate (Hz)')
%%
sb = reshape(1:(num_units-1)^2,num_units-1,num_units-1)';
t = -(2*twin-1):(2*twin-1);
figure
for i = 1:num_units
    for j = 1:num_units
        if j > i
            all_xc = NaN(size(saccade_locked_firing,1),4*twin-1);
            count = 0;
            for eye = 1:size(saccade_locked_firing{i},1);
                if ~all(isnan(saccade_locked_firing{i}(eye,:))) &&  ~all(isnan(saccade_locked_firing{j}(eye,:)))
                    xc = xcorr(saccade_locked_firing{i}(eye,:),saccade_locked_firing{j}(eye,:));
                    all_xc(eye,:) = xc;
                    count = count+1;
                end
            end
            all_xc = laundry(all_xc);
            all_xc = all_xc/count;
            sm = nandens(all_xc,10,'gauss',Fs,'nanflt');
            subplot(num_units-1,num_units-1,sb(i,j-1))
            plot(t,sm,'k')
            yl = ylim;
            if yl(1) < 0 
                ylim([0 yl(2)]);
            end
            xlim([-250 250])
            xlabel('Time (ms)')
            ylabel('Probability')
        end
    end
end