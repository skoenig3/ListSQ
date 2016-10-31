data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
% load([data_dir 'TO151208-Eyemovement_Locked_List_results.mat'])
% unit = 7;
load([data_dir 'TO160105-Eyemovement_Locked_List_results.mat'])
% for unit = 1:length(fixation_information)
unit = 1;


t = -twin1:twin2-1;

%ignore anything within 500 ms of image onset
info = fixation_information{unit}((fixation_information{unit}(:,4) > image_on_twin),:);
fixation_firing = fixation_locked_firing{unit}(fixation_information{unit}(:,4) > image_on_twin,:);

%remove trials without spikes, used for shuffling anlaysis above
fixation_firing_plus = fixation_firing;
fixation_firing_plus(nansum(fixation_firing,2) == 0,:) = [];

figure

%plot firing rate curve over time
subplot(2,2,1)
hold on
[~,~,~,y,~] = dofill(t,fixation_firing_plus,'black',1,smval); %all trials with spikes
yl = ylim;
if yl(1) < 0;
    yl(1) =0;
end
hold on
plot([0 0],[yl(1) yl(2)],'k--')
hold off
ylim(yl);
ylabel('Firing Rate (Hz)')
xlabel('Time from fixation Start')
set(gca,'Xtick',[-200 -100 0 100 200 300 400])
set(gca,'XMinorTick','on','YMinorTick','on')
grid on
ylim(yl);
title(['Bit ' num2str(temporal_info.fixation.shuffled_rate_prctile(unit),3) '% '...
    '\rho_{1/2} = ' num2str(temporal_info.fixation.temporalstability(1,unit),2) ...
    ' (' num2str(temporal_info.fixation.shuffled_temporalstability_prctile(1,unit),3) ...
    '%) \rho_{e/o} = ' num2str(temporal_info.fixation.temporalstability(2,unit),2) ...
    ' (' num2str(temporal_info.fixation.shuffled_temporalstability_prctile(2,unit),3) '%)'])

%plot raster over time by fixation occurence with session
subplot(2,2,3)
[trial,time] = find(fixation_firing_plus == 1);
plot(time-twin1,(trial),'.k')
ylim([0 max(trial)])
ylabel('Occurence #')
set(gca,'Xtick',[-100 0 100 200 300 400])
xlabel('Time from fixation Start')
title('Raster whole period')


figure
%for inside field
all_xc = NaN(size(fixation_firing_plus,1),2*size(fixation_firing_plus,2)-1);
for f = 1:size(fixation_firing_plus,1)
    xc = xcorr(fixation_firing_plus(f,:),fixation_firing_plus(f,:));
    all_xc(f,:) = xc;
end
all_xc(:,(size(all_xc,2)+1)/2) = all_xc(:,(size(all_xc,2)+1)/2+1); %remove 0 lag it will always be 1
all_xc = all_xc/f;

t = -((size(all_xc,2)+1)/2-1):((size(all_xc,2)+1)/2-1);

subplot(1,2,1)
dofill(t,all_xc,'black',1,10)
% hold on
% plot(t,1000*mean(all_xc),'r')
% hold off
xlim([-250 250])
xlabel('lag (ms)')
ylabel('Spike Probability')
title('5 ms Gaussian smoothing')

subplot(1,2,2)
[~,~,~,ysmth,y]=dofill(t,all_xc,'black',1,30);
% hold on
%plot(t,1000*mean(all_xc),'r')
% hold off
xlim([-400 400])
xlabel('lag (ms)')
ylabel('Spike Probability')
title('15 ms Gaussian smoothing')

% 
% end
