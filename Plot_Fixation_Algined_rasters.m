%plots fixation aligned rasters so that can export to eps

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';

task_file = 'TO160311_3';
unit_name = 'sig001a';

load([data_dir  task_file(1:8) '-Place_Cell_Analysis.mat']);

this_unit = [];
for unit = 1:size(unit_stats,2)
    if strcmpi(unit_stats{1,unit},unit_name)
        this_unit = unit;
        break
    end
end

fix_locked_firing = list_fixation_locked_firing{this_unit};
fix_in_out = in_out{this_unit};
t = -twin1:twin2-1;

num_in = sum(fix_in_out == 1);
num_out = sum(fix_in_out == 4);
downsample = floor(num_out/num_in)/2; %many more fixations outside so downsample to show ~equal number
figure
%---Fixations in->out vs out->out---%
subplot(2,2,1)
out_matrix = fix_locked_firing(fix_in_out == 4,:);
out_matrix = out_matrix(1:downsample:end,:);
[trial,time] = find(out_matrix == 1);
plot(time-twin1,(trial),'.b')
hold on
if ~isempty(trial)
    b4 = max(trial);
else
    b4 = 0;
end
[trial,time] = find(fix_locked_firing(fix_in_out == 1,:) == 1);
trial = trial+b4;
plot(time-twin1,(trial),'.r')
if ~isempty(trial)
    ylim([0 max(trial)])
else
    ylim([0 b4])
end
box off
plot([0 0],[0 max(trial)+1],'k--')
ylabel('Occurence')
set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
xlabel('Time from Fixation Start (ms)')
title('Fixation Aligned Rasters')
axis square

subplot(2,2,2)
hold on
[~,~,~,y_list,~] = dofill(t,fix_locked_firing(fix_in_out == 1,:),'red',1,smval);%out-> in
dofill(t,fix_locked_firing(fix_in_out == 4,:),'blue',1,smval);%out->out
%     plot(t,list_95_curve{1,unit},'k','linewidth',2);%95% confidence interval
[pks,locs] = findpeaks(y_list,'MinPeakWidth',40);
locs(pks < 0.66*max(y_list)) = [];
pks(pks < 0.66*max(y_list)) = [];
plot(locs-twin1,pks,'*k')
yl = ylim;
if yl(1) < 0
    yl(1) = 0;
    ylim(yl);
end
plot([0 0],[yl(1) yl(2)],'k--')
% gaps = findgaps(sig_ind);
% if ~isempty(gaps)
%     for g = 1:size(gaps,1)
%         gp = gaps(g,:);
%         gp(gp == 0) = [];
%         if length(gp) > 40
%             h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
%                 [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
%             uistack(h,'down')
%             set(h,'facealpha',.25,'EdgeColor','None')
%         end
%     end
% end
xlim([-twin1 twin2]);
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Firing Rate (Hz)')
legend('out->in','out->out','Location','NorthWest')
set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
% title(sprintf(['n_{out->in} = ' num2str(sum(fix_in_out == 1))...
%     ', n_{out->out} = ' num2str(sum(fix_in_out == 4))]));
title(['Peak Firing Rate of ' num2str(pks,3) ' Hz at ' num2str(locs-twin1) ' ms'])
axis square

