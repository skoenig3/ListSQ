unit = 1;
event = 13;
epoch_name = 'First Reward';
figure
hold on
time_matrix1 = time_lock_firing{event,unit}(which_sequence == 1,:);
time_matrix2 = time_lock_firing{event,unit}(which_sequence == 2,:);
t = 1:size(time_matrix1,2);
dofill(t,time_matrix1,'blue',1,smval);
dofill(t,time_matrix2,'red',1,smval);
set(gca,'Xtick',0:250:size(time_matrix1,2))
set(gca,'XtickLabel',num2cell((0:250:size(time_matrix1,2))-500));
ylabel('Firing Rate (Hz)')
xlabel(['Time from ' epoch_name ' (ms)'])
ylimit = ylim;
ylim([0 ylimit(2)])
hold off
export_fig C:\Users\seth.koenig\Desktop\File.pdf