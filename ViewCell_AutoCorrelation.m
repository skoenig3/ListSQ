data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
%data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
screen_size = get(0, 'ScreenSize');

w = warning ('off','all');

a = what(data_dir);
a = a.mat;

all_lags = [];
min_all_lags = []; 
n_cells = 0;
n_seq_in_out = 0; 
p_seq_in_out = 0;
for file = 1:length(a)
    if ~isempty(strfind(a{file},'-View_Cell_Analysis.mat'))
         load([data_dir a{file}],'list_fixation_locked_firing','in_out','twin',...
             'smval','multiunit','task_file','seq_p_in_out','sequences_inside')
    else
        continue
    end
    
    for unit = 1:length(list_fixation_locked_firing)
        if multiunit(unit) 
            continue
        end
        if ~isempty(list_fixation_locked_firing{unit})
           
            in_fixloc = list_fixation_locked_firing{unit}(in_out{unit} == 1,:);
            out_fixloc = list_fixation_locked_firing{unit}(in_out{unit} == 0,:);

            all_in_xc = [];
            for f = 1:size(in_fixloc,1);
                xc = xcorr(in_fixloc(f,:),in_fixloc(f,:));
                all_in_xc = [all_in_xc; xc];
            end
            all_in_xc(:,1000) = 0;
            
            all_out_xc = [];
            for f = 1:size(out_fixloc,1);
                xc = xcorr(out_fixloc(f,:),out_fixloc(f,:));
                all_out_xc = [all_out_xc; xc];
            end
            all_out_xc(:,1000) = 0;
            
            all_xc = [];
            for f = 1:size(list_fixation_locked_firing{unit},1)
                xc = xcorr(list_fixation_locked_firing{unit}(f,:),list_fixation_locked_firing{unit}(f,:));
                all_xc = [all_xc; xc];
            end
            all_xc(:,1000) = 0;
             
            t = (1:size(all_xc,2))-1000;
            figure
            subplot(1,2,1)
            dofill(t,all_xc,'black',1,20);
            hold on
            dofill(t,all_out_xc,'red',1,20);
            dofill(t,all_in_xc,'blue',1,20);
            hold off
            yl = ylim;
            yl(1) = 0;
            ylim(yl);
            xlim([-400 400])
            ylabel('Correlation (a.u.)')
            legend('All Fix','Fix out','Fix in')
            xlabel('lag')
            
            t = -499:500;
            subplot(1,2,2)
            dofill(t,list_fixation_locked_firing{unit},'black',1,20);
            hold on
            dofill(t,out_fixloc,'red',1,20);
            dofill(t,in_fixloc,'blue',1,20);
            hold off
            xlabel('Time from Fixation Onset (ms)') 
            ylabel('Firing Rate')
            set(gcf, 'Position', [0 0 screen_size(3) screen_size(4)]);
           
            close all

        end
    end
    
end