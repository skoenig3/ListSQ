% Calculate the distribution of time to fixation
% written by Seth Koenig December 5, 2014

% 
% dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\';
% alltime2fixation = NaN(1,4);
% 
% cd(dir);
% a = what;
% m = a.mat;
% for file = 1:size(m,1);
%     if ~isempty(strfind(m{file},'Eyemovement_Locked_Sequence_results.mat'))
%         load(m{file},'time_to_fixation');
%         if size(time_to_fixation,1) > 100%so got through at least some of the task
%             alltime2fixation = [alltime2fixation; time_to_fixation(21:end,:)];%ignore familiarization block
%         end
%     end
% end
% 
% a = alltime2fixation(:,2:end);
% hist(a(1:end),200)
% xlim([-450 400])
% xlabel('Time to fixation (ms)')
% ylabel('Count')

%%

dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\';
alltime2fixation = NaN(1,4);

cd(dir);
a = what;
m = a.mat;
predict_count = 0;
item_count = 0;
for file = 1:size(m,1);
    if ~isempty(strfind(m{file},'-preprocessed.mat'))
        item_set = [];
        load(m{file},'cfg','item_set');
        if ~isempty(strfind(item_set,'ckwenzr')) || isempty(item_set)%%ckwenzr or cvt task
            continue
        end
        for t = 1:length(cfg.trl)
            if  any(cfg.trl(t).allval == 3)%rewarded sequence trial
                events = cfg.trl(t).allval;
                item_count = item_count  + 4;
                predict_count = predict_count + sum(events == 20)-1;
            end
        end
    end
end

