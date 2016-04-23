% Written to determine if fixation durations are increasing in ListSQ task
% Written October 21, 2014 by Seth  Konig

dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\';
fixation_durations = {};

cd(dir);
a = what;
m = a.mat;
set = 0;
for file = 1:size(m,1);
    if ~isempty(strfind(m{file},'preprocessed'))
        fixationstats = [];
        load(m{file},'fixationstats','cfg','item_set');
        if ~isempty(strfind(item_set,'ckwenzr'))
           continue 
        end
        if ~isempty(fixationstats) %so is ListsQ
            set = set +1;
            fixation_durations{1,set} = NaN(96,40);
            fixation_durations{2,set} = NaN(96,40);
            [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
            [which_img,novel_vs_repeat] = get_image_numbers(cfg,itmlist,sequence_items,23);
            imgnum = 0;
            for t = 1:length(fixationstats);
                if any(cfg.trl(t).allval == 23) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                    imgnum = imgnum+1;
                    trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == 15);
                    imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == 23)-trial_start+500;%to ignore 1st fixation
                    fixationtimes = fixationstats{t}.fixationtimes;
                    if ~isempty(fixationtimes)
                        tooearly = find(fixationtimes(1,:) <= imgon);
                        fixationtimes(:,tooearly) = [];
                        fixdurs = diff(fixationtimes)+1;
                        if novel_vs_repeat(imgnum) == 1 %then novel
                            fixation_durations{1,set}(which_img(imgnum),1:length(fixdurs))=fixdurs;
                        else
                            fixation_durations{2,set}(which_img(imgnum),1:length(fixdurs))=fixdurs;
                        end
                    end
                end
            end
        end
    end
end
%%
all_durs = cell(1,2);
for set = 1:size(fixation_durations,2);
    all_durs{1} = [all_durs{1}; fixation_durations{1,set}(:,1:40)];
    all_durs{2} = [all_durs{2}; fixation_durations{2,set}(:,1:40)];
end

figure
hold on
plot(nanmean(all_durs{1}(:,1:20)))
errorb(1:20,nanmean(all_durs{1}(:,1:20)),nanstd(all_durs{1}(:,1:20))./sqrt(sum(~isnan(all_durs{1}(:,1:20)))),'color','b')
plot(nanmean(all_durs{2}(:,1:20)),'r')
errorb(1:20,nanmean(all_durs{2}(:,1:20)),nanstd(all_durs{2}(:,1:20))./sqrt(sum(~isnan(all_durs{2}(:,1:20)))),'color','r')
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
xlim([0 21])
legend('Novel','Repeat')
%%
%looking at power for combining fixations 4-10; January 18
meannov = mean(all_durs{1}(:,4:10)');
meannov(isnan(meannov)) = [];

meanrep = mean(all_durs{2}(:,4:10)');
meanrep(isnan(meanrep)) = [];

n = sampsizepwr('t',[mean(meannov) std(meannov)],[mean(meanrep)],0.8);
%% 
%looking at power by analyzing the mean fixation duration by ordinal
%fixation by averaging over all images within a set to create a single data
%point. Only going to take fully completed sets. 
set_data = cell(1,2); % to store average fixations for novel (row 1) and
% repeat (row 2) presentations

for set = 1:size(fixation_durations,2)
    if sum(~isnan(fixation_durations{2,set}(:,1))) >= 64 %then completed all images
        set_data{1} = [set_data{1}; nanmean(fixation_durations{1,set}(:,1:40))];
        set_data{2} = [set_data{2}; nanmean(fixation_durations{2,set}(:,1:40))];
    end
end

figure
hold on
errorbar(mean(set_data{1}(:,1:20)),std(set_data{1}(:,1:20))...
    ./sqrt(size(set_data{1},1)),'b')
errorbar(mean(set_data{2}(:,1:20)),std(set_data{2}(:,1:20))...
    ./sqrt(size(set_data{2},1)),'r')
hold off

xlabel('Fixation Number')
ylabel('Fixation Duration (ms)')

%looking at fixations 1-20 roughly the median number of fixations (I think
%18 is the median)
num_samps = NaN(2,20);
for f = 1:20
    num_samps(1,f) = sampsizepwr('t',[mean(set_data{1}(:,f)),std(set_data{1}(:,f))],...
        mean(set_data{2}(:,f)),0.8);
    num_samps(2,f) = sampsizepwr('t',[mean(set_data{1}(:,f)),std(set_data{1}(:,f))],...
        mean(set_data{2}(:,f)),0.9);
end
    
%% looking only at grouping fixations 2-12 or 4-10 
num_samps_combined = NaN(2,2);
group = {2:12,4:10};
power = [0.8, 0.9];
for p = 1:2
    for g = 1:2
        data1 = set_data{1}(:,group{g});
        data2 = set_data{2}(:,group{g});
        data1=data1(1:end);
        data2=data2(1:end);
    num_samps_combined(p,g) = sampsizepwr('t',[mean(data1),std(data1)],...
        mean(data2),power(p));
    end
end