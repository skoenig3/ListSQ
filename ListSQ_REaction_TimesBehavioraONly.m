%ListSQ Sequence Reaction times Analysis

% data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\';
% %
% listsq_files = {'PW140725_3-sorted.nex','PW140728_2_3-sorted.nex','PW140729_3-sorted.nex',...
%     'PW140730_3-sorted.nex','PW140801_3-sorted.nex','PW140805_3-sorted.nex','PW140806_3-sorted.nex',...
%     'PW140825_3-sorted.nex','PW140826_3-sorted.nex','PW140829_3-sorted.nex','PW140908_3-sorted.nex',...
%     'PW140910_3-sorted.nex','PW140915_3-sorted.nex','PW140101_3-sorted.nex',...
%     'PW1floor(424/block_size)006_3-sorted.nex','PW1floor(424/block_size)007_3-sorted.nex','PW1floor(424/block_size)013_3-sorted.nex','PW1floor(424/block_size)014_3-sorted.nex',...
%     'PW1floor(424/block_size)009_3-sorted.nex','PW1floor(424/block_size)015_3-sorted.nex','PW1floor(424/block_size)008_3-sorted.nex','PW1floor(424/block_size)010_3-sorted.nex',...
%     'PW1floor(424/block_size)024_3-sorted.nex','PW1floor(424/block_size)027_3-sorted.nex','PW1floor(424/block_size)028_3-sorted.nex','PW1floor(424/block_size)029_3-sorted.nex',...
%     'PW1floor(424/block_size)023_3-sorted.nex','PW1floor(424/block_size)031_3-sorted.nex','PW1floor(424/block_size)103_3-sorted.nex','PW1floor(424/block_size)105_3-sorted.nex',...
%     'PW1floor(424/block_size)106_3-sorted.nex','PW1floor(424/block_size)110_3-sorted.nex','PW150121_3-sorted.nex','PW150122_3-sorted.nex',...
%     'PW150123_3-sorted.nex','PW150126_3-sorted.nex','PW150127_3-sorted.nex','PW150204_3-sorted.nex',...
%     'PW150205_3-sorted.nex'};

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';

listsq_files = {'TO151203_3-sorted.nex','TO151204_3-sorted.nex','TO151211_3-sorted.nex',...
                'TO151214_3-sorted.nex','TO151217_3-sorted.nex','TO151218_3-sorted.nex',...
                'TO151222_3-sorted.nex','TO151228_3-sorted.nex','TO160105_3-sorted.nex',...
                'TO160107_3-sorted.nex','TO160108_3-sorted.nex','TO160114_3-sorted.nex',...
                'TO160115_3-sorted.nex','TO160121_3-sorted.nex','TO160122_3-sorted.nex',...
                'TO160126_3-sorted.nex','TO160127_3-sorted.nex','TO160201_3-sorted.nex',...
                'TO160202_3-sorted.nex','TO151208_3-sorted.nex','TO151216_3-sorted.nex',...
                'TO160113_3-sorted.nex'};

min_rt = 138;%ms; 138 ms for tobii, 156 ms for vivian

fixwin = 5;
all_time_to_fixation = cell(1,length(listsq_files));
all_extrafixations = cell(1,length(listsq_files));
all_fixation_accuracy = cell(1,length(listsq_files));
all_which_sequence = cell(1,length(listsq_files));
for ls = 1:length(listsq_files)
    load([data_dir,listsq_files{ls}(1:end-11) '-preprocessed'],'cfg','item_set','fixationstats');%for vivian
    
    if strcmpi(listsq_files{ls}(1:2),'TO')
        if strcmpi('ListSQ04.itm',item_set)
            cnd_file = 'ListSQ53.cnd';
        elseif strcmpi('ListSQ05.itm',item_set)
            cnd_file = 'ListSQ10.cnd';
        elseif strcmpi('ListSQ06.itm',item_set)
            cnd_file = 'ListSQ41.cnd';
        else
            cnd_file = [item_set(1:end-3) 'cnd'];
        end
    end
    
    if cfg.trl(end).cnd > 110
        
        %get important task specific information
        if strcmpi(listsq_files{ls}(1:2),'TO')
            [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set,cnd_file);
        else
            [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
            
        end
        overlap = find((sequence_locations{1}(1,:) ==  sequence_locations{2}(1,:)) & ...
            (sequence_locations{1}(2,:) ==  sequence_locations{2}(2,:)));
        
        num_trials = length(cfg.trl);
        
        %preallocate space and parallel structure of cfg
        successful_sequence_trials = NaN(1,length(cfg.trl));
        which_sequence = NaN(1,length(cfg.trl));
        for t = 1:length(cfg.trl);
            if sum(cfg.trl(t).allval == 3) >= 6; %in which sequence trials were rewarded
                which_sequence(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
                successful_sequence_trials(t) = t;
            end
        end
        successful_sequence_trials = laundry(successful_sequence_trials);
        which_sequence = laundry(which_sequence);
        
        fixationstats = fixationstats(successful_sequence_trials);
        cfg.trl = cfg.trl(successful_sequence_trials);
        
        
        time_to_fixation = NaN(length(fixationstats),4); %how fast did they get to it the first time around
        fixation_accuracy = NaN(length(fixationstats),4); %how far off
        extrafixations = NaN(length(fixationstats),4); %how many noncross hair fixations do they mak
        
        num_trials = length(cfg.trl);
        for trial = 1:num_trials
            locs = sequence_locations{which_sequence(trial)};
            
            %convert to DVA for this analysis
            locs(1,:) = (locs(1,:)-400)/24;
            locs(2,:) = (locs(2,:)-300)/24;
            fixationstats{trial}.fixations(1,:) = (fixationstats{trial}.fixations(1,:)-400)/24;
            fixationstats{trial}.fixations(2,:) = (fixationstats{trial}.fixations(2,:)-300)/24;
            
            event_codes = cfg.trl(trial).allval;
            event_codes(event_codes == 100)= 0;
            event_codes(1) = 100;%eye data starts for recording right away
            event_times = cfg.trl(trial).alltim;
            
            trialdata = analyze_sequence_trial(fixationstats{trial},locs,fixwin,...
                event_codes,event_times);
            
            time_to_fixation(trial,:) = trialdata.t2f;
            extrafixations(trial,:) = trialdata.extrafixations;
            fixation_accuracy(trial,:) =  trialdata.accuracy;
            
        end
        
        all_time_to_fixation{ls} = time_to_fixation;
        all_extrafixations{ls} = extrafixations;
        all_fixation_accuracy{ls} = fixation_accuracy;
        all_which_sequence{ls} = which_sequence;
    end
end
%%
block_size = 5;
avg_rt = NaN(length(listsq_files),floor(424/block_size));
predict_pct = NaN(length(listsq_files),floor(424/block_size));
extra = NaN(length(listsq_files),floor(424/block_size));
accuracy = NaN(length(listsq_files),floor(424/block_size));

for ls = 1:length(listsq_files)
    time_to_fixation = all_time_to_fixation{ls}(:,2:4);
    extrafixations = all_extrafixations{ls}(:,2:4);
    fixation_accuracy = all_fixation_accuracy{ls}(:,2:4);
    for blk = 1:floor(size(time_to_fixation)/block_size)
        ind = block_size*(blk-1)+1:block_size*blk;
        rt = time_to_fixation(ind,:);
        avg_rt(ls,blk) = nanmean(rt(:));
        
        predict_pct(ls,blk) = 100*sum(rt(:) < min_rt)/sum(~isnan(rt(:)));
        
        ex = extrafixations(ind,:);
        extra(ls,blk) = nanmean(ex(:));
        
        fa = fixation_accuracy(ind,:);
        accuracy(ls,blk) = nanmean(fa(:));
        
    end
end

figure
subplot(2,2,1)
hold on
bar(nanmean(accuracy))
errorb(nanmean(accuracy),nanstd(accuracy)./sqrt(sum(~isnan(accuracy))))
hold off
xlim([0 floor(424/block_size)])
xlabel(['Block # (' num2str(block_size) ' trials)'])
ylabel('Fixation Accuarcy (dva)')
title('Fixation accuracy by Block')

subplot(2,2,2)
hold on
bar(nanmean(extra))
errorb(nanmean(extra),nanstd(extra)./sqrt(sum(~isnan(extra))))
hold off
xlim([0 floor(424/block_size)])
xlabel(['Block # (' num2str(block_size) ' trials)'])
ylabel('# of Extra Fixations)')
title('Extra Fixations by Block')

subplot(2,2,3)
hold on
bar(nanmean(avg_rt))
errorb(nanmean(avg_rt),nanstd(avg_rt)./sqrt(sum(~isnan(avg_rt))))
hold off
xlim([0 floor(424/block_size)])
xlabel(['Block # (' num2str(block_size) ' trials)'])
ylabel('Average Reaction time (ms)')
title('Reaction time by Block')
ylim([150 250])

subplot(2,2,3)
hold on
bar(nanmean(avg_rt))
errorb(nanmean(avg_rt),nanstd(avg_rt)./sqrt(sum(~isnan(avg_rt))))
hold off
xlim([0 floor(424/block_size)])
xlabel(['Block # (' num2str(block_size) ' trials)'])
ylabel('Average Reaction time (ms)')
title('Reaction time by Block')

subplot(2,2,4)
hold on
bar(nanmean(predict_pct))
errorb(nanmean(predict_pct),nanstd(predict_pct)./sqrt(sum(~isnan(predict_pct))))
hold off
xlim([0 floor(424/block_size)])
xlabel(['Block # (' num2str(block_size) ' trials)'])
ylabel('% Predicted')
title('Prediction Rate by Block')

subtitle(['Monkey: ' listsq_files{1}(1:2) ' n = ' num2str(length(listsq_files))])

for ls = 1:length(listsq_files)
    100*sum(all_time_to_fixation{ls} < 50)./sum(~isnan(all_time_to_fixation{ls}))
end

avg_rt_1 = NaN(length(listsq_files),floor(424/block_size));
predict_pct_1 = NaN(length(listsq_files),floor(424/block_size));
extra_1 = NaN(length(listsq_files),floor(424/block_size));
accuracy_1= NaN(length(listsq_files),floor(424/block_size));

avg_rt_2 = NaN(length(listsq_files),floor(424/block_size));
predict_pct_2 = NaN(length(listsq_files),floor(424/block_size));
extra_2 = NaN(length(listsq_files),floor(424/block_size));
accuracy_2= NaN(length(listsq_files),floor(424/block_size));

for ls = 1:length(listsq_files)
    time_to_fixation = all_time_to_fixation{ls}(all_which_sequence{ls} == 1,2:4);
    extrafixations = all_extrafixations{ls}(all_which_sequence{ls} == 1,2:4);
    fixation_accuracy = all_fixation_accuracy{ls}(all_which_sequence{ls} == 1,2:4);
    for blk = 1:floor(size(time_to_fixation)/block_size)
        ind = block_size*(blk-1)+1:block_size*blk;
        rt = time_to_fixation(ind,:);
        avg_rt_1(ls,blk) = nanmean(rt(:));
        
        predict_pct_1(ls,blk) = 100*sum(rt(:) < min_rt)/sum(~isnan(rt(:)));
        
        ex = extrafixations(ind,:);
        extra_1(ls,blk) = nanmean(ex(:));
        
        fa = fixation_accuracy(ind,:);
        accuracy_1(ls,blk) = nanmean(fa(:));
        
    end
end

for ls = 1:length(listsq_files)
    time_to_fixation = all_time_to_fixation{ls}(all_which_sequence{ls} == 2,2:4);
    extrafixations = all_extrafixations{ls}(all_which_sequence{ls} == 2,2:4);
    fixation_accuracy = all_fixation_accuracy{ls}(all_which_sequence{ls} == 2,2:4);
    for blk = 1:floor(size(time_to_fixation)/block_size)
        ind = block_size*(blk-1)+1:block_size*blk;
        rt = time_to_fixation(ind,:);
        avg_rt_2(ls,blk) = nanmean(rt(:));
        
        predict_pct_2(ls,blk) = 100*sum(rt(:) < min_rt)/sum(~isnan(rt(:)));
        
        ex = extrafixations(ind,:);
        extra_2(ls,blk) = nanmean(ex(:));
        
        fa = fixation_accuracy(ind,:);
        accuracy_2(ls,blk) = nanmean(fa(:));
        
    end
end


figure
subplot(2,2,1)
hold on
errorbar(nanmean(accuracy_1),nanstd(accuracy_1)./sqrt(sum(~isnan(accuracy_1))))
errorbar(nanmean(accuracy_2),nanstd(accuracy_2)./sqrt(sum(~isnan(accuracy_2))),'r')
hold off
xlim([0 floor(424/block_size)/2])
xlabel(['Block # (' num2str(block_size) ' trials)'])
ylabel('Fixation Accuarcy (dva)')
title('Fixation accuracy by Block')

subplot(2,2,2)
hold on
errorbar(nanmean(extra_1),nanstd(extra_1)./sqrt(sum(~isnan(extra_1))))
errorbar(nanmean(extra_2),nanstd(extra_2)./sqrt(sum(~isnan(extra_2))),'r')
hold off
xlim([0 floor(424/block_size)/2])
xlabel(['Block # (' num2str(block_size) ' trials)'])
ylabel('# of Extra Fixations)')
title('Extra Fixations by Block')

subplot(2,2,3)
hold on
errorbar(nanmean(avg_rt_1),nanstd(avg_rt_1)./sqrt(sum(~isnan(avg_rt_1))))
errorbar(nanmean(avg_rt_2),nanstd(avg_rt_2)./sqrt(sum(~isnan(avg_rt_2))),'r')
hold off
xlim([0 floor(424/block_size)/2])
xlabel(['Block # (' num2str(block_size) ' trials)'])
ylabel('Average Reaction time (ms)')
title('Reaction time by Block')
ylim([125 250])

subplot(2,2,4)
hold on
errorbar(nanmean(predict_pct_1),nanstd(predict_pct_1)./sqrt(sum(~isnan(predict_pct_1))))
errorbar(nanmean(predict_pct_2),nanstd(predict_pct_2)./sqrt(sum(~isnan(predict_pct_2))),'r')
hold off
xlim([0 floor(424/block_size)/2])
xlabel(['Block # (' num2str(block_size) ' trials)'])
ylabel('% Predicted')
title('Prediction Rate by Block')

subtitle(['Monkey: ' listsq_files{1}(1:2) ' n = ' num2str(length(listsq_files))])
%% Only look at the "learend sequence"
thresh = 0.1;%10% or > of rts are predictive
subset = cell(2,10);
%figure
hold all
index = [1 1];
for ls = 1:length(listsq_files)
    rts = all_time_to_fixation{ls};
    seq = all_which_sequence{ls};
    
    rt1 = rts(seq == 1,2:4);
    predicted = find(sum(rt1 < min_rt)./sum(~isnan(rt1)) > thresh);
    
    for p = 1:length(predicted)
        %         plot(rt1(:,predicted(p)))
        subset{1,index(1)} = rt1(:,predicted(p));
        index(1) = index(1)+1;
    end
    
    
    rt1 = rts(seq == 2,2:4);
    predicted = find(sum(rt1 < min_rt)./sum(~isnan(rt1)) > thresh);
    
    for p = 1:length(predicted)
        %plot(rt1(:,predicted(p)))
        subset{2,index(2)} = rt1(:,predicted(p));
        index(2) = index(2)+1;
    end
end

avg_rt_1 = NaN(length(subset),20);
predict_pct_1 = NaN(length(subset),20);

avg_rt_2 = NaN(length(subset),21);
predict_pct_2 = NaN(length(subset),20);


for ls = 1:size(subset,2)
    time_to_fixation = subset{1,ls};
    if ~isempty(time_to_fixation)
        for blk = 1:floor(size(time_to_fixation)/block_size)
            ind = block_size*(blk-1)+1:block_size*blk;
            rt = time_to_fixation(ind,:);
            avg_rt_1(ls,blk) = nanmean(rt(:));
            predict_pct_1(ls,blk) = 100*sum(rt(:) < min_rt)/sum(~isnan(rt(:)));
        end
    end
    
    time_to_fixation = subset{2,ls};
    if ~isempty(time_to_fixation)
        for blk = 1:floor(size(time_to_fixation)/block_size)
            ind = block_size*(blk-1)+1:block_size*blk;
            rt = time_to_fixation(ind,:);
            avg_rt_2(ls,blk) = nanmean(rt(:));
            predict_pct_2(ls,blk) = 100*sum(rt(:) < min_rt)/sum(~isnan(rt(:)));
        end
    end
end

figure
subplot(1,2,1)
hold on
errorbar(nanmean(avg_rt_1),nanstd(avg_rt_1)./sqrt(sum(~isnan(avg_rt_1))))
errorbar(nanmean(avg_rt_2),nanstd(avg_rt_2)./sqrt(sum(~isnan(avg_rt_2))),'r')
hold off
xlim([0 floor(424/block_size)/2])
xlabel(['Block # (' num2str(block_size) ' trials)'])
ylabel('Average Reaction time (ms)')
title('Reaction time by Block')
legend('Seq 1','Seq 2')
ylim([125 250])

subplot(1,2,2)
hold on
errorbar(nanmean(predict_pct_1),nanstd(predict_pct_1)./sqrt(sum(~isnan(predict_pct_1))))
errorbar(nanmean(predict_pct_2),nanstd(predict_pct_2)./sqrt(sum(~isnan(predict_pct_2))),'r')
hold off
xlim([0 floor(424/block_size)/2])
xlabel(['Block # (' num2str(block_size) ' trials)'])
ylabel('% Predicted')
title('Prediction Rate by Block')
legend('Seq 1','Seq 2')

