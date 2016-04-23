% test to see if ListSQ was running correctly without a monkey present
cd('C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Old and Misc Code\');
[time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata('TT141020.2');

ITMFile = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Item and Conditions Files\ListSQ07.itm';

CNDFile = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Item and Conditions Files\ListSQ07.cnd';


itmfil=[];
h =fopen(ITMFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(itmfil)
        if length(tline)>size(itmfil,2)
            tline=tline(1:size(itmfil,2));
        end
    end
    tline = [tline ones(1,(size(itmfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        itmfil=[itmfil; tline];
    else
        break
    end
end
fclose(h);

cndfil=[];
h=fopen(CNDFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(cndfil)
        if length(tline)>size(cndfil,2)
            tline=tline(1:size(cndfil,2));
        end
    end
    tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        cndfil=[cndfil; tline];
    else
        break
    end
end
fclose(h);

itmlist = zeros(size(cndfil,1)-1,1);
for i = 2:size(cndfil,1);
    str = textscan(cndfil(i,:),'%d');
    itmlist(i-1) = str{1}(end);
end

all_trials = [];
numrpt = size(event_arr,2);
new_eog_arr=[];
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if event_arr((find(event_arr(:,rptlop)>500,1,'last')),rptlop)-500 == 2 %1st block is color change always
        if size(find(event_arr(:,rptlop) == 151)) ~=0 %fine if trial was rewarded
            perbegind = find(event_arr(:,rptlop) == 100); % eye data starts at 100; might have predictive looking
            perendind = find(event_arr(:,rptlop) == 101); % eye data stops collecting after rewards so can stop here
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
            begtimdum = time_arr(perbegind,rptlop);
            endtimdum = time_arr(perendind(end),rptlop);
            
            cnd = event_arr(cndnumind,rptlop)-1000;
            cndline = textscan(cndfil(cnd+1,:),'%d');
            trial_itm = cndline{1}(end);
            
            if trial_itm < 20
                
                test0_on = find(event_arr(:,rptlop) == 23);
                test0_off = find(event_arr(:,rptlop) == 24);
                test1_on = find(event_arr(:,rptlop) == 25);
                test1_off = find(event_arr(:,rptlop) == 26);
                test2_on = find(event_arr(:,rptlop) == 27);
                test2_off = find(event_arr(:,rptlop) == 28);
                test3_on = find(event_arr(:,rptlop) == 29);
                test3_off = find(event_arr(:,rptlop) == 30);
                
                %for looking at trials with break fixations
                test0_on = time_arr(test0_on,rptlop)-begtimdum;
                if isempty(test0_off);
                    test0_off = NaN;
                else
                    test0_off = time_arr(test0_off,rptlop)-begtimdum;
                end
                if isempty(test1_off);
                    test1_off = NaN;
                else
                    test1_off = time_arr(test1_off,rptlop)-begtimdum;
                end
                if isempty(test2_off);
                    test2_off = NaN;
                else
                    test2_off = time_arr(test2_off,rptlop)-begtimdum;
                end
                if isempty(test3_off);
                    test3_off = NaN;
                else
                    test3_off = time_arr(test3_off,rptlop)-begtimdum;
                end
                if isempty(test1_on);
                    test1_on = NaN;
                else
                    test1_on = time_arr(test1_on,rptlop)-begtimdum;
                end
                if isempty(test2_on);
                    test2_on = NaN;
                else
                    test2_on = time_arr(test2_on,rptlop)-begtimdum;
                end
                if isempty(test3_on);
                    test3_on = NaN;
                else
                    test3_on = time_arr(test3_on,rptlop)-begtimdum;
                end
                if endtimdum > begtimdum
                    valrptcnt = valrptcnt + 1;
                    per(valrptcnt).begsmpind = begtimdum;
                    per(valrptcnt).endsmpind = endtimdum(1);
                    per(valrptcnt).begpos = 1;
                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop)-1000;
                    per(valrptcnt).blk = event_arr(blknumind,rptlop)-500;
                    per(valrptcnt).allval = event_arr(:,rptlop);
                    per(valrptcnt).alltim = time_arr(:,rptlop);
                    per(valrptcnt).event = rptlop;
                    per(valrptcnt).test= ...
                        [test0_on test0_off;...
                        test1_on test1_off;...
                        test2_on test2_off;...
                        test3_on test3_off];
                    new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
                    all_trials = [all_trials [event_arr(cndnumind,rptlop)-1000; 1]];
                end
            else
                valrptcnt = valrptcnt + 1;
                cross_on = find(event_arr(:,rptlop) == 9);
                cross_off = find(event_arr(:,rptlop) == 10);
                img_on = find(event_arr(:,rptlop) == 23);
                img_off = find(event_arr(:,rptlop) == 24);
                
                per(valrptcnt).crossfixdur = time_arr(cross_off,rptlop)-time_arr(cross_on,rptlop);
                per(valrptcnt).cnd = event_arr(cndnumind,rptlop)-1000;
                per(valrptcnt).imgdur = time_arr(img_off,rptlop)-time_arr(img_on,rptlop);
            end
        end
    end
end

imgdur = [];
crossdur = [];
for trial = 1:length(per);
    imgdur = [imgdur per(trial).imgdur];
    crossdur = [crossdur per(trial).crossfixdur];
end
  %%
crosshair_locations = cell(2,size(itmlist,1)-1); %x in row 1 y in row 2
for cnd = 1:size(itmlist,1);
    for itm = 1:size(itmlist,2);
        itmline = itmlist(cnd,itm)+6;
        str = textscan(itmfil(itmline,:),'%d');
        crosshair_locations{1,cnd}(1,itm) = str{1}(4); %xpos
        crosshair_locations{1,cnd}(2,itm) = str{1}(5); %ypos
    end
    if itmlist(cnd,itm) < 20 && itmlist(cnd,itm) > 3
        crosshair_locations{2,cnd} =  itmlist(cnd,itm);
    end
end
crosscnd = cell2mat(crosshair_locations(2,:));

%%
imgs = zeros(2,96);
img_trial = 1;
for trial = 1:length(per);
    cndline = textscan(cndfil(per(trial).cnd+1,:),'%d');
    imgnum = cndline{1}(end)-19;
    if imgnum >= 1;
        if imgs(1,imgnum) == 0;
        imgs(1,imgnum) = img_trial;
        else
             imgs(2,imgnum) = img_trial;
        end
        img_trial = img_trial+1;
    end
end
%%
cnd = zeros(1,501);
for trial = 1:length(per);
   cnd(per(trial).cnd) =   cnd(per(trial).cnd)+1;
end

%%
type = [];
for trial =sum 1:length(per);
    type(trial) = per(trial).cnd;
end
