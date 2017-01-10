%public domain images for examples in paper
%from Manfred

cortexfile = 'MF161214.3';
itemfile = 'xtrImg99.itm';
img_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\FlickrImages\Public Domain Images\LSQ99\';

samprate = 5;%number of ms between samples for ISCAN i.e. 200 Hz
imageX = 800; %horizontal image size
imageY = 600; %vertical image size

%---Import Cortex data file---%
if strcmpi(cortexfile(1:2),'PW')
    cortexfile = ['\\research.wanprc.org\research\Buffalo Lab\Cortex Data\Vivian\' cortexfile];
elseif strcmpi(cortexfile(1:2),'TT')
    cortexfile = ['\\research.wanprc.org\research\Buffalo Lab\Cortex Data\Timmy\' cortexfile];
elseif strcmpi(cortexfile(1:2),'RR')
    cortexfile = ['\\research.wanprc.org\research\Buffalo Lab\Cortex Data\Red\' cortexfile];
elseif strcmpi(cortexfile(1:2),'TO')
    cortexfile = ['\\research.wanprc.org\research\Buffalo Lab\Cortex Data\Tobii\' cortexfile];
elseif strcmpi(cortexfile(1:2),'MF')
    cortexfile = ['\\research.wanprc.org\research\Buffalo Lab\Cortex Data\Manfred\' cortexfile];
end


[time_arr,event_arr,eog_arr,epp_arr,~,~]  = get_ALLdata(cortexfile);

%---Import Item/Condition file informoation---%
% [itmlist,clrchng_locations,first_img_item,imgs,image_names,trialtype] = read_SCMRM_itm_and_cnd_files(itemfile);
ITMFile = ['\\research.wanprc.org\research\Buffalo Lab\eblab\Cortex Programs\SCM Relational Memory Project\Item Files\' itemfile];
CNDFile = '\\research.wanprc.org\research\Buffalo Lab\eblab\Cortex Programs\SCM Relational Memory Project\SCMRM.cnd';

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

first_img_item = [];
for i = 6:size(itmfil,1); %first 5 are hearder, background, etc
    str = textscan(itmfil(i,:),'%s');
    if ~isempty(strfind(str{1}{end},'.bmp')) && isempty(first_img_item)
        first_img_item  = str2num((str{1}{1}));
        break
    end
end

itmlist = zeros(size(cndfil,1)-1,1);
for i = 2:size(cndfil,1);
    str = textscan(cndfil(i,:),'%d');
    itmlist(i-1) = str{1}(end);
end

%probably a more efficient way of doing this since there's only 25 color
%change locations, but it may be easier for processing later
% also get which images were displayed. Should be the same across all items
% and condition sets but what the hey.
imgs = zeros(1,72);
image_names = cell(1,72);
clrchng_locations = cell(1,length(itmlist));
for cnd = 1:length(itmlist)
    if itmlist(cnd) < first_img_item %then it is a clrchng trial
        str = textscan(itmfil(itmlist(cnd)+6,:),'%d');
        clrchng_locations{cnd} =  double(str{1}(4:5)); %location in dva
    else
        str = textscan(itmfil(itmlist(cnd)+6,:),'%s');
        imgnum = str2double(str{1}{end}(14:15));
        if imgs(1,imgnum) == 0; %novel image
            imgs(1,imgnum) = cnd;
            image_names{1,imgnum} = str{1}{end}(10:end);
        else
            imgs(2,imgnum) = cnd;
            image_names{2,imgnum} = str{1}{end}(10:end);
        end
    end
end


%---Get the calibration for the eye data---%
%essentially the same as all other task's we use
numrpt = size(event_arr,2);
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if itmlist(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop)-1000) < first_img_item
        if size(find(event_arr(:,rptlop) == 3)) ~=0 %if clrchng trial was rewarded
            perbegind = find(event_arr(:,rptlop) == 24);%was originally 23, changed this and begtimdum line below to optimize
            perendind = find(event_arr(:,rptlop) == 24);
            if length( perbegind) > 1
                perbegind = perbegind(2);
                perendind = perendind(2);
            end
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
            begtimdum = time_arr(perbegind,rptlop)-100;
            endtimdum = time_arr(perendind,rptlop);
            if endtimdum > begtimdum
                valrptcnt = valrptcnt + 1;
                clrchgind(valrptcnt)=rptlop;
                per(valrptcnt).begsmpind = begtimdum;
                per(valrptcnt).endsmpind = endtimdum;
                per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                per(valrptcnt).blk = event_arr(blknumind,rptlop);
                per(valrptcnt).allval = event_arr(:,rptlop);
                per(valrptcnt).alltim = time_arr(:,rptlop);
                per(valrptcnt).event = rptlop;
            end
        end
    end
end

%Don't keep first 6 successful trials. These trials are all central
%cailbration trials for calibration offset correct
per(1:5) = [];

clear cnd
numrpt = size(per,2);
cnd = zeros(1,numrpt);
for rptlop = 1:numrpt
    cnd(rptlop)=per(rptlop).cnd-1000;
end

% Create structures x and y of the corresponding average eye data for each trial
% instance (l) of each condition (k)
spacex = [-12,-6,0,6,12];%what actually gets displayed
spacey = [-8,-4,0,4,8];%what actually gets displayed
x = cell(length(spacex),length(spacey));%---For Calibration with Eye tracking data with cp2tform---%
y = cell(length(spacex),length(spacey));
control = NaN(length(cnd),2);
clr = 'rgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmk';
figure
hold on
for k = 1:length(cnd)
    control(k,:) = clrchng_locations{cnd(k)}';
    
    xi = find(clrchng_locations{cnd(k)}(1) == spacex);
    yi = find(clrchng_locations{cnd(k)}(2) == spacey);
    eyeind = floor(((per(k).begsmpind-1000)/samprate)*2):(floor((per(k).endsmpind-1000)/samprate))*2;
    evenind = eyeind(logical(~rem(eyeind,2)));
    oddind =  eyeind(logical(rem(eyeind,2)));
    x{xi,yi} = [x{xi,yi} mean(eog_arr(oddind,per(k).event))];
    y{xi,yi} = [y{xi,yi} mean(eog_arr(evenind,per(k).event))];
    plot(mean(eog_arr(oddind,per(k).event)),mean(eog_arr(evenind,per(k).event)),[clr(xi*yi) '+'])
end
title(['Calibration transformation for ' cortexfile(end-9:end)])
% save_and_close_fig(figure_dir,['Calibration_'  cortexfile(end-9:end-2)])

%Test for errors%
count = zeros(length(spacey),length(spacex));
for xi = 1:length(spacex);
    for yi = 1:length(spacey);
        count(yi,xi) = sum(control(:,1) == spacex(xi) & control(:,2) == spacey(yi));
    end
end
if any(count < 10);
    disp('Calibration trial analysis incomplete or error')
    disp('Check number of calibration pionts or task not finished')
end

clear meanx meany
for k=1:numel(x)
    xss = x{k};
    low = mean(xss)-std(xss);
    high = mean(xss)+std(xss);
    xss(xss < low) = [];
    xss(xss > high) = [];
    meanx(k)=median(xss);
end
for k=1:numel(y)
    yss = y{k};
    low = mean(yss)-std(yss);
    high = mean(yss)+std(yss);
    yss(yss < low) = [];
    yss(yss > high) = [];
    meany(k)=median(y{k});
end

controlx = [];
controly = [];
for i = 1:length(spacex);
    for ii = 1:length(spacey);
        controly = [controly spacey(i)];
        controlx = [controlx spacex(ii)];
    end
end

%want to change this to MSE estimate with different calibration functions,
%probably affine and polynomial 3 or 4
tform = cp2tform([controlx' controly'], [meanx' meany'],'affine');
tform.forward_fcn = tform.inverse_fcn;

newx = [];
newy = [];
figure
hold on
for i = 1:length(controlx);
    plot(controlx(i),controly(i),'r+')
    [x,y] = tformfwd(tform,meanx(i),meany(i));
    plot(x,y,'*b')
    newx(i) = x;
    newy(i) = y;
end
title(['Calibration transformation for ' cortexfile(end-9:end)])
xlim([-17.5 17.5])
ylim([-12.5 12.5])
% save_and_close_fig(figure_dir,['Estimated_Calibration_'  cortexfile(end-9:end-2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Import Image Viewing Data---%

new_eog_arr = [];
if ~isempty(epp_arr);%if there is pupil data
    new_epp_arr = [];
end

numrpt = size(event_arr,2);
new_eog_arr=[];
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if itmlist(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop)-1000) >= first_img_item
        if size(find(event_arr(:,rptlop) == 23)) ~=0 %if image actually turned on
            perbegind = find(event_arr(:,rptlop) == 100,1,'first'); %eye data on
            perendind = find(event_arr(:,rptlop) == 101,1,'first'); %eye data off
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
            begtimdum = time_arr(perbegind,rptlop);
            endtimdum = time_arr(perendind,rptlop);
            if endtimdum > begtimdum
                valrptcnt = valrptcnt + 1;
                vpcind(valrptcnt)=rptlop;
                per(valrptcnt).begsmpind = begtimdum;
                per(valrptcnt).endsmpind = endtimdum;
                per(valrptcnt).begpos = 1;
                per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                per(valrptcnt).blk = event_arr(blknumind,rptlop);
                per(valrptcnt).allval = event_arr(:,rptlop);
                per(valrptcnt).alltim = time_arr(:,rptlop);
                new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
                new_epp_arr = cat(2,new_epp_arr,epp_arr(:,rptlop));
            end
        end
    end
end

%---get eye data for only when fixation cross or picture is displayed---%
eyedat = cell(1,length(per));
if ~isempty(epp_arr);
    pupildata = cell(1,length(per)); %horizontal pupil diameter
end
for trlop=1:size(per,2)
    trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
    horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
    vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
    
    picstart=1*samprate;
    picend=per(trlop).endsmpind-per(trlop).begsmpind;%added in case sometimes get weird
    %indexing artifacts that are off by 1 index due to cortex having a
    %clock speed with a 1 ms resoultion and the eye data collected at a 5 ms
    %resoultion
    
    try
        eyedat{trlop}(1,:) = horeog(round(picstart/samprate):floor(picend/samprate));
        eyedat{trlop}(2,:) = vrteog(round(picstart/samprate):floor(picend/samprate));
    catch
        if picend > 20000
            disp(['EOG overflow. Image Trial ' num2str(trlop) '/' num2str(size(per,2))])
            eyedat{trlop}(1,:) = horeog;
            eyedat{trlop}(2,:) = vrteog;
        else
            try
                eyedat{trlop}(1,:) = horeog(round(picstart/samprate):floor(picend/samprate)-1);
                eyedat{trlop}(2,:) = vrteog(round(picstart/samprate):floor(picend/samprate)-1);
            catch
                error('Unknown reason EOG and pic duration do not line up')
            end
        end
    end
    try
        trlepp=new_epp_arr(~isnan(new_epp_arr(:,trlop)),trlop); % epp for this trial
        pupildata{trlop} = trlepp(2:2:2*floor(picend/samprate));%odd indexes contains nothing but noise
    catch
        if picend > 20000
            pupildata{trlop} = trlepp(2:2:end); %epp for this trial
        else
            try
                pupildata{trlop} = trlepp(2:2:(2*floor(picend/samprate)-1));%odd indexes contains nothing but noise
            catch
                %error('Unknown reason EOG and pic duration do not line up')
            end
        end
    end
end

%---Recalibrate and automatically scale eye data---%
fixationstats = {};
raw_eyedat = eyedat;
for eye = 1:length(eyedat)
    x = eyedat{eye}(1,:);
    y = eyedat{eye}(2,:);
    [x,y] = tformfwd(tform,x,y);
    
    x = 24*x; %convert from cortex dva to pixels
    y = 24*y; %convert from cortex dva to pixels
    x = x+imageX/2;
    y = y+imageY/2;
    
    raw_eyedat{eye} = [x;y]; 
    y(x < -24) = NaN;
    x(x < -24) = NaN;
    y(x > imageX+24) = NaN;
    x(x > imageX+24)= NaN;
    x(y < -24) = NaN;
    y(y < -24) = NaN;
    x(y > imageY+24) = NaN;
    y(y > imageY+24) = NaN;
    
    eyedat{eye} = [x;y];
end
%%
for img = 1:72
if img >= 10
im = imread([img_dir 'S99I' num2str(img) '.bmp']);
else
im = imread([img_dir 'S99I0' num2str(img) '.bmp']);

end

eye_start = per(img).begsmpind;
imgon = per(img).alltim(per(img).allval == 23);
imgoff = per(img).alltim(per(img).allval == 24); 
eye_img_start = ceil(imgon-eye_start)/5;
eye_img_end = ceil(imgoff-eye_start)/5;

if eye_img_end-eye_img_start < 1000 %5000 ms 
disp('what')
end 
x = eyedat{img}(1,eye_img_start:eye_img_end);
y = 600-eyedat{img}(2,eye_img_start:eye_img_end);
figure
imshow(im)
hold on
plot(x,y,'y')


end