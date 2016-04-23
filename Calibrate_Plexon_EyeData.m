function [tform] = Calibrate_Plexon_EyeData(data_dir,cch25_file,figure_dir)
% written July 26, 2014 by Seth Konig

% Import and process cch25f data through plexon and then outputs
% transformation function (tform). Function will also plot calibration data
% to visual confirm that the data is good.

% Plexon eye tracking values (i.e. voltage values)
% are different from what cortex recieves though the they seem to map
% nearly linearly but not 100% sure...likely due to differences in cable
% length or other differences in inherent resistnaces throughout the
% system.
%
% INPUTS:
%   1) data_dir: the directory containing the 25 point calibration plexon data file
%   2) cch25_file: the 25 point calibration plexon data file
%
% OUPUTS:
%   1) tform: calibration function

% define the trials solely based on the marker events
cfg=[];
cfg.dataset       = [data_dir cch25_file];
cfg.trialfun      = 'trialfuncch25';
cfg = ft_definetrial(cfg);

% read the data from file and preprocess them
% only need x and y data for this task unless there's interesting LFP and
% spike data
cfg.channel       = {'X','Y'};
cfg.datatype      = 'continuous';
data = getPlexonTrialData(cfg);

item_cnd_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Item and Conditions Files\';

%%----Color Change Calibration----%%

ITMFile = [item_cnd_dir 'cch25.itm'];
CNDFile =   [item_cnd_dir 'cch25.cnd'];
% this is different becasue the spacing is different and I don't have
% a new item file on the network for the new spacing
ind_spacex = [-6,-3,0,3,6]; %whats currently in the file
ind_spacey = [-6,-3,0,3,6];%whats currently in the file
spacex = [-12,-6,0,6,12];%what actually gets displayed
spacey = [-8,-4,0,4,8];%what actually gets displayed

%---read in itm file---%
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

%---read in cnd file---%
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

%---Determine which cnds go with which items---%
itmlist = zeros(size(cndfil,1)-1,1);
for i = 2:size(cndfil,1);
    str = textscan(cndfil(i,:),'%d');
    itmlist(i-1) = str{1}(end);
end

%sometimes 1st trial gets cutoff depending on when recording sessions
%starts
if ~isempty(cfg.trl(1).allval == 15)
    cfg.trl(1) = [];
    data(1).values(1) = [];
    data(2).values(1) = [];
end

% which conditions were displayed
clear cnd
numrpt = length(data(1).values);
cnd = zeros(1,numrpt);
for rptlop = 1:numrpt
    cnd(rptlop)=cfg.trl(rptlop).cnd-1000;
end

%---For Calibration with Eye tracking data with cp2tform---%
% Create structures x and y of the corresponding average eye data for each trial for each condition (k)
x = cell(length(spacex),length(spacey));
y = cell(length(spacex),length(spacey));
control = NaN(length(cnd),2);
clr = ['rgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmk'];
figure
hold on
for k = 1:length(cnd)
    C = textscan(itmfil(itmlist(cnd(k))+5,:),'%d');
    control(k,:) = C{1}(9:10)';
    
    xi = find(C{1}(9) == ind_spacex);
    yi = find(C{1}(10) == ind_spacey);
    eye_x = data(1).values{k};
    eye_y = data(2).values{k};
    x{xi,yi} = [x{xi,yi} mean(eye_x(end-100:end))];
    y{xi,yi} = [y{xi,yi} mean(eye_y(end-100:end))];
    plot(x{xi,yi},y{xi,yi},[clr(xi*yi) '+']) %plot to make sure it looks ok
end
save_and_close_fig(figure_dir,[cch25_file(1:10) '-calibration data from all trials'])


%---Test to make sure enough trials where displayed---%
count = cellfun(@numel,x);
if any(count < 7);
    disp('Calibration trial analysis incomplete or error')
    disp('Check number of calibration pionts or task not finished')
end

%---remove outliers and get median or mena values for each calibration point---%
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

%---setup the control structure to be equivalent to the collected structure---%
controlx = [];
controly = [];
for i = 1:length(spacex);
    for ii = 1:length(spacey);
        controly = [controly spacey(i)];
        controlx = [controlx spacex(ii)];
    end
end

tform = get_calibration_fcn([controlx; controly],[meanx; meany]);

tform.forward_fcn = tform.inverse_fcn;

%---check to make sure the transforemed calibrtion looks reasonable---%
% if not try a different transformation above
figure
hold on
for i = 1:length(controlx);
    plot(controlx(i),controly(i),'r+')
    [x,y] = tformfwd(tform,meanx(i),meany(i));
    plot(x,y,'*b')
end
title(['Calibration transformation for '  cch25_file(1:10)])
xlim([-17.5 17.5])
ylim([-12.5 12.5])
save_and_close_fig(figure_dir,[cch25_file(1:10) '-accuracy of  calibration'])
