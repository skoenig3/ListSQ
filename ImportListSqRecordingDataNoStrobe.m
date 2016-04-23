function ImportListSQRecordingDataNoStrobe(data_dir,cch25_file,listsq_file,item_set,multiunit,figure_dir)
% modified from ImportListSQRecordingData on 3/30/16 to handle data from
% 3/30/16 when strobe signal was unavailable. Using eye signal to correlate
% with eye signal and then importing time stamps from cortex signal. 
% written July 26, 2014 by Seth Konig updated August 19 for automation
% import and process cch25f data through plexon. Plexon eye tracking values (i.e. voltage values)
% are different from what cortex recieves though the they seem to map
% nearly linearly but not 100% sure...likely due to differences in cable
% length or other differences in inherent resistnaces throughout the system

% Inputs:
%   1) data_dir: directory containing cvtnew_file
%   2) cch25_file: plexon_nex file containing 25 point calibration trials
%       -Here is for calibration only for the purposes of the listsq task
%   3) listsq_file: plexon_nex file containing data from the ListSQ task
%   4) item_set: listsq item file name e.g. ListSQ04.itm
%   5  multiunit: 1 for mulitunit 0 for single unit
% Outputs:
%   1) saves preprocessed data to data_dir

%import cortex data
cortexfile = ['\\research.wanprc.org\research\Buffalo Lab\Cortex Data\Tobii\'...
    listsq_file(1:8) '.' listsq_file(10)];

[time_arr,event_arr,eog_arr,epp_arr,~,~]  = get_ALLdata(cortexfile);

samprate = 5;
numrpt = size(event_arr,2);
new_eog_arr=[];
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
    cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
    cnd = event_arr(cndnumind,rptlop)-1000;
    
    perbegind = find(event_arr(:,rptlop) == 100); % eye data starts at 100; might have predictive looking
    perendind = find(event_arr(:,rptlop) == 101); % eye data stops collecting after rewards so can stop here
    cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
    blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
    begtimdum = time_arr(perbegind,rptlop);
    endtimdum = time_arr(perendind(end),rptlop);
    
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
        new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
    end
end

%---get eye data for only when fixation cross or picture is displayed---%
eyedat = cell(1,length(per));
if ~isempty(epp_arr);
    pupildata = cell(1,length(per)); %horizontal pupil diameter
end
cnd=[];
teststart = [];
for trlop=1:size(per,2)
    trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
    horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
    vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
    picstart=1*samprate;
    picend=per(trlop).endsmpind-per(trlop).begsmpind;
    picend(picend > samprate*length(horeog)) = length(horeog)*samprate; 
    
    eyedat{trlop}(1,:) = horeog(round(picstart/samprate):floor(picend/samprate));
    eyedat{trlop}(2,:) = vrteog(round(picstart/samprate):floor(picend/samprate));
    cnd(trlop)=per(trlop).cnd;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Import Plexon_nex File Data---%%%
% define the trials solely based on the marker events
cfg=[];
cfg.channel = {};
cfg.datafile      = [data_dir listsq_file];
cfg.trialfun      = 'trialfunListSQ';
cfg.datatype      = 'continuous';
cfg.headerformat  = 'plexon_nex';
cfg.dataformat    = 'plexon_nex';
cfg.headerfile = cfg.datafile;
% cfg = ft_definetrial(cfg);

% read the data from file and preprocess them
% Want spike data (valid units only), LFPs, and X and Y eye data
fieldtripdefs
hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
num_units = 0;
waveform_channels = [];
for l = 1:length(hdr.label);
    if ~isempty(strfind(hdr.label{l},'sig')); %all recorded "units"
        if isempty(strfind(hdr.label{l},'wf')); %all recorded "units" waveforms, dont want for now
            if isempty(strfind(hdr.label{l}(end),'i'));%i is for invalidated waveforms?, rest are valid waveforms
                cfg.channel = [cfg.channel hdr.label(l)];
                num_units = num_units+1;
            end
        else
            if isempty(strfind(hdr.label{l}(7),'i'));%i is for invalidated waveforms, rest are valid waveforms
                waveform_channels= [waveform_channels l];
            end
        end
    elseif length(hdr.label{l}) == 1 && ~isempty(strfind(hdr.label{l},'X')); %horiziontal eye data
        cfg.channel = [cfg.channel hdr.label(l)];
    elseif length(hdr.label{l}) == 1 && ~isempty(strfind(hdr.label{l},'Y')); %vertical eye data,somehow LED3_Y made it into my recordings
        cfg.channel = [cfg.channel hdr.label(l)];
    elseif ~isempty(strfind(lower(hdr.label{l}),'pupil')) %pupil data does not exist for all recordings
        cfg.channel = [cfg.channel hdr.label(l)];
    elseif ~isempty(strfind(hdr.label{l},'AD'));  %LFP data by channel
        cfg.channel = [cfg.channel hdr.label(l)];
    end
end
cfg.trl(1).begsmpind = 1;
cfg.trl(1).endsmpind = hdr.nSamples;
cfg.channel = {'X','Y'};
data = getPlexonTrialData(cfg);

[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_set);
img_trials = find(itmlist > sequence_items(end));

%look at the first 5 minutes of eye data one of the trials should show good
%corrleations 
x = data(1).values{1}(1:60*1000*5)/2;
figure
for i = 1:length(img_trials)
    
    xx = resample(eyedat{img_trials(i)}(1,200:end),5,1); %ignore the first second due to cross
    xx = resample(eyedat{27}(1,200:end),5,1); %ignore the first second due to cross
    [acor,lag] = xcorr(x,xx);
    [~,I] = max(abs(acor));
    lagDiff = lag(I);
    
    plot(xx)
    hold on
    plot(x(lagDiff:lagDiff+length(xx)),'r');
    hold off

    reply = input('Does this look aligned? yes or no?');
    if strcmpi(reply,'y')
        break
    end
        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---calibrate eye data---%%%y
imageX = 800;
imageY = 600;

[tform] = Calibrate_Plexon_EyeData(data_dir,cch25_file,figure_dir); %get the calibration function

eyechans = find_desired_channels(cfg,'eye');
% replaces outliers (more than 1 dva outside image with NaNs not sure how Cluster Fix will process that
for t = 1:length(data(eyechans(1)).values);
    x = data(eyechans(1)).values{t}; %second to last index in structure array is horizontal eye data;
    y = data(eyechans(2)).values{t}; %last index in structure array is vertical eye data;
    
    [x,y] = tformfwd(tform,x,y); %calibrate: transform votlages into "dva", 24 pixels/dva
    
    x = 24*x; %convert from cortex dva to pixels
    y = 24*y; %convert from cortex dva to pixels
    x = x+imageX/2; 
    y = y+imageY/2;
    y(x < -24) = NaN;
    x(x < -24) = NaN;
    y(x > imageX+24) = NaN;
    x(x > imageX+24)= NaN;
    x(y < -24) = NaN;
    y(y < -24) = NaN;
    x(y > imageY+24) = NaN;
    y(y > imageY+24) = NaN;
    
    %store calibrated data
    data(eyechans(1)).values{t} = x;
    data(eyechans(2)).values{t} = y;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Run Cluster Fix---%%%
% use cluster Fix to Detect fixations and saccades in XY eye data
% since timing is very important NaNs must remain. Trials with NaNs or with eye
% data less lasting less than 200 ms (less than 1 saccade and fixation) are ignored
% and trials are processed in chunks of valid eye data. 
% reconcatendated together
disp('Running Cluster Fix')
fixationstats = cell(1,length(data(end-1).values));
for t = 1:length(data(eyechans(1)).values);
    if ~isempty(data(end-1).values{t})
        fixationstats{t} = ClusterFix_Plexon([data(eyechans(1)).values{t};data(eyechans(2)).values{t}]);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Save the Data---%%%
save([data_dir listsq_file(1:end-11) '-preprocessed.mat'],'cch25_file','cfg','data','multiunit',...
    'fixationstats','hdr','item_set','listsq_file','num_units','waveforms','tform','item_set');
disp(['Data  from ' listsq_file ' Successfully Preprocessed and Saved'])
end