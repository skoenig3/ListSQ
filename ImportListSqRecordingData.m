function ImportListSQRecordingData(data_dir,cch25_file,listsq_file,item_set,multiunit,figure_dir)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Import Plexon_nex File Data---%%%
% define the trials solely based on the marker events
cfg=[];
cfg.channel = {};
cfg.dataset       = [data_dir listsq_file];
cfg.trialfun      = 'trialfunListSQ';
cfg.datatype      = 'continuous';
cfg = ft_definetrial(cfg);

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
data = getPlexonTrialData(cfg);

% import waveforms seperately for some reason couldn't get it to work above
waveforms = cell(1,length(waveform_channels));
% for wv = 1:length(waveform_channels);
%     spikewaves = read_plexon_nex(cfg.dataset , 'channel',waveform_channels(wv));
%     waveforms{wv} = spikewaves.dat;
% end

disp('Data imported Successfully')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---calibrate eye data---%%%
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