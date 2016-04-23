function ImportCVTNEWRecordingDataV2(data_dir,session_data)
% Version 2.0 of ImportCVTNEWRecordingData. It was re-written slighlty to
% handle importing sesion data directly from an excel sheet. Seth Konig 1/19/2016
%
% Function imports and preprocesses CVTNEW data from plexon recordings.
% NOTE Plexon eye tracking values (i.e. voltage values) are different from what cortex
% recieves though the they seem to map likely due to differences in cable
% length or other differences in inherent resistnaces throughout the
% system. Currently this function does not apply Cluster Fix to the XY data.
%
% Inputs:
%   1) data_dir: directory containing task_file
%   2) session_data: data containing info about the session data
%
% Outputs:
%   1) saves preprocessed data to data_dir

%grab important file information from session_data
task = 'cvtnew';
 [task_file,~,~,multiunit,unit_names,unit_confidence,sorting_quality,waveform_count]...
    = get_task_data(session_data,task);
if isempty(task_file)
    warning('No CVTNEW file could be found. Exiting function...')
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Import Plexon_nex File Data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define the trials solely based on the marker events
cfg=[];
cfg.channel = {};
cfg.dataset       = [data_dir task_file];
cfg.trialfun      = 'trialfuncvtnew';
cfg.datatype      = 'continuous';
cfg = ft_definetrial(cfg);

%get position of the dot over time by trial
meta = cvtnew_dot_pos(cfg);

% read the data from file and preprocess them
% Want spike data (valid units only), LFPs, and X and Y eye data
fieldtripdefs
hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
num_units = 0;
waveform_channels = [];
for l = 1:length(hdr.label);
    if ~isempty(strfind(hdr.label{l},'sig'));
        if isempty(strfind(hdr.label{l},'wf')); %all recorded unit spike times
            if isempty(strfind(hdr.label{l}(7),'i'));%i is for invalidated waveforms rest are valid waveforms
                cfg.channel = [cfg.channel hdr.label(l)]; %all recorded "units"
                num_units = num_units+1;
            end
        else
            if isempty(strfind(hdr.label{l}(7),'i'));%i is for invalidated waveforms, rest are valid waveforms
                waveform_channels= [waveform_channels l];
            end
        end
    elseif ~isempty(strfind(hdr.label{l},'X')); %horiziontal eye data
        cfg.channel = [cfg.channel hdr.label(l)];
    elseif ~isempty(strfind(hdr.label{l},'Y')); %vertical eye data
        cfg.channel = [cfg.channel hdr.label(l)];
    elseif ~isempty(strfind(hdr.label{l},'AD'));  %LFP data by channel
        cfg.channel = [cfg.channel hdr.label(l)];
    end
end
data = getPlexonTrialData(cfg);

% import waveforms seperately for some reason couldn't get it to work above
% import waveforms seperately for some reason couldn't get it to work above
waveforms = cell(2,length(waveform_channels));
for wv = 1:length(waveform_channels);
    spikewaves = read_plexon_nex(cfg.dataset , 'channel',waveform_channels(wv));
    waveforms{1,wv} = spikewaves.dat;%waveorm shapes
    waveforms{2,wv} = spikewaves.ts;%waveform timestamps
    
    if length(spikewaves.dat) ~= waveform_count(wv)
        disp(['Warning number of waveforms imported differnet than expected: ' ...
            'Imported ' num2str(length(spikewaves.dat)) ' but expected ' num2str(waveform_count(wv))]);
    end
end

disp('Data imported Successfully')


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Save the Data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
save([data_dir task_file(1:end-11) '-preprocessed.mat'],'cfg','data','hdr',...
    'num_units','meta','task_file','num_units','waveforms',...
    'multiunit','unit_confidence','sorting_quality');
disp(['Data  from ' task_file ' Successfully Preprocessed and Saved'])