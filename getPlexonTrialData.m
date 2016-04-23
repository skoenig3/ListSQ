function [dataout,skipped_trials] = getPlexonTrialData(cfg)
% written 7/26/14. Simplified and streamlined version of ft_preprocessing
% sufficient to grab data from invidiaul trials

% Known bug. Found that several recent nex files cch25f.sav have not had
% the end of the session recorded/coded properly? I'm not sure why this keeps
% happending. These trials will just have to be cutoff and ignored. SDK 1/6/2016


fieldtripdefs

cfg.continuous = 'yes'; % get data as if it were continous

% read the header
hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);

ntrl = size(cfg.trl,2);
if ntrl<1
    error('no trials were selected for preprocessing, see FT_DEFINETRIAL for help');
end

% determine the channel numbers of interest for preprocessing
[~,rawindexes] = match_str(cfg.channel, hdr.label);

dataout.values = []; %raw data for each trial
dataout.hdr = []; %header info
dataout.allchannels = []; %which channels worth of data
dataout.fsample= [];%sampling rate

for j=1:length(rawindexes)
    % read one channel group at a time, this speeds up combined datasets
    % a multiplexed dataformat is faster if you read all channels, one trial at a time
    rawindx = rawindexes(j);
    trial_dat = cell(1,ntrl);
    fprintf('processing channel { %s}\n', sprintf('''%s'' ', hdr.label{rawindexes(j)}));

    skipped_trials = []; 
    for i=1:ntrl
        nsamples = cfg.trl(i).endsmpind-cfg.trl(i).begsmpind;
        begsample = cfg.trl(i).begsmpind; %trail start sample
        endsample = cfg.trl(i).endsmpind; %trial end sample
        if begsample<1
            warning('cannot apply enough padding at begin of file');
            begsample  = 1;
        end
        if endsample > hdr.nSamples %for whaterver reason there isn't data always at the end of the file so have to skip
            skipped_trials = [skipped_trials i];
            disp(['Skipped trial #s ' num2str(i) '-' num2str(ntrl) ' Not enough data???'])
            break
        end
        %read in the data 1 trial at a time
        dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', rawindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
        trial_dat{i} = dat;
    end % for all trials
    trial_dat(skipped_trials:end)=[];
    
    dataout(j).values = trial_dat;
    dataout(j).hdr = hdr; % header details of the datafile
    dataout(j).whichchannel = rawindx; %which channel. Channel is is label
    dataout(j).fsample= hdr.Fs;%sampling rate
    
end % for all channel groups
