% ft_defaults;

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
data_file = 'TO151214_3-sorted.nex';
load([data_dir data_file(1:end-11) '-preprocessed.mat'],'cfg');
all_cfg = cfg;

spike = ft_read_spike([data_dir data_file]); 

cfg = [];
cfg.spikechannel = {'sig001b_wf','sig003a_wf','sig004a_wf','sig004b_wf'};
spike = ft_spike_select(cfg, spike); 

cfg             = [];
cfg.fsample     = 40000;
cfg.interpolate = 1; % keep the density of samples as is
[wave, spikeCleaned] = ft_spike_waveform(cfg,spike);
%%
for k = 1:length(spike.label)
  figure, 
  sl = squeeze(wave.dof(k,:,:))>1000; % only keep samples with enough spikes
  plot(wave.time(sl), squeeze(wave.avg(k,:,sl)),'k') % factor 10^6 to get microseconds
  hold on
 
  % plot the standard deviation
  plot(wave.time(sl), squeeze(wave.avg(k,:,sl))+sqrt(squeeze(wave.var(k,:,sl))),'k--') 
  plot(wave.time(sl), squeeze(wave.avg(k,:,sl))-sqrt(squeeze(wave.var(k,:,sl))),'k--')
 
  axis tight
  set(gca,'Box', 'off') 
  xlabel('time')
  ylabel('normalized voltage')
end

%%
cfg          = []; 
cfg.dataset  = [data_dir data_file];
cfg.trialfun = 'trialfun_stimon';
cfg = ft_definetrial(cfg);
cfg.timestampspersecond =  spike.hdr.FileHeader.Frequency; % 40000
spikeTrials = ft_spike_maketrials(cfg,spike); 
%%
cfg                     = [];
hdr                     = ft_read_header([data_dir data_file]);
cfg.trl                 = [0 hdr.nSamples*hdr.TimeStampPerSample 0];
cfg.timestampspersecond =  spike.hdr.FileHeader.Frequency; % 40000
spike_notrials   = ft_spike_maketrials(cfg,spike); 

dat = ft_checkdata(spikeTrials,'datatype', 'raw', 'fsample', 1000);
%%
cfg       = [];
cfg.bins  = [0:0.0005:0.1]; % use bins of 0.5 milliseconds
cfg.param = 'coeffvar'; % compute the coefficient of variation (sd/mn of isis)
isih = ft_spike_isi(cfg,spikeTrials);
%%
for k = 1:length(spike.label)
  cfg              = [];
  cfg.spikechannel = isih.label{k};
  cfg.interpolate  = 5; % interpolate at 5 times the original density
  cfg.window       = 'gausswin'; % use a gaussian window to smooth
  cfg.winlen       = 0.004; % the window by which we smooth has size 4 by 4 ms
  cfg.colormap     = jet(300); % colormap
  cfg.scatter      = 'no'; % do not plot the individual isis per spike as scatters
  figure, ft_spike_plot_isireturn(cfg,isih)
end
%%
cfg             = []; 
cfg.binsize     =  0.1; % if cfgPsth.binsize = 'scott' or 'sqrt', we estimate the optimal bin size from the data itself
cfg.outputunit  = 'rate'; % give as an output the firing rate
cfg.latency     = [-1 3]; % between -1 and 3 sec.
cfg.vartriallen = 'yes'; % variable trial lengths are accepted
cfg.keeptrials  = 'yes'; % keep the psth per trial in the output
psth = ft_spike_psth(cfg,spikeTrials);
%%
cfg         = [];
cfg.binsize =  [0.05];
cfg.latency = [-1 3];
psth = ft_spike_psth(cfg,spikeTrials);
 
cfg              = [];
cfg.topplotfunc  = 'line'; % plot as a line
cfg.spikechannel = spikeTrials.label([1 2]);
cfg.latency      = [-1 3];
cfg.errorbars    = 'std'; % plot with the standard deviation
cfg.interactive  = 'no'; % toggle off interactive mode
figure, ft_spike_plot_raster(cfg,spikeTrials, psth) 
%%
% read in the data, select the channels and define the trials
spike = ft_read_spike([data_dir data_file]); 
 
cfg              = [];
cfg.spikechannel = {'sig001b_wf'};
spike = ft_spike_select(cfg, spike);
 
cfg          = []; 
cfg.dataset  = [data_dir data_file];
cfg.trialfun = 'trialfun_stimon';
cfg = ft_definetrial(cfg);
cfg.timestampspersecond =  spike.hdr.FileHeader.Frequency; % 40000
spikeTrials = ft_spike_maketrials(cfg,spike);  

cfg             = [];
cfg.maxlag      = 0.5; % maximum 200 ms
cfg.binsize     = 0.001; % bins of 1 ms
cfg.outputunit  = 'proportion'; % make unit area
cfg.latency     = [-2.5 0];
cfg.vartriallen = 'yes'; % do not allow variable trial lengths
cfg.method      = 'xcorr'; % compute the normal cross-correlogram
Xc = ft_spike_xcorr(cfg,spikeTrials);  
 
% compute the shuffled correlogram
cfg.method      = 'shiftpredictor'; % compute the shift predictor
Xshuff = ft_spike_xcorr(cfg,spikeTrials);
%%
iCmb = 1;
jCmb = 1;
figure
xcSmoothed = conv(squeeze(Xc.xcorr(iCmb,jCmb,:)),ones(1,5)./5,'same'); % do some smoothing
hd = plot(Xc.time(3:end-2),xcSmoothed(3:end-2),'k'); % leave out borders (because of smoothing)
hold on
xcSmoothed = conv(squeeze(Xshuff.shiftpredictor(iCmb,jCmb,:)),ones(1,5)./5,'same');    
plot(Xc.time(3:end-2),xcSmoothed(3:end-2),'r')
hold on
xlabel('delay')
ylabel('proportion of coincidences')        
title([Xc.label{iCmb} Xc.label{jCmb}])
axis tight