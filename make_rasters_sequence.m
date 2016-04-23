function make_rasters_sequence(file,cfg,data,figure_dir,item_set,multiunit)
% written by Seth Konig August 2014
% plots spike trains aligned to a key event in the task
% also plots waveform divided into 1/4 of the task at time
% plots will help visually indicate whether changes in firing rate are due
% to changes in spike waveforms i.e. are the recordings stable over time
% or are changes in firing rate indicative of "remapping" because waveforms
% are stable across quarters of the task
reward_code = 3;
num_units = length(cfg.channel)-6;
for unit = 1:num_units
    type = NaN(1,length(cfg.trl));
    allspikes = zeros(length(cfg.trl),5000);
    for t = 1:length(cfg.trl);
        if any(cfg.trl(t).allval == reward_code)
            trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == 15);
            event = cfg.trl(t).alltim(cfg.trl(t).allval == 23)-trial_start;
            
            spikes = find(data(unit).values{t});
            spikes = spikes-event;
            spikes(spikes < 1) = [];
            spikes(spikes > 7000) = [];
            allspikes(t,spikes) = 1;
        end
    end
    
    figure
    [trial,time] = find(allspikes);
    plot(time,trial,'.k')
    xlabel('Time from Stimulus Onset (ms)')
    ylabel('Trial #')
    
    if multiunit(unit)
        subtitle(['Multiunit ' cfg.channel{unit}]);
    else
        subtitle(cfg.channel{unit});
    end
    save_and_close_fig(figure_dir,[file(1:end-11) '_' cfg.channel{unit} '_raster']);
end
end