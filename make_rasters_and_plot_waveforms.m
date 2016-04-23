function make_rasters_and_plot_waveforms(file,cfg,data,task,figure_dir,item_set,multiunit)
% written by Seth Konig August 2014. Updated January 2016 by SDK
% 1) Plots spike trains aligned to trial start. Below is the PSTH.
% 2) Plots waveform divided into 1/4 of the task at time
% plots will help visually indicate whether changes in firing rate are due
% to changes in spike waveforms i.e. are the recordings stable over time
% or are changes in firing rate indicative of "remapping" because waveforms
% are stable across quarters of the task

%aligns data to item 1 on for sequence task
%aligns data to image onset for image task
%aligns data to dot onset for cvtnew task

reward_code = 3;
binsize = 35;%ms probably want 25-100 ms

switch task
    case 'ListSQ'
        [itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_set);
        num_units = length(find_desired_channels(cfg,'sig'));
        for unit = 1:num_units
            type = NaN(1,length(cfg.trl));
            allspikes = zeros(length(cfg.trl),5000);
            for t = 1:length(cfg.trl);
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == 15);
                event = cfg.trl(t).alltim(cfg.trl(t).allval == 23)-trial_start;
                
                if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(1) % for the sequence trials
                    if  any(cfg.trl(t).allval == reward_code)
                        type(t) = 1;
                    else
                        continue;
                    end
                elseif  itmlist(cfg.trl(t).cnd-1000) <= sequence_items(2)
                    if  any(cfg.trl(t).allval == reward_code)
                        type(t) = 2;
                    else
                        continue;
                    end
                else
                    if  any(cfg.trl(t).allval == 23) %if image was shown
                        type(t) = 3;
                    else
                        continue;
                    end
                end
                
                spikes = find(data(unit).values{t});
                spikes = spikes-event;
                spikes(spikes < 1) = [];
                spikes(spikes > 7000) = [];
                allspikes(t,spikes) = 1;
            end
            
            allspikes(isnan(type),:) = [];
            type(isnan(type)) = [];
            
            maxtime = zeros(1,3);
            figure
            for i = 1:3;
                subplot(2,3,i)
                [trial,time] = find(allspikes(type == i,:));
                plot(time,trial,'.k')
                maxtrial = ceil(max(trial)/10)*10;
                ylim([0 maxtrial])
                ylabel('Trial #')
                maxtime(i) = ceil(max(time)/1000)*1000;
            end
            
            for i = 1:3;
                subplot(2,3,3+i)
                asp = bin1(allspikes(type == i,:),binsize,'lower','sum');
                bar(binsize:binsize:binsize*length(asp),asp,'k');
                xlim([0 maxtime(i)])
                ylabel('Count')
                if i < 3
                    xlabel('Time from Trial Start (ms)')
                else
                    xlabel('Time from Image On (ms)')
                end
                title('PSTH')
            end
            
            if multiunit(unit)
                subtitle(['Multiunit ' cfg.channel{unit}]);
            else
                subtitle(cfg.channel{unit});
            end
            save_and_close_fig(figure_dir,[file(1:end-11) '_' cfg.channel{unit} '_raster']);
        end
    case 'CVTNEW'
        num_units = length(find_desired_channels(cfg,'sig'));
        for unit = 1:num_units
            allspikes = zeros(length(cfg.trl),5000);
            for t = 1:length(cfg.trl);
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == 15);
                event = cfg.trl(t).alltim(cfg.trl(t).allval == 25)-trial_start;
                event = event(1);
                spikes = find(data(unit).values{t});
                spikes = spikes-event;
                spikes(spikes < 1) = [];
                spikes(spikes > 3500) = [];
                allspikes(t,spikes) = 1;
            end
            
            figure
            subplot(1,2,1)
            [trial,time] = find(allspikes);
            plot(time,trial,'.k')
            xlabel('Time from Dot On (ms)')
            ylabel('Trial #')
            maxtrial = ceil(max(trial)/10)*10;
            ylim([0 maxtrial])
            maxtime = ceil(max(time)/1000)*1000;
            
            subplot(1,2,2)
            asp = bin1(allspikes,binsize,'lower','sum');
            bar(binsize:binsize:binsize*length(asp),asp,'k');
            xlim([0 maxtime])
            xlabel('Time from Dot On (ms)')
            ylabel('Count')
            title('PSTH')

            
            if multiunit(unit)
                subtitle(['Multiunit ' cfg.channel{unit}]);
            else
                subtitle(cfg.channel{unit});
            end
            save_and_close_fig(figure_dir,[file(1:end-11) '_' cfg.channel{unit} '_raster'])
            
        end
end
end