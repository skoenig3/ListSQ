function [shuffled_95_direction_STA,observed_STA,observed_direction] = CVTNEW_Direction_Analysis(spike_times,dotdirection,numshuffs,correct_trials)

twin = 1000;


spk = spike_times;
these_nans = isnan(spike_times);%get nanind
spk(correct_trials == 0,:) = 0; %remove spikes from incorrect trials
spk(these_nans) = NaN; %restor nanind
spk = spk';
spk = spk(:);
spk(isnan(spk)) = [];%remove nans
spikes = find(spk == 1);


dd = dotdirection';
dd = dd(:);
dd(isnan(dd)) = [];%remove nans


directions = NaN(length(spikes),twin*2+1);
for s = 1:length(spikes);
    tind = spikes(s);
    if tind <= twin
        continue
    elseif tind >= length(dd)-twin
        continue
    end
    directions(s,:) = dd(tind-twin:tind+twin);
end


observed_direction = NaN(1,size(directions,2));
observed_STA = NaN(1,size(directions,2));
for t = 1:size(directions,2)
    vals =  directions(:,t);
    vals(isnan(vals)) = [];
    if ~isempty(vals)
        observed_direction(t) = circ_mean(vals);
        observed_STA(t) = circ_r(vals);
    end
end

%%
shuffled_95_direction_STA = NaN(numshuffs,twin*2+1);
parfor shuff = 1:numshuffs
    shuffled_spike_times = circshift_cvtnew(spike_times);
    these_nans = isnan(shuffled_spike_times);
    shuffled_spike_times(correct_trials == 0,:) = 0;
    shuffled_spike_times(these_nans) = NaN;
    shuffled_spike_times = shuffled_spike_times';
    shuffled_spike_times = shuffled_spike_times(:);
    shuffled_spike_times(isnan(shuffled_spike_times)) = [];
    
    spikes = find(shuffled_spike_times == 1);
    directions = NaN(length(spikes),twin*2+1);
    for s = 1:length(spikes);
        tind = spikes(s);
        if tind <= twin
            continue
        elseif tind >= length(dd)-twin
            continue
        end
        directions(s,:) = dd(tind-twin:tind+twin);
    end
   
    shuffled_circ_vars = NaN(1,2*twin+1);
    for t = 1:size(directions,2)
        vals =  directions(:,t);
        vals(isnan(vals)) = [];
        if ~isempty(vals)
            shuffled_circ_vars(t) = circ_r(vals);
        end
    end
     shuffled_95_direction_STA(shuff,:) = shuffled_circ_vars;
end
shuffled_95_direction_STA = prctile(shuffled_95_direction_STA,95);

end