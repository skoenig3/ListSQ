function [filtered_space] = filter_space(eyepos,spike_times,imageX,imageY,binsize,H)
%caluclate total spikes over space
[firing_location] = pixel_spike_location(eyepos,spike_times,imageX,imageY);
filtered_space = bin2(firing_location,binsize,binsize);
filtered_space = imfilter(filtered_space,H);
% filtered_space(filtered_space == 0) = NaN;
filtered_space = filtered_space(end:-1:1,:); %flip so rightside up
end