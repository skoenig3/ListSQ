function [filtered_time] = filter_time(eyepos,imageX,imageY,Fs,binsize,H)
%calculate the total time spent at any locaitons in binned pixels
spatial_time = time_per_pixel(eyepos,imageX,imageY,Fs);
filtered_time = bin2(spatial_time,binsize,binsize);
filtered_time = imfilter(filtered_time,H);
filtered_time(filtered_time == 0) = NaN; %can cause aribitrarily high firing rates
filtered_time = filtered_time(end:-1:1,:); %flip so rightside up
end

