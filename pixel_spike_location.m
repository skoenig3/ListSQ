function [firing_location] = pixel_spike_location(eyepos,spike_times,imageX,imageY)
%written by Seth Konig August, 2014
firing_location = zeros(imageY,imageX);

[trial,time] = find(spike_times == 1);
spikeind = sub2ind(size(spike_times),trial,time);

x = eyepos(1:2:end,:);
y = eyepos(2:2:end,:);
x = stretch(x,imageX);
y = stretch(y,imageY);

xs = x(spikeind);
ys = y(spikeind);
[xs,ys] = remove_nans(xs,ys);

%probably can't use logical indexing fast here since need to add spikes to pixels
for i = 1:length(xs);
    firing_location(ys(i),xs(i)) =  firing_location(ys(i),xs(i)) +1;
end

end

function vector = stretch(vector,max)
vector = vector(1:end);
vector(vector < 1) = 1;
vector(vector > max) = max;
end