function [spatial_time] = time_per_pixel(eyepos,imageX,imageY,Fs)
% written by Seth Konig August, 2014
spatial_time = zeros(imageY,imageX);
x = eyepos(1:2:end,:);
y = eyepos(2:2:end,:);
x = stretch(x,imageX);
y = stretch(y,imageY);

%probably can't use logical indexing fast here since need to add spikes to pixels
xy_ind = sub2ind(size(spatial_time),y,x);
for i = 1:length(xy_ind);
    spatial_time(xy_ind(i)) =  spatial_time(xy_ind(i)) +1;
end
spatial_time = spatial_time/Fs; %convert from ms to sec
end

function vector = stretch(vector,max)
    vector = vector(1:end);
    vector(isnan(vector)) = [];
    vector(vector < 1) = 1;
    vector(vector > max) = max;
end