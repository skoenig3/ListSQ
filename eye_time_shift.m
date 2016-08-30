function eyeposition = eye_time_shift(eyepos,shift)
%written by Seth Konig August 23, 2016
%code shifts eye position relative to the spike trains so that can see if
%place fields differ by various time points in the past and future


if shift == 0
    eyeposition = eyepos;
elseif shift > 0
    len = size(eyepos,2);
    eyeposition = [NaN(size(eyepos,1),shift-1) eyepos(:,1:(len-shift))];
elseif shift < 0
    shift = -shift;
    eyeposition = [eyepos(:,shift+1:end) NaN(size(eyepos,1),shift)];
end
end