function [selected] = select_eyepos(eyepos,select_rows)
%x eye position is in odd rows and y eye position is in even rows
%select rows should be a logical index

% Code rechecked for bugs October 17, 2016 SDK

%split x and y data,then select desired trials
x = eyepos(1:2:end,:); %horizontal eye data
y = eyepos(2:2:end,:); %vertical eye data
x = x(select_rows,:);
y = y(select_rows,:);
selected = NaN(sum(select_rows)*2,size(x,2));

%put back in original format
selected(1:2:end,:) = x;
selected(2:2:end,:) = y;
end