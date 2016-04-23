function meta = cvtnew_dot_pos(cfg)
% written by Niklas Wilming, adapted by Seth Konig on August 12, 2014
%% Returns meta information for each trial defined in trial_times
%
% Reads the events and parses the data structure to extract where the
% moving dot was at which time point.

task_code = 666;
meta = struct();
k=0;
for trial = 1:length(cfg.trl);
    
    e = cfg.trl(trial).allval;%event codes
    et = cfg.trl(trial).alltim;%event times
    
    tflip = et(e==task_code);
    
    condnum = cfg.trl(trial).cnd-1000;
    condnum = condnum(1);%in case there are multiple ones for whatever reason??
    
    %[trial, condnum, trl_start,trial_times(trial, 1), trl_end,trial_times(trial, 2), length(e)]
    
    A = dlmread(fullfile('C:\Users\seth.koenig\Documents\MATLAB\ListSQ\cvtnew', [num2str(condnum) 'hm.trj']));
    X = A(2:2:end)*24 + 400;
    Y = A(3:2:end)*24 + 300;

    % find x,y pos
    idmod = 2000<=e & e<3000;
    iddiv = e>=3000 & e<4000;
    cmod = e(idmod)-2000;
    cdiv = (e(iddiv)-3000)*1000;
    counter = cdiv+cmod+1;
    if max(counter) > length(X)
        % In some rare cases trials did not end at the end of the
        % trajectory. This happened when the counter was advanced and a
        % break error occured before the counter position was checked.
        % Stupid programming error. I'm throwing away the trials where this
        % happened. Keeping the trajectory though.
        x = nan*counter;
        y = nan*counter;
    else
        x = X(counter);
        y = Y(counter);
    end
    assert(length(tflip)==length(x));
    meta(trial).x = x;
    meta(trial).y = y;
    meta(trial).sample = tflip;
    meta(trial).counter = counter;
    meta(trial).cond = condnum;
end

end
