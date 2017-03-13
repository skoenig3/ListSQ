function [meta_fixed,cfg_fixed] = CVTNEW_remove_bad_condition(meta,cfg)
%function removes a few trials in which condition is confused because of
%the way cortex timing file reports it. This confusion occurs at the 
%transition between blocks/paths and should be set as the condition following 
%this. I visually and manually confirmed this to be true in the code but 
%it's just easier to remove this trials completely. This is important for
%spatial analysis but not for timing analysis. The bug in the timing file 
%was fixed introduced and then fixed. This bug only applies to
%PW141007-PW141024.

%---Grab All Condition Numbers---%
all_cnds = NaN(1,length(cfg.trl));
for trl = 1:length(cfg.trl)
    all_cnds(trl) = meta(trl).cond;
end

%---Based on Condition Numbers Determine When Block/Path Changed---%
change_blocks = find(abs(diff(all_cnds)) > 0);
if rem(length(change_blocks),2) ~=0 
    error('Should come in pairs or two')
end

%even pairs are the error trials
error_trials = change_blocks(2:2:end);

%---Remove These Bad Trials From Data---%
%could also set to next condition number
for et = 1:length(error_trials)
    cfg.trl(error_trials(et)).cnd = NaN;
    cfg.trl(error_trials(et)).allval = NaN;
    cfg.trl(error_trials(et)).alltim = NaN;
    meta(error_trials(et)).cond = NaN;
    meta(error_trials(et)).x = NaN;
    meta(error_trials(et)).y = NaN;
end

%---check fixed results---%
% all_cnds = NaN(1,length(cfg.trl));
% for trl = 1:length(cfg.trl)
%     all_cnds(trl) = meta(trl).cond;
% end
% plot(all_cnds)


meta_fixed = meta;
cfg_fixed = cfg;