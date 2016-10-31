function [namevars] = extractvarargin(inputargs,lowercase)
%EXTRACTVARARGIN -- Extracts dyads defining variables in caller 
%   varargin - NX1 cell array of variable definitions occuring in dyads.
%       The first element of each dyad should be a variable name (string).
%       The second element is the value assigned to that variable.
%   lowercase - binary. converts all variable names to lowercase. 
%       Default: false
%
%Jon Rueckemann 2014

if isempty(inputargs)
    namevars={};
    return
end

if nargin==1
    lowercase=false;
end

assert(iscell(inputargs),'Input must be a cell array.');
assert(rem(numel(inputargs),2)==0,'Input must occur in dyads.');
assert(iscellstr(inputargs(1:2:numel(inputargs))),['First element in '...
    'each pair must be a variable name (string).'])%Variable names=strings?

if lowercase
    inputargs(1:2:numel(inputargs))=lower(inputargs(1:2:numel(inputargs)));
end

for m = 1:2:numel(inputargs)
  assignin('caller',inputargs{m},inputargs{m+1});
end

namevars=inputargs(1:2:numel(inputargs));%Return names of variables defined
end