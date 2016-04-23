function [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,multiunit,unit_confidence,sorting_quality)
%function just grabs unit names. Session data contains sorting quality,unit
%names, multiunit seperability, etc. Some units have 0 waveforms and are
%not imported (and therefore not present in data or cfg structures) so need to match 
%the unit names in session)_data to recording data structure. Multiunit is 
% outputed seperatly since this is often used seperately. Written 1/12/16 Seth Koenig
%
% Inputs:
%   1) cfg: configuration file and trial data
%   2) hdr: header file
%   3) data: recording data 
%   4) unit_names: unit names from session_data (see listsq_read_excel)
%   5) multiunit: mulitunit score 0-5 in session data (see listsq_read_excel)
%   6) unit_confidence: confidence that unit is really a neuron 0-100% (see listsq_read_excel)
%   7) sorting_quality: how seperable unit is from noise cluster (see listsq_read_excel)
%
% Outputs:
%   1) multiunit: 1s for multiunit and 0's for single units
%   2) unit_stats: column by unit
%       a) row 1: unit name
%       b) row 2: unit confidence
%       c) row 3: sorting_quality
%   3) num_units: number of units in data file

unit_channels = find_desired_channels(cfg,'sig');
num_units = length(unit_channels);
unit_stats = cell(3,num_units);
mu = NaN(1,num_units);
for unit = 1:num_units
    unit_stats{1,unit} = hdr.label{data(unit).whichchannel};
    for n = 1:length(unit_names)
        if strcmpi(unit_stats{1,unit},strtrim(unit_names{n}))
            if  multiunit(n) > 3 %likley single unit
                mu(unit) = 0;
            else %likely multi unit
               mu(unit) = 1; 
            end
            unit_stats{2,unit} = unit_confidence(n);
            unit_stats{3,unit} = sorting_quality(n);
            break
        end
    end
end
if any(any(cellfun(@isempty,unit_stats)))
    warning('No units found!!! Unit data may not have been not imported correctly!!!')
end
multiunit = mu;