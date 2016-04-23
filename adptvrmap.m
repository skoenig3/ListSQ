function [adptvRateMap shortind] = adptvrmap(spikemap, timemap)
% by MM July 2014.
% Compute rate map using adaptive smoothing method.

%% ADAPTIVE SMOOTHING METHOD
% This method of smoothing was originally described by Skaggs et al. (1996)
% and also used by Zhang...Moser (2013) for the purpose of smoothing the rate map
% before computing Skaggs' spatial information measure (bits per second or spike).

% Interestingly, Zhang et al. used this method of producing a rate map only
% before computing the Skaggs' Info Measure, and did not use this rate map for
% the grid cell computation (used a gaussian for that rate map). The
% gaussian- and adaptive-smoothed rate maps do look similar, except that the
% adaptive smoothing causes a more noisy, pixelated map than the gaussian,
% which brings some associated value changes for the firing rate of spatial
% bins. Also this adaptive smoothing method does not exclude spatial bins
% without occupancy, although those spatial bins could of course be removed post-hoc from the
% rate map, and also are not "overlooked" in in the smoothing method, which
% takes spatial bin occupancy into account.

% A good quality of using this adaptive smoothing for the rate map is that
% the smoothing is "blind" for each spatial bin, whereas with a gaussian filter,
% the smoothing for the entire rate map is *determined arbitrarily*
% (like sigma of 2 or 3) before computing the rate map.

% This smoothing is "used to optimize the trade-off between blurring error
% and sampling error" (Zhang 2014 and Skaggs 1996). Firing rate at each point
% in the environment was estimated by expanding a circle around the point until

% Zhang:
% radius of the circle in bins >= constant / (# of occupancy samples in the circle * square root of the total # spikes in those occupancy samples)

% Skaggs (and me):
% # of spikes within the circle >  constant / [ (# of occupancy samples in the circle)^2 * (radius of circle in bins)^2 ]


% Then:

% Firing rate at that point = position sampling rate * (total # spikes in the circle / # of occupancy samples in the circle).



%%
constant = 1000000 ; % a scaling parameter (10,000 in Zhang and 1,000,000 in Skaggs 1996.
rmvEyelessBins = 0 ; % If "1", then spatial bins where monkey did not look get a NaN value in the rate map.
sampRate = 1000   ; % I-SCAN eye tracker sampling rate in Hz.
maxBinRadius =  min(abs(floor(size(spikemap,2)/2)), abs(floor(size(spikemap,1)/2))) ;

% For adaptive smoothing: bins where the hypothetical "circle" expands
% outside rate map edges will just have "no occupancy" for that part of the
% circle, as this seems rationally accounted for by the inclusion of number
% of occupancy bins in the circle in the smoothing procedure.

% expand matrix for adaptive smoothing circles that go beyond rate map
% edges.
topAdd = zeros(maxBinRadius,2*maxBinRadius + size(spikemap,2)) ;
botAdd = topAdd;
sideAdd = zeros(size(spikemap,1),maxBinRadius);
spikemap2  = spikemap ;
timemap2   = timemap ;
spikemap2 = [sideAdd,  spikemap2, sideAdd] ;
spikemap2 = [topAdd;  spikemap2; botAdd];
timemap2 = [sideAdd,  timemap2, sideAdd] ;
timemap2 = [topAdd;  timemap2; botAdd];
binCentered_rows =  [-maxBinRadius:maxBinRadius] ;
binCentered_cols =  binCentered_rows ; % bin centered map is symmetric, so these values are the same.
binCenteredMap_size = [size(binCentered_rows,2), size(binCentered_cols,2)] ;
allradii = zeros(binCenteredMap_size);
for yind = 1:binCenteredMap_size(1)
    for xind = 1:binCenteredMap_size(2)
        allradii(yind,xind) = norm([binCentered_cols(xind) binCentered_rows(yind)]);
    end
end
% The "allradii" matrix can be applied to *any* bin for pinpointing the
% "growing circle" bins around the bin of interest (in process of adaptive smoothing).
for radius = 1:maxBinRadius
    circind{radius,1} = (  allradii <= radius )  ;
end

%%  Now iterate through *only the rate map bins* within the expanded matrix to
% determine each bin's smoothed firing rate.

FR = NaN(size(spikemap));
shortind =  NaN(size(spikemap));
bigenoughRadius = NaN(size(spikemap)); % just to see what were the FR circle radii (in bin units) for FR computation at each bin.
spi =  NaN(size(spikemap)); % just to see what individual bin spi_scores were.
totalTime = sum(sum(timemap)) ; % total ms looking at images (outside of first 500 ms of every image presentation to avoid onset transients)

for bin = 1:size(spikemap,2)*size(spikemap,1) % for every spatial bin in the rate map.
    
    temp = zeros(size(spikemap)); temp(bin)=1;
    
    [binRow binCol] = find(temp); % real coordinates of bin in rate map.
    
    % If the monkey spent no time looking at this bin and we want to remove those
    %  spatial bins from the rate map (as it will afffect grid score and border score)
    if timemap(binRow, binCol) == 0 ; shortind(binRow, binCol) = 1; else shortind(binRow,binCol) = 0; end
    
    if timemap(binRow, binCol) == 0 && rmvEyelessBins == 1;
      
        FR(binRow, binCol) = NaN;
        
    else % Get a FR for the bin:
        
        % convert the bin of interest coordinates to expanded-map coordinates:
        binRow2 = binRow + maxBinRadius; binCol2 = binCol + maxBinRadius;
        binCenteredMap_spikes = spikemap2(binRow2 - maxBinRadius: binRow2 + maxBinRadius, binCol2 - maxBinRadius: binCol2 + maxBinRadius) ;
        binCenteredMap_time   = timemap2(binRow2 - maxBinRadius: binRow2 + maxBinRadius, binCol2 - maxBinRadius: binCol2 + maxBinRadius) ;
        
        bigenoughCirc = zeros(maxBinRadius,1);
        
        for radius = 1:maxBinRadius; % for every radius from 1 bin to the maximum circle radius possible in this rate map.
            
            circ_spikes(radius,1) = sum(sum(binCenteredMap_spikes(circind{radius,1}))) ;
            circ_times(radius,1) = sum(sum(binCenteredMap_time(circind{radius,1}))) ;
            
            %     # of spikes within the circle >  constant / [ (# of occupancy samples in the circle)^2 * (radius of circle in bins)^2 ]
            
            num_occSamples(radius,1) = circ_times(radius,1) .* (sampRate / 1000) ;
            
            if circ_spikes(radius,1) > constant ./  ( (num_occSamples(radius,1)^2).*(radius^2) ) ;
                bigenoughCirc(radius,1) = 1;
            else
                bigenoughCirc(radius,1) = 0;
            end
            
        end % end loop of possible bin radii of circle.
        if sum(bigenoughCirc) == 0; % if there are no circles around the bin that meet criterion.
            FR(binRow, binCol) = 0; % (A cell with time but no spikes around gets an FR of 0-- as circles expands until criterion spikes are collected within).
        else
            bigenoughI =  find(bigenoughCirc);
            bigenoughRadius(binRow, binCol) = bigenoughI(1) ; % First radius big enough where equation is true.
            
            % Firing rate at that point = position sampling rate * (total # spikes in the circle / # of occupancy samples in the circle).
            FR(binRow, binCol) = sampRate .* ( circ_spikes(bigenoughRadius(binRow, binCol)) ./ num_occSamples(bigenoughRadius(binRow, binCol)) ) ;
            % "sampRate" here can be left in Hz, bc we want firing rate in Hz, and samples cancels out leaving spikes divieded by seconds for firing rate.
        end
    end
end % end loop of spatial bins


% % To evaluate results of above adaptive smoothing:
%  Also see variables: "FR" and see "bigenoughRadius" to know how many bins
%  were in the radius of the circle used for computation of that bin's FR.
%  Should see more bins for areas not looked at much since there needs to
%  be at least 1 spike for that equation determining circle radius to become true.

%%  TO PLOT RATE MAP PRODUCED:
% minFR = min(FR(:)); medFR = median(FR(:)); diff = medFR - minFR;
% figure(100); imsc(FR,[minFR  medFR+diff],'jet'); axis xy ; axis image % Need to flip axis to compare to rate maps plotted

adptvRateMap = FR;