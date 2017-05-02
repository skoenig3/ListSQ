function xc = calculate_shuffled_xc(spike_times,twin)
    spike_prime = shake2(spike_times);
    spike_prime = spike_prime(1:end);
    xc = xcorr(spike_prime,spike_prime,twin);
    xc(:,twin+1) = 0;%set center value to zero for now
    xc = xc/max(xc);%normalize to 1
    xc(twin+1) = 1;%set zeroth lag to 1
    xc = xc-mean(xc); %subtract mean
end

% function [Y] = shake2(X)
% % from SHAKE(X,DIM) randomizes along the dimension DIM.
% % for Matlab R13
% % version 4.1 (may 2008)
% % (c) Jos van der Geest
% % email: jos@jasen.nl
% 
% % we are shaking the indices
% I = reshape(1:numel(X),size(X)) ;
% sz = size(I) ;
% % reshape it into a 2D matrix
% % we'll randomize along rows
% [~,ri] = sort(rand(size(I)),1) ;  % get new row indices
% ci = repmat([1:size(I,2)],size(I,1),1) ; % but keep old column indices
% I = I(sub2ind(size(I),ri,ci)) ; % retrieve values
% % re-index
% Y = X(I) ;
% end