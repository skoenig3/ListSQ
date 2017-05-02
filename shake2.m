function [Y] = shake2(X)
% modified from SHAKE(X,DIM) randomizes along the dimension DIM.
% for Matlab R13
% version 4.1 (may 2008)
% (c) Jos van der Geest
% email: jos@jasen.nl

% we are shaking the indices
I = reshape(1:numel(X),size(X)) ;
sz = size(I) ;
% reshape it into a 2D matrix
% we'll randomize along rows
[~,ri] = sort(rand(size(I)),1) ;  % get new row indices
ci = repmat([1:size(I,2)],size(I,1),1) ; % but keep old column indices
I = I(sub2ind(size(I),ri,ci)) ; % retrieve values
% re-index
Y = X(I) ;
end