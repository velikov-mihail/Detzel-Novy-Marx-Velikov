function [ind] = ptfindex_from_inds(x1,x2)
% PURPOSE: assigns indexes on double sorted portfolios
%---------------------------------------------------
% USAGE: [ind]=ptfindex_from_inds(ind1,ind2)      
%        Function that spits out a matrix of which combined bin each firm falls
%        under as well as a matrix that shows how the individual sorts map into
%        the double one. Note that B/M is available only in June of each year!
%---------------------------------------------------
% Inputs
%        -x1 - index of assigned firms to bins 
%        -x2 - index of assigned firms to bins 
% Output
%        -ind - matrix with assigned firms to n1*n2 bins

n1 = max(max(x1));                           % Number of bins in first sort
n2 = max(max(x2));                           % Number of bins in second sort
ind  = zeros(size(x1));                      % Output matrix

for i = 1:n1   
    for j = 1:n2                             % For each common bin
        ind(x1 == i & x2 == j) = (i-1)*n2 + j;   % Give it a number, in the end, all will be from 1 to n1*n2
    end
end