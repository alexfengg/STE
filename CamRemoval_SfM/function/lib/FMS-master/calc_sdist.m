function [ dist ] = calc_sdist( S1 , S2 )
%% CALC_SDIST Summary of this function goes here
%   Calculate distance between subspaces using principal angles
%   S1,S2 column bases for subspaces
%
% Author: Tyler Maunu, 2014

    A = S1'*S2;
    [u,s,v] = svd(A);
    s = diag(s);
    for i=1:size(s,1)
        s(i,1) = acos(s(i,1));
    end
    dist = s'*s;
    dist = dist^.5;
end

