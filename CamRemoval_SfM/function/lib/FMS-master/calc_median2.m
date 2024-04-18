function [ median ] = calc_median2( X )
%% CALC_MEDIAN calculate the median of the NxD matrix X
%   Input set of data N: number of points, D: dimension
%   Output geometric median of data set using iterative weiszfeld
%   algorithm, see wikipedia
%   Note: this is unregularized
%
% Author: Tyler Maunu, 2014

[N,D]=size(X);
median = zeros(1,D);
delta=10^-9;
convergence = bitmax;
iter=1;

while convergence>delta && iter<500
    
    %Weiszfeld iterations
    prev_median=median;
    
    scale=sum((X-repmat(prev_median,N,1)).^2,2).^.5;
    median=sum(X./repmat(scale,1,D),1)/sum(scale.^-1);
    
    convergence = norm(prev_median-median);
    
    iter=iter+1;
    iter;
end 
end

