function [ Lf,Lc ] = fms( X , d , options  )
%% FMS: Robust subspace recovery algorithm
%   X: NxD data set with N points dim D
%   d: dim of the supspace to find
%   options: struct containing any of the following
%       epsilon: the regularization paramter, epsilon=p*\delta
%       maxiter: the maximum number of desired iterations
%       svdopt: 'randomized' for fast randomized svd, 'normal' for full svd
%       scaleopt: 'normal' for scaling by dist^{2-p}, 'log' for log scaling (p=1 by
%           default)
%       initopt: 'random' for completely random initialization, 'pca' for pca
%           intialization
%       lambda: regularization term for log scaling
%   
%   Lf: column basis for the final estimate of the subspace
%
% Reference: “Fast, Robust and Non-convex Subspace Recovery” 
%	https://arxiv.org/abs/1406.6145
%
% Author: Tyler Maunu, 2014


%% Set settings, default setting included
if isfield(options,'scaleopt')
    scaleopt=options.scaleopt;
else
    scaleopt='normal';
end
if isfield(options,'p')
    p=options.p;
else
    p=1;
end
if isfield(options,'initopt')
    initopt=options.initopt;
else
    initopt='random'; 
end
if isfield(options,'svdopt')
    svdopt=options.svdopt;
else
    svdopt='normal';
end
if isfield(options,'maxiter')
    maxiter=options.maxiter;
else
    maxiter=100;
end
if isfield(options,'epsilon')
    epsilon=options.epsilon;
else
    epsilon=10^-10;
end
if isfield(options,'lambda')
    lambda=options.lambda;
else
    lambda=0.0001;
end
% addpath('random pca');

[N,D] = size(X);
Lf=zeros(D,d);


%Spherize the data
X = X./(repmat(sum(X.*X,2).^.5,1,D));

%% Set convergence criteria, initialize iterations
c = 1e7;
iter = 0;

%Initialize subspace
if strcmp(initopt,'random')
    %Completely random
    Vi=rand(D,d);
    Vi=orth(Vi);
elseif strcmp(initopt,'pca')
    %svd initialization
    [Ui,Si,Vi] = svd(X);
    Vi = Vi(:, 1:d);
else
    error('Bad initopt tag')
end

Vi_prev = Vi;


%% Iterative procedure
while c>10^-9 && iter<maxiter
    
    %Project datapoints onto the orthogonal complement
    C = X' - Vi*(Vi'*X');
    
    %Create the scaled Y matrix
    if strcmp(scaleopt,'log')
    % weights = dist/logdist, lambda is a regularization parameter
        scale = ((sum(C.*C,1)).^.5);
        logscale=log(scale+ones(size(scale,1),size(scale,2))*lambda);
        scale=scale.*(logscale.^(-1/2));
        Y=X.*repmat(min(scale'.^-1,1/epsilon),1,D);
    elseif strcmp(scaleopt,'normal')
    %Powers of distance: weights = 1/dist^((2-p)/2)
        scale = ((sum(C.*C,1)).^.5).^(2-p);
        Y=X.*repmat(min(scale'.^-0.5,1/epsilon),1,D);
    else
        error('Bad scaleopt tag')
    end
    
    %calculate new subspace via svd
    if strcmp(svdopt,'normal')
        [Ui,Si,Vi] = svd(Y,'econ');Vi = Vi(:,1:d);
    elseif strcmp(svdopt,'randomized')
        [Ui,Si,Vi] = randpca(Y,d);
    else
        error('Bad svdopt tag')
    end

    %update convergence criteria
    c = calc_sdist(Vi(:,1:d),Vi_prev(:,1:d));
    
    %update subspace and iteration count
    Vi_prev = Vi(:,1:d);
    iter = iter+1;
    
end

%% Final subspace

[Ui,Si,Vi] = svd(Y,'econ');
Lf = Vi(:,1:d);
Lc = Vi(:,d+1:end);
end
