%%*************************************************************************  
%% This function implements the heuristic IRLS method used to solve problem
%% (12) in [1]. In essence, it takes as input noisy 2D subspace samples 
%% nu_ij^k's, and at each iteration solves for the PCA problem (where the 
%% contribution of each subspace sample is weighted according to how well 
%% it was approximated in the subspace of the previous iteration) to find 
%% the "robust" 2D subspace of the samples nu_ij^k's (or equivalently, the 
%% line, represented by the output gamVec, orthogonal to this 2D subspace). 
%%
%% [1] O. Ozyesil, A. Singer, ``Robust Camera Location Estimation by Convex
%% Programming'', CVPR 2015.
%%
%% Author: Onur Ozyesil --- Written using MATLAB R2012a (7.14.0.739)
%%*************************************************************************  
%% Input:
%% nuMat    : 3-by-m_ij matrix of noisy 2D subspace samples, computed by 
%%            using (11) in [1]. 
%% Output:
%% gamVec   : 3-by-1 vector, representing the (undirected) line orthogonal 
%%            to the computed "robust" 2D subspace (which is used as the 
%%            line between cameras i and j to estimate the direction vector 
%%            between them).
%%*************************************************************************

function [gamVec,costVal] = Reaper2DIRLS(nuMat)

%% Algortihmic parameters for IRLS:
epslon = 1e-5;                % Level of accuracy
kmax = 100;                   % Max number of iterations
max_in_iter = 10;             % Max number of stagnated iterations
p_w = 1;                      % Penalty exponent, set to 1 to solve (12) in [1]
w_k = ones(1,size(nuMat,2));  % Initial PCA weights

% Trivial parameters
in_iter = 0;
k = 0;
costVal = 0;
gamVec = normc(randn(3,1));

%% IRLS iterations:

while (in_iter < max_in_iter)&&(k < kmax)
    Amat = nuMat*diag(1./max(w_k.^p_w,10^(-7)))*nuMat'; % Matrix of PCA
    [Uk,~,~] = svd(Amat); % Solve PCA
    gamVecNew = Uk(:,3); % Update optimal point
    w_k = abs(gamVecNew'*nuMat); % Update PCA weights
    costNew = sum(w_k); % Update cost value
    
    % Stagnation check:
    delt = max(abs(costVal-costNew),1-gamVecNew'*gamVec);
    if (delt <= epslon)
        in_iter = in_iter + 1;
    else
        in_iter = 0;
    end
    
    % Trivial updates:
    costVal = costNew; gamVec = gamVecNew;
    k = k+1;
end

return