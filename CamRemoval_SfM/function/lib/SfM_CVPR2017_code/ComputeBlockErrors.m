%%*************************************************************************  
%% This function computes the consistency errors of pairwise rotation 
%% measurements (estimates of Rij = Ri^(-1)*Rj) with respect to a given set
%% of estimates of rotations Ri 
%%
%% Author: Onur Ozyesil
%%*************************************************************************  
%% Input: 
%% H          : 3n-by-3n matrix of estimates of ratios Rij, ij'th 3-by-3 
%%              block Hij of H is equal to the measurement of Rij (if there
%%              is no measurement for ij, Hij = 0_{3x3})
%% SO3Mats_est : 3-by-3-by-n matrix of rotation estimates 
%% AdjMat      : n-by-n adjacency matrix of the measurement graph
%%
%% Output:
%% BlockErrMat : n-by-n matrix of consistency errors 
%%*************************************************************************

function BlockErrMat = ComputeBlockErrors(H,SO3Mats_est,AdjMat)

n = size(SO3Mats_est,3);
BlockErrMat = zeros(n,n);
for i = 1:n
    for j = i+1:n
        if (AdjMat(i,j) == 1)
            BlockErrMat(i,j) = norm(SO3Mats_est(:,:,i)'*SO3Mats_est(:,:,j) - H(3*(i-1)+1:3*i,3*(j-1)+1:3*j),'fro');
        end
    end
end
BlockErrMat = BlockErrMat + BlockErrMat';

return