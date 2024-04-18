%%*************************************************************************  
%% This function estimates the directed line \gamma_ij between the camera 
%% locations ti and tj. Each direction is computed by solving the 
%% non-convex program (12) in our paper using a "heuristic" IRLS approach.
%%
%% Author: Onur Ozyesil
%%*************************************************************************  
%% Input:
%% R_Mats   : 3-by-3-by-n matrix of estimated camera orientations
%% CorrPnts : n-by-n cell of corresponding points between images (ij'th 
%%            entry of CorrPnts is a 2-by-mij-by-2 matrix qij, where 
%%            qij(:,:,1) is the matrix of points in the i'th image and 
%%            qij(:,:,2) is the matrix of points in the j'th image)
%% invK_Mats: 3-by-3-by-n matrices of inverses of camera calibration 
%%            matrices
%% 
%% Output:
%% t_ij_Mat : 3-by-(# of corresponding image pairs) matrix of 
%%            direction vectors


function [tijMat,CostMat,TijCell] = RobustSignedLineEstimInvK(R_Mats,CorrPnts,invK_Mats)

n = size(R_Mats,3);
AdjMat = zeros(n,n);
for i = 1:n
    for j = i+1:n
        AdjMat(i,j) = size(CorrPnts{i,j},2);
    end
end
AdjMat = AdjMat + AdjMat';
[Ind_j, Ind_i] = find(tril(AdjMat,-1));
ss_num = length(Ind_i);
tijMat = zeros(3,ss_num);
CostMat=nan(size(CorrPnts));
TijCell=cell(size(CorrPnts));

i_k_print = 1;
for k = 1:ss_num
    i_k = Ind_i(k); j_k = Ind_j(k);
    if (i_k ~= i_k_print)
        %fprintf('Completed subspace estimations for camera i = %d \n ----------- \n',i_k_print);
    end
    mij = AdjMat(i_k,j_k);
    
    % Compute noisy subspace samples
    R_pi = R_Mats(:,:,i_k)*invK_Mats(:,:,i_k)*[CorrPnts{i_k,j_k}(:,:,1);ones(1,mij)];
    R_pj = R_Mats(:,:,j_k)*invK_Mats(:,:,j_k)*[CorrPnts{i_k,j_k}(:,:,2);ones(1,mij)];
    nuijMat = normc(cross(R_pi,R_pj,1));
    
    [tijInit,costval] = Reaper2DIRLS(double(nuijMat));
    xi_normed = normc(R_pi);
    xj_normed = normc(R_pj);
    
    % determine the signs considering the inliers!!!
    
    votePlus = 0;
    voteMinus = 0;
    voteWeights = abs(tijInit'*nuijMat); voteWeights = voteWeights./sum(voteWeights);
    for kk = 1:mij
        xi_k = -xi_normed(:,kk);
        xj_k = xj_normed(:,kk);
        Xij = [xj_k xi_k];
        [Uij,~,~] = svd((Xij*Xij'));
        Pij = Uij(:,1:2)*Uij(:,1:2)';
        tijProj = normc(Pij*tijInit);
        ang_tx_i = acos(tijProj'*xi_k);
        ang_tx_j = acos(tijProj'*xj_k);
        if ((ang_tx_i + ang_tx_j) <= pi)
            votePlus = votePlus + 1/(min(voteWeights(kk),100));
        else
            voteMinus = voteMinus + 1/(min(voteWeights(kk),100));
        end 
    end
    tijMat(:,k) = ((votePlus >= voteMinus) - (votePlus < voteMinus))*tijInit;
    CostMat(i_k,j_k)=costval;
    TijCell{i_k,j_k}=tijMat(:,k);
    
    i_k_print = i_k;
    
    Eij=crossProdMat(R_Mats(:,:,i_k)'*tijMat(:,k))*R_Mats(:,:,i_k)'*R_Mats(:,:,j_k);
    P1=[CorrPnts{i_k,j_k}(:,:,1);ones(1,mij)];
    P2=[CorrPnts{i_k,j_k}(:,:,2);ones(1,mij)];
    e=sum(P1.*(Eij*P2));
    d=abs(e);
end

return