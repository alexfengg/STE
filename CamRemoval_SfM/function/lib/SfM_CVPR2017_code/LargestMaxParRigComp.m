%%*************************************************************************  
%% This function computes the largest maximally parallel rigid (in R^d)  
%% component of a graph with adjacency matrix AdjMat using the algorithm
%% in [1]. 
%%
%% [1] R. Kennedy, K. Daniilidis, O. Naroditsky, and C. J. Taylor, 
%%     Identifying maximal rigid components in bearing-based localization, 
%%     in IEEE/RSJ International Conference on Intelligent Robots and 
%%     Systems, Vilamoura, Algarve, Portugal, October 2012
%%
%% Author: Onur Ozyesil
%%*************************************************************************  
%% Input: 
%% AdjMat : n-by-n adjacency matrix of the graph
%% d      : dimension of the problem
%% 
%% Output:
%% LarMaxRigCompInds : Indices of nodes in the largest maximally 
%%                     parallel rigid component of the graph
%%*************************************************************************


function [LarMaxRigCompInds,MaxRigidComps] = LargestMaxParRigComp(AdjMat,d)

n = size(AdjMat,1);
[Ind_j, Ind_i] = find(tril(AdjMat,-1));
ss_num = length(Ind_i);

%% Generate random locations and the corresponding lines
t_rand = randn(d,n); t_rand = t_rand - sum(t_rand,2)*ones(1,n);
Q_ij_Mats = zeros(d,d,ss_num);
for k = 1:ss_num
    t_ij = normc(t_rand(:,Ind_i(k)) - t_rand(:,Ind_j(k)));
    Q_ij_Mats(:,:,k) = eye(d) - (t_ij*t_ij');
end

%% Compute R'*R, where R is the parallel rigidity matrix 
j_Vec_Lmat = [Ind_i Ind_j]'; j_Vec_Lmat = j_Vec_Lmat(:);
i_Vec_Lmat = kron([1:ss_num]',ones(2,1));
val_Vec_Lmat = kron(ones(ss_num,1),[1;-1]);
Lmat = kron(sparse(i_Vec_Lmat,j_Vec_Lmat,val_Vec_Lmat,ss_num,n,2*ss_num),speye(d));

i_Vec_RTRmat = reshape(kron(reshape(1:d*ss_num,d,ss_num),ones(1,d)),d^2*ss_num,1);
j_Vec_RTRmat = kron([1:d*ss_num]',ones(d,1));
val_Vec_RTRmat = reshape(Q_ij_Mats,d^2*ss_num,1);
RTRmat = Lmat'*sparse(i_Vec_RTRmat,j_Vec_RTRmat,val_Vec_RTRmat,d*ss_num,d*ss_num,d^2*ss_num)*Lmat;

%% Compute the maximally parallel rigid components
display('Starting computation of max-par-rig components');
Nmat = null(full(RTRmat));
totSizeComps = 0;
CostDistMats = zeros(n-1,n-1,d);
MaxRigidComps = [];
Nodes_in_max_rig_comps = [];
i = 1;
while (totSizeComps < n)&&(i<=n)
    indVec_i = [1:i-1,i+1:n];
    Nmat_i = [Nmat(1:d*(i-1),:);Nmat(d*i+1:end,:)]-kron(ones(n-1,1),Nmat(d*(i-1)+1:d*i,:));
    for kk = 1:d
        CostDistMats(:,:,kk) = abs(abs(squareform(pdist(Nmat_i(kk:d:end,:),'cosine'))-1)-1);
    end
    SameCompInds = (max(CostDistMats,[],3) < 1e-8) - eye(n-1);
    Nodes_i = sum(SameCompInds) > 0;
    indVec_i = indVec_i(Nodes_i);
    % [S_Comps, C_Comps] = graphconncomp(sparse(SameCompInds(Nodes_i,Nodes_i)));
    Graph = graph(sparse(SameCompInds(Nodes_i,Nodes_i)));
    [C_Comps,S_Comps] = conncomp(Graph);
    
    for i_k = 1:S_Comps
        MaxRigidComps{i_k} = sort([indVec_i(C_Comps == i_k) i]);
        Nodes_in_max_rig_comps = union(Nodes_in_max_rig_comps,MaxRigidComps{i_k});
    end
    totSizeComps = length(Nodes_in_max_rig_comps);
    i = i + 1;
end

CompSizes = zeros(length(MaxRigidComps),1);
if ~isempty(CompSizes)
    for kk = 1:length(MaxRigidComps)
        CompSizes(kk) = length(MaxRigidComps{kk});
    end
    [~,i_LarMaxComp] = max(CompSizes);
    LarMaxRigCompInds = MaxRigidComps{i_LarMaxComp};
else
    fprintf('\n No maximally rigid component of size > 2 ! \n');
    LarMaxRigCompInds = [];
    MaxRigidComps = [];
end

return