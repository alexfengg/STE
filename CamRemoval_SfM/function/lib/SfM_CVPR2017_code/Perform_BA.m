function [res]=Perform_BA(result,data_mod,norm_lim,max_th,av_th,id,rid)
addpath('../../sba-1.6/matlab/');
ext=['_' num2str(rid)];
%% Create neccessary stuff for BA
[X_init,R,t,K,W_x,W_y,W_mask,cams,pts3D]=create_3D_initial(data_mod,result,norm_lim,max_th,av_th,ext);
%[X_init,R,t,K,W_x,W_y,W_mask,cams,pts3D]=create_3D_iter(data_mod,result,norm_lim,max_th,av_th);

%% convert to rotation and t
Ta=cams(:,5:7)'; Ta=reshape(Ta,[3,1,size(Ta,2)]);
for i=1:size(cams,1)
    Ra(:,:,i)=quat2rotm(cams(i,1:4));
end

%%

SO3Mats_estimates = zeros(size(Ra));
t_est = zeros(size(Ta));
nn = size(SO3Mats_estimates,3);
for i = 1:nn
    SO3Mats_estimates(:,:,i) = Ra(:,:,i)';
    t_est(:,:,i) = -Ra(:,:,i)'*Ta(:,:,i);
end
% Res.t = t_est; Res.R = Results_0.R; Res.K = Results_0.K;
t_est = squeeze(t_est);

%% calculate reprojection error after BA
[W_err_PBA, mse_err_PBA, mean_err_PBA, num_reproj_PBA] = ComputeReprojectionErr_Fast(pts3D', Ra, Ta, K, W_x, W_y, W_mask, 1);

S_pba_origin = median(pts3D',2);
S_pba_ctrd = bsxfun(@minus, pts3D', S_pba_origin);


%% Compare PBA output translations and rotations to Snavely's:

ComparisonInds = data_mod.ComparisonInds;
Jmat = diag([1 1 -1]);
SO3Mats_orig = zeros(3,3,data_mod.n);
t_orig2 = zeros(3,data_mod.n);

for k = 1:data_mod.n
    SO3Mats_orig(:,:,k) = Jmat*data_mod.R(:,:,k)'*Jmat;
    t_orig2(:,k) = -SO3Mats_orig(:,:,k)*Jmat*data_mod.t(:,k);
end
t_orig = t_orig2(:,ComparisonInds) - mean(t_orig2(:,ComparisonInds),2)*ones(1,size(t_orig2(:,ComparisonInds),2));


[SO3Mats_corr2, MSE_rots2, R_global2] = GlobalSOdCorrectLeft(SO3Mats_estimates, SO3Mats_orig);
[E_res,e,normE]=CompareRotations(SO3Mats_corr2,SO3Mats_orig);
res.Rot_err_mean=E_res(1); res.Rot_err_median=E_res(2);

[t_fit2,t_opt2,c_opt2,NRMSE2,~] = SimpleTransScaleRemove(R_global2*t_est(:,ComparisonInds), t_orig,'L1');

res.Tran_err_mean=mean(sqrt(sum((t_fit2 - t_orig).^2)));
res.Tran_err_median=median(sqrt(sum((t_fit2 - t_orig).^2)));