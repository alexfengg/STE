% This function computes a simple global rotation, translation and scale
% factor given (generic assumed) estimated rotations and translations and
% the ground truth rotations and translations \hat{S}_i's, also it assumes a
% ratio measurement of the form S_i^{-1}S_j, where S_i = (R_i,t_i)

function [R_Mat_GloballyCorrected, t_Mat_GloballyCorrected, MSE_rots, ...
    NRMSE_trans, R_global, R2_global, t_global, alf_global] = RotTransScaleRemove(R_Mat_est, R_Mat_orig, t_Mat_est, t_Mat_orig)

n = size(R_Mat_orig,3);

%%% Estimate the global rotation
Q = zeros(3,3);

for i = 1:n
    Q = Q + R_Mat_est(:,:,i)*R_Mat_orig(:,:,i)';
end

[Uq,Sq,Vq] = svd(Q);

R_global = Vq*Uq';
%R_global = Vq*[eye(2) zeros(2,1); zeros(1,2) det(Vq*Uq')]*Uq';

gam_est = R_global*t_Mat_est;
[R2_global t_global alf_global] = computeOrientation(t_Mat_orig,gam_est,'absolute');
t_Mat_GloballyCorrected = alf_global*R2_global*gam_est+repmat(t_global,[1 n]);

R_Mat_GloballyCorrected = zeros(size(R_Mat_est));

t_err = t_Mat_GloballyCorrected - t_Mat_orig;
ti0_mean = mean(t_Mat_orig,2);
ti0_ctrd = t_Mat_orig - ti0_mean*ones(1,n);
NRMSE_trans = sqrt(sum(bsxfun(@dot,t_err,t_err))/sum(bsxfun(@dot,ti0_ctrd,ti0_ctrd)));

MSE_rots = 0;
for i = 1:n
    R_Mat_GloballyCorrected(:,:,i) = R_global*R_Mat_est(:,:,i);
    MSE_rots = MSE_rots + norm(R_Mat_GloballyCorrected(:,:,i) - R_Mat_orig(:,:,i),'fro')^2;
end
MSE_rots = (1/n)*MSE_rots;

return