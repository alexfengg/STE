% This file computes "robust" transform between point clouds using RANSAC

function [R_opt, t_opt, c_opt, t_fit, MeanError, MedianError, SqErrVec, InlierInds] = RobustRANSACFit(t, t_0, iterNums, threshR)

n = size(t,2);
HypoSize = 3;
maxInliers = -1;

for k = 1:iterNums
    Ind_i = randsample(n,HypoSize);
    [R_i, t_i, c_i] = computeOrientation(t_0(:,Ind_i),t(:,Ind_i),'absolute');
    DiffMat = t_0 - (c_i*R_i*t + t_i*ones(1,n));
    SqErrVec = bsxfun(@dot,DiffMat,DiffMat);
    InlierInds = SqErrVec < threshR^2;
    numInliers = sum(InlierInds);
    if (numInliers > maxInliers)
        R_opt = R_i;
        t_opt = t_i;
        c_opt = c_i;
        maxInliers = numInliers;
    end
    if (mod(k,1000) == 0)
        fprintf('RANSAC iteration %d done\n', k);
    end
end

DiffMat = t_0 - (c_opt*R_opt*t + t_opt*ones(1,n));
SqErrVec = bsxfun(@dot,DiffMat,DiffMat);
InlierInds = SqErrVec < threshR^2;
[R_opt, t_opt, c_opt] = computeOrientation(t_0(:,InlierInds),t(:,InlierInds),'absolute'); 

t_fit = c_opt*R_opt*t + t_opt*ones(1,n);
DiffMat = t_0 - t_fit;
SqErrVec = bsxfun(@dot,DiffMat,DiffMat);
MeanError = mean(sqrt(SqErrVec));
MedianError = median(sqrt(SqErrVec));

return