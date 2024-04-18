function [cov,res]=STE(X,d,gam)
% STE with fixed gamma

[D,N]=size(X);
initcov=eye(D);
eps=10^-10;

oldcov=initcov-1;
cov=initcov;
iter=0;
res = {};
while norm((oldcov-cov),'fro')>10^-10 & iter<1000
    oldcov=cov;
    temp = (X'/(cov+eps*eye(D))).*X';
    w = 1./(sum(temp,2)+eps);
    cov = X * (w .* (X'))/N*D;

    [U,S]=svd(cov);
    S1=diag(real(S));
    S1((dd+1):end)=mean(S1((d+1):end))/gam;

    cov = U*diag(S1)*U';
    cov = cov/trace(cov);
    iter=iter+1;
    res{iter} = cov;
end

end