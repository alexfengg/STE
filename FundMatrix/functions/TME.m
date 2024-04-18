function [cov,res]=TME(X)
% Tyler's M-estimator (TME)
[D,N]=size(X);

initcov=eye(D);

oldcov=initcov-1;
cov=initcov;
eps=10^-10;%regulirization parameter

iter=0;
res = {};
while norm((oldcov-cov),'fro')>10^-10 & iter<1000

    oldcov=cov;

    temp = (X'/(cov+eps*eye(D))).*X';
    w = 1./(sum(temp,2)+eps);
    cov = X * (w .* (X'));
    cov = cov/trace(cov);

    iter=iter+1;
    res{iter} = cov;
end



end