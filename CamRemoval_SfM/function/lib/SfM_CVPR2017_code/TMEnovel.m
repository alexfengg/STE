function [cov,iter]=TMEnovel(X,dd,lam)
[D,N]=size(X);

initcov=eye(D);
oldcov=initcov-1;
cov=initcov;
iter=1;
eps=10^-10;


while norm((oldcov-cov),'fro')>10^-10 & iter<1000

    oldcov=cov;
    temp = (X'/(cov+eps*eye(D))).*X';
    w = 1./(sum(temp,2)+eps);

    % cov = (X./repmat(w',[D,1])) * X'/N*D;
    
    cov = X * (w .* (X'))/N*D;

    [U,S]=svd(cov);
    U=real(U);

    S1=diag(real(S));

    S1((dd+1):end)=mean(S1((dd+1):end))/lam;

    % S1(end) = 0;

    cov=U*diag(S1)*U';

    cov = cov/trace(cov);

    iter=iter+1;
end



end