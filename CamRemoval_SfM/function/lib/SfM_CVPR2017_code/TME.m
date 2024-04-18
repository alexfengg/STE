function [cov,iter]=TME(X)
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

    cov = cov/trace(cov);

    iter=iter+1;
end



end