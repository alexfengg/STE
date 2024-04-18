function [cov,obj]=TME(X)
[D,N]=size(X);

initcov=eye(D);
oldcov=initcov-1;
cov=initcov;
iter=1;
eps=10^-10;

obj=[];
while norm((oldcov-cov),'fro')>10^-10 & iter<100
    % iter
    oldcov=cov;
    temp = (X'/(cov+eps*eye(D))).*X';
    w = 1./(sum(temp,2)+eps);

    % cov = (X./repmat(w',[D,1])) * X'/N*D;
    
    cov = X * (w .* (X'))/N*D;

    cov = cov/trace(cov);

    iter=iter+1;
    
    
    % obj(iter) = sum(log(svd(cov)))+mean(log(1./w));
end



end