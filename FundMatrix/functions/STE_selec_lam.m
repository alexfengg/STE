function cov=STE_selec_lam(X,d,Gam)
% STE with selection of best gamma out of a set of gammas

[D,N]=size(X);
initcov=eye(D);
eps=10^-10;

for i = 1:length(Gam)
    gam = Gam(i);
    oldcov=initcov-1;
    cov=initcov;
    iter=1;
    while norm((oldcov-cov),'fro')>10^-10 & iter<1000
        oldcov=cov;
        temp = (X'/(cov+eps*eye(D))).*X';
        w = 1./(sum(temp,2)+eps);
        cov = X * (w .* (X'))/N*D;

        [U,S]=svd(cov);
        S1=diag(real(S));
        S1((d+1):end)=mean(S1((d+1):end))/gam;

        cov = U*diag(S1)*U';
        cov = cov/trace(cov);
        iter=iter+1;
    end
    COV{i} = cov;
    C = U(:,end)*U(:,end)'*X;
    dist = ((sum(C.*C,1)).^.5);
    Dis(i,:) = dist;
end

threshold = median(Dis(:));
inliers = sum(Dis<threshold,2);
[B,ind] = max(inliers);
cov = COV{ind};

end