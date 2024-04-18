function [F1]=STE_FundMat(X)
% estimate fundamental matrix by applying STE

[D,N]=size(X);
dd = 8;
initcov=eye(D);
eps=10^-10;
Lam = [2,4,6,7,8,9,10];

for i = 1:length(Lam)

    lam = Lam(i);
    oldcov=initcov-1;
    cov=initcov;
    iter=1;

    while norm((oldcov-cov),'fro')>10^-10 & iter<100
        oldcov=cov;
        temp = (X'/(cov+eps*eye(D))).*X';
        w = 1./(sum(temp,2)+eps);
        cov = X * (w .* (X'))/N*D;

        [U,S]=svd(cov);
        S1=diag(real(S));
        S1((dd+1):end)=mean(S1((dd+1):end))/lam;

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

[U,S,V] = svd(cov);
F = reshape(V(:,end),[3,3])';
[U,S,V] = svd(F);
S(3,3)=0;
F1 = U*S*V';

end