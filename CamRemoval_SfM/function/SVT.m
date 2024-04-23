function X = SVT(M,Omega,delta)
% A Singular Value Thresholding Algorithm for Matrix Completion
% omega1 = Omega(1:3:end,1:3:end);

maxiter = 1000;
[m,n] = size(M);
n=n/3;
p = sum(Omega(:));

if nargin<3
    delta = 1*m*n/p;
end

if isempty(delta)
    delta = 1*m*n/p;
end

tau = 10*m;
OmegaMnorm = norm(Omega.*M,'fro');
k0 = tau/delta/OmegaMnorm;

Y = k0*delta*Omega.*M;

stringNum = 0;
relerr = 0.9;
for k=1:maxiter
    X = svtopt(Y,tau);
    relerr1 = relerr;
    relerr = norm(Omega.*(M-X),'fro')/OmegaMnorm;
    fprintf(repmat('\b', 1, stringNum)); % Erase previous status
    fprintf('Relative Error: %.4f, iter: %d',relerr,k)
    stringNum = length(sprintf('Relative Error: %.4f, iter: %d',relerr,k));
    % norm(Omega.*(M-X),'fro')/norm(Omega.*M,'fro')
    if relerr<0.1
        break
    end
    if relerr1<relerr && k>50
        break
    end
    Y = Y+delta*Omega.*(M-X);
end

end

function Y = svtopt(X,tau)
% [U,S,V] = svd(X);
[U,S,V] = randpca(X,100);
Y = U*diag(max(diag(S)-tau,0))*V';
end