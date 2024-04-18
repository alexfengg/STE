% Each test is eceuted twice. The error should be exactly the same.
% Errors should be at most of the order of sigma.

clear;

n = 10000;
n2 = 1000;
k = 10;
its = 3;
sigma = 1.e-04;

s=zeros(n,1);
for j = 1:k
  s(j) = sigma^((j-1)/(k-1));
end

for j = k+1:n
  s(j) = sigma * (n-j)/(n-k-1);
end

%
% Provide A explicitly
%
A = sparse(1:n,1:n,s,n,n);
rand('seed',42)
[U,S,V] = pca(A,k,its);
randn('state',0)
adiff = diffsnorm(A,U,S,V);
disp(adiff);

A = sparse(1:n,1:n,s,n,n);
rand('seed',42)
[U,S,V] = pca(A,k,its);
randn('state',0)
adiff = diffsnorm(A,U,S,V);
disp(adiff);

B = sparse(1:n,1:n,s,n,2*n);
rand('seed',42)
[U,S,V] = pca(B,k,its);
randn('state',0)
bdiff = diffsnorm(B,U,S,V);
disp(bdiff);

B = sparse(1:n,1:n,s,n,2*n);
rand('seed',42)
[U,S,V] = pca(B,k,its);
randn('state',0)
bdiff = diffsnorm(B,U,S,V);
disp(bdiff);

C = sparse(1:n,1:n,s,2*n,n);
rand('seed',42)
[U,S,V] = pca(C,k,its);
randn('state',0)
cdiff = diffsnorm(C,U,S,V);
disp(cdiff);

C = sparse(1:n,1:n,s,2*n,n);
rand('seed',42)
[U,S,V] = pca(C,k,its);
randn('state',0)
cdiff = diffsnorm(C,U,S,V);
disp(cdiff);

D = diag(s(1:n2));
rand('seed',42)
[U,S,V] = pca(D,k,its);
randn('state',0)
ddiff = diffsnorm(D,U,S,V);
disp(ddiff);

D = diag(s(1:n2));
rand('seed',42)
[U,S,V] = pca(D,k,its);
randn('state',0)
ddiff = diffsnorm(D,U,S,V);
disp(ddiff);

E = [diag(s(1:n2)) zeros(n2)];
rand('seed',42)
[U,S,V] = pca(E,k,its);
randn('state',0)
ediff = diffsnorm(E,U,S,V);
disp(ediff);

E = [diag(s(1:n2)) zeros(n2)];
rand('seed',42)
[U,S,V] = pca(E,k,its);
randn('state',0)
ediff = diffsnorm(E,U,S,V);
disp(ediff);

F = [diag(s(1:n2)); zeros(n2)];
rand('seed',42)
[U,S,V] = pca(F,k,its);
randn('state',0)
fdiff = diffsnorm(F,U,S,V);
disp(fdiff);

F = [diag(s(1:n2)); zeros(n2)];
rand('seed',42)
[U,S,V] = pca(F,k,its);
randn('state',0)
fdiff = diffsnorm(F,U,S,V);
disp(fdiff);


%
% Provide function handle
%

params=struct('m',n,'n',n,'k',k);
rand('seed',42)
[U,S,V] = pca('T','Tt',params,k,its);
randn('state',0)
adiff = diffsnorm('T','Tt',params,U,S,V);
disp(adiff);

params=struct('m',n,'n',n,'k',k);
rand('seed',42)
[U,S,V] = pca('T','Tt',params,k,its);
randn('state',0)
adiff = diffsnorm('T','Tt',params,U,S,V);
disp(adiff);

params=struct('m',n,'n',2*n,'k',k);
rand('seed',42)
[U,S,V] = pca('T','Tt',params,k,its);
randn('state',0)
bdiff = diffsnorm('T','Tt',params,U,S,V);
disp(bdiff);

params=struct('m',n,'n',2*n,'k',k);
rand('seed',42)
[U,S,V] = pca('T','Tt',params,k,its);
randn('state',0)
bdiff = diffsnorm('T','Tt',params,U,S,V);
disp(bdiff);

params=struct('m',2*n,'n',n,'k',k);
rand('seed',42)
[U,S,V] = pca('T','Tt',params,k,its);
randn('state',0)
cdiff = diffsnorm('T','Tt',params,U,S,V);
disp(cdiff);


params=struct('m',2*n,'n',n,'k',k);
rand('seed',42)
[U,S,V] = pca('T','Tt',params,k,its);
randn('state',0)
cdiff = diffsnorm('T','Tt',params,U,S,V);
disp(cdiff);

params=struct('m',n2,'n',n2,'k',k);
rand('seed',42)
[U,S,V] = pca('T','Tt',params,k,its);
randn('state',0)
ddiff = diffsnorm('T','Tt',params,U,S,V);
disp(ddiff);

params=struct('m',n2,'n',n2,'k',k);
rand('seed',42)
[U,S,V] = pca('T','Tt',params,k,its);
randn('state',0)
ddiff = diffsnorm('T','Tt',params,U,S,V);
disp(ddiff);

