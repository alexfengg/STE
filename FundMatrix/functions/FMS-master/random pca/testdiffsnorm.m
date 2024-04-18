% Each pair of outputs must be identical.

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

params=struct('m',n,'n',n,'k',k);
rand('seed',42)
[U,S,V] = pca('T','Tt',params,k,its);
randn('state',0)
adiff = diffsnorm_old(T(params,eye(params.n,params.n)),U,S,V);
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
bdiff = diffsnorm_old(T(params,eye(params.n,params.n)),U,S,V);
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
cdiff = diffsnorm_old(T(params,eye(params.n,params.n)),U,S,V);
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
ddiff = diffsnorm_old(T(params,eye(params.n,params.n)),U,S,V);
disp(ddiff);

params=struct('m',n2,'n',n2,'k',k);
rand('seed',42)
[U,S,V] = pca('T','Tt',params,k,its);
randn('state',0)
ddiff = diffsnorm('T','Tt',params,U,S,V);
disp(ddiff);

