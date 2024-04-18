function [m,n,cflag]=T(params,v)

if nargin==1
    n=params.n;
    m=params.m;
    cflag=0;
    return;
end

n=params.n;
m=params.m;
k=params.k;

sigma = 1.e-04;

ll=min(m,n);

s=zeros(ll,1);
for j = 1:k
  s(j) = sigma^((j-1)/(k-1));
end

for j = k+1:ll
  s(j) = sigma * (n-j)/(n-k-1);
end

A = sparse(1:ll,1:ll,s,m,n);
w=A*v;
m=w;