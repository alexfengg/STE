function E_pro = projectToEssential(E)
if sum(sum(~isnan(E)))==0
    E_pro=E;
    return;
end
[U,W,V]=svd(E);
w=diag(W);
v=[(w(1)+w(2))/2;(w(1)+w(2))/2;0];
E_pro=U*diag(v)*V';
end