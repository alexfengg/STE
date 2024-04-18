function f=getFFerr(E_est,E_gt,K1,K2)

F_est=inv(K1')*E_est*inv(K2); F_est=F_est/norm(F_est,'fro');
F_gt=inv(K1')*E_gt*inv(K2); F_gt=F_gt/norm(F_gt,'fro');

f1=mean(abs(F_est(:)-F_gt(:))./abs(F_gt(:)));
f2=mean(abs(F_est(:)+F_gt(:))./abs(F_gt(:)));
f=min(f1,f2);