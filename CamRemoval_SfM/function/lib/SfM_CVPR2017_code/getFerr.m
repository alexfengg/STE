function f=getFerr(E_est,E_gt)

% E_est1=E_est*(norm(E_gt,'fro')/norm(E_est,'fro'));
% %same scale
% 
% f=norm(E_est1-E_gt,'fro');
E_est= projectToEssential(E_est);
E_gt= projectToEssential(E_gt);

f1=norm(E_est/norm(E_est,'fro') - E_gt/norm(E_gt,'fro'),'fro');
f2=norm(E_est/norm(E_est,'fro') + E_gt/norm(E_gt,'fro'),'fro');
f=min(f1,f2);