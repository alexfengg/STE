function [W,err,errGT,errGTB]=estimate_weights(var,p)

delta=10^-3;
W=zeros(var.n);
for i=1:var.n
    for j=i+1:var.n
        if (var.keep(i,j)~=0)
         W(i,j)=(1/max(delta,(norm(s3(var.E_est,i,j)-s3(var.E,i,j)*var.lam(i,j),'fro')^(2-p))));
        end
        
    end
end
W=W+W';

W(~var.keep)=0;
W(logical(eye(size(W,1))))=0;

W=W/norm(W,'fro');

end