function [objL2,objLp]=eval_obj_B_admm(var,p)
N=var.n;
e=var.E_est-var.E.*kron(var.lam,ones(3,3));
w=kron(var.W,ones(3,3));
objL2=0.5*sum(sum(w.*(e.^2)));
%reduction
residual=zeros(N);
    for i=1:N
        for j=i+1:N
            if var.keep(i,j)==0
                continue;
            end
           residual(i,j)=norm(s3(e,i,j),'fro');
           residual(j,i)=residual(i,j);
        end
    end
    
    objLp=sum((residual(:)).^p);

end

