function [result]=iterative_tukey(var,c1)

for iter=1:3
%find s
var.S=find_s(var);

%update weights - can make it before 1st iteration also
c=4.685;
%c=1.345;
var.W=update_tukey(var,c);

% do ADMM
[result] = iterative_ADMM_ndset(var,c1);
var=result.var;

end



end


function S=find_s(var)

err=nan(var.n);
for i=1:var.n
    for j=i+1:var.n
        if (var.keep(i,j)==0)
            continue;
        end
        err(i,j)=getFerr(s3(var.E_est,i,j),s3(var.E,i,j)*var.lam(i,j));
        err(j,i)=err(i,j);

    end
end

err=err(:); err(isnan(err))=[];
MAD=mad(err);
S=MAD*1.4826;

end

function w=update_tukey(var,c)
err=nan(var.n); w=nan(var.n);
for i=1:var.n
    for j=i+1:var.n
        if (var.keep(i,j)==0)
            continue;
        end
        err(i,j)=getFerr(s3(var.E_est,i,j),s3(var.E,i,j)*var.lam(i,j))/var.S;

        % tukey update weights
        if err(i,j)<=c
            w(i,j)=(1-(err(i,j)/c)^2)^2;
        else
            w(i,j)=0;
        end

        % huber update weights
%         if err(i,j)<=c
%             w(i,j)=1;
%         else
%             w(i,j)=c/abs(err(i,j));
%         end
    end
end
w(isnan(w))=0;
w=w+w';


end
