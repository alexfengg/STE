function [var,data] = initialize_admm_ndset(res,data)
% data is the output of the baseline
n=length(res.R_out);

var.n=n;

U=[]; V=[];
for i=1:n
    U=[U;res.R_out{i}];
    V=[V;res.R_out{i}*crossProdMat(res.T_out(:,i))'];
end
var.X=U*V';
var.E=var.X+var.X';


%normalize E
var.E_est=zeros(3*data.n);
for i=1:n
    for j=i+1:n
        if data.keep(i,j)
            E_tmp{i,j}=data.E_est{i,j}/norm(data.E_est{i,j},'fro');


            if isempty(E_tmp{i,j})
                E_tmp{i,j}=zeros(3);
            end
            E_tmp{j,i}=E_tmp{i,j}';

            var.E_est(s3(i),s3(j))=E_tmp{i,j};
            var.E_est(s3(j),s3(i))=E_tmp{j,i};
        end
    end
end


var.W=logical(data.keep);


var.W=full(var.W);

%initialize lam
lam=zeros(data.n);
for i=1:n
    for j=i+1:n
        lam(i,j)=sum(sum(s3(var.E_est,i,j).*s3(var.E,i,j)))/(norm(s3(var.E,i,j),'fro')^2);
        lam(j,i)=lam(i,j);
    end
end

var.lam=lam;

var.keep=data.keep;

var.gamma=zeros(size(var.X));
var.Y=var.X;
var.W=var.W/norm(double(var.W),'fro');

end

