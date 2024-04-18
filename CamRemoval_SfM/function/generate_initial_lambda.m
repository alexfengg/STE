function lam1 = generate_initial_lambda(res,data_lud,data,RemComp)
% data_lud is the output of the LUD pipline, some cameras are removed to make
% the graph is parallel rigid.
% RemComp: the location of nodes are remained.
% data_lud: the output of the data (some cameras are removed from LUD)
% data: whole dataset


n=length(res.R_out);

U=[]; V=[];
for i=1:n
    U=[U;res.R_out{i}];
    V=[V;res.R_out{i}*crossProdMat(res.T_out(:,i))'];
end
var.X=U*V';
var.E=var.X+var.X';

var.E_est=zeros(3*data_lud.n);
E_tmp={};

for i=1:n
    for j=i+1:n
        if data_lud.keep(i,j)
            Ftmp = data_lud.E_est{i,j};
            E_tmp{i,j}=Ftmp/norm(Ftmp,'fro');
            
            if isempty(E_tmp{i,j})
                E_tmp{i,j}=zeros(3);
            end
            E_tmp{j,i}=E_tmp{i,j}';

            var.E_est(s3(i),s3(j))=E_tmp{i,j};
            var.E_est(s3(j),s3(i))=E_tmp{j,i};
        end
    end
end

%initialize lam for data_lud
lam=zeros(data_lud.n);
for i=1:n
    for j=i+1:n
        lam(i,j)=sum(sum(s3(var.E_est,i,j).*s3(var.E,i,j)))/(norm(s3(var.E,i,j),'fro')^2);
        lam(j,i)=lam(i,j);
    end
end


lam1 = zeros(data.n);
Idx = [1:data.n];
Idx1 = Idx'.*RemComp;
Idx1(Idx1==0)=[];
Idx2 = Idx'.*(~RemComp);
Idx2(Idx2==0)=[];

lam1(Idx1,Idx1) = lam;

for u = 1:length(Idx2)
    for v = 1:length(Idx2)
        i=Idx2(u);
        j=Idx2(v);
        E_tmp1 = data.E_est{i,j}/norm(data.E_est{i,j},'fro');
        E_tmp2 = data.E(s3(i),s3(j));
        if ~isempty(E_tmp1)
            lam1(i,j) = sum(E_tmp1.*E_tmp2,'all')/norm(E_tmp2,'fro')^2;
        end
    end
end

end