function bigF = generateFstacked(data,lam)

n = data.n;
E_est=zeros(3*data.n);
E_tmp={};
bigF=zeros(3*data.n);

for i=1:n
    for j=i+1:n
        if data.keep(i,j)
            % K1 = data.K(:,:,i);
            % K2 = data.K(:,:,j);
            Ftmp = data.E_est{i,j};
            E_tmp{i,j}=Ftmp/norm(Ftmp,'fro');

            if isempty(E_tmp{i,j})
                E_tmp{i,j}=zeros(3);
            end
            E_tmp{j,i}=E_tmp{i,j}';

            E_est(s3(i),s3(j))=E_tmp{i,j};
            E_est(s3(j),s3(i))=E_tmp{j,i};
        end
    end
end

for i=1:n
    for j=i+1:n
        if data.keep(i,j)

            if lam(i,j)~=0
                bigF(s3(i),s3(j))=E_tmp{i,j}/lam(i,j);
                bigF(s3(j),s3(i))=E_tmp{j,i}/lam(j,i);
            end
        end
    end
end

bigF(isnan(bigF))=0;


end