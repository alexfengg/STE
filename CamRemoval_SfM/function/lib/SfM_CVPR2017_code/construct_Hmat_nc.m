function [Hmat,tij] = construct_Hmat_nc(data,sym)
Hmat=eye(3*data.n); tij=cell(data.n);
%% To do : sign for Hmat will it be a problem?
for i=1:data.n
    for j=1:data.n%i+1:data.n
%         if isempty(data.corr{i,j})
%             Hmat(s3(i),s3(j))=eye(3); Hmat(s3(j),s3(i))=eye(3);
%             continue;
% 
%         end
 if i==j
     continue;
 end
        if sym==1 %comes from data
            E=data.E_est{i,j};
            
        elseif sym==2 %comes from res, converting to essential format
            E=s3(data.E_est,i,j); 
            
            E=projectToEssential(E); E=E/norm(E,'fro');
        end
        
        if (sum(sum(~isnan(E)))==0)||(sum(E(:))==0)
                continue;
        end
        
        
        
        [Tc,R1,R2] = getCameraParametersFromE(E);
        
%         err1 = norm(R1 - s3(data.Hmat,i,j),'fro');
%         err2 = norm(R2 - s3(data.Hmat,i,j),'fro');
%          err10 = norm(R1 + s3(data.G_gt,i,j),'fro');
%         err20 = norm(R2 + s3(data.G_gt,i,j),'fro');

        Hgt=data.R(:,:,i)*data.R(:,:,j)';

        err1 = norm(R1 - Hgt,'fro');
        err2 = norm(R2 - Hgt,'fro');
         err10 = norm(R1 + Hgt,'fro');
        err20 = norm(R2 + Hgt,'fro');
        
        if err10<err1
            R1=-R1; err1=err10;
        end
        if err20<err2
            R2=-R2; err2=err20;
        end
        
        if err1<err2
            Hmat(s3(i),s3(j))=R1; Hmat(s3(j),s3(i))=R1';
        else
            Hmat(s3(i),s3(j))=R2; Hmat(s3(j),s3(i))=R2';
        end
        
        if det(s3(Hmat,i,j))<0
            Hmat(s3(i),s3(j))=-s3(Hmat,i,j);
        end
   
        T=[Tc(3,2),Tc(1,3),Tc(2,1)]';
 tijgt=data.R(:,:,i)'*(data.t(:,j)-data.t(:,i));
 Et=crossProdMat(T)*s3(Hmat,i,j); Et=Et/norm(Et,'fro');
 
 a=norm(E-Et); b=norm(E+Et);
 if a<b
     tij{i,j}=T/norm(T);
 else
     tij{i,j}=-T/norm(T);
 end
        
        %% correct
%         if isempty(R)
%            Hmat(s3(i),s3(j))=eye(3); Hmat(s3(j),s3(i))=eye(3);
%        end
%         Hmat(s3(i),s3(j))=R; Hmat(s3(j),s3(i))=R';
%         err=min(norm(data.Hmat(s3(i),s3(j))-R,'fro'),norm(data.Hmat(s3(i),s3(j))+R,'fro'));
%         fprintf('Err : %.4e\n',err);
     end
    %disp(num2str(i));
end

