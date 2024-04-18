function [tij] = construct_Tmat_nc(data,res,sym)
tij=cell(data.n);
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
            
            E=projectToEssential(E); 
            E=E/norm(E,'fro');
            
        elseif sym==3
            E=s3(data.E_est,i,j); 
            
            E=projectToEssential(E); 
            
        end
        
        if (sum(sum(~isnan(E)))==0)||(sum(E(:))==0)
                continue;
        end
        
        

        Hgt=res.R_out{i}*res.R_out{j}';

        Ta=E*pinv(Hgt);
        Tc=0.5*(Ta-Ta');
   
        T=[Tc(3,2),Tc(1,3),Tc(2,1)]';
 
     
     
     if sym==3
         tij{i,j}=T;
     else
         tij{i,j}=T/norm(T);
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

