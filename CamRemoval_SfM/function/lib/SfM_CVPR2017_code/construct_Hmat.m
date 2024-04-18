function Hmat = construct_Hmat(data,sym)
Hmat=eye(3*data.n);
%% To do : sign for Hmat will it be a problem?
for i=1:data.n
    for j=i+1:data.n
        if isempty(data.corr{i,j})
            Hmat(s3(i),s3(j))=eye(3); Hmat(s3(j),s3(i))=eye(3);
            continue;

        end
        if sym==1 %comes from data
            E=data.E_est{i,j};
        elseif sym==2 %comes from res, converting to essential format
            E=s3(data.E_est,i,j); E=projectToEssential(E); E=E/norm(E,'fro');
        end

        %pt=data.corr{1,2}; pt1=pt(:,:,1); pt2=pt(:,:,2);
        %idx=randperm(size(pt1,2)); idx=idx(1:10);
        %[T,R,R1,R2] = getCameraParametersFromE_new(E,pt1(:,idx),pt2(:,idx),data.K(:,:,i),data.K(:,:,j));
        [T,R1,R2] = getCameraParametersFromE(E);

        err1 = norm(R1 - s3(data.G_gt,i,j),'fro');
        err2 = norm(R2 - s3(data.G_gt,i,j),'fro');
        err10 = norm(R1 + s3(data.G_gt,i,j),'fro');
        err20 = norm(R2 + s3(data.G_gt,i,j),'fro');

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
            Hmat(s3(j),s3(i))=-s3(Hmat,j,i);
        end



        %% correct
        %         if isempty(R)
        %            Hmat(s3(i),s3(j))=eye(3); Hmat(s3(j),s3(i))=eye(3);
        %            continue;
        %            disp('Empty');
        %        end
        %         Hmat(s3(i),s3(j))=R; Hmat(s3(j),s3(i))=R';
        %         err=min(norm(data.Hmat(s3(i),s3(j))-R,'fro'),norm(data.Hmat(s3(i),s3(j))+R,'fro'));
        %         fprintf('Err : %.4e\n',err);
        %         if err >1
        %             1;
        %         end
    end
    %disp(num2str(i));
end

