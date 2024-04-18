function [T,R,R1,R2] = getCameraParametersFromE_new(E,pt1,pt2,k1,k2)
E=projectToEssential(E); E=E/norm(E,'fro');
s=svd(E); E=E/s(1); %makes E have singular values of 1

W = [0 -1 0; 1 0 0; 0 0 1];
Z = [0 1 0; -1 0 0; 0 0 0];

[U,S,V] = svd(E);

T1=U(:,3); T2=-U(:,3);
R1 = U*W*V';
R2 = U*W'*V';

rot(:,:,1)=R1; rot(:,:,2)=R1; rot(:,:,3)=R2; rot(:,:,4)=R2;
t(:,:,1)=T1; t(:,:,2)=T2;t(:,:,3)=T1;t(:,:,4)=T2;
[R,T,correct] = SelectCorrectEssentialCameraMatrix(rot,t,pt1,pt2,k1,k2);
if correct==0
    R=[]; T=[];
end
% 
% R1 = R1 / det(R1);
% R2 = R2 / det(R2);
end