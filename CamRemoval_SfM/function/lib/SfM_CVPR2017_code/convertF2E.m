function [var]=convertF2E(var,data)
var.F=var.E;
Etmp=zeros(size(var.E));
for i=1:var.n
    for j=1:var.n
        Etmp(s3(i),s3(j))=data.K(:,:,i)'*var.T(:,:,i)'*s3(var.E,i,j)*var.T(:,:,j)*data.K(:,:,j);
    end
end
var.E=Etmp;
        