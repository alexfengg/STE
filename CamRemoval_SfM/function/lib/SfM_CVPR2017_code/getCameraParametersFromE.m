function [T R1 R2] = getCameraParametersFromE(E)
W = [0 -1 0; 1 0 0; 0 0 1];
Z = [0 1 0; -1 0 0; 0 0 0];

[U,S,V] = svd(E);

T  = U*Z*U'*trace(S)/2;
R1 = U*W*V';
R2 = U*W'*V';

R1 = R1 / det(R1);
R2 = R2 / det(R2);
end