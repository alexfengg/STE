function [xu1,T1] = normalize_points(xu)
% Input: xuu: 2 by N matrix
% Output: xu1: normalized points
%          T1: Projection matrix

mean_u1 = mean(xu(1,:));
mean_u2 = mean(xu(2,:));
var_u1 = var(xu(1,:)-mean_u1);
var_u2 = var(xu(2,:)-mean_u2);

T1 = [1/sqrt(var_u1), 0,0; 0,1/sqrt(var_u2), 0; 0,0,1]*[1,0,-mean_u1;0,1,-mean_u2;0,0,1];

if size(xu,1)==2
    xuu = [xu;ones(1,size(xu,2))];
else
    xuu = xu;
end

xu1 = T1*xuu;

end