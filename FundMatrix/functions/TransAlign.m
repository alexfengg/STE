function [T_out,R_align,err] = TransAlign(T_est, T_gt)

% Translation Alignment
% minimize sum_i ||tij - R*tij^*||_F^2
% Input: T_est: 3 by n matrix
%
% Output: T_out

% Remove columns with NaN values, some of T_gt are NaN since some rotation
% matrices from the datasets are not available

T_gt1 = T_gt(:, ~any(isnan(T_gt)));
T_est1 = T_est(:, ~any(isnan(T_gt)));

A = T_gt1*T_est1';
[U,~,V] = svd(A);

if det(U*V')<0
    V(:,end) = -V(:,end);
end

R_align = U*V';

T_out = R_align*T_est1;
err = acos(abs(sum(T_out.*T_gt1)))/pi*180;

end