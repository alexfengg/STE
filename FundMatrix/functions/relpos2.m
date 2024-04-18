function [Ruv,tuv] = relpos2(F,x_u,x_v,K1,K2)
% extract relative pose for given F from the ambiguity
Vec2Skew = @(v) [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];

K1(2,2) = -K1(2,2);
K2(2,2) = -K2(2,2);

E = K1'*F*K2;
[U,~,V] = svd(E);

t1 = U(:,3);R1 = U*[0,-1,0;1,0,0;0,0,1]*V';R1 = R1/det(R1);t1 = t1/det(R1);
t2 = -U(:,3);R2 = U*[0,-1,0;1,0,0;0,0,1]'*V';R2 = R2/det(R2);t2 = t2/det(R2);
count = zeros(2,2);t = [t1,t2]; R = zeros(3,3,2); R(:,:,1) = R1; R(:,:,2) = R2;R0 = eye(3);t0=zeros(3,1);

n_tmp = min(size(x_u,2),50);

for i1=1:n_tmp
    p = x_u(:,i1);
    q = x_v(:,i1);

    for k1 = 1:2
        for k2=1:2
            P0 = R0*[eye(3),-t0];
            Pi = (R(:,:,k2))'*[eye(3),-t(:,k1)];
            [~,~,x] = svd([Vec2Skew([p;1]) * K1 * P0; Vec2Skew([q;1]) * K2 * Pi]);x = x(:,4);
            x = x/x(4); x = x(1:3,1);
            if(Pi(3,1:3)*(x-t(:,k1)) > 0) && (P0(3,1:3)*(x-P0(:,4)) > 0)
                count(k1,k2) = count(k1,k2)+1;
            end
        end
    end

end
[k1,k2] = find(count == max(max(count))); k1=k1(1); k2=k2(1);
J=diag([1,1,-1]);
Ruv = J*R(:,:,k2)*J;
tuv = t(:,k1);

end