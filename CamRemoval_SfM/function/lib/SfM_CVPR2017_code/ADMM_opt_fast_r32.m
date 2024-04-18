function result = ADMM_opt_fast_r32(var)

max_ite = 500;
threshold = 10^-8;
count = 1;
diff = inf;

[obj] = objective_fun(var); 
%fprintf('count: %d objective: %.5e \n',1,obj);

out = struct;
out.residual = zeros(max_ite,1);
out.obj = zeros(max_ite,1);

%profile on
actv=1;
while count < max_ite && actv && diff > threshold %
        
    

    out.obj(count) = obj;
    obj_old = obj; old_var=var;
    
    var = opt_X(var);
    var = opt_Y(var);
    
    var = opt_gamma(var);

    
    
    count = count + 1;

%    
%      [obj] = objective_fun(var,count);
%      fprintf('count: %d objective: %.5e \n', count,obj);
%     diff = abs((obj - obj_old));
%     out.X_diff(count)=norm(var.X-old_var.X,'fro')/max(norm(var.X),norm(old_var.X));
%     out.Y_diff(count)=norm(var.Y-old_var.Y,'fro')/max(norm(var.Y),norm(old_var.Y));
%     out.gamma_diff(count)=norm(var.gamma-old_var.gamma,'fro')/max(norm(var.X),norm(var.Y));
%     out.XYdiff(count)=norm(var.X-var.Y,'fro')/max(norm(var.Y),norm(var.X));
%     out.Xnorm(count)=norm(var.X,'fro');
%     out.Ynorm(count)=norm(var.Y,'fro');
%      
    

     if mod(count,20)==1
        obj=objective_fun(var);
         %fprintf('count: %d objective: %.5e tau = %.4e \n', count,obj,var.tau);
    end
    
     if mod(count,20)==10
         obj1=objective_fun(var);
%         %fprintf('count: %d objective: %.5e\n', count,obj1);
%         
         diff=abs((obj - obj1))/max(abs(obj),abs(obj1));
%     
         if diff < threshold
             actv=0;
         end
% 
     end
      
     
    
end

result = struct;
result.var = var;
result.out = out;

end

function [obj] = objective_fun(var)

obj1=eval_obj_t(var);

V=svds(var.Y,size(var.Y,1));

obj3 = norm(var.Y - var.X + var.gamma, 'fro')^2;

obj =obj1 + var.tau/2.0*obj3;

end

function [var] = opt_X(var)
n=var.n;
G= var.Y + var.gamma;
GS=(G+G'); GN=(G-G');

%% for all rest
iter=0; f2_old=0; diff1=Inf; MAX=10;
    
while iter < MAX && diff1 > 10^-6
        iter=iter+1;
 % update for AS
 mat1=kron(var.W.*var.lam,ones(3,3)).*var.E_est+var.tau*GS/4;
 mat2=kron(var.W.*var.lam.^2,ones(3,3))+var.tau/4;
 AS=mat1./mat2;
 
 %make block diagonals 0
 idx_block=logical(kron(eye(n),ones(3,3)));
 AS(idx_block)=0;
 
 %update for AN
 AN=GN;
 
 var.X=(AS+AN)/2;
        
var.E=var.X+var.X';
       
       % update for lambda
       m1=var.E_est.*var.E; m2=var.E.*var.E;
       bl1=conv2(m1,ones(3,3),'valid'); bl1=bl1(1:3:end,1:3:end);
       bl2=conv2(m2,ones(3,3),'valid'); bl2=bl2(1:3:end,1:3:end);
       var.lam=bl1./bl2;
       var.lam(logical(eye(n)))=0;

       f2=eval_obj_t(var);

       diff1=abs(f2-f2_old);
       f2_old=f2;
       
end


end

function var = opt_Y(var)

tmpY = var.X - var.gamma;
 [U, S, V] = svds(tmpY, 3);
 var.Y = U*S*V';

end

function var = opt_gamma(var)
var.gamma = var.gamma + (var.Y - var.X);
end

function f=eval_obj_t(var)
e=var.E_est-var.E.*kron(var.lam,ones(3,3));
w=kron(var.W,ones(3,3));
f=0.5*sum(sum(w.*e.^2));
end


       
