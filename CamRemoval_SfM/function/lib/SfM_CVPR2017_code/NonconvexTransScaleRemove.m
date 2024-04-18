% This file removes the scale and translation between two point clouds
% using non-convex programming

function [t_fit, t_opt, c_opt, MeanError, MedianError] = NonconvexTransScaleRemove(c0, t0, ti, ti0, p)

n = size(ti,2); d = size(ti,1);

options = optimset('MaxFunEvals',1e4,'Display','iter-detailed');

x_opt = fminunc(@RegCostFnc,[c0;t0],options);

c_opt = x_opt(1); t_opt = x_opt(2:d+1);    

%% Evaluate the fit and the error
t_fit = c_opt*ti+t_opt*ones(1,n);
t_err = ti0 - t_fit;
MeanError = mean(norms(t_err,2,1));
MedianError = median(norms(t_err,2,1));

    function CostVal = RegCostFnc(x)
        CostVal = sum(norms(ti0 - (x(1)*ti+x(2:d+1)*ones(1,n)),2,1).^p);
    end

end