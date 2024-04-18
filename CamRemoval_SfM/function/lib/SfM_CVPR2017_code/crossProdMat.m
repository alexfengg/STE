function out = crossProdMat(v, method)

if ~exist('method', 'var')
    method='col';
end

n = size(v,2);

if n ==1
    out = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
elseif strcmp(method, 'row')
    out = zeros(3,3*n);
    for i=0:n-1
        out(:, 3*i+1 : 3*i+3) = crossProdMat(v(:,i+1));
    end
elseif strcmp(method, 'col')
    out = -crossProdMat(v, 'row')';
elseif strcmp(method, 'diag')
    out = zeros(3*n);
    
    for i=0:n-1
        out(3*i+1 : 3*i+3, 3*i+1 : 3*i+3) = crossProdMat(v(:,i+1));
    end  
end

end
