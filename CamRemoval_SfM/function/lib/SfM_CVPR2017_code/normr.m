function n = normr(x)
[~, col] = size(x);
n = x./repmat(sqrt(sum(x.^2,2)), [1, col]);
end