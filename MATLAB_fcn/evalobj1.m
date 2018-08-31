function obj = evalobj1(X, y, ftype, beta, beta0)
n = size(X,1);
F = find(beta ~= 0); 
if strcmp(ftype, 'LS')
    obj = 0.5/n * norm(ones(n,1)*beta0 + X(:, F)*beta(F) - y, 2)^ 2;
elseif strcmp(ftype, 'Logistic')
    yX = repmat(y, 1, size(X, 2)) .* X;
    obj = 1/n*sum(log(1+ exp(-y*beta0-yX(:, F) * beta(F))));
else
    fprintf('function type error!\n');
end
end