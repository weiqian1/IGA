function vgrad = evalgrad(X, y, ftype, beta, beta0, lambda)
n = size(X,1);
F = find(beta~=0); 
if strcmp(ftype,'LS')
    if isempty(F)
        vgrad = X'*(beta0-y)/n+lambda*beta;
    else
        vgrad = X'*(X(:,F)*beta(F)+beta0-y)/n+lambda*beta;
    end    
elseif strcmp(ftype,'Logistic')
    if isempty(F)
        eta = repmat(beta0,n,1);
    else
        eta = X(:,F)*beta(F)+beta0;
    end
    vgrad = X'*(-y./(1+exp(y.*eta)))/n+lambda*beta;    
else
    fprintf('function type error!\n');    
end

end