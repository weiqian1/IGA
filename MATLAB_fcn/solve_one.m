function objval = solve_one(X, y, ftype, beta, beta0, g, G, lambda, maj)
n = length(y);
if strcmp(ftype,'LS')
    F = find(beta~=0);
    res = y-repmat(beta0,n,1)-X(:,F)*beta(F);  
    idx = find(G==g);
    beta(idx) = (X(:,idx)'*X(:,idx)+lambda*eye(length(idx)))\X(:,idx)'*res;  
    F = find(beta~=0);
    objval = 0.5/n * norm(ones(n,1)*beta0 + X(:, F)*beta(F) - y)^2 + 0.5*lambda*norm(beta(F))^2;       
elseif strcmp(ftype,'Logistic')
    F = find(beta~=0);
    off_set = repmat(beta0,n,1)+X(:,F)*beta(F);
    idx = find(G==g);
    beta(idx) = glmfit_one(X(:,idx),y,off_set,lambda,maj);
    objval = evalobj(X,y,ftype,beta,beta0,lambda);   
end
end