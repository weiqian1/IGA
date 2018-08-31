function [objval, beta, beta0] = solve_set(X, y, ftype, FG, G, lambda)
n = length(y);
F = FG_map(FG,G);
beta = zeros(length(G),1);
if strcmp(ftype,'LS')
    if isempty(FG)
        beta0 = mean(y);
        objval = 0.5/n*norm(ones(n,1)*beta0-y)^2;
    else
        ymean = mean(y);
        Xtemp = X(:,F);
        Xmean = mean(Xtemp);
        ytemp = y-ymean;
        Xtemp = Xtemp-repmat(Xmean,[n,1]);
        beta(F) = (Xtemp'*Xtemp+lambda*eye(length(F)))\Xtemp'*ytemp;
        beta0 = ymean-Xmean*beta(F);                    
        objval = 0.5/n * norm(ones(n,1)*beta0 + X(:, F)*beta(F) - y)^2+0.5*lambda*norm(beta(F))^2;
    end
elseif strcmp(ftype,'Logistic')
    if isempty(FG)
        beta0 = glmfit(ones(n,1),(y+1)/2,'binomial','constant','off');
        objval = evalobj(X,y,ftype,beta,beta0,lambda);
    else
        [beta0,beta(F)] = glmfit_set(X(:,F),y,lambda);
        objval = evalobj(X,y,ftype,beta,beta0,lambda);
    end   
end
end
