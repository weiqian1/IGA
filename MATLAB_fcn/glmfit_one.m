function res = glmfit_one(x, y, off_set, lambda, maj)  
n = length(y);
d = size(x,2);
maxit = 1000;
tol = 1e-10;
res = zeros(d,1);
nyx = -repmat(y,1,d).*x; 
yax = (off_set+x*res).*y; 
for i = 1:maxit    
    res0 = res;
    Gj = mean(nyx./repmat(ones(n,1)+exp(yax),1,d))';
    res = maj/(maj+lambda)*(res0-1/maj*Gj);
    diff = res-res0;
    if max(diff.^2) < tol
       break; 
    end    
    yax = yax+x*diff.*y;
end     
end




