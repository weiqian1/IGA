function [res0,res] = glmfit_set(x, y, lambda) 
n = length(y);
d = size(x,2);
maxit = 1000;
tol = 1e-10;
res0 = 0;
res = zeros(d,1);
Xtemp = [ones(n,1),x];
Hm = 1/4/n*(Xtemp'*Xtemp);
maj = eigs(Hm,1);
nyx = -repmat(y,1,d+1).*Xtemp; 
yax = (repmat(res0,n,1)+x*res).*y; 
for i = 1:maxit   
    res0_t = res0;
    res_t = res;
    Gg = mean(nyx./repmat(ones(n,1)+exp(yax),1,d+1))';
    res0 = res0_t-1/maj*Gg(1);
    res = maj/(maj+lambda)*(res_t-1/maj*Gg(2:d+1));
    diff = [res0-res0_t;res-res_t];
    if max(diff.^2) < tol
       break; 
    end
    yax = yax+Xtemp*diff.*y;  
end
end