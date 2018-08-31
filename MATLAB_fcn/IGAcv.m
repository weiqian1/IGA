 function cvres = IGAcv(X, y, ftype, altype, Kfold, standard, G, eps, flag, nu_cand, wt, lambda)
% This function implements IGA algorithm with tuning parameter selection by cross validation

% X: predictor
% y: response
% ftype: function types, 'LS' or 'Logistic'
% altype: algorithm type, choose 'obj' for plain-vanilla IGA or 'gdt' for
% Gradiant-based IGA
% Kfold: number of folds for cross validation
% standard: whether to standadize design matrix first
% G: a p-dim vector of group numbers, each element from 1 to m
% eps: stopping criteria, currently indicating maximal step size
% flag: 1 is forward-backward, 0 is forward
% nu_cand: candidate values for interactive parameter
% lambda: a ridge factor
% wt: specify expert opinion list: if a group g is in the list, set
% wt(g)>1; otherwise, set wt(g)=1

[n,p] = size(X);
% standadize X to have mean 0 and sd 1
if (standard)
    xbar0 = mean(X);
    sd0 = std(X)*sqrt((n-1)/n);
    X = (X-repmat(xbar0,n,1))./repmat(sd0,n,1);
else
    sd0 = ones(p,1);
    xbar0 = zeros(p,1);
end

% make the splits and start CV
ind = randsample(n,n);
flen = floor(n/Kfold);
fstart = zeros(1,Kfold);
fend = zeros(1,Kfold);
for k = 1:(Kfold-1)
    fstart(k) = (k-1)*flen+1;
    fend(k) = k*flen;
end
fstart(Kfold) = (Kfold-1)*flen+1;
fend(Kfold) = n;
klen = eps+1;
nulen = length(nu_cand);
cvloss = Inf(klen*nulen,Kfold);

parfor k = 1:Kfold 
    [~, id] = lastwarn;
    warning('off', id);
    indval = ind(fstart(k):fend(k));
    indtr = setdiff(1:n,indval);   
    Xtr = X(indtr,:);
    Xval = X(indval,:);
    Ytr = y(indtr);
    Yval = y(indval);
    nval = length(Ytr);
    cvlossK = Inf(klen*nulen,1);    
    for d = 1:nulen       
        restemp = IGA(Xtr,Ytr,ftype,altype,G,eps,flag,nu_cand(d),wt,lambda,1);        
        for s = 1:(restemp.maxs+1)
            cvlossK(klen*(d-1)+s) = nval*evalobj1(Xval,Yval,ftype,restemp.beta{s},restemp.beta0(s)); 
        end    
    end
    cvloss(:,k) = cvlossK;    
end

cverr = sum(cvloss,2);
cverr = reshape(cverr,klen,nulen);
if p >= n
for i = 1:nulen
    cverr(:,i) = movingmean(cverr(:,i),5);
end
end
cverr = cverr(:);
[~,cvix] = min(cverr);
dix = floor((cvix-1)/klen)+1;
nu = nu_cand(dix);
sbest = mod((cvix-1),klen)+1;
fprintf('regression type = %s; the best nu is %d, and the best k (including intercept) is %d.\n', ftype, nu, sbest);
res = IGA(X,y,ftype,altype,G,eps,flag,nu,wt,lambda,0);
sbestm = find(~isnan(res.beta0),1,'last');
sbest = min(sbestm,sbest);
beta0 = res.beta0(sbest);
beta = res.beta{sbest};
if (standard)
    beta = beta./sd0;
    beta0 = beta0-xbar0*beta;
end

% recording
cvres.nu = nu;
cvres.sbest = sbest;
cvres.res = res;
cvres.sd0 = sd0;
cvres.xbar0 = xbar0;
cvres.ftype = ftype;
cvres.altype = altype;
cvres.Kfold = Kfold;
cvres.standard = standard;
cvres.nu_cand = nu_cand;
cvres.cverr = reshape(cverr,[klen,nulen]);
cvres.wt = wt;
cvres.lambda = lambda;
cvres.beta = beta;
cvres.beta0 = beta0;

end




