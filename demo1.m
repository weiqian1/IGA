%%%%% This is an example to run IGA for sparse linear regression

%%% simulate a high-dimensional data
rng(999); 
n = 300; 
p = 1000;
si = 5;
K = p/si;
s = 11; 
Sigma = eye(p);
rho = 0.5;
for i = 1:(p-1)
    for j = (i+1):p
        Sigma(i,j) = rho^(j-i);
        Sigma(j,i) = Sigma(i,j);
    end
end
X = mvnrnd(zeros(p,1),Sigma,n); % generate feature covariates
G = reshape(repmat((1:K)', 1, si)', [], 1);
F = zeros(p,1);
for i = 1:s
    F(10*(i-1)+1 : 10*(i-1)+si) = 1;
end
FG_true = unique(G(F==1));
beta = zeros(p,1);
for i = 1:length(FG_true)
    j = FG_true(i);
    jlen = sum(G==j);
    beta(G==j) = ((rand(jlen,1)>0.5)*2-1);
end
y = X*beta + normrnd(0,2,n,1); % generate response
    
%%% running IGA (without expert opinion set)
ftype = 'LS';
altype = 'obj';
Kfold = 3;
standard = 0;
flag = 1;
eps = 20;
cvres = IGAcv(X,y,ftype,altype,Kfold,standard,G,eps,flag,1,ones(K,1),0);
% estimated coefficients
path = cvres.res.solpath;
beta0 = cvres.beta0;
beta = cvres.beta;
    
%%% running IGA (with expert opinion set)
wt = ones(K,1); % incorporate expert opinion on variable importance
fntemp = ceil(length(FG_true)*3/5);
wt(FG_true(1:fntemp)) = 3;  % wt > 1 corresponds to expert set
nu_cand = (6-(1:5))*.2;
cvres = IGAcv(X,y,ftype,altype,Kfold,standard,G,eps,flag,nu_cand,wt,0);
% estimated coefficients
path = cvres.res.solpath;
beta0 = cvres.beta0;
beta = cvres.beta;    
    
%%% running Gradient-based IGA to speed up computation
ftype = 'LS';
altype = 'gdt';
Kfold = 3;
standard = 0;
flag = 1;
eps = 20;
cvres = IGAcv(X,y,ftype,altype,Kfold,standard,G,eps,flag,1,ones(K,1),0);
% estimated coefficients
path = cvres.res.solpath;
beta0 = cvres.beta0;
beta = cvres.beta;

















