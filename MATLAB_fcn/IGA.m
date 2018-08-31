function path = IGA(X, y, ftype, altype, G, eps, flag, nu, wt, lambda, QUIET)
% This function implements IGA algorithm to find path and coefficients

% X: predictor
% y: response
% ftype: function types, choose 'LS' or 'Logistic'
% altype: algorithm type, choose 'obj' for plain-vanilla IGA or 'gdt' for
% Gradiant-based IGA
% G: a p-dim vector of group numbers, each element from 1 to m
% eps: stopping criteria, currently indicating maximal step size
% flag: 1 is forward-backward, 0 is forward
% nu: interactive parameter
% lambda: a ridge factor
% wt: specify expert opinion list: if a group g is in the list, set
% wt(g)>1; otherwise, set wt(g)=1
% QUIET: prints IGA path if 0 

maxit = max(G)+1;
if ~QUIET
    fprintf('%3s\t%3s\t%10s\t%10s\t%10s\n','Iter','k','objval','forwardNum','backwardNum');
end
path.forwardNum = 0;
path.backwardNum = 0;
p = length(G);
n = length(y);
beta = zeros(p,1);  
if strcmp(ftype, 'LS')
    beta0 = mean(y);
elseif strcmp(ftype, 'Logistic')
    beta0 = glmfit(ones(n,1),(y+1)/2,'binomial','constant','off');
else
    error('function type error!');
end

FG = []; % selected variable groups
k = 1; 
objval = evalobj(X, y, ftype, beta, beta0, lambda);
top = [];

% some initialization
path.beta = cell(maxit,1); % coefficients
path.beta0 = NaN(maxit,1); % beta0's
path.objval = NaN(maxit,1); % objective values
path.FG = cell(maxit,1); % variable selection results
path.top = cell(maxit,1);
path.delta =  NaN(maxit,1); % delta_k values
path.solpath = []; % record solution path
path.maxs = 0; % maximum sparsity
path.sparsity = NaN(maxit,1);

path.beta{k} = beta;
path.beta0(k) = beta0;
path.objval(k) = objval;
path.FG{k} = FG;
path.delta(k) = 0;
path.top{k} = top;
path.sparsity(k) = 0;

iter = 0;
major = zeros(max(G),1);
if strcmp(ftype,'Logistic')
    for g = 1:max(G)
        idx = find(G==g);
        Hg = 1/4/n*X(:,idx)'*X(:,idx);
        major(g) = eigs(Hg,1);
    end
end

while iter <= maxit
    iter = iter+1;    
    grpval = zeros(max(G),1);       
    if strcmp(altype,'obj')
        for g = 1:max(G)
            if ~ismember(g, path.FG{k})
                grpval(g) = path.objval(k) - solve_one(X,y,ftype,beta,beta0,g,G,lambda,major(g));
            end
        end 
    elseif strcmp(altype,'gdt')
        vgrad = evalgrad(X,y,ftype,beta,beta0,lambda);
        for g = 1:max(G)
            if ~ismember(g, path.FG{k})
                grpval(g) = norm(vgrad(G==g),2);
            end
        end 
    else 
        error('algorithm type error!');
    end        
    [a, I] = sort(grpval,'descend');
    
    if length(path.solpath) >= min(eps,n) 
        break;
    end
    
    % generate candidate set
    if (abs(nu-1)<1e-8) % no variable preference
        selg = I(1);
    else % incorporate expert opinion
        top = find(grpval>=a(1)*nu);
        if size(top)<=1
            selg = I(1);
        else            
            selg = get_one(wt(top),top,grpval(top));
        end
    end
    FG = sort([FG,selg]);   
    [objval,beta,beta0] = solve_set(X,y,ftype,FG,G,lambda);
    delta = path.objval(k)-objval;
    
    % record
    path.beta{k+1} = beta;
    path.beta0(k+1) = beta0;
    path.objval(k+1) = objval;
    path.FG{k+1} = FG;
    path.delta(k+1) = delta;
    path.top{iter} = top;
    path.solpath = [path.solpath,selg];
    path.forwardNum = path.forwardNum+1;
    path.sparsity(k+1) = k;
    path.maxs = max(path.maxs,k);   
    k = k+1;
    
    % backward step    
    while flag
        deleted = NaN(length(FG),1);
        for j = 1:length(FG)
            temp = beta;
            temp(G==FG(j)) = 0;
            deleted(j) = evalobj(X,y,ftype,temp,beta0,lambda)-objval;
        end        
        [a,I] = sort(deleted,'ascend');
        if a(1) >= path.delta(k)/2
            break;
        end       
        selg = FG(I(1));
        path.solpath = [path.solpath,-selg];
        FG = FG(FG~=selg);        
        [objval,beta,beta0] = solve_set(X,y,ftype,FG,G,lambda);       
        % recording
        k = k-1;
        path.beta{k} = beta;
        path.beta0(k) = beta0;
        path.objval(k) = objval;
        path.FG{k} = FG;
        path.backwardNum = path.backwardNum+1;        
        if k==1
            break;
        end        
    end    
    if ~QUIET
        fprintf('%3d\t%3d\t%10.4f\t%10d\t%10d\n',iter,k-1,path.objval(k),...
            path.forwardNum,path.backwardNum);
    end    
end

end


