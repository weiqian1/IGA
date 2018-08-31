function a = movingmean(b,k)
    k0 = (k-1)/2;
    n = length(b);
    a = zeros(n,1);
    for i = 1:n
        ix = (i-k0):(i+k0);
        new = (ix >= 1) & (ix <= n);
        ix = ix(new);
        a(i) = mean(b(ix));
    end
end