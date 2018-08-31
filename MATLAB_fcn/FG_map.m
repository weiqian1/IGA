function F = FG_map(FG, G)
F = zeros(length(G), 1);
if iscell(G)
    for g = 1 : length(G)
        F(g) = sum(ismember(G{g},  FG)) > 0;
    end
else
    for i = 1 : length(FG)
        idx = find(G == FG(i));
        F(idx) = 1;
    end
end
F = find(F == 1);
end
