function selg = get_one(wtt,toplist,chgval)
    if sum(wtt>1) == 0
        [~,tidx] = max(chgval);
        selg = toplist(tidx);
    else
        toptemp = toplist(wtt>1);
        chgvalt = chgval(wtt>1);
        [~,tidx] = max(chgvalt);
        selg = toptemp(tidx);
    end
end