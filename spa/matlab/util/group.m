function [X, ID, TXT, COL, LOC, GRP] = group(X, ID, TXT, COL, LOC, GRP)

[d, m] = size(X);

UTXT = unique(TXT);
pm = zeros(m, m);
pr = 0;

for i = 1:length(UTXT)
    pop = UTXT(i);
    
    pidx = find(strcmp(TXT, pop));
    
    subp = length(pidx);
    
    for j = 1:subp
        pm(pidx(j), pr + j) = 1;
    end
        
    pr = pr + subp;
end

X = X * pm;

[row, col] = find(pm);
ID = ID(row);
TXT = TXT(row);
COL = COL(row, :);
LOC = LOC(row, :);
GRP = GRP(row);
