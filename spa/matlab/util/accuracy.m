function [al, ag, d50, d90, c50, c90] = accuracy(loc, ploc, grp, pgrp, ppos)


N = size(loc, 1);

al = sum((loc(:, 1) == ploc(:, 1)) & (loc(:, 2) == ploc(:, 2))) / N;
ag = sum(grp == pgrp) / N;


ddist = zeros(N, 1);
cdist = zeros(N, 1);
for i = 1:N
    ddist(i) = pos2dist(loc(i, 2), loc(i, 1), ploc(i, 2), ploc(i, 1));
    cdist(i) = pos2dist(loc(i, 2), loc(i, 1), ppos(i, 2), ppos(i, 1));
end

ddist = sort(ddist, 'ascend');
cdist = sort(cdist, 'ascend');

i50 = round(0.5 * N);
i90 = round(0.9 * N);
if i50 < 1; i50 = 1; end;
if i90 < 1; i90 = 1; end;
if i50 > N; i50 = N; end;
if i50 > N; i90 = N; end;

d50 = ddist(i50);
d90 = ddist(i90);

c50 = cdist(i50);
c90 = cdist(i90);
