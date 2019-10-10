function [pPos, Y] = evaluation(name, loc, grp, Y, opt)

%%
% loc  - N * 2 matrix indicating the latitude and longitude of the sample
% grp  - N * 1 matrix indicating the locational group of the sample
% Y    - N * L matrix resulted from dimensional reduction
% opt  - options
%        opt.rotate  - 1 for rotation, 0 otherwise
%        opt.dim     - no. fo used dimensions for regression

%%

% get available locations and groups
[uLoc, cm, cn]= unique(loc, 'rows');
uGrp = grp(cm);
uNam = name(cm);

params.uLoc = uLoc;
params.uGrp = uGrp;



[N, L] = size(Y);
Y = Y(:, 1:opt.dim);

if opt.rotate == 1
    [direction, theta] = find_angle(loc, Y(:, 1:2));

    Y2 = rotate_angle(theta) * direction * Y(:, 1:2)';    
    Y(:, 1:2) = Y2';    
end



%% leave-one-out prediction

pPos = zeros(N, 2);
pLoc = zeros(N, 2);
pGrp = zeros(N, 1);


for i = 1:N    
    trLoc = loc;
    trY = Y;    
    trLoc(i, :) = [];
    trY(i, :) = [];
    
    tsY = Y(i, :);
    
    [tsPos, tsLoc, tsGrp] = predict(trLoc, trY, tsY, params);

    pPos(i, :) = tsPos;
    pLoc(i, :) = tsLoc;
    pGrp(i, :) = tsGrp;
end


%% 



U = size(uLoc, 1);
RE = zeros(U, 9);

accLoc = zeros(U, 1);
accGrp = zeros(U, 1);
dq50 = zeros(U, 1);
dq90 = zeros(U, 1);
cq50 = zeros(U, 1);
cq90 = zeros(U, 1);
ns = zeros(U, 1);
for i = 1:U
    ind = (loc(:, 1) == uLoc(i, 1)) & (loc(:, 2) == uLoc(i, 2));
    
    subLoc = loc(ind, :);
    subpLoc = pLoc(ind, :);
    subGrp = grp(ind, :);
    subpGrp = pGrp(ind, :);
    subPos = pPos(ind, :);
    
    [al, ag, d50, d90, c50, c90] = accuracy(subLoc, subpLoc, subGrp, subpGrp, subPos);
    
    ns(i) = sum(ind);
    accLoc(i) = al;
    accGrp(i) = ag;    
    dq50(i) = d50;
    dq90(i) = d90; 
    cq50(i) = c50;
    cq90(i) = c90;
    
    fprintf(1, '%s (%d), locacc=%f, grpacc=%f, dq50=%f, dq90=%f, cq50=%f, cq90=%f\n', ...
        uNam{i}, sum(ind), al, ag, d50, d90, c50, c90);
    RE(i, :) = [uLoc(i, 1), uLoc(i, 2), sum(ind), al, ag, d50, d90, c50, c90];
end

fprintf(1, 'AVERAGE\n');
fprintf(1, '(%d), locacc=%f, grpacc=%f, dq50=%f, dq90=%f, cq50=%f, cq90=%f\n', ...
    sum(ns), mean(accLoc), mean(accGrp), mean(dq50), mean(dq90), mean(cq50), mean(cq90));

ind = (ns >= 6);

fprintf(1, 'AVERAGE for n > 6\n');
fprintf(1, '(%d), locacc=%f, grpacc=%f, dq50=%f, dq90=%f, cq50=%f, cq90=%f\n', ...
    sum(ns(ind)), mean(accLoc(ind)), mean(accGrp(ind)), mean(dq50(ind)), ...
    mean(dq90(ind)), mean(cq50(ind)), mean(cq90(ind)));



[C, I] = sort(-RE(:, 3));
RE = RE(I, :);
uNam = uNam(I);

for i = 1:U
    fprintf(1, '%s & $%d$ & $%.2f (%.2f)$ & $%.2f$ & $%d$ & $%d$ & $%d$ & $%d$ \\\\ \n', ...
        uNam{i}, RE(i, 3), RE(i, 4), sqrt(RE(i, 4) * (1 - RE(i, 4)) / RE(i, 3)) ,...
        RE(i, 5), round(RE(i, 6)), round(RE(i, 7)), round(RE(i, 8)), round(RE(i, 9)));
end
fprintf(1, 'Mean & $%.1f$ & $%.2f$ & $%.2f$ & $%d$ & $%d$ & $%d$ & $%d$ \\\\ \n', ...
    sum(ns)/U, mean(accLoc), mean(accGrp), round(mean(dq50)), round(mean(dq90)),...
    round(mean(cq50)), round(mean(cq90)));
fprintf(1, 'Mean when $n \\ge 6$ & $%.1f$ & $%.2f$ & $%.2f$ & $%d$ & $%d$ & $%d$ & $%d$ \\\\ \n', ...
    sum(ns(ind)) / sum(ind), mean(accLoc(ind)), mean(accGrp(ind)), round(mean(dq50(ind))), round(mean(dq90(ind))),...
    round(mean(cq50(ind))), round(mean(cq90(ind))));


for i = 1:U
    fprintf(1, '%s & %d & %.2f & %.2f \\ \n', ...
        uNam{i}, RE(i, 3), RE(i, 1), RE(i, 2));
end














