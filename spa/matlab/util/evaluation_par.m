function [predX1Pos, predX2Pos, X1Pos, X2Pos] = evaluation_par(name, loc, grp, Y, loc1, loc2, grp1, grp2, X1, X2, predX1, predX2, opt)


[uLoc, cm, cn]= unique(loc, 'rows');
uGrp = grp(cm);
uNam = name(cm);

params.uLoc = uLoc;
params.uGrp = uGrp;
params.uNam = uNam;



[N, K] = size(X1);

Y = Y(:, 1:opt.dim);

if opt.rotate == 1
    [direction, theta] = find_angle(loc, Y(:, 1:2));

    Y2 = rotate_angle(theta) * direction * Y(:, 1:2)';    
    Y(:, 1:2) = Y2';
    
    X1 = rotate_angle(theta) * direction * X1(:, 1:2)';
    X2 = rotate_angle(theta) * direction * X2(:, 1:2)';
    predX1 = rotate_angle(theta) * direction * predX1(:, 1:2)';
    predX2 = rotate_angle(theta) * direction * predX2(:, 1:2)';
    
    X1 = X1'; X2 = X2'; predX1 = predX1'; predX2 = predX2';
end


%% location prediction for the parents


[X1Pos, X1Loc, X1Grp] = predict(loc, Y, X1, params);
[X2Pos, X2Loc, X2Grp] = predict(loc, Y, X2, params);
[predX1Pos, predX1Loc, predX1Grp] = predict(loc, Y, predX1, params);
[predX2Pos, predX2Loc, predX2Grp] = predict(loc, Y, predX2, params);
    

% revese the parents if better, for symmetric problem
for i = 1:N
    
    dist_now = pos2dist(predX1Pos(i, 1), predX1Pos(i, 2), loc1(i, 1), loc1(i, 2)) ...
        + pos2dist(predX2Pos(i, 1), predX2Pos(i, 2), loc2(i, 1), loc2(i, 2));
    dist_rev = pos2dist(predX2Pos(i, 1), predX2Pos(i, 2), loc1(i, 1), loc1(i, 2)) ...
        + pos2dist(predX1Pos(i, 1), predX1Pos(i, 2), loc2(i, 1), loc2(i, 2));
    
    if dist_now > dist_rev
        
        tmppos = predX1Pos(i, :);
        tmploc = predX1Loc(i, :);
        tmpgrp = predX1Grp(i, :);
        tmpx   = predX1(i, :);
        
        predX1Pos(i, :) = predX2Pos(i, :);
        predX1Loc(i, :) = predX2Loc(i, :);
        predX1Grp(i, :) = predX2Grp(i, :);
        predX1(i, :)    = predX2(i, :);
        
        predX2Pos(i, :) = tmppos;
        predX2Loc(i, :) = tmploc;
        predX2Grp(i, :) = tmpgrp;
        predX2(i, :)    = tmpx;        
    end
end 


%% measure accuracy


U = size(uLoc, 1);
RE = zeros(U*(U+1)/2, 9);

accLoc = zeros(U, 1);
accGrp = zeros(U, 1);
dq50 = zeros(U, 1);
dq90 = zeros(U, 1);
cq50 = zeros(U, 1);
cq90 = zeros(U, 1);
ns = zeros(U, 1);

cnt = 0;
for i = 1:U
    for j = (i+1):U
        ind_now = sum(abs(loc1(:, 1) - uLoc(i, 1))) < 0.1 & sum(abs(loc1(:, 2) - uLoc(i, 2))) < 0.1  ...
            & sum(abs(loc2(:, 1) - uLoc(j, 1))) < 0.1 & sum(abs(loc2(:, 2) - uLoc(j, 2))) < 0.1;

        ind_rev = sum(abs(loc2(:, 1) - uLoc(i, 1))) < 0.1 & sum(abs(loc2(:, 2) - uLoc(i, 2))) < 0.1  ...
            & sum(abs(loc1(:, 1) - uLoc(j, 1))) < 0.1 & sum(abs(loc1(:, 2) - uLoc(j, 2))) < 0.1;
        
        ind = (ind_now | ind_rev);
        
        if sum(ind) > 0

            subLoc = [loc1(ind, :); loc2(ind, :)];
            subpLoc = [predX1Loc(ind, :); predX2Loc(ind, :)];
            subGrp = [grp1(ind, :); grp2(ind, :)];
            subpGrp = [predX1Grp(ind, :); predX2Grp(ind, :)];
            subpPos = [predX1Pos(ind, :); predX2Pos(ind, :)];

            [al, ag, d50, d90, c50, c90] = accuracy(subLoc, subpLoc, subGrp, subpGrp, subpPos);

            ns(i) = sum(ind);
            accLoc(i) = al;
            accGrp(i) = ag;    
            dq50(i) = d50;
            dq90(i) = d90; 
            cq50(i) = c50;
            cq90(i) = c90;

            fprintf(1, 'location: %s, %s (%d), locacc=%f, grpacc=%f, dq50=%f, dq90=%f, cq50=%f, cq90=%f\n', ...
                uNam{i}, uNam{j}, sum(ind), al, ag, d50, d90, c50, c90);
            
            cnt = cnt + 1;
            RE(cnt, :) = [i, j, sum(ind), al, ag, d50, d90, c50, c90];
            
        end
    end
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

for i = 1:cnt
    fprintf(1, '%s & %s & $%d$ & $%.2f (%.2f)$ & $%.2f$ & $%d$ & $%d$ & $%d$ & $%d$ \\\\ \n', ...
        uNam{RE(i, 1)}, uNam{RE(i, 2)}, RE(i, 3), RE(i, 4), sqrt(RE(i, 4) * (1 - RE(i, 4)) / RE(i, 3)), ...
        RE(i, 5), round(RE(i, 6)), round(RE(i, 7)), round(RE(i, 8)), round(RE(i, 9)));
end


fprintf(1, 'Mean & $%.1f$ & $%.2f$ & $%.2f$ & $%d$ & $%d$ & $%d$ & $%d$ \\\\ \n', ...
    sum(ns)/U, mean(accLoc), mean(accGrp), round(mean(dq50)), round(mean(dq90)),...
    round(mean(cq50)), round(mean(cq90)));
fprintf(1, 'Mean when $n \\ge 6$ & $%.1f$ & $%.2f$ & $%.2f$ & $%d$ & $%d$ & $%d$ & $%d$ \\\\ \n', ...
    sum(ns(ind)) / sum(ind), mean(accLoc(ind)), mean(accGrp(ind)), round(mean(dq50(ind))), round(mean(dq90(ind))),...
    round(mean(cq50(ind))), round(mean(cq90(ind))));














