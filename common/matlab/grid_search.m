function [x, recx, recobj]= grid_search(objective, xlim, ngrid, eps)

% this function finds the minimum point by grid search

% objective - function handler for objective function
% xlim - L * 2 matrix, specifies the lower and upper bound in each dimension
% ngrid - number of regions in each division
% eps -  tolerance of error

fprintf(1, 'Searching ');
for i = 1:size(xlim, 1)
  fprintf(1, '[%.5f %.5f] ', xlim(i, 1), xlim(i, 2));
end
fprintf(1, '\n');

L = size(xlim, 1);

if norm(xlim(:, 2) - xlim(:, 1)) < eps
    x = mean(xlim, 2);
    return;
end

step = (xlim(:, 2) - xlim(:, 1)) / (ngrid + 1);

obj = Inf;
x = 0;
recx = [];
recobj = [];
for i = 1:ngrid^L
    [xt, index] = get_point(xlim, L, ngrid, step, i-1);
    objt = objective(xt);        
    grid_output(xt, objt);
    recx = [recx; xt'];
    recobj = [recobj; objt];
    if objt < obj
        obj = objt;
        x = xt;
    end
end

x = grid_search(objective, [x - step, x + step], ngrid, eps);

function [xt, index] = get_point(xlim, L, ngrid, step, idx)

xt = zeros(L, 1);
index = zeros(L, 1);
for i = 1:L
    j = rem(idx, ngrid);    
    index(i) = j + 1;
    xt(i) = xlim(i, 1) + step(i) * (j + 1);
    idx = (idx - j) / ngrid;
end

function grid_output(xt, objt)

fprintf(1, '(');
for i = 1:length(xt)
  fprintf(1, '%.5f ', xt(i));
end
fprintf(1, ') = %.5f\n', objt);
