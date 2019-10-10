function [Y, V] = eigmap(X, p, params)


[d, m] = size(X);

cov = snp_distance(X, X, 'corr');


A = adjacency(cov, params);

L = laplacian(A);

options.disp = 0;
[Y, eigenvals] = eigs(L, p+1, 0, options);
eigenvals = diag(eigenvals);
[eigenvals, index] = sort(eigenvals);
V = eigenvals(2:p+1);
Y = Y(:,index(2:p+1))';


function L = laplacian(A)

D = diag(sum(A));

L = D - A;

L = L ./ sqrt(diag(D) * diag(D)');


function A = adjacency(cov, params)

[m, m] = size(cov);

A = spconvert([m, m, 0]);

if strcmp(params.type, 'eps') == 1
    
    [row, col] = find(cov > params.eps);
    
    for i = 1:length(row)
        val = cov(row(i), col(i));
        A(row(i), col(i)) = exp(-(1 - val)^2 / params.t);
    end
    
elseif strcmp(params.type, 'nn') == 1
    
    [sorted, index] = sort(- cov);
    neighbors = index(2:(1 + params.K), :);
    values = - sorted(2:(1 + params.K), :);
    
    for i = 1:m        
        for j = 1:params.K            
            A(i, neighbors(j, i)) = exp(-(1 - values(j, i))^2 / params.t);
            A(neighbors(j, i), i) = A(i, neighbors(j, i));
        end
    end              
end