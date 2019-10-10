function [Comps, Coefs, Mean, Eigens] = mypca(X, p)
% This function is for performing principle component analysis (PCA) on
% given data matrix and get the eigenfaces
% 
% Input: X is d * m matrix where d is the data dimension and m is the
% number of sample. p is the required number of principle components.
%
% Output: Comps is a d * p matrix with principle components in each row.
% Coefs is a p * m matrix with the projected coordinates for each sample
%
% (c) Wen-Yun Yang 2009  wenyun@ucla.edu


[d, m] = size(X);

if length(p) ~= 1
    error('p should be an integer');
end

if p > d
    error('p can not be larger then the data dimension');
end



mX = mean(X, 2);
sX = X - repmat(mX, 1, m);

[u, e, v] = svd(sX'*sX, 0);   % compute m * m matrix's eigenvectors, instead of d*d
[a,b]=sort(-diag(e));
% e = diag(e(b(1:m),b(1:m)));
e = diag(e);
e = e(b(1:m));
U = u(:, b(1:m));   % re-order the eigenvectors by eigenvalues

Comps = sX * U(:, 1:p);  % recover the d*d matrix's eigenvector
Comps = Comps ./ repmat(((sum(Comps.^2)).^0.5), d, 1);  % normalize
Coefs = Comps' * sX;  % projection weight
Mean = mX;
Eigens = e(1:p);


