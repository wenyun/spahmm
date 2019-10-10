function [Y, distance] = my_mvu(X, eps, K)


X = normalize(X);  
distance = dissimilar(X, eps, K, 0);


[Y, details] = fastmvu(distance, 2, 'leigsdim', 10, 'eta', 1e-10, 'maxiter', 2000);

