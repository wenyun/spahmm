function [Y,eigenvals,neighbors] = lle(X,K,d,tol,distance)

% LLE ALGORITHM
%
% function [Y,eigenvals,neighbors] = lle(X,K,d,tol)
%
% X = data as D x N matrix (D = dimensionality, N = #points)
% K = number of neighbors
% d = embedding dimensionality
% tol = regularizer (defaults to 1e-4)
% Y = embedding as d x N matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PAIRWISE DISTANCES
[D,N] = size(X);

% NEIGHBORS
[sorted,index] = sort(distance);
neighbors = index(2:(1+K),:);


% RECONSTRUCTION WEIGHTS
if (nargin<4), tol=1e-4; end;
W = zeros(K,N);
for i=1:N
  z = X(:,neighbors(:,i))-repmat(X(:,i),1,K);
  C = z'*z;
  C = C + tol*trace(C)*eye(K)/K; % REGULARIZATION
  uW =  C \ ones(K, 1);
  W(:,i) = uW ./ sum(uW);
end;

% COST MATRIX
M = eye(N);
for i=1:N 
  w = W(:,i);
  j = neighbors(:,i);
  M(i,j) = M(i,j) - w';
  M(j,i) = M(j,i) - w;
  M(j,j) = M(j,j) + w*w';
end;

% CALCULATION OF EMBEDDING
options.disp = 0;
[Y,eigenvals] = eigs(M,d+1,0,options);
eigenvals = diag(eigenvals);
[eigenvals,indx] = sort(eigenvals);
eigenvals = eigenvals(2:d+1);
Y = Y(:,indx(2:d+1))'*sqrt(N);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
