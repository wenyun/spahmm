function [Y,eigenvals] = lle2(X,eps,d,distance,tol)

% LLE ALGORITHM
%
% function [Y,eigenvals,neighbors] = lle(X,K,d,tol)
%
% X = data as D x N matrix (D = dimensionality, N = #points)
% eps = epsilon ball
% d = embedding dimensionality
% tol = regularizer (defaults to 1e-4)
% Y = embedding as d x N matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PAIRWISE DISTANCES
[D,N] = size(X);

% NEIGHBORS

warning off    %% Next line causes an unnecessary warning, so turn it off
distance = distance ./(distance <= eps); 
distance = min(distance, Inf); 
distance = distance + diag(ones(N, 1)*Inf);
warning on

[rows, cols] = find(distance < Inf);
N2 = length(rows);
fprintf(1, 'average no. of neighbors: %f\n', N2 / N);

% RECONSTRUCTION WEIGHTS
if (nargin<5), tol=1e-4; end;
W = zeros(N2, 1);
for i=1:N
  ind = (rows == i);
  subN = sum(ind);
  subngb = cols(ind);
  z = X(:,subngb)-repmat(X(:,i),1,subN);
  C = z'*z;
  C = C + tol*trace(C)*eye(subN)/subN; % REGULARIZATION
  uW =  C \ ones(subN, 1);
  W(ind) = uW ./ sum(uW);
end;

% COST MATRIX
M = eye(N);
for i=1:N 
  ind = (rows == i);
  subngb = cols(ind);
  
  w = W(ind);
  j = subngb;
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
