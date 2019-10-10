function [xmin, fmin] = newton_descent(objective, gradient, x, param)

% this function performs the gradient descent algorithm
% objective(x) -  function handler for f(x)
% gradient(x) - function handler for computing gradient f'(x), and Hessian f''(x)
%               must return a postive definite hessian matrix
% x0 - initial point of x
% param - parameters for optimization

% xmin - the x minimizing f(x) locally
% fmin - the minimal f(x)

if nargin < 3
  error('objective and gradient function are required');
end

if nargin < 4  % default parameters
  param.maxiter = 1e3;
  param.alpha = 0.01;
  param.beta = 0.5;
  param.eps = 1e-3;
  param.verbose = 2;
end

if nargin > 4
  error('at most 4 parameters');
end

for iter = 1:param.maxiter
  
  % compute direction by gradient
  obj = objective(x);
  [grad, hess] = gradient(x);
  xd = - hess \ grad;
  
  % test termination
  lambda = (- grad' * xd) ^ 0.5;
  if param.verbose >= 2; fprintf(1, 'iter = %d : obj = %f, lambda = %.10f\n', iter, obj, lambda); end
  if lambda < param.eps; break; end;
  
  % line search
  t = param.tinit;
  while objective(x + t * xd) > obj + param.alpha * t * grad' * xd
    t = param.beta * t;
  end
  x = x + t * xd;
  
  if t < 1e-8 && lambda < sqrt(param.eps)
    %fprintf('numerical problem or gradient/objective wrong\n');
    break;
  end
  
end

xmin = x;
fmin = obj;