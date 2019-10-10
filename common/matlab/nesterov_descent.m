function [xmin, fmin] = nesterov_descent(objective, gradient, x0, param)

% this function performs the Nesterov's fast gradient descent algorithm
% objective(x) -  function handler for f(x)
% gradient(x) - function handler for computing gradient f'(x)
% x0 - initial point of x
% param - parameters for optimization

% xmin - the x minimizing f(x) locally
% fmin - the minimal f(x)

if nargin < 3
  error('objective and gradient function are required');
end

if nargin < 4  % default parameters
  param.maxiter = 1e3;  
  param.eps = 1e-3;
  param.verbose = 2;
  param.t = 1;
end

if nargin > 4
  error('at most 4 parameters');
end

x = x0;
y = x;

for k = 1:param.maxiter
  
  % compute direction by gradient  
  grad = gradient(y); 
  xn = y - param.t * grad;      
  y = xn + (k - 1) / (k + 2) * (xn - x);
      
  % test termination
  gnorm = norm(grad);
  
  if param.verbose >= 2
    obj = objective(y);
    fprintf(1, 'iter = %d : obj = %f, ||g|| = %.10f\n', k, obj, gnorm);
  elseif param.verbose >= 1  
    fprintf(1, 'iter = %d : ||g|| = %.10f\n', k, gnorm);
  end
  if gnorm < param.eps; break; end;
  
end

xmin = x;
fmin = obj;