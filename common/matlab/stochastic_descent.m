function xmin = stochastic_descent(gradient, x, param)

% this function performs the gradient descent algorithm
% gradient(x) - function handler for computing gradient f'(x)
% x0 - initial point of x
% param - parameters for optimization

% xmin - the x minimizing f(x) locally

if nargin < 3  % default parameters
    param.maxiter = 1e3;    
    param.verbose = 2;
    param.decay = 1;
end

if nargin > 3
    error('at most 3 parameters');
end

for k = 1:param.maxiter
    
    % compute projected gradient 
    grad = gradient(x);    
    
    % update    
    if param.verbose >= 1 
      fprintf(1, 'iter = %d : ||g|| = %.10f', k, norm(grad));
    end
    if param.verbose >= 2
      fprintf(1, ', x = (%.5f %.5f)', x(1), x(2));
    end
    fprintf(1, '\n');
    
    x = x - (1 / (k + param.decay)) * grad;
end

xmin = x;
fprintf(1, 'final : ||g|| = %.10f\n', norm(grad));