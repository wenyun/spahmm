function [xmin, fmin] = gradient_descent(objective, gradient, x, param, projection)

% this function performs the gradient descent algorithm
% objective(x) -  function handler for f(x)
% gradient(x) - function handler for computing gradient f'(x)
% x0 - initial point of x
% param - parameters for optimization
% projection - function handler for gradient projection descent

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
    param.tinit = 1;
    param.t = 1;
    param.linesearch = 1;
end

if nargin < 5
    mode = 'no_projection';
elseif isa(projection, 'function_handle')
    mode = 'projection';
end 

if nargin > 5
    error('at most 5 parameters');
end

for iter = 1:param.maxiter
    
    % compute projected gradient 
    grad = gradient(x);            
    t = param.tinit;
    if strcmp(mode, 'projection')
        gp = 1 / t * (x - projection(x - t * grad));
    elseif strcmp(mode, 'no_projection')
        gp = grad;
    end
    
    % test termination
    gpnorm = norm(gp);
    lambda = gpnorm / length(x(:));        
    if lambda < param.eps; break; end;
    
    % update 
    if param.linesearch > 0
      % line search
      obj = objective(x);
      while objective(x - t * gp) > obj - param.alpha * t * grad' * gp
          t = param.beta * t;
          if strcmp(mode, 'projection')
              gp = 1 / t * (x - projection(x - t * grad));        
          end
      end
    else
      % constant step size
      t = param.t;      
      
      if param.verbose >= 2
        obj = objective(x);
      else
        obj = 0.0;
      end
    end
    
    if param.verbose >= 2; 
      fprintf(1, 'iter = %d : obj = %.5f, lambda = %.10f, forward t = %.5f\n', ...
                 iter, obj, lambda, t);
    end
    
    x = x - t * gp;
end

xmin = x;
fmin = objective(x);

fprintf(1, 'final: obj = %.5f, lambda = %.10f\n', fmin, lambda);