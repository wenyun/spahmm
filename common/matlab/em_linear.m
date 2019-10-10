function x = em_linear(expectation, maximization, x0, Param)

%em_linear    EM algorithm for linear chain
%  [prob, vSum] = expectation(x) is function handler for expectation
%  x = maximization(prob, x) is function handler for maximization step
%  x0 is starting location for the optimization
%  Param is parameters
%    @required Param.maxIter
%    @required Param.eps

x = x0;
for iter = 1:Param.maxIter
  % E step (forward backward algorithm)
  [prob, vSum] = expectation(x);
  if iter > 1
    fprintf(1, 'EM iter %d : objective = %.5f, step = %.5f\n\n', iter-1, vSum, xstep);    
  end
  oldX = x;
  
  % M step
  x = maximization(prob, x);
  
  xstep = norm(x(:) - oldX(:)) / norm(oldX(:));
  
  if xstep < Param.eps
    [prob, vSum] = expectation(x);
    fprintf(1, 'EM final: objective = %.5f\n\n', vSum);
    break;
  end
end
