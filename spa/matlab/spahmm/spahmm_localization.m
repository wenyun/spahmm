function [loc, recLoc, recObj] = spahmm_localization(readProb1, readProb2, Hapref, hapIdx, Param)

%spahmm_localization    localize individual by EM algorithm
%  hap is the input haplotype, can be one or two 
%  Hapref is the compressed structure for haplotype reference
%  index is the including reference haplotype
%  Param is parameters for the algorithm


% grid search with forward-backward algorithm

if strcmp(Param.method, 'grid')
  [x, recLoc, recObj] = ...
    grid_search(@(loc) -spahmm_genotype_objective(loc, readProb1, readProb2, Hapref, hapIdx, Param), ...
                [-1.5, 2.5; -1.6, 1.3], ...
                Param.nGrid, ...
                Param.eps);
  loc = x';

elseif strcmp(Param.method, 'sgd')
  recObj = [];
  recLoc = [];
  loc = [0, 0];
  for k = 1:Param.maxIter
    % compute stochastic gradient 
    [obj, grad] = gradient_adapter(loc, readProb1, readProb2, Hapref, 1, hapIdx, Param);

    recLoc = [recLoc; loc];
    recObj = [recObj; obj];
    fprintf(1, 'iter = %d : loc = (%.5f %.5f), likelihood = %.5f, g = (%.5f %.5f)\n', ...
            k, loc(1), loc(2), obj, grad(1), grad(2));

    % update
    loc = loc + (1 / (k + Param.decay)) * grad';
  end
  [C, I] = max(recObj);
  loc = recLoc(I, :);
  obj = recObj(I);

end

function [obj, grad] = gradient_adapter(loc, readProb1, readProb2, Hapref, snpIdx, hapIdx, Param)

[lamGrad1, locGrad1, vsum1] = spahmm_grad_mex(loc, readProb1, Hapref, snpIdx, hapIdx, Param);
[lamGrad2, locGrad2, vsum2] = spahmm_grad_mex(loc, readProb2, Hapref, snpIdx, hapIdx, Param);

grad = locGrad1 + locGrad2;
obj = vsum1 + vsum2;