function [lambda, obj, obj0] = spahmm_lambda(loc, readProb, Hapref, snpIdx, hapIdx, Param)

%spahmm_lambda    estimate the optimal lambda by grid search
%  loc is the location for donor haplotype
%  readProb is nSnp * 2 matrix containing likelihood for allele 0 or 1
%  Hapref is compressed reference panel
%  snpIdx is the starting snp index, default is 1
%  hapIdx is the indices of reference haplotype
%  Param is parameters
%    nSnp is the number of impute snps
%    nHap is the number of used haplotypes
%    omega is the copying error
%    lambda is the spatial effect


lambda = grid_search(@(lambda) -likelihood_adapter(loc, readProb, Hapref, snpIdx, hapIdx, Param, lambda), ...
                     [0, 5], ...
                     5, ...
                     0.01);

obj = likelihood_adapter(loc, readProb, Hapref, snpIdx, hapIdx, Param, lambda);
obj0 = likelihood_adapter(loc, readProb, Hapref, snpIdx, hapIdx, Param, 0);
% recObj = [];
% recLambda = [];
% lambda = 0;
% for k = 1:Param.maxIter
%   % compute stochastic gradient 
%   [obj, grad] = gradient_adapter(loc, readProb, Hapref, snpIdx, hapIdx, Param, lambda);
%   
%   recLambda = [recLambda, lambda];
%   recObj = [recObj, obj];
%   fprintf(1, 'iter = %d : lambda = %.5f, likelihood = %.5f, ||g|| = %.10f\n', ...
%           k, lambda, obj, grad);
%   
%   % update
%   lambda = lambda + (1 / (k + Param.decay)) * grad;
% end
% 
% [C, I] = max(recObj);
% lambda = recLambda(I);
% obj = recObj(I);

function vsum = likelihood_adapter(loc, readProb, Hapref, snpIdx, hapIdx, Param, lambda)

Param.lambda = lambda;
[vsum, dist, tranBase] = spahmm_fb_mex(loc, readProb, Hapref, snpIdx, hapIdx, Param);


function [obj, grad] = gradient_adapter(loc, readProb, Hapref, snpIdx, hapIdx, Param, lambda)

Param.lambda = lambda;
[lamGrad, locGrad, vsum] = spahmm_grad_mex(loc, readProb, Hapref, snpIdx, hapIdx, Param);

grad = lamGrad;
obj = vsum;