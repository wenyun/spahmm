function [pred, prob, copyDist, copyRef, vsum] = spahmm_imputation(loc, readProb, Hapref, snpIdx, hapIdx, Param)

%spahmm_imputation    perform imputation using read likelihood and reference
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
%    

[vsum, posterior, dist, tranBase] = ...
  spahmm_fb_mex(loc, readProb, Hapref, snpIdx, hapIdx, Param);
   
prob = zeros(Param.nSnp, 2);
copyDist = zeros(Param.nSnp, 1);
copyRef = zeros(Param.nSnp, 1);
for i = 1:Param.nSnp
  snp = get_binary(Hapref.haps, snpIdx + i - 1, hapIdx);
  
  posteriorI = exp(posterior(:, i));
  prob(i, 1) = sum(posteriorI(snp == 0));
  prob(i, 2) = sum(posteriorI(snp == 1));
  
  copyDist(i) = sum(dist .* posteriorI);
  [C, I] = max(posteriorI);
  copyRef(i) = I;
end

[C, I] = max(prob, [], 2);
pred = I - 1;