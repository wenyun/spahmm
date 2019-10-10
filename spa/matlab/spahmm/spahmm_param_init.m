function Param = spahmm_param_init(snp, Param)

%spahmm_param_init    initialize the parameters used in SPAHMM
%  snp is a nSnp * 3 matrix loaded from data, columns are chromosome, offset, rs
%  Param is a starting parameter, nSnp, recombRate are set.

% compute the recombination rate based on distance
nonRecomb = exp(- Param.rho * (snp(2:Param.nSnp, 2) - snp(1:(Param.nSnp - 1), 2)));
nonRecomb(snp(2:Param.nSnp, 1) ~= snp(1:(Param.nSnp - 1), 1)) = 0;
recomb = 1 - nonRecomb;

Param.recombProb = reallog(recomb);
Param.nonRecombProb = reallog(nonRecomb);