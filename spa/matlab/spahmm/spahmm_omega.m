function omega = spahmm_omega(Param)

theta = sum(1 ./ (1:(Param.nHap - 1)));
omega = 0.5 * theta / (Param.nHap + theta);