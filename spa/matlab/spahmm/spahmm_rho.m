function rho = spahmm_rho(Param)

Ne = 2e4;  % european effective size
rho = 4 * Ne * 1e-8 / Param.nHap;