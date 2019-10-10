function obj = spahmm_genotype_objective(loc, readProb1, readProb2, Hapref, index, Param)

[vsum1, dist1, tranBase1] = spahmm_fb_mex(loc, readProb1, Hapref, 1, index, Param);
[vsum2, dist2, tranBase2] = spahmm_fb_mex(loc, readProb2, Hapref, 1, index, Param);

obj = vsum1 + vsum2;