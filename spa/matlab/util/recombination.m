function recomb = recombination(SNP, param)

% SNP - L * 3 matrix, columns: chrom, position, rs number

% recomb - (L-1) * 1 vector, recombination rate between every nearby SNPs
% recomb_rate is default to be 1.19 * 1e-8

recomb = param.recomb_rate * (SNP(2:param.L, 2) - SNP(1:(param.L - 1), 2));
recomb(SNP(2:param.L, 1) ~= SNP(1:(param.L - 1), 1)) = 0.5;

if param.window > 0
    recomb = [0.5; recomb];
    for i = 1:22
        idx = find(SNP(:, 1) == i);
        offset = floor((SNP(idx, 2) - SNP(idx(1), 2)) / param.window);
        in_idx = offset(2:end) == offset(1:end-1);
        idx = idx(2:end);
        recomb(idx(in_idx)) = 0;
    end
    recomb = recomb(2:end);
end
