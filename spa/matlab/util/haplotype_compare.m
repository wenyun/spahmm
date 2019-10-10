function lens = haplotype_compare(hap1, hap2, snp)

%haplotype_compare    compute the number and average length of shared haplotype

comp = [1 == 0; (hap1 == hap2); 1 == 0] ;
snp = [snp(1, :); snp; snp(end, :)];

first = (comp(1:end-1) ~= comp(2:end)) & (comp(1:end-1) == 0);
last = (comp(1:end-1) ~= comp(2:end)) & (comp(1:end-1) == 1);

firstInd = logical(first);
lastInd = logical(last);
firstPos = snp(firstInd, 2);
lastPos = snp(lastInd, 2);

lens = lastPos - firstPos;