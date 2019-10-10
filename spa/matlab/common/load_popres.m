function [txt, loc, col, Genome, pruneList] = load_popres(pruneListFile)

%load_popres    load POPRES dataset
%  txt is country name
%  loc is country center for each individual
%  col is color assignemtn for each individual
%  Genome is a structure containing all haplotypes for individuals

ylab_init;

[Genome.haps, Genome.nSnp, Genome.nHap] = ...
  loadmatrix_binary(strcat(YLAB, ...
                           '/data/snp/popres/full/popres_phased_mat1.txt'));
if nargin == 0
  isPruned = 0;
else
  isPruned = 1;
end

if isPruned
  pruneList = logical(load(pruneListFile));  
  Genome.nSnp = sum(pruneList);
  Genome.haps = Genome.haps(pruneList, :);
end

txt = textread(strcat(YLAB, '/data/snp/popres/john/europe_txt.txt'), '%s');
col = load(strcat(YLAB, '/data/snp/popres/john/europe_col.txt'));
loc = load(strcat(YLAB, '/data/snp/popres/john/europe_loc.txt'));

idx = 1:1387;
% these two individuals are missing from POPRES data to European data
idx([1017, 1023]) = [];

txt = txt(idx, :);
col = col(idx, :);
loc = loc(idx, :);

loc = loc ./ 100; % make it within [0, 1] for numerical reason