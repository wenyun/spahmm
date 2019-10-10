function [txt, txt2, loc, col, col2, Genome, pruneList] = load_thousand(pruneListFile)

%load_popres    load 1000 Genome dataset
%  txt is population name
%  loc is country center for each individual
%  col is color assignemtn for each individual
%  Genome is a structure containing all haplotypes for individuals

ylab_init;

txt = textread(strcat(YLAB, '/data/snp/thousand/chr22_txt.txt'), '%s');
txt2 = textread(strcat(YLAB, '/data/snp/thousand/chr22_txt2.txt'), '%s');
col = load(strcat(YLAB, '/data/snp/thousand/chr22_col.txt'));
col2 = load(strcat(YLAB, '/data/snp/thousand/chr22_col2.txt'));
loc = load(strcat(YLAB, '/data/snp/thousand/chr22_loc.txt'));

if nargout > 5
  [Genome.haps, Genome.nSnp, Genome.nHap] = ...
    loadmatrix_binary(strcat(YLAB, '/data/snp/thousand/chr22_mat.txt'));

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
end