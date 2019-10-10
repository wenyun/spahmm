function [txt, loc, col, genome, pruneList] = load_indian(pruneListFile)

%load_indian    load indian dataset
%  txt is population name
%  loc is country center for each individual
%  col is color assignemtn for each individual
%  genome is a matrix of genotype

ylab_init;

genome = load(strcat(YLAB, '/data/snp/indian/clean_india.geno'));
if nargin == 0
  isPruned = 0;
else
  isPruned = 1;
end

if isPruned
  pruneList = logical(load(pruneListFile));  
  genome = genome(pruneList, :);
end

fid = fopen(strcat(YLAB, '/data/snp/indian/clean_india.ind'));
info = textscan(fid, '%s %s %s %d %d %d');
fclose(fid);

nInd = length(info{1});
loc = zeros(nInd, 2);

txt = info{3};
col = double([info{4}, info{5}, info{6}]);