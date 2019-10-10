function [txt, txt2, loc, col, genome, pruneList] = load_hgdp_stanford(population, pruneListFile)

%load_hgdp_stanford    load hgdp_stanford data
%  txt is population name
%  txt2 is group name
%  loc is country center for each individual
%  col is color assignemtn for each individual
%  genome is a matrix of genotype

ylab_init;

fid = fopen(strcat(YLAB, '/data/snp/hgdp/hgdp_stanford.ind'));
info = textscan(fid, '%s %s %s %s %f %f %d %d %d');
fclose(fid);

txt = info{3};
txt2 = info{4};
loc = [info{5}, info{6}];
loc = loc / 100;
col = double([info{7}, info{8}, info{9}]);

if nargin < 1
  ind = ones(size(loc, 1), 1);
else
  ind = zeros(size(loc, 1), 1);
  for i = 1:size(loc, 1)
    ind(i) = sum(strcmp(txt2(i), population));
  end
end
ind = logical(ind);

txt = txt(ind);
txt2 = txt2(ind);
col = col(ind, :);
loc = loc(ind, :);

genome = load(strcat(YLAB, '/data/snp/indian/clean_hgdp_stanford.geno'));
if nargin == 1
  isPruned = 0;
else
  isPruned = 1;
end

if isPruned
  pruneList = logical(load(pruneListFile));  
  genome = genome(pruneList, :);
end

genome = genome(:, ind);