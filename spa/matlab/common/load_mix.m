function [txt, loc, col, pid, Genome, pruneList] = load_mix(pruneListFile)

%load_mix    load POPRES dataset's admixed individual from mix_* files
%  txt is country name
%  loc is country center for each individual
%  col is color assignemtn for each individual
%  Genome is a nInd * nSnp genotype matrix

ylab_init;

Genome = load(strcat(YLAB, '/data/snp/popres/full/mix_mat.txt'));

if nargin == 0
  isPruned = 0;
else
  isPruned = 1;
end

if isPruned
  pruneList = logical(load(pruneListFile));  
  Genome = Genome(:, pruneList);
else
  pruneList = logical(ones(size(Genome, 2), 1) == 1);
end

nInd = size(Genome, 1);
txt = cell(nInd, 1);
col = cell(nInd, 1);
loc = cell(nInd, 1);
pid = zeros(nInd, 1);

fid = fopen(strcat(YLAB, '/data/snp/popres/full/mix_txt.txt'), 'r');
for i = 1:nInd
  line = fgetl(fid);
  token = regexp(line, ' ', 'split');    
  txt{i} = token;
end
fclose(fid);

fid = fopen(strcat(YLAB, '/data/snp/popres/full/mix_col.txt'), 'r');
for i = 1:nInd
  line = fgetl(fid);
  token = regexp(line, ' ', 'split');
  token = str2double(token);
  col{i} = reshape(token, 3, [])';
end
fclose(fid);

fid = fopen(strcat(YLAB, '/data/snp/popres/full/mix_loc.txt'), 'r');
for i = 1:nInd
  line = fgetl(fid);
  token = regexp(line, ' ', 'split'); 
  token = str2double(token) ./ 100; % make it within [0, 1] for numerical reason
  loc{i} = reshape(token, 2, [])';
end
fclose(fid);

fid = fopen(strcat(YLAB, '/data/snp/popres/full/mixlist.txt'), 'r');
for i = 1:nInd
  line = fgetl(fid);
  token = regexp(line, ' ', 'split'); 
  token = str2double(token);
  pid(i) = token(1);
end
fclose(fid);
