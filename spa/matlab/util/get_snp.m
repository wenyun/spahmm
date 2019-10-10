function g = get_snp(id)


fid = floor((id - 1) / 1e4);
offset = rem(id - 1, 1e4) + 1;

fname = strcat('europe_mat_', int2str(fid), '.txt');

Gsub = load(fname);

ID = load('europe_idv.txt');

TXT = textread('europe_txt.txt', '%s');
COL = load('europe_col.txt');
LOC = load('europe_loc.txt');
GRP = load('europe_grp.txt');

[Gsub, ID, TXT, COL, LOC, GRP] = group(Gsub', ID, TXT, COL, LOC, GRP);
Gsub = Gsub';

g = Gsub(:, offset);
