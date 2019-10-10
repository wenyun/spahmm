function snp = load_indian_snp(file)

fid = fopen(file);
info = textscan(fid, '%s %f %f %d %c %c');
fclose(fid);

snp = [info{2}, info{3}];