function savematrix_matlab(dmat, fname, fmt, flg, fmt2)

if nargin < 4
    flg = 0;
end

[m, n] = size(dmat);
str = [];
for i = 1:n
    str = [str, fmt];
end
str = [str, '\n'];

if (flg == 1)
    fid = fopen(fname, 'a');
else
    fid = fopen(fname, 'w');    
end

for i = 1:m
    fprintf(fid, str, dmat(i,:));
end
fclose(fid);