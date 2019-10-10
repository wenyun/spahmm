function r = selection_corr(a, b)

ind = (a ~= 0) & (b ~= 0);

a = a(ind);
b = b(ind);

R = corrcoef(a, b);
r = R(1, 2);
