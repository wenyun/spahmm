function likelihood = haplotype_likelihood(count, eps)

%haplotype_likelihood    compute the haplotype likelihood from read counts
%  count is L * 2 matrix
%  eps is the sequencing error

correct = (1 - eps) .^ count;
error = eps .^ count;
total = sum(count, 2);

likelihood = zeros(size(count));
for i = 1:size(count, 1)
  likelihood(i, 1) = correct(i, 1) * error(i, 2) * nchoosek(total(i), count(i, 1));
  likelihood(i, 2) = correct(i, 2) * error(i, 1) * nchoosek(total(i), count(i, 2));
end
