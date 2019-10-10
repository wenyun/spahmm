function [H, n_snp, n_individual] = filter_matrix_binary(H, n_snp, n_individual, freq)

ind = ones(1, n_snp);
for i = 1:n_snp
    G = get_binary(H, i, 1:n_individual);
    
    if sum(G) < freq * n_individual || sum(G) > (1 - freq) * n_individual        
        ind(i) = 0;
    end
end

H = H(ind == 1, :);
n_snp = sum(ind);
