function ctable = compute_ctable(ind1, ind2)

ctable = zeros(3, 3);

ctable(1, 1) = sum((ind1 == 0) & (ind2 == 0));
ctable(1, 2) = sum((ind1 == 0) & (ind2 == 1));
ctable(1, 3) = sum((ind1 == 0) & (ind2 == 2));
ctable(2, 1) = sum((ind1 == 1) & (ind2 == 0));
ctable(2, 2) = sum((ind1 == 1) & (ind2 == 1));
ctable(2, 3) = sum((ind1 == 1) & (ind2 == 2));
ctable(3, 1) = sum((ind1 == 2) & (ind2 == 0));
ctable(3, 2) = sum((ind1 == 2) & (ind2 == 1));
ctable(3, 3) = sum((ind1 == 2) & (ind2 == 2));


