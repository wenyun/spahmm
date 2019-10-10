function n = bisec(A, b)

% this function returns the number of elements in a sorted A 
% less than or equal to b

L = length(A);

i = 1;
j = L;


if A(j) < b
    n = L;
elseif A(i) > b
    n = 0;
else
    while i + 1 < j 
        k = floor((i + j) / 2);

        if A(k) > b
            j = k;
        else
            i = k;
        end
    end

    n = i;
end