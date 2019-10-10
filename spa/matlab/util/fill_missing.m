function X = fill_missing(X, mv)


[d, m] = size(X);


if nargin == 2
    for i = 1:d
        ind = (X(i, :) == -1);
        X(i, ind) = mv(i);
    end
    
else
    for i = 1:d
        ind = (X(i, :) == -1);
        rind = (X(i, :) ~= -1);
        avg = mean(X(i, rind));    
        X(i, ind) = avg;
    end
end
   
