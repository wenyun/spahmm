function X = standardize(X)

d = size(X, 1);

%% this is for mean 0 and variance 1

m = mean(X, 2);
s = std(X, 0, 2);
for i = 1:d
  if s(i) > 0
    X(i, :) = (X(i, :) - m(i)) / s(i);
  end
end

%% this is for price's method

% f = mean(X, 2) / 2;
% 
% for i = 1:d
%     X(i, :) = (X(i, :) - 2 * f(i)) / sqrt(2 * f(i) * (1 - f(i)));
% end


