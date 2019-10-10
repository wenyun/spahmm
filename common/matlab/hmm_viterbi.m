function [X, Prob] = hmm_viterbi(emission, transition, N, L, isLog)

% this function performs the Viterbi algorithm
% emission(n, l) - function handler for P(x_l = n)
% transition(n1, n2, l) - function handler for P(x_l = n2| x_l-1 = n1), prior is transition(-, n2, 1)
% N - a vector L * 1 indicating the number of possible values for each node
% L - a scalar indicating the length of the HMM

% X - a vector L * 1 of the most likely variable sequences
% Prob - a scalar, the log of associating probability

if nargin == 4
    isLog = 0;
end

% forward
F = zeros(max(N), L);
U = zeros(max(N), L);
X = zeros(L, 1);

if isLog
    F(1:N(1), 1) = emission(1:N(1), 1) + transition(1, 1:N(1), 1);
else
    F(1:N(1), 1) = reallog(emission(1:N(1), 1)) + reallog(transition(1, 1:N(1), 1));
end

for i = 2:L
    for n2 = 1:N(i)
        if isLog
            [F(n2, i), U(n2, i)]= max(transition(1:N(i-1), n2, i) + F(:, i-1));            
        else
            [F(n2, i), U(n2, i)]= max(reallog(transition(1:N(i-1), n2, i)) + F(:, i-1));
        end
    end
    if isLog
        F(1:N(i), i) = F(1:N(i), i) + emission(1:N(i), i);
    else
        F(1:N(i), i) = F(1:N(i), i) + reallog(emission(1:N(i), i));
    end
end

[Prob, X(L)] = max(F(:, L));

for i = L-1:-1:1
    X(i) = U(X(i+1), i+1);
end

