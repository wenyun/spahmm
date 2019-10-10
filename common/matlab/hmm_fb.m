function [Prob, Vsum] = hmm_fb(emission, transition, N, L, isLog)

% this function performs the forward-backward algorithm
% emission(n, l) - function handler for P(x_l = n)
% transition(n1, n2, l) - function handler for P(x_l = n2| x_l-1 = n1), prior is transition(-, n2, 1)
% N - a vector L * 1 indicating the number of possible values for each node
% L - a scalar indicating the length of the HMM

% Prob - a max(N) * L matrix of posterior probablity for each node
% Vsum - the total probability by summing over forward algorithm

if nargin == 4
    isLog = 0;
end

% forward
F = zeros(max(N), L);
B = zeros(max(N), L);
Prob = zeros(max(N), L);

if isLog
    F(1:N(1), 1) = emission(1:N(1), 1) + transition(1, 1:N(1), 1);
else
    F(1:N(1), 1) = reallog(emission(1:N(1), 1)) + reallog(transition(1, 1:N(1), 1));
end

for i = 2:L
    if isLog
        F(1:N(i), i) = emission(1:N(i), i);
    else
        F(1:N(i), i) = reallog(emission(1:N(i), i));
    end
    for n2 = 1:N(i)
        if isLog
            F(n2, i) = F(n2, i) + ...
                       logsum_mex(transition(1:N(i-1), n2, i) + F(:, i-1));
        else
            F(n2, i) = F(n2, i) + ...
                       logsum_mex(reallog(transition(1:N(i-1), n2, i)) + F(:, i-1));
        end
    end
end

B(:, L) = 0;
for i = (L-1):-1:1
    for n1 = 1:N(i)
        if isLog
            B(n1, i) = logsum_mex(transition(n1, 1:N(i+1), i+1) + ...
                                  emission(1:N(i+1), i+1) + ...
                                  B(:, i+1));
        else
            B(n1, i) = logsum_mex(reallog(transition(n1, 1:N(i+1), i+1)) + ...
                                  reallog(emission(1:N(i+1), i+1)) + ...
                                  B(:, i+1));
        end
    end
end

Vsum = logsum_mex(F(:, L));
for i = 1:L
    Prob(:, i) = exp(F(:, i) + B(:, i) - Vsum);
    Prob(:, i) = Prob(:, i) / sum(Prob(:, i));  % normalize for numerical reason
end