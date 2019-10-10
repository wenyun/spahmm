function [direction, rotation] = find_angle(loc, Y)

%%
% loc  - N * 2 matrix indicating the latitude and longitude of the sample
% grp  - N * 1 matrix indicating the locational group of the sample
% Y    - N * 2 matrix resulted from dimensional reduction

%%


reverse = zeros(2, 2, 4);
reverse(:, :, 1) = eye(2);
reverse(:, :, 2) = [-1, 0; 0 1];
reverse(:, :, 3) = [1, 0; 0 -1];
reverse(:, :, 4) = [-1, 0; 0 -1];


maxrho = -Inf;
maxr = 0;
maxtheta = 0;
for r = 1:4
    for theta = -30:50
        YY = rotate_angle(theta) * reverse(:, :, r) * Y';
        YY = YY';
                
        rho = sum_correlation(loc, YY);
        if rho > maxrho
            maxrho = rho;
            maxr = r;
            maxtheta = theta;
        end
    end
end

direction = reverse(:, :, maxr);
rotation = maxtheta;




function rho = sum_correlation(loc, YY)

rho = corrcoef(loc(:, 1), YY(:, 1)) + corrcoef(loc(:, 2), YY(:, 2));
