function mapping(x, tex, col, xrev, yrev, angle, id1, id2, p1, p2)



% otxt = {'SP'; 'UK'; 'SZ'};
% ntxt = {'ES'; 'GB'; 'CH'};

% for i = 1:size(otxt, 1)
%    ind = strcmp(otxt(i, :), tex);
%    tex(ind, :) = ntxt(i, :);
% end


[d, m] = size(x);

figure('Position',[100, 100, 1000, 900]);


x = x([2, 1], :);

if nargin == 10
    p1 = p1([2, 1], :);
    p2 = p2([2, 1], :);
end

if xrev == 1
    x(1, :) = - x(1, :);
    if nargin == 10
        p1(1, :) = - p1(1, :);
        p2(1, :) = - p2(1, :);
    end
end

if yrev == 1
    x(2, :) = - x(2, :);
    if nargin == 10
        p1(2, :) = - p1(2, :);
        p2(2, :) = - p2(2, :);
    end
end

rt = rotation(angle);
x = rt * x;

if nargin == 10
    p1 = rt * p1;
    p2 = rt * p2;
end


xmin = min(x(1, :));
xmax = max(x(1, :));
ymin = min(x(2, :));
ymax = max(x(2, :));
xlen = xmax - xmin;
ylen = ymax - ymin;
if xlen > ylen
    ymin = ymin - (xlen - ylen) / 2;
    ymax = ymax + (xlen - ylen) / 2;
else
    xmin = xmin - (ylen - xlen) / 2;
    xmax = xmax + (ylen - xlen) / 2;    
end
xmean = (xmax + xmin) / 2;
ymean = (ymax + ymin) / 2;
center = [xmean; ymean];


pc1 = rt * [0; 1];
pc2 = rt * [1; 0];

hold on;
% plot(center(1) + [100 * pc1(1); -100 * pc1(1)], center(2) + [100 * pc1(2); -100 * pc1(2)], 'k-');
% plot(center(1) + [100 * pc2(1); -100 * pc2(1)], center(2) + [100 * pc2(2); -100 * pc2(2)], 'k-');
% 



for i = 1:m
    ctex = strcat(['\color[rgb]{' num2str(col(1, i)/255) ' ' num2str(col(2, i)/255) ' '...
        num2str(col(3, i)/255) '}'], tex(i));
    text(x(1, i), x(2, i), ctex, 'FontSize', 10, 'HorizontalAlignment','center');         
end 


utex = unique(tex);
p = length(utex);

for i = 1:p
    
    pop = utex(i);    
    pind = strcmp(tex, pop);
    pidx = find(pind);
    
    xm = mean(x(:, pind), 2);
    cm = col(:, pidx(1)) / 255;
    
    
    hold on;
    plot(xm(1), xm(2), 'o', 'MarkerSize',20, 'MarkerFaceColor', cm, 'MarkerEdgeColor', cm);
    text(xm(1), xm(2), pop, 'FontSize', 10, 'HorizontalAlignment','center');
end



if nargin == 10
    hold on;
    plot(p1(1), p1(2), 'd', 'MarkerEdgeColor', 'r', 'LineWidth', 5, 'MarkerSize', 20);
    hold on;
    plot(p2(1), p2(2), 'd', 'MarkerEdgeColor', 'r', 'LineWidth', 5, 'MarkerSize', 20);

    ctex1 = strcat(['\color[rgb]{' '0' ' ' '0' ' '...
        '0' '}'], tex(id1));   
    ctex2 = strcat(['\color[rgb]{' '0' ' ' '0' ' '...
        '0' '}'], tex(id2));     
    text(x(1, id1), x(2, id1), ctex1, 'FontSize', 30, 'HorizontalAlignment','center');
    text(x(1, id2), x(2, id2), ctex2, 'FontSize', 30, 'HorizontalAlignment','center');
end

axis([xmin xmax ymin ymax]);



hold on;


set(gca, 'Visible', 'off');



function rt = rotation(angle)

rt = [cos(pi*angle/180), -sin(pi*angle/180); sin(pi*angle/180), cos(pi*angle/180)];

