function spahmm_mapping(x, tex, col)

[d, m] = size(x);

figure('Position',[100, 100, 500, 450]);
hold on;

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
  
  plot(xm(1), xm(2), 'o', 'MarkerSize', 30, 'MarkerFaceColor', cm, 'MarkerEdgeColor', cm);
  text(xm(1), xm(2), pop, 'FontSize', 15, 'HorizontalAlignment','center');
end

xmin = min(x(1, :));
xmax = max(x(1, :));
ymin = min(x(2, :));
ymax = max(x(2, :));
axis([xmin - 0.2 xmax + 0.2 ymin ymax]);

xlabel('Mapping Distance', 'FontSize', 30);
ylabel('Shared Haplotype Length (bp)', 'FontSize', 30);
set(gca, 'FontSize', 30);