function mappingm(x, tex, col, y, ycol, latlim, lonlim, fontSize)

if nargin == 7
  fontSize = 10;
end

position = 'center';

otxt = {'SP'; 'UK'; 'SZ'};
ntxt = {'ES'; 'GB'; 'CH'};

for i = 1:size(otxt, 1)
    ind = strcmp(otxt(i, :), tex);
    tex(ind, :) = ntxt(i, :);
end

axesm('mercator','MapLatLimit', latlim, 'MapLonLimit', lonlim, 'Frame', 'on');
geoshow('cntry02.shp', 'FaceColor', [0.98 1.0 0.98]);

[d, m] = size(x);
for i = 1:m
    ctex = strcat(['\color[rgb]{' num2str(col(1, i)/255) ' ' num2str(col(2, i)/255) ' '...
        num2str(col(3, i)/255) '}'], tex(i));
    hold on;
    if x(1, i) > latlim(1) + 2 && x(1, i) < latlim(2) - 2 && ...
       x(2, i) > lonlim(1) + 2 && x(2, i) < lonlim(2) - 2
      textm(x(1, i), x(2, i), ctex, 'FontSize', fontSize, 'HorizontalAlignment', position);     
    end
end 

[d, p] = size(y);
for i = 1:p
    hold on;
    plotm(y(1, i), y(2, i), 'o', 'MarkerSize',10, 'MarkerFaceColor', ycol(:, i) / 255, 'MarkerEdgeColor', 'k');
end

set(gcf,'Renderer','zbuffer');
set(gca, 'Visible', 'off');
set(gcf, 'InvertHardCopy', 'on');