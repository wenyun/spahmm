function mhtplot(chrom, offset, score, left, right)


len = [247199719, 242751149, 199446827, 191263063, 180837866, ...
       170896993, 158821424, 146274826, 140442298, 135374737, ...
       134452384, 132289534, 114127980, 106360585, 100338915, ...
       88822254, 78654742, 76117153, 63806651, 62435965, ...
       46944323, 49528953];



   
plot(offset, score, 'k.');

mins = min(score);
maxs = max(score);
lens = maxs - mins;
if lens > 8
    step = 4;
elseif lens > 4
    step = 1;
elseif lens > 1
    step = 0.5;
elseif lens > 0.4
    step = 0.2;
elseif lens > 0.2
    step = 0.1;
elseif lens > 0.08
    step = 0.04;
elseif lens > 0.04
    step = 0.02;
elseif lens > 0.02
    step = 0.01;
end

if nargin == 3
    axis([0, len(chrom), mins, maxs]);
elseif nargin == 5
    axis([left, right, mins, maxs]);
end

set(gca,'YTick', 0:step:100);
