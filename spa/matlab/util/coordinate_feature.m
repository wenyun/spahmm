function feature = coordinate_feature(coordiate)

%get_coordiante_feature    get features from coordiate

[nSample, nDim] = size(coordiate);

feature1 = coordiate;
feature2 = zeros(nSample, nDim*(nDim+1)/2);

for n = 1:nSample
    c = 1;
    for i = 1:nDim
        for j = 1:i
            feature2(n, c) = coordiate(n, i) * coordiate(n, j);
            c = c + 1;
        end
    end
end

feature = [feature1, feature2, ones(nSample, 1)];