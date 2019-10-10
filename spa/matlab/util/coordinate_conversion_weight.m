function [latitudeWeight, longtitudeWeight] = ...
  coordinate_conversion_weight(trainLocation, trainCoordinate)

%coordinate_conversion    convert unsupervised coordinates to geography
% trainCoordinate - N * L matrix from dimension reduction
% trainLocation - N * 2 matrix indicating the latitude and longitude
% latitudeWeight - column vector to convert to latitude
% longitudeWeight - column vector to convert to longitude

%% regression

% regression variables
trainFeature = coordinate_feature(trainCoordinate);
trainLatitude = trainLocation(:, 1);
trainLongitude = trainLocation(:, 2);

latitudeWeight = lin_regression(trainFeature, trainLatitude);
longtitudeWeight = lin_regression(trainFeature, trainLongitude);

function beta = lin_regression(x, y)
beta = (x' * x) \ x' * y;