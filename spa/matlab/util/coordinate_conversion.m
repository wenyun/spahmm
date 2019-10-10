function testLocation = ...
  coordinate_conversion(latitudeWeight, longtitudeWeight, testCoordinate)

%coordinate_conversion    convert unsupervised coordinates to geography
% latitudeWeight - column vector to convert to latitude
% longitudeWeight - column vector to convert to longitude
% testCoordinate    - M * L matrix
% testLocation  - M * 2 matrix indicating prediction of latitude and longitude

%% prediction

testFeature = coordinate_feature(testCoordinate);
testLocation = [testFeature * latitudeWeight, testFeature * longtitudeWeight];
testLocation = location_bound(testLocation);