function location = location_bound(location)

%location_bound    validate the latitude and longitude
%  location - M * 2 matrix, each row is (latitude, longitude)

ind = location(:, 1) > 90;
location(ind, 1) = 90;
ind = location(:, 1) < -90;
location(ind, 1) = -90;

ind = location(:, 2) > 360;
location(ind, 2) = 360;
ind = location(:, 2) < -360;
location(ind, 2) = -360;