function dist = pos2dist(lag1,lon1,lag2,lon2)
% calculates sphereic geodesic distance for points farther apart,
% but ignores flattening of the earth:
% d =
% R_aver * acos(cos(lag1)cos(lag2)cos(lon1-lon2)+sin(lag1)sin(lag2))
% Output dist is in km.
% Returns -99999 if input argument(s) is/are incorrect.

if nargin < 4
    dist = -99999;
    disp('Number of input arguments error! distance = -99999');
    return;
end
if abs(lag1)>90 || abs(lag2)>90 || abs(lon1)>360 || abs(lon2)>360
    dist = -99999;
    disp('Degree(s) illegal! distance = -99999');
    return;
end
if lon1 < 0
    lon1 = lon1 + 360;
end
if lon2 < 0
    lon2 = lon2 + 360;
end

R_aver = 6374;
deg2rad = pi/180;
lag1 = lag1 * deg2rad;
lon1 = lon1 * deg2rad;
lag2 = lag2 * deg2rad;
lon2 = lon2 * deg2rad;
dist = R_aver * acos(cos(lag1)*cos(lag2)*cos(lon1-lon2) + sin(lag1)*sin(lag2));
dist = real(dist);