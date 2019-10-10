
function rt = rotate_angle(angle)

rt = [cos(pi*angle/180), -sin(pi*angle/180); sin(pi*angle/180), cos(pi*angle/180)];