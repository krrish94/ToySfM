function Rx = rotx(theta)
% ROTX  Creates a rotation matrix for a rotation about the X-axis by an
% angle theta (in radians).

Rx = [1, 0, 0; ...
    0, cos(theta), -sin(theta); ...
    0, sin(theta), cos(theta)];

end