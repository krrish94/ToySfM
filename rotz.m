function Rz = rotz(theta)
% ROTZ  Creates a rotation matrix for a rotation about the Z-axis by an
% angle theta (in radians).

Rz = [cos(theta), -sin(theta), 0; ...
    sin(theta), cos(theta), 0; ...
    0, 0, 1];

end