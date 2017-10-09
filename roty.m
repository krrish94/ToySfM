function Ry = roty(theta)
% ROTY  Creates a rotation matrix for a rotation about the Y-axis by an
% angle theta (in radians).

Ry = [cos(theta), 0, sin(theta); ...
    0, 1, 0; ...
    -sin(theta), 0, cos(theta), 0];

end