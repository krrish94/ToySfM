function g = logMapSE3(G)
% LOGMAPSE3  Takes in the 4 x 4 SE(3) matrix G and returns its 6 x 1
% exponential coordinate vector g.

% Extract the rotation matrix and translation vector from G
R = G(1:3,1:3);
t = G(1:3,4);

% First compute the log map of the SO(3) part
theta = acos(0.5 * (trace(R)-1));
omega_cross = (theta / (2*sin(theta))) * (R - R');
omega = [-omega_cross(2,3); omega_cross(1,3); -omega_cross(1,2)];

% If theta is small, use an approximation of only the first 4 terms in the
% Taylor series expansion for sin and cos.
if theta < 1e-5
    fprintf('Using small angle approximation to compute the logarithm map.\n');
    V = eye(3) + ((0.5*(1-((thetasq/12) * ((1 - thetasq/30) * (1-thetasq/56))))) * omega_cross) + ...
        (((1/6)*((1-thetasq/20)*((1-thetasq/42)*(1-thetasq/72)))) * (omega_cross * omega_cross));
else
    V = eye(3) + ( ((1-cos(theta))/(theta*theta)) * omega_cross) + ...
        ( ((theta-sin(theta))/(theta*theta*theta)) * (omega_cross * omega_cross) );
end

% Logarithm map for translational part
u = inv(V) * t;

g = [omega; u];

end
