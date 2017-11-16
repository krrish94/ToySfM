function G = expMapSE3(g)
% EXPMAPSE3  Takes in the 6 x 1 exponential coordinates g of an se(3)
% element. Returns the 4 x 4 SE(3) matrix G.

% We assume that our se(3) exponential coordinates are omega and u, i.e.,
% the axis-angle rotation vector and the translation vector respectively.
omega = g(1:3);
u = g(4:6);

% Magnitude of omega (theta: the angle of rotation)
theta = norm(omega);

% Cross-product skew-symmetric matrix corresponding to omega
omega_cross = skew3(omega);

% If theta is small, use an approximation of only the first 4 terms in the
% Taylor series expansion for sin and cos.
if theta < 1e-5
    
    fprintf('Using small angle approximation to compute the exponential map.\n');
    
    % Exponential map for the rotation part
    thetasq = theta^2;
    R = eye(3) + ((1-((thetasq/6) * ((1-thetasq/20) * (1-thetasq/42)))) * omega_cross) + ...
        ((0.5*(1-((thetasq/12) * ((1 - thetasq/30) * (1-thetasq/56))))) * (omega_cross * omega_cross));
    
    % Exponential map for the translation part
    V = eye(3) + ((0.5*(1-((thetasq/12) * ((1 - thetasq/30) * (1-thetasq/56))))) * omega_cross) + ...
        (((1/6)*((1-thetasq/20)*((1-thetasq/42)*(1-thetasq/72)))) * (omega_cross * omega_cross));
    
    
else
    
    % Exponential map for the rotation part (use Rodrigues formula)
    R = eye(3) + ((sin(theta) / theta) * omega_cross) + ...
        ( ((1-cos(theta))/(theta*theta)) * (omega_cross * omega_cross) );
    
    % Exponential map for the translation part
    V = eye(3) + ( ((1-cos(theta))/(theta*theta)) * omega_cross) + ...
        ( ((theta-sin(theta))/(theta*theta*theta)) * (omega_cross * omega_cross) );
end


G = [R, V*u; zeros(1,3), 1];

end
