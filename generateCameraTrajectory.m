function [Rs, ts] = generateCameraTrajectory(numViews, radius)
% GENERATECAMERATRAJECTORY  Generates a circular camera trajectory about
% the origin, assuming the radius of the trajectory is r, and that all
% views (a total of numViews) are equally spaced along the circle.
% The first point on the trajectory is assumed to lie on the X-axis.

% Angle increment at each time step
angleIncrement = 2*pi / numViews;

% Variables to hold the rotation matrices and translation vectors
Rs = zeros(numViews,3,3);
ts = zeros(numViews, 3);

% Generate each trajectory point
for i = 1:numViews
    
    % Current camera yaw angle (about Z-axis in world frame)
    curYawAngle = angleIncrement * (i-1);
    
    % Rotation matrix (first, rotate about the world Z-axis by pi/2 and
    % then about the obtained frame's X-axis by -pi/2, and then about the
    % thus obtained frame's Y-axis by -1*curYawAngle
    Rs(i,:,:) = rotz(pi/2) * rotx(-pi/2) * roty(-1*curYawAngle);
    
    % Translation vector
    % ts(i,:) = (-1 * squeeze(Rs(i,:,:))' * [radius*cos(curYawAngle); radius*sin(curYawAngle); 0])';
    ts(i,:) = [radius*cos(curYawAngle); radius*sin(curYawAngle); 0]';
    
    % Note that the R and t here can be used to transform a point Xw from
    % world coordinates to camera coordinates as X_cam = R^T*(Xw - t)
    
end


end
