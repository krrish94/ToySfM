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
    
    % Rotation matrix (first, rotate about the world Y-axis by -pi/2 and
    % then about the obtained frame's Z-axis by pi/2, and then about the
    % thus obtained frame's Y-axis by -1*curYawAngle
    Rs(i,:,:) = rotx(-pi/2)*roty(-pi/2-curYawAngle);    % KM's version
    % Rs(i,:,:) = rotz(curYawAngle) * roty(-pi/2);      % SKT's version
    
    % Translation vector
    % ts(i,:) = (rotz(curYawAngle) * [radius;0;0])';
    ts(i,:) = [radius * cos(curYawAngle); radius * sin(curYawAngle); 0]';
    
    % Note that the R and t here can be used to transform a point Xw from
    % world coordinates to camera coordinates as X_cam = R' * (Xw - t)
    
end


end
